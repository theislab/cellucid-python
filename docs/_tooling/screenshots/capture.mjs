#!/usr/bin/env node
// Deterministic documentation screenshots for the Cellucid web app.
//
// USAGE
//   node capture.mjs list
//   node capture.mjs shoot --all              --out DIR
//   node capture.mjs shoot <scenario-id> ...  --out DIR
//   node capture.mjs shoot <scenario-id> ...  --emit
//   node capture.mjs check  <scenario-id> ...
//
//   --out DIR   write the PNGs under DIR (required unless --emit)
//   --emit      write into docs/_static/screenshots/<topic>/<name>.png
//   --check     recapture and compare against the recorded digest; write
//               nothing and exit non-zero on any difference
//
// Ports: the application is served by the web repository's own browser-test
// server, so this honours the same two variables it does —
// CELLUCID_BROWSER_TEST_PORT and CELLUCID_BROWSER_TEST_SAMPLE_PORT — and adds
// CELLUCID_DOCS_EXPORTS_PORT (default: main port + 2) for the dataset host.
// Set them when another run already owns the default 4173/4174.
//
// ---------------------------------------------------------------------------
// THE FIVE THINGS THIS EXISTS TO GUARANTEE
//
// 1. DETERMINISM. The same scenario produces the same bytes. Fixed viewport,
//    fixed device scale factor, fixed colour scheme, fixed locale and time
//    zone, reduced motion, a seeded `Math.random` (the renderer's particle
//    system and parts of the analysis layer draw from it), Playwright's own
//    `animations: 'disabled'`, and — the part that actually does the work — a
//    settle loop that recaptures until two consecutive frames are byte-equal
//    before anything is written. A WebGL canvas is not settled because a
//    promise resolved; it is settled when it stops changing.
//
// 2. RESOLUTION AND SCALING. Capture happens at deviceScaleFactor 2, so every
//    image is supersampled relative to the CSS geometry the application lays
//    out. The emitted file is that native capture, downscaled only if it is
//    wider than MAX_EMIT_WIDTH.
//
//      emitted width = min(css width of the captured region x 2, 1440)
//
//    1440 is not arbitrary. The documentation theme lays the article column out
//    at `max-width: 60em` with `1rem` padding — 928 CSS px — and
//    `docs/_static/custom.css` gives every figure `max-width: 100%`, so an
//    image wider than 928 px is scaled down by the browser and gains sharpness
//    on a high-DPI display instead of overflowing. 1440 leaves 1.55x of
//    oversampling there while keeping files to hundreds of kilobytes, and it is
//    the width the eighteen existing full-window captures already use, so
//    full-window shots stay drop-in compatible with the published set.
//
//    The floor of that rule matters more than the ceiling. A 246 CSS px sidebar
//    panel used to be emitted at 246 px and rendered 1:1, which is why the
//    existing panel captures are unreadable on a high-DPI screen. The same
//    panel now emits at 492 px and renders at 492 px: twice the size on the
//    page and sharp. That is a deliberate, visible change to how panel figures
//    look, and it is the whole point.
//
//    `:width:` is not free to choose. `tests/test_docs_current_api_contract.py`
//    requires every `{figure}` to declare `:width:` equal to the PNG's
//    intrinsic pixel width, so the emitted width *is* the width the page
//    declares. `--emit` prints the exact directive block to paste.
//
// 3. CROPPING BY NAME. A crop is a DOM selector, never a pixel box. The
//    bounding rectangle is read at capture time, so a control that moves takes
//    its crop with it; the failure mode of a hand-tuned box — silently framing
//    the wrong thing after a layout change — cannot happen. `region` unions
//    several selectors for a control group that spans siblings.
//
//    The rectangle is read *again* after the settle loop and must still be the
//    same one. A rectangle read once is only true at the instant it was read,
//    and a subject that disappears while the frames are settling leaves the
//    crop framing the wallpaper — at exactly the right size, which is what
//    makes it convincing. See `assertCropStillFramed`.
//
// 4. THE POINTER. See `overlay.mjs`. A real `page.hover()` gives the control
//    its true hover styling; the drawn arrow tells the reader where the pointer
//    is. Both, or neither.
//
// 5. PROVENANCE. Every capture appends to `captures.json` beside this file:
//    scenario, dataset, viewport, device scale factor, theme, emitted size,
//    SHA-256, the web repository commit the application was served from, the
//    browser build, and whether the image settled. An image whose provenance is
//    not recorded cannot be re-verified later, and `check` is what re-verifies
//    it.
//
//    `web_repository_revision` describes the *served tree*, not `HEAD`. The
//    server reads files from disk, so a capture taken over uncommitted changes
//    shows an application no commit contains; that record carries a `-dirty`
//    suffix and a warning rather than a bare SHA that would read as clean. Only
//    a suffix-free value identifies a rebuildable application.
// ---------------------------------------------------------------------------

import { createHash } from 'node:crypto';
import { spawn } from 'node:child_process';
import { execFile } from 'node:child_process';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createRequire } from 'node:module';
import path from 'node:path';
import process from 'node:process';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { promisify } from 'node:util';

import { clearOverlay, OVERLAY_ROOT_ID, paintOverlay } from './overlay.mjs';
import { SCENARIOS } from './scenarios.mjs';
import { startExportsServer } from './serve-exports.mjs';

const run = promisify(execFile);

const HERE = path.dirname(fileURLToPath(import.meta.url));
const DOCS_ROOT = path.resolve(HERE, '..', '..');
const PUBLISHED_SCREENSHOTS = path.join(DOCS_ROOT, '_static', 'screenshots');
const MANIFEST_PATH = path.join(HERE, 'captures.json');

// The workspace layout is a sibling checkout per repository. Both are
// overridable so the tool is not welded to one machine's directory names.
const WORKSPACE = path.resolve(DOCS_ROOT, '..', '..');
const WEB_REPO =
  process.env.CELLUCID_WEB_REPO ?? path.join(WORKSPACE, 'cellucid');
const DATASETS_REPO =
  process.env.CELLUCID_DATASETS_REPO ??
  path.join(WORKSPACE, 'cellucid-datasets');

/** Every capture uses this CSS viewport; it is the app's documentation size. */
const VIEWPORT = Object.freeze({ width: 1440, height: 1000 });

/** Supersampling factor. Two is enough; three triples the settle cost. */
const DEVICE_SCALE_FACTOR = 2;

/** Widest file the documentation ever needs. See note 2 in the header. */
const MAX_EMIT_WIDTH = 1440;

/** Fixed seed for the page-side `Math.random` replacement. */
const RANDOM_SEED = 0x5eed1e5;

const SETTLE_DEFAULT = Object.freeze({ attempts: 14, gapMs: 220 });

// ---------------------------------------------------------------------------
// Environment
// ---------------------------------------------------------------------------

const require = createRequire(pathToFileURL(path.join(WEB_REPO, 'package.json')));

async function resolvePorts() {
  const module = await import(
    pathToFileURL(path.join(WEB_REPO, 'scripts', 'browser-test-ports.mjs'))
  );
  const { host, port, samplePort, origin } = module.resolveBrowserTestPorts();
  const exportsPort = Number(
    process.env.CELLUCID_DOCS_EXPORTS_PORT ?? String(port + 2)
  );
  if (
    !Number.isInteger(exportsPort) ||
    exportsPort < 1024 ||
    exportsPort > 65535 ||
    exportsPort === port ||
    exportsPort === samplePort
  ) {
    throw new RangeError(
      'CELLUCID_DOCS_EXPORTS_PORT must be a free port between 1024 and 65535 ' +
        `that differs from ${port} and ${samplePort}; received ${exportsPort}.`
    );
  }
  return { host, port, samplePort, exportsPort, origin };
}

/**
 * Identify the application that was actually served.
 *
 * `HEAD` alone is a lie whenever the working tree carries uncommitted changes:
 * the server reads files from disk, not from the object database, so a capture
 * taken over an edited tree would be stamped with a commit that does not
 * contain what the image shows — and `check` would later re-serve a *different*
 * application while believing it had reproduced the same one. The suffix makes
 * that visible instead. It is deliberately the shape `git describe --dirty`
 * uses, so the value still reads as a revision to a human and still starts with
 * the resolvable SHA.
 *
 * `git status --porcelain` is the whole tree, which is the right question:
 * `index.html`, every module under `assets/`, and the fixtures the scenarios
 * load are all served, so there is no subset that is safely ignorable.
 */
async function webRepositoryRevision() {
  let head;
  try {
    const { stdout } = await run('git', ['rev-parse', 'HEAD'], {
      cwd: WEB_REPO,
    });
    head = stdout.trim();
  } catch {
    return 'unknown';
  }
  try {
    const { stdout } = await run(
      'git',
      ['status', '--porcelain', '--untracked-files=normal'],
      { cwd: WEB_REPO }
    );
    return stdout.trim() === '' ? head : `${head}-dirty`;
  } catch {
    // The revision resolved but its cleanliness did not, so the one thing that
    // must not be recorded is a bare SHA implying a clean tree.
    return `${head}-unknown-worktree`;
  }
}

/**
 * True when a revision cannot be trusted to identify the application. A bare
 * commit SHA contains no `-`, so every qualifier this tool appends — and an
 * absent or malformed value in a hand-edited manifest — reads as provisional.
 */
function revisionIsProvisional(revision) {
  return typeof revision !== 'string' || revision === 'unknown' || revision.includes('-');
}

async function startApplicationServer(ports) {
  const child = spawn('node', ['scripts/serve-browser-tests.mjs'], {
    cwd: WEB_REPO,
    env: {
      ...process.env,
      CELLUCID_BROWSER_TEST_PORT: String(ports.port),
      CELLUCID_BROWSER_TEST_SAMPLE_PORT: String(ports.samplePort),
    },
    stdio: ['ignore', 'pipe', 'pipe'],
  });

  const failures = [];
  child.stderr.on('data', chunk => failures.push(String(chunk)));

  const deadline = Date.now() + 30_000;
  for (;;) {
    if (child.exitCode !== null) {
      throw new Error(
        `Application server exited before it was ready.\n${failures.join('')}`
      );
    }
    try {
      const response = await fetch(`${ports.origin}/index.html`, {
        method: 'HEAD',
      });
      if (response.ok) break;
    } catch {
      /* not listening yet */
    }
    if (Date.now() > deadline) {
      child.kill('SIGINT');
      throw new Error(
        `Application server did not answer on ${ports.origin} within 30s.\n` +
          failures.join('')
      );
    }
    await new Promise(resolve => setTimeout(resolve, 200));
  }

  return {
    close: () =>
      new Promise(resolve => {
        child.once('exit', () => resolve());
        child.kill('SIGINT');
        setTimeout(() => {
          child.kill('SIGKILL');
          resolve();
        }, 4000).unref();
      }),
  };
}

// ---------------------------------------------------------------------------
// Page-side helpers (serialised into the browser; they close over nothing)
// ---------------------------------------------------------------------------

/**
 * Replace `Math.random` with a seeded generator before any application code
 * runs. The renderer's particle overlay and parts of the analysis layer draw
 * from it, so an unseeded run produces a different picture every time.
 *
 * @param {number} seed
 */
function installSeededRandom(seed) {
  let state = seed >>> 0;
  Math.random = function seededRandom() {
    state = (state + 0x6d2b79f5) >>> 0;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Read the bounding rectangles of a list of selectors, in CSS pixels relative
 * to the viewport.
 *
 * @param {string[]} selectors
 * @returns {Array<{x: number, y: number, width: number, height: number}>}
 */
/**
 * Viewport-relative CSS rectangles for a list of selectors.
 *
 * Resolved through Playwright locators rather than `document.querySelector`, so
 * a crop, a ring and a cursor accept exactly the same selector vocabulary the
 * step list already accepts — including `:has-text(…)`, which is the only way
 * to name several of this application's controls, because many of its buttons
 * carry a shared utility class and no id.
 *
 * `boundingBox()` reports the same units and origin as
 * `getBoundingClientRect()` for the main frame, so the clip arithmetic below is
 * unchanged.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string[]} selectors
 */
async function readRects(page, selectors) {
  const rects = [];
  for (const selector of selectors) {
    const locator = page.locator(selector).first();
    if ((await locator.count()) === 0) {
      throw new Error(`Capture selector matched nothing: ${selector}`);
    }
    const rect = await locator.boundingBox();
    if (rect === null || rect.width <= 0 || rect.height <= 0) {
      throw new Error(
        `Capture selector has no painted area: ${selector} ` +
          `(${rect === null ? 'not visible' : `${rect.width}x${rect.height}`})`
      );
    }
    rects.push({ x: rect.x, y: rect.y, width: rect.width, height: rect.height });
  }
  return rects;
}

/**
 * Downscale a PNG with the browser's own high-quality resampler.
 *
 * Both directions travel as base64 because Playwright's evaluate boundary does
 * not carry typed arrays.
 *
 * @param {{data: string, width: number}} request
 * @returns {Promise<string>}
 */
async function resizePngInPage(request) {
  const binary = atob(request.data);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) {
    bytes[index] = binary.charCodeAt(index);
  }
  const source = await createImageBitmap(new Blob([bytes], { type: 'image/png' }));
  const height = Math.round((source.height * request.width) / source.width);
  const canvas = new OffscreenCanvas(request.width, height);
  const context = canvas.getContext('2d');
  context.imageSmoothingEnabled = true;
  context.imageSmoothingQuality = 'high';
  context.drawImage(source, 0, 0, request.width, height);
  source.close();
  const blob = await canvas.convertToBlob({ type: 'image/png' });
  const buffer = new Uint8Array(await blob.arrayBuffer());
  let out = '';
  for (let index = 0; index < buffer.length; index += 1) {
    out += String.fromCharCode(buffer[index]);
  }
  return btoa(out);
}

// ---------------------------------------------------------------------------
// Capture
// ---------------------------------------------------------------------------

function digest(bytes) {
  return createHash('sha256').update(bytes).digest('hex');
}

function pngDimensions(bytes) {
  if (
    bytes.length < 24 ||
    bytes.readUInt32BE(0) !== 0x89504e47 ||
    bytes.readUInt32BE(4) !== 0x0d0a1a0a
  ) {
    throw new Error('Capture did not produce a PNG.');
  }
  return { width: bytes.readUInt32BE(16), height: bytes.readUInt32BE(20) };
}

function anchorPoint(rect, at, dx, dy) {
  const points = {
    center: [rect.x + rect.width / 2, rect.y + rect.height / 2],
    left: [rect.x + Math.min(18, rect.width / 2), rect.y + rect.height / 2],
    right: [rect.x + rect.width - Math.min(18, rect.width / 2), rect.y + rect.height / 2],
    top: [rect.x + rect.width / 2, rect.y + Math.min(12, rect.height / 2)],
    bottom: [rect.x + rect.width / 2, rect.y + rect.height - Math.min(12, rect.height / 2)],
  };
  const point = points[at];
  if (point === undefined) {
    throw new Error(
      `Unknown cursor anchor "${at}"; expected one of ${Object.keys(points).join(', ')}.`
    );
  }
  return { x: point[0] + dx, y: point[1] + dy };
}

async function applyStep(page, step) {
  if (typeof step.run === 'function') return step.run(page);
  if ('waitFor' in step) {
    return page.waitForSelector(step.waitFor, { state: 'visible' });
  }
  if ('waitForText' in step) {
    return page.waitForFunction(
      ([selector, text]) =>
        document.querySelector(selector)?.textContent?.includes(text) === true,
      [step.waitForText, step.text]
    );
  }
  if ('waitForHidden' in step) {
    return page.waitForSelector(step.waitForHidden, { state: 'hidden' });
  }
  if ('click' in step) return page.click(step.click);
  if ('hover' in step) return page.hover(step.hover);
  if ('press' in step) return page.keyboard.press(step.press);
  if ('waitForStable' in step) {
    // The application rebuilds the field selectors when the observation
    // manifest arrives, which is *after* the point count settles. A selection
    // made in that window is silently discarded — `selectOption` reports
    // success and the control comes back holding a different value. Waiting
    // for the element's markup to stop changing is the honest fix; retrying the
    // selection until it sticks would only hide the race.
    const forMs = step.forMs ?? 700;
    return page.waitForFunction(
      async ([selector, quietMs]) => {
        const read = () => document.querySelector(selector)?.outerHTML ?? null;
        const first = read();
        if (first === null) return false;
        await new Promise(resolve => setTimeout(resolve, quietMs));
        return read() === first;
      },
      [step.waitForStable, forMs],
      { polling: forMs }
    );
  }
  if ('select' in step) {
    // `label` selects by the visible option text, which is what the
    // documentation names; the underlying option values are field indices.
    const chosen = await page.selectOption(
      step.select,
      step.label !== undefined ? { label: step.label } : step.value
    );
    // Assert the effect, not the call. A screenshot captioned "coloured by
    // clusters" that shows "Field: None" is the precise documentation defect
    // this tool exists to prevent, so a discarded selection must fail loudly.
    const settledValue = await page.inputValue(step.select);
    if (chosen.length !== 1 || settledValue !== chosen[0]) {
      throw new Error(
        `${step.select} did not keep the requested option: asked for ` +
          `${JSON.stringify(step.label ?? step.value)} (value ` +
          `${JSON.stringify(chosen)}), control now holds ` +
          `${JSON.stringify(settledValue)}. Add a ` +
          `{ waitForStable: '${step.select}' } step before this one.`
      );
    }
    return undefined;
  }
  if ('fill' in step) return page.fill(step.fill, step.value);
  if ('scrollIntoView' in step) {
    return page.locator(step.scrollIntoView).scrollIntoViewIfNeeded();
  }
  if ('openSection' in step || 'closeSection' in step) {
    // The sidebar accordions are `<details class="accordion-section">`.
    // Setting `.open` is idempotent; clicking the `<summary>` would toggle,
    // so a scenario that ran twice would leave the panel in the wrong state.
    const selector = step.openSection ?? step.closeSection;
    return page.evaluate(
      ([target, open]) => {
        const section = document.querySelector(target);
        if (section === null) {
          throw new Error(`Accordion section not found: ${target}`);
        }
        if (section instanceof HTMLDetailsElement) {
          section.open = open;
          return;
        }
        throw new Error(
          `${target} is not a <details> accordion; open it with an explicit ` +
            'click step instead.'
        );
      },
      [selector, 'openSection' in step]
    );
  }
  if ('wait' in step) {
    // Present for scenes with no observable ready signal. Every use is recorded
    // in the manifest, because a timeout is an admission that the scenario has
    // no deterministic anchor yet.
    return page.waitForTimeout(step.wait);
  }
  throw new Error(`Unrecognised step: ${JSON.stringify(step)}`);
}

function clipFromRects(rects, padding, selectors) {
  const left = Math.min(...rects.map(rect => rect.x)) - padding;
  const top = Math.min(...rects.map(rect => rect.y)) - padding;
  const right = Math.max(...rects.map(rect => rect.x + rect.width)) + padding;
  const bottom = Math.max(...rects.map(rect => rect.y + rect.height)) + padding;

  // Clamp into the viewport: a clip that leaves it captures transparent pixels.
  const x = Math.max(0, Math.round(left));
  const y = Math.max(0, Math.round(top));
  const clip = {
    x,
    y,
    width: Math.min(VIEWPORT.width - x, Math.round(right) - x),
    height: Math.min(VIEWPORT.height - y, Math.round(bottom) - y),
  };
  if (clip.width <= 0 || clip.height <= 0) {
    throw new Error(
      `Crop for ${JSON.stringify(selectors)} resolved to an empty region.`
    );
  }
  return clip;
}

function cropSelectors(shot) {
  return shot.mode === 'element' ? [shot.selector] : [...shot.selectors];
}

async function resolveClip(page, shot) {
  if (shot.mode === 'window') {
    return { x: 0, y: 0, width: VIEWPORT.width, height: VIEWPORT.height };
  }
  const selectors = cropSelectors(shot);
  for (const selector of selectors) {
    await page.waitForSelector(selector, { state: 'visible' });
  }
  return clipFromRects(
    await readRects(page, selectors),
    shot.padding ?? 0,
    selectors
  );
}

/**
 * The subject has to still be in frame when the shutter closes.
 *
 * The clip is a rectangle, read once, before the settle loop. The loop then
 * runs for as long as it takes two frames to agree — and a crop named after an
 * element that disappears in that window keeps its rectangle and loses its
 * subject, so what gets written is a correctly-sized photograph of whatever was
 * behind it. Worse, the loop *converges on that state*: a transient element is
 * the thing that was moving, so the first pair of identical frames is usually
 * the pair taken after it went away, and the image is recorded as `settled`.
 *
 * Two documentation images shipped that way — `data_loading/fail-h5ad-wrong-file`
 * and `data_loading/fail-missing-umap-embedding`, both crops of
 * `#notification-center` around a failure the notification center auto-dismisses
 * four seconds later. Each was 776 px of background scatter with no text in it,
 * carrying an `:alt:` describing a message that was never in the frame.
 *
 * Re-reading the rectangle afterwards is what closes that hole. It costs one
 * `boundingBox()` per crop and turns a plausible-looking wrong image into a
 * failed run.
 *
 * @param {import('@playwright/test').Page} page
 * @param {Object} shot
 * @param {{x: number, y: number, width: number, height: number}} clip
 */
async function assertCropStillFramed(page, shot, clip) {
  if (shot.mode === 'window') return;
  const selectors = cropSelectors(shot);
  let after;
  try {
    after = clipFromRects(
      await readRects(page, selectors),
      shot.padding ?? 0,
      selectors
    );
  } catch (error) {
    throw new Error(
      `${JSON.stringify(selectors)} left the frame while the image settled, ` +
        'so the capture is of whatever was behind it. If the subject is ' +
        'transient — a notification, a toast, a progress row — the scenario ' +
        'has to settle inside its lifetime. ' +
        error.message,
      { cause: error }
    );
  }
  // Two pixels of tolerance: sub-pixel reflow is normal, a collapsing or
  // relocating subject is not.
  const drift = ['x', 'y', 'width', 'height'].filter(
    key => Math.abs(after[key] - clip[key]) > 2
  );
  if (drift.length > 0) {
    throw new Error(
      `${JSON.stringify(selectors)} moved while the image settled ` +
        `(${drift.join(', ')} changed: ` +
        `${clip.width}x${clip.height}+${clip.x}+${clip.y} became ` +
        `${after.width}x${after.height}+${after.x}+${after.y}), so the crop ` +
        'no longer frames what it was measured against.'
    );
  }
}

async function paint(page, scenario, clip) {
  const rings = [];
  for (const ring of scenario.rings ?? []) {
    const [rect] = await readRects(page, [ring.on]);
    rings.push({ ...rect, label: ring.label ?? null });
  }

  let cursor = null;
  if (scenario.cursor != null) {
    const { on, at = 'center', dx = 0, dy = 0, state = 'hover' } =
      scenario.cursor;
    const [rect] = await readRects(page, [on]);
    const point = anchorPoint(rect, at, dx, dy);
    // The real hover is what gives the control its true styling; the drawing
    // only says where the pointer is.
    await page.mouse.move(point.x, point.y);
    if (state === 'press') await page.mouse.down();
    cursor = { x: point.x, y: point.y, state };
  }

  if (cursor === null && rings.length === 0) return () => Promise.resolve();

  await page.evaluate(paintOverlay, {
    rootId: OVERLAY_ROOT_ID,
    cursor,
    rings,
  });

  return async () => {
    if (scenario.cursor?.state === 'press') await page.mouse.up();
    await page.evaluate(clearOverlay, OVERLAY_ROOT_ID);
  };
}

async function settledScreenshot(page, clip, settle) {
  const options = {
    clip,
    scale: 'device',
    animations: 'disabled',
    caret: 'hide',
    type: 'png',
  };
  if (settle?.frames !== undefined) {
    // An animated scene never stops changing. Advancing a fixed number of
    // frames from a seeded start is the only reproducible anchor available, and
    // the manifest records that this image is frame-pinned rather than settled.
    await page.evaluate(async frames => {
      for (let index = 0; index < frames; index += 1) {
        await new Promise(resolve => requestAnimationFrame(() => resolve()));
      }
    }, settle.frames);
    return { bytes: await page.screenshot(options), settled: false };
  }

  const attempts = settle?.attempts ?? SETTLE_DEFAULT.attempts;
  const gapMs = settle?.gapMs ?? SETTLE_DEFAULT.gapMs;
  let previous = await page.screenshot(options);
  let previousDigest = digest(previous);
  for (let attempt = 1; attempt < attempts; attempt += 1) {
    await page.waitForTimeout(gapMs);
    const current = await page.screenshot(options);
    const currentDigest = digest(current);
    if (currentDigest === previousDigest) {
      return { bytes: current, settled: true };
    }
    previous = current;
    previousDigest = currentDigest;
  }
  return { bytes: previous, settled: false };
}

async function emitWidthFor(scenario, clip) {
  const native = clip.width * DEVICE_SCALE_FACTOR;
  if (scenario.emitWidth !== undefined) {
    if (scenario.emitWidth > native) {
      throw new Error(
        `${scenario.id}: emitWidth ${scenario.emitWidth} exceeds the native ` +
          `capture width ${native}. Upscaling a screenshot is never correct.`
      );
    }
    return scenario.emitWidth;
  }
  return Math.min(native, MAX_EMIT_WIDTH);
}

async function captureScenario(browser, scenario, urls) {
  // A fresh context per scenario, not a fresh page. The application persists
  // theme, sidebar geometry and onboarding state, so a shared context makes the
  // result depend on what ran before it — which is the exact property this tool
  // exists to eliminate.
  const context = await browser.newContext({
    viewport: VIEWPORT,
    deviceScaleFactor: DEVICE_SCALE_FACTOR,
    colorScheme: scenario.theme === 'dark' ? 'dark' : 'light',
    reducedMotion: 'reduce',
    locale: 'en-US',
    timezoneId: 'UTC',
  });
  context.setDefaultTimeout(45_000);
  await context.addInitScript(installSeededRandom, RANDOM_SEED);

  const page = await context.newPage();
  const consoleErrors = [];
  page.on('pageerror', error => consoleErrors.push(String(error.message)));

  try {
    // Anything that has to exist before the first byte of application code
    // runs — `page.route` interceptors, `addInitScript` seeds — belongs here.
    // A step cannot do it: steps run after navigation, by which time the
    // application has already made the requests a mock exists to answer. The
    // community annotation scenarios use this to stand up a synthetic GitHub
    // Worker origin, so the whole feature is reachable without an account and
    // without touching a real repository.
    if (typeof scenario.setup === 'function') {
      await scenario.setup({ context, page, urls });
    }

    const url = new URL('/index.html', urls.appOrigin);
    url.searchParams.set('exportsBaseUrl', `${urls.exportsOrigin}/`);
    if (scenario.dataset != null) url.searchParams.set('dataset', scenario.dataset);
    for (const [key, value] of Object.entries(scenario.query ?? {})) {
      url.searchParams.set(key, value);
    }
    await page.goto(url.toString(), { waitUntil: 'domcontentloaded' });

    // The welcome overlay is part of first-run, and a scenario that wants it
    // says so; every other scenario dismisses it the way the browser suite
    // does. It must be *waited for* before Escape is sent: the overlay paints
    // before its key handler is wired, and an early Escape is swallowed, which
    // leaves the overlay covering every later capture.
    const welcome = page.locator('#welcome-modal');
    await welcome.waitFor({ state: 'visible' });
    if (scenario.keepWelcome !== true) {
      await page.keyboard.press('Escape');
      await welcome.waitFor({ state: 'hidden' });
    }

    const steps = scenario.steps ?? [];
    for (let index = 0; index < steps.length; index += 1) {
      try {
        await applyStep(page, steps[index]);
      } catch (error) {
        // Without this the failure reads as an anonymous 45-second timeout and
        // says nothing about which of a dozen steps never came true.
        throw new Error(
          `${scenario.id}: step ${index + 1} of ${steps.length} ` +
            `${JSON.stringify(steps[index], (key, value) =>
              typeof value === 'function' ? '[run]' : value
            )} failed: ${error.message}`,
          { cause: error }
        );
      }
    }

    const clip = await resolveClip(page, scenario.shot);
    const undoPaint = await paint(page, scenario, clip);
    const { bytes, settled } = await settledScreenshot(
      page,
      clip,
      scenario.settle
    );
    await assertCropStillFramed(page, scenario.shot, clip);
    await undoPaint();

    const targetWidth = await emitWidthFor(scenario, clip);
    const native = pngDimensions(bytes);
    let finalBytes = bytes;
    if (native.width !== targetWidth) {
      const encoded = await page.evaluate(resizePngInPage, {
        data: bytes.toString('base64'),
        width: targetWidth,
      });
      finalBytes = Buffer.from(encoded, 'base64');
    }
    const emitted = pngDimensions(finalBytes);

    return {
      bytes: finalBytes,
      record: {
        id: scenario.id,
        topic: scenario.topic,
        name: scenario.name,
        describes: scenario.describes,
        dataset: scenario.dataset ?? null,
        theme: scenario.theme ?? 'light',
        viewport: `${VIEWPORT.width}x${VIEWPORT.height}`,
        device_scale_factor: DEVICE_SCALE_FACTOR,
        crop: scenario.shot.mode === 'window' ? 'full window' : describeCrop(scenario.shot),
        captured_css: `${clip.width}x${clip.height}`,
        native_px: `${native.width}x${native.height}`,
        emitted_px: `${emitted.width}x${emitted.height}`,
        cursor: scenario.cursor ? `${scenario.cursor.state ?? 'hover'} on ${scenario.cursor.on}` : null,
        rings: (scenario.rings ?? []).map(ring => ring.on),
        settled,
        uses_timeout: (scenario.steps ?? []).some(step => 'wait' in step),
        page_errors: consoleErrors,
        sha256: digest(finalBytes),
      },
    };
  } finally {
    await context.close();
  }
}

function describeCrop(shot) {
  const selectors =
    shot.mode === 'element' ? [shot.selector] : [...shot.selectors];
  const padding = shot.padding ? ` +${shot.padding}px` : '';
  return `${selectors.join(' ∪ ')}${padding}`;
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

function figureBlock(record, relativePrefix) {
  return [
    '```{figure} ' +
      `${relativePrefix}_static/screenshots/${record.topic}/${record.name}.png`,
    `:alt: ${record.describes}`,
    `:width: ${record.emitted_px.split('x')[0]}px`,
    '',
    'CAPTION — replace this line with one sentence saying what the reader',
    'should take from the image.',
    '```',
  ].join('\n');
}

async function loadManifest() {
  try {
    return JSON.parse(await readFile(MANIFEST_PATH, 'utf8'));
  } catch {
    return { version: 1, captures: {} };
  }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

function parseArguments(argv) {
  const command = argv[0] ?? 'list';
  const ids = [];
  const flags = { emit: false, check: false, all: false, out: null };
  for (let index = 1; index < argv.length; index += 1) {
    const token = argv[index];
    if (token === '--emit') flags.emit = true;
    else if (token === '--check') flags.check = true;
    else if (token === '--all') flags.all = true;
    else if (token === '--out') {
      index += 1;
      flags.out = argv[index];
    } else if (token.startsWith('--')) {
      throw new Error(`Unknown flag ${token}`);
    } else ids.push(token);
  }
  return { command, ids, flags };
}

async function main() {
  const { command, ids, flags } = parseArguments(process.argv.slice(2));

  if (command === 'list') {
    for (const scenario of SCENARIOS) {
      process.stdout.write(
        `${scenario.id.padEnd(34)} ${scenario.topic}/${scenario.name}.png\n` +
          `${' '.repeat(35)}${scenario.describes}\n`
      );
    }
    return;
  }
  if (command !== 'shoot' && command !== 'check') {
    throw new Error(`Unknown command "${command}"; expected list, shoot, check.`);
  }
  const checking = command === 'check' || flags.check;

  const selected = flags.all
    ? SCENARIOS
    : ids.map(id => {
        const scenario = SCENARIOS.find(entry => entry.id === id);
        if (scenario === undefined) {
          throw new Error(
            `No scenario "${id}". Run \`node capture.mjs list\` for the registry.`
          );
        }
        return scenario;
      });
  if (selected.length === 0) {
    throw new Error('Name at least one scenario, or pass --all.');
  }
  if (!checking && !flags.emit && flags.out === null) {
    throw new Error(
      'Pass --out DIR to write somewhere safe, or --emit to write into ' +
        'docs/_static/screenshots/. An unreferenced PNG under ' +
        '_static/screenshots/ fails the Python contract suite, so --emit is ' +
        'only correct when the {figure} block lands in the same change.'
    );
  }

  const ports = await resolvePorts();
  const { chromium } = require('playwright');

  const exportsRoot = path.join(DATASETS_REPO, 'exports');
  const exportsServer = await startExportsServer(exportsRoot, ports.exportsPort);
  const appServer = await startApplicationServer(ports);
  const browser = await chromium.launch();

  const manifest = await loadManifest();
  const revision = await webRepositoryRevision();
  const capturedAt = new Date().toISOString().replace(/\.\d+Z$/, 'Z');
  const failures = [];

  try {
    // The qualifier must survive the abbreviation: `slice(0, 12)` on
    // `<sha>-dirty` would print exactly the string a clean tree prints.
    const [revisionSha, ...revisionQualifier] = revision.split('-');
    const shortRevision = [revisionSha.slice(0, 12), ...revisionQualifier].join(
      '-'
    );
    process.stdout.write(
      `app ${ports.origin}  exports ${exportsServer.origin}  ` +
        `chromium ${browser.version()}  web@${shortRevision}\n\n`
    );
    if (revisionIsProvisional(revision)) {
      process.stdout.write(
        `  WARNING: the web repository working tree is not clean, so ` +
          `${checking ? 'this check re-served' : 'these images show'} an ` +
          'application that no commit contains. Commit or stash the web ' +
          'repository before a capture whose provenance has to mean ' +
          'something.\n\n'
      );
    }

    for (const scenario of selected) {
      const started = Date.now();
      let captured;
      try {
        captured = await captureScenario(browser, scenario, {
          appOrigin: ports.origin,
          exportsOrigin: exportsServer.origin,
        });
      } catch (error) {
        // One broken scenario must not cost the whole sweep. The failure is
        // reported at the end and the exit status is non-zero.
        failures.push(`${scenario.id}: ${error.message}`);
        process.stdout.write(`  ${scenario.id}  FAILED: ${error.message}\n`);
        continue;
      }
      const { bytes, record } = captured;
      record.web_repository_revision = revision;
      record.browser = `chromium ${browser.version()}`;
      record.captured_at = capturedAt;

      const previous = manifest.captures[scenario.id];
      const changed = previous !== undefined && previous.sha256 !== record.sha256;

      if (checking) {
        if (previous === undefined) {
          failures.push(`${scenario.id}: no recorded capture to check against`);
        } else if (changed) {
          failures.push(
            `${scenario.id}: recaptured bytes differ from the recorded image ` +
              `(${previous.sha256.slice(0, 12)} -> ${record.sha256.slice(0, 12)})`
          );
        } else if (revisionIsProvisional(previous.web_repository_revision)) {
          // Bytes agreeing proves the image is reproducible; it does not
          // recover which application produced it. Say so rather than letting
          // a clean run imply provenance the record never had.
          process.stdout.write(
            `  ${scenario.id}  NOTE: recorded against ` +
              `${previous.web_repository_revision}, so its provenance names no ` +
              'commit. Recapture from a clean tree to fix the record.\n'
          );
        }
      } else {
        const directory = flags.emit
          ? path.join(PUBLISHED_SCREENSHOTS, record.topic)
          : path.join(path.resolve(flags.out), record.topic);
        await mkdir(directory, { recursive: true });
        const file = path.join(directory, `${record.name}.png`);
        await writeFile(file, bytes);
        if (flags.emit) {
          // Provenance describes the *published* bytes, so only a run that
          // writes them may record it. A scratch iteration into `--out` used to
          // overwrite the record too, leaving `captures.json` describing a file
          // nobody can see and `check` comparing the published image against a
          // digest it was never supposed to have.
          manifest.captures[scenario.id] = record;
        }
        process.stdout.write(`${file}\n`);
        if (flags.emit) {
          // `test_every_documented_screenshot_uses_its_intrinsic_responsive_width`
          // asserts `referenced == screenshots`: every PNG under
          // `_static/screenshots/` must be named by a `{figure}` block. An
          // emitted image with no reference fails `python -m pytest` for the
          // whole repository until the block lands, so the two are one change.
          process.stdout.write(
            '  NOTE: this PNG now fails python -m pytest until the {figure} ' +
              'block below is pasted into a page.\n'
          );
        }
      }

      const seconds = ((Date.now() - started) / 1000).toFixed(1);
      const flagText = [
        record.settled ? 'settled' : 'NOT SETTLED',
        record.page_errors.length > 0
          ? `${record.page_errors.length} page error(s)`
          : null,
        changed && !checking ? 'bytes changed' : null,
      ]
        .filter(Boolean)
        .join(', ');
      process.stdout.write(
        `  ${scenario.id}  ${record.captured_css} css -> ` +
          `${record.native_px} native -> ${record.emitted_px} emitted  ` +
          `[${flagText}]  ${seconds}s\n`
      );
      if (record.page_errors.length > 0) {
        for (const error of record.page_errors) {
          process.stdout.write(`      page error: ${error}\n`);
        }
      }
      if (!checking) {
        process.stdout.write(
          `\n${figureBlock(record, '../../../')}\n\n`
        );
      }
    }

    // Only an `--emit` run changed anything worth persisting; `check` records
    // nothing and a `--out` iteration is scratch work. Rewriting the file
    // anyway would churn shared state that several units edit at once.
    if (flags.emit && !checking) {
      manifest.version = 1;
      await writeFile(
        MANIFEST_PATH,
        `${JSON.stringify(manifest, null, 2)}\n`,
        'utf8'
      );
    }
  } finally {
    await browser.close();
    await appServer.close();
    await exportsServer.close();
  }

  if (failures.length > 0) {
    for (const failure of failures) process.stderr.write(`${failure}\n`);
    process.exitCode = 1;
  }
}

await main();
