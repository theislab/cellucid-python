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
//   --out DIR   write the images under DIR (required unless --emit)
//   --emit      write into docs/_static/screenshots/<topic>/<name>.<png|webp>
//   --check     recapture and compare against the recorded digest; write
//               nothing and exit non-zero on any difference
//
// Ports: the application is served by the web repository's own browser-test
// server, so this honours the same two variables it does —
// CELLUCID_BROWSER_TEST_PORT and CELLUCID_BROWSER_TEST_SAMPLE_PORT — and adds
// CELLUCID_DOCS_CATALOG_PORT (default: main port + 2) for the sample catalog
// host (`cellucid-datasets/anndata/`, read through `?sampleCatalog=`).
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
import { execFile, execFileSync } from 'node:child_process';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createRequire } from 'node:module';
import path from 'node:path';
import process from 'node:process';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { promisify } from 'node:util';

import { clearOverlay, OVERLAY_ROOT_ID, paintOverlay } from './overlay.mjs';
import { SCENARIOS } from './scenarios.mjs';
import { startCatalogServer } from './serve-catalog.mjs';

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

/**
 * The sample catalog folder the scenarios load: `datasets.json` (schema 2) and
 * the AnnData stores it lists, as published at
 * https://theislab.github.io/cellucid-datasets/anndata/. Refused by name when
 * the datasets checkout does not hold it, rather than letting every scenario
 * time out on an empty dataset list.
 *
 * @returns {Promise<string>}
 */
async function requireSampleCatalog() {
  const root = path.join(DATASETS_REPO, 'anndata');
  let catalog;
  try {
    catalog = JSON.parse(await readFile(path.join(root, 'datasets.json'), 'utf8'));
  } catch (error) {
    throw new Error(
      `No sample catalog at ${root}/datasets.json (${error.code ?? error.message}). ` +
        'The scenarios load the AnnData sample datasets of cellucid-datasets; ' +
        'set CELLUCID_DATASETS_REPO to a checkout that holds anndata/.'
    );
  }
  if (catalog?.schema !== 2 || !Array.isArray(catalog.datasets)) {
    throw new Error(
      `${root}/datasets.json is not a schema 2 dataset catalog.`
    );
  }
  return root;
}

/** The CSS viewport of every capture that names no other; the app's documentation size. */
const VIEWPORT = Object.freeze({ width: 1440, height: 1000 });

/**
 * The browser context a scenario runs in. `viewport: { width, height,
 * touch }` is for a page about another screen — the phone layout, which the
 * app switches to at 720 CSS px and below. `touch` makes the context a
 * touch device (`isMobile`, `hasTouch`), so the app sees the coarse pointer
 * a phone has, not a narrow desktop window.
 *
 * @param {Object} scenario
 * @returns {{width: number, height: number, touch: boolean}}
 */
function viewportOf(scenario) {
  const { width = VIEWPORT.width, height = VIEWPORT.height, touch = false } =
    scenario.viewport ?? {};
  if (!Number.isInteger(width) || !Number.isInteger(height) || width < 200 || height < 200) {
    throw new Error(`${scenario.id}: viewport must be whole CSS pixels, at least 200x200.`);
  }
  return { width, height, touch: touch === true };
}

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
  const catalogPort = Number(
    process.env.CELLUCID_DOCS_CATALOG_PORT ?? String(port + 2)
  );
  if (
    !Number.isInteger(catalogPort) ||
    catalogPort < 1024 ||
    catalogPort > 65535 ||
    catalogPort === port ||
    catalogPort === samplePort
  ) {
    throw new RangeError(
      'CELLUCID_DOCS_CATALOG_PORT must be a free port between 1024 and 65535 ' +
        `that differs from ${port} and ${samplePort}; received ${catalogPort}.`
    );
  }
  return { host, port, samplePort, catalogPort, origin };
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
 * Downscale the settled PNG with the browser's own high-quality resampler and
 * encode it as the scenario's format (`image/png` or `image/webp`).
 *
 * Both directions travel as base64 because Playwright's evaluate boundary does
 * not carry typed arrays.
 *
 * @param {{data: string, width: number, type: string, quality?: number}} request
 * @returns {Promise<{data: string, width: number, height: number}>}
 */
async function encodeInPage(request) {
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
  const blob = await canvas.convertToBlob(
    request.quality === undefined
      ? { type: request.type }
      : { type: request.type, quality: request.quality }
  );
  if (blob.type !== request.type) {
    throw new Error(`This browser cannot encode ${request.type} (got ${blob.type}).`);
  }
  const buffer = new Uint8Array(await blob.arrayBuffer());
  let out = '';
  for (let index = 0; index < buffer.length; index += 1) {
    out += String.fromCharCode(buffer[index]);
  }
  return { data: btoa(out), width: request.width, height };
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
    // The application rebuilds the field selectors when the obs field list
    // arrives, which is *after* the point count settles. A selection
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
    // Some controls are rebuilt while the choice is applied (the Shows picker
    // re-renders its options once the view holds the new embedding), so the
    // value is given a few seconds to come back before it counts as discarded.
    await page.waitForFunction(
      ([selector, expected]) => document.querySelector(selector)?.value === expected,
      [step.select, chosen[0]],
      { timeout: 5_000 }
    ).catch(() => {});
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

function clipFromRects(rects, padding, selectors, viewport) {
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
    width: Math.min(viewport.width - x, Math.round(right) - x),
    height: Math.min(viewport.height - y, Math.round(bottom) - y),
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
  const viewport = page.viewportSize();
  if (shot.mode === 'window') {
    return { x: 0, y: 0, width: viewport.width, height: viewport.height };
  }
  const selectors = cropSelectors(shot);
  for (const selector of selectors) {
    await page.waitForSelector(selector, { state: 'visible' });
  }
  return clipFromRects(
    await readRects(page, selectors),
    shot.padding ?? 0,
    selectors,
    viewport
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
 * and a since-retired refusal crop, both of `#notification-center` around a
 * failure the notification center auto-dismisses four seconds later. Each was
 * 776 px of background scatter with no text in it, carrying an `:alt:`
 * describing a message that was never in the frame.
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
      selectors,
      page.viewportSize()
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

  if (cursor === null && rings.length === 0) {
    return { undo: () => Promise.resolve(), clip: clip };
  }

  await page.evaluate(paintOverlay, {
    rootId: OVERLAY_ROOT_ID,
    cursor,
    rings,
  });

  return {
    undo: async () => {
      if (scenario.cursor?.state === 'press') await page.mouse.up();
      await page.evaluate(clearOverlay, OVERLAY_ROOT_ID);
    },
    clip: clipWithPointer(clip, cursor, page.viewportSize()),
  };
}

/**
 * The crop, grown to hold the whole pointer drawn on it.
 *
 * The pointer is drawn with its tip on a control, and the arrow extends
 * 22x32 CSS px down and to the right of the tip (a press adds a 30 px ring
 * centred on it). A control at the edge of its crop — a button at the bottom
 * of a panel — used to leave half an arrow in the picture. The crop is still
 * named by its selectors; only the pointer can widen it, and never past the
 * viewport.
 *
 * @param {{x: number, y: number, width: number, height: number}} clip
 * @param {{x: number, y: number, state: string}|null} cursor
 * @param {{width: number, height: number}} viewport
 */
function clipWithPointer(clip, cursor, viewport) {
  if (cursor === null) return clip;
  const margin = 4;
  const ring = cursor.state === 'press' ? 15 : 0;
  const left = Math.max(0, Math.min(clip.x, Math.floor(cursor.x - ring - margin)));
  const top = Math.max(0, Math.min(clip.y, Math.floor(cursor.y - ring - margin)));
  const right = Math.min(
    viewport.width,
    Math.max(clip.x + clip.width, Math.ceil(cursor.x + Math.max(22, ring) + margin))
  );
  const bottom = Math.min(
    viewport.height,
    Math.max(clip.y + clip.height, Math.ceil(cursor.y + Math.max(32, ring) + margin))
  );
  return { x: left, y: top, width: right - left, height: bottom - top };
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

/**
 * The file format a scenario emits. PNG is lossless and is what every
 * `{figure}` in the guides uses. WebP is for images a page places itself with
 * `<img>`/`srcset` (the documentation home page), where a full-window
 * capture of a dense point cloud would otherwise weigh megabytes; its
 * `quality` is part of the scenario, so the bytes stay reproducible.
 *
 * @param {Object} scenario
 * @returns {{type: string, extension: string, quality?: number}}
 */
function formatOf(scenario) {
  const format = scenario.format ?? 'png';
  if (format === 'png') {
    if (scenario.quality !== undefined) {
      throw new Error(`${scenario.id}: quality applies to webp only; PNG is lossless.`);
    }
    return { type: 'image/png', extension: 'png' };
  }
  if (format === 'webp') {
    const quality = scenario.quality;
    if (typeof quality !== 'number' || !(quality > 0 && quality <= 1)) {
      throw new Error(`${scenario.id}: a webp scenario names its quality in (0, 1].`);
    }
    return { type: 'image/webp', extension: 'webp', quality };
  }
  throw new Error(`${scenario.id}: unknown format "${format}"; expected png or webp.`);
}

async function captureScenario(browser, scenario, urls) {
  // A fresh context per scenario, not a fresh page. The application persists
  // theme, sidebar geometry and onboarding state, so a shared context makes the
  // result depend on what ran before it — which is the exact property this tool
  // exists to eliminate.
  const viewport = viewportOf(scenario);
  const context = await browser.newContext({
    viewport: { width: viewport.width, height: viewport.height },
    isMobile: viewport.touch,
    hasTouch: viewport.touch,
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
  // Errors the application logs (a failed request, a refused reading) are
  // printed with the scenario's result so a stray notification in a frame
  // comes with its cause; they are not part of the provenance record.
  const loggedErrors = [];
  page.on('console', message => {
    if (message.type() === 'error') loggedErrors.push(message.text().slice(0, 400));
  });

  try {
    // Anything that has to exist before the first byte of application code
    // runs — `page.route` interceptors, `addInitScript` seeds — belongs here.
    // A step cannot do it: steps run after navigation, by which time the
    // application has already made the requests a mock exists to answer. The
    // community annotation scenarios use this to stand up a synthetic GitHub
    // Worker origin, so the whole feature is reachable without an account and
    // without touching a real repository.
    // `servedAs` shows the application under another origin, for features
    // that exist only on origins a provider registered (Google and Microsoft
    // sign-in answer `http://127.0.0.1:8765`, not this run's port). Every
    // request to that origin is answered by this run's own server through a
    // context route, so no second server binds the port — which may belong to
    // a real `cellucid serve` — and nothing leaves the machine.
    // A page answered by a route has no network address, so Chromium's local
    // network access check refuses its requests to this run's other loopback
    // hosts (the sample catalog): those are answered through the route too.
    const appOrigin = scenario.servedAs ?? urls.appOrigin;
    if (scenario.servedAs !== undefined) {
      const served = new URL(scenario.servedAs).origin;
      await context.route(
        url => url.hostname === '127.0.0.1' || url.hostname === 'localhost',
        async route => {
          const requested = new URL(route.request().url());
          const response = await route.fetch(requested.origin === served
            ? { url: new URL(requested.pathname + requested.search, urls.appOrigin).toString() }
            : {});
          await route.fulfill({ response });
        }
      );
    }

    if (typeof scenario.setup === 'function') {
      await scenario.setup({ context, page, urls });
    }

    // `entry` opens an address a product printed — the Viewer URL of
    // `cellucid serve`, of `cellucid_view()` in R — instead of the run's own
    // application with its sample catalog: that address is what a reader
    // opens, and the server behind it serves its own copy of the app.
    let url;
    if (typeof scenario.entry === 'function') {
      url = new URL(await scenario.entry(urls));
    } else {
      url = new URL('/index.html', appOrigin);
      url.searchParams.set('sampleCatalog', `${urls.catalogOrigin}/`);
      if (scenario.dataset != null) url.searchParams.set('dataset', scenario.dataset);
      // A query may depend on the run's origins (a `?data=` link to a loopback
      // host, or to a host the scenario starts), so it can be a function of
      // `urls`, and an async one.
      const query = typeof scenario.query === 'function'
        ? await scenario.query(urls)
        : scenario.query ?? {};
      for (const [key, value] of Object.entries(query)) {
        url.searchParams.set(key, value);
      }
    }
    await page.goto(url.toString(), { waitUntil: 'domcontentloaded' });

    // The welcome overlay is part of first-run, and a scenario that wants it
    // says so; every other scenario dismisses it the way the browser suite
    // does. It must be *waited for* before Escape is sent: the overlay paints
    // before its key handler is wired, and an early Escape is swallowed, which
    // leaves the overlay covering every later capture.
    // A startup that names its data exactly (a `?data=` link, a server)
    // keeps onboarding hidden, so the overlay is waited for only until
    // startup is either showing it or done without it. A scenario whose
    // subject is a startup that has not finished — the opening progress over
    // the empty viewer — says `awaitStartup: false` and waits for its own
    // subject instead.
    if (scenario.awaitStartup !== false) {
      const welcome = page.locator('#welcome-modal');
      await page.waitForFunction(
        () => document.getElementById('welcome-modal')?.classList.contains('hidden') === false
          || document.documentElement.dataset.startupPhase === 'ready',
        undefined,
        { timeout: 120_000 }
      );
      if (await welcome.isVisible()) {
        if (scenario.keepWelcome !== true) {
          await page.keyboard.press('Escape');
          await welcome.waitFor({ state: 'hidden' });
        }
      } else if (scenario.keepWelcome === true) {
        throw new Error(`${scenario.id}: the welcome overlay never opened.`);
      }
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

    const subject = await resolveClip(page, scenario.shot);
    const painted = await paint(page, scenario, subject);
    const clip = painted.clip;
    const { bytes, settled } = await settledScreenshot(
      page,
      clip,
      scenario.settle
    );
    await assertCropStillFramed(page, scenario.shot, subject);
    await painted.undo();

    const targetWidth = await emitWidthFor(scenario, clip);
    const native = pngDimensions(bytes);
    const format = formatOf(scenario);
    let finalBytes = bytes;
    let emitted = native;
    if (native.width !== targetWidth || format.type !== 'image/png') {
      const encoded = await page.evaluate(encodeInPage, {
        data: bytes.toString('base64'),
        width: targetWidth,
        type: format.type,
        quality: format.quality,
      });
      finalBytes = Buffer.from(encoded.data, 'base64');
      emitted = { width: encoded.width, height: encoded.height };
    }

    return {
      bytes: finalBytes,
      loggedErrors,
      record: {
        id: scenario.id,
        topic: scenario.topic,
        name: scenario.name,
        describes: scenario.describes,
        dataset: scenario.dataset ?? null,
        theme: scenario.theme ?? 'light',
        viewport: `${viewport.width}x${viewport.height}${viewport.touch ? ' touch' : ''}`,
        device_scale_factor: DEVICE_SCALE_FACTOR,
        crop: scenario.shot.mode === 'window' ? 'full window' : describeCrop(scenario.shot),
        captured_css: `${clip.width}x${clip.height}`,
        native_px: `${native.width}x${native.height}`,
        emitted_px: `${emitted.width}x${emitted.height}`,
        format: format.extension,
        ...(format.quality === undefined ? {} : { quality: format.quality }),
        cursor: scenario.cursor ? `${scenario.cursor.state ?? 'hover'} on ${scenario.cursor.on}` : null,
        rings: (scenario.rings ?? []).map(ring => ring.on),
        settled,
        uses_timeout: (scenario.steps ?? []).some(step => 'wait' in step),
        page_errors: consoleErrors,
        sha256: digest(finalBytes),
      },
    };
  } catch (error) {
    // What the page logged is usually the reason a step never came true.
    error.loggedErrors = [...consoleErrors.map(text => `page: ${text}`), ...loggedErrors];
    throw error;
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
  if (record.format === 'webp') {
    // Placed by the page's own markup (see `formatOf`), so the snippet is an
    // `<img>` with the intrinsic size that reserves its box before it loads.
    const [width, height] = record.emitted_px.split('x');
    return (
      `<img src="${relativePrefix}_static/screenshots/${record.topic}/` +
      `${record.name}.webp" width="${width}" height="${height}" ` +
      `alt="${record.describes}" loading="lazy" decoding="async">`
    );
  }
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
        `${scenario.id.padEnd(34)} ${scenario.topic}/${scenario.name}.` +
          `${formatOf(scenario).extension}\n` +
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

  const restoreLocale = await pinNumberLocale();
  let catalogServer = null;
  let appServer = null;
  let browser = null;
  try {
    catalogServer = await startCatalogServer(
      await requireSampleCatalog(),
      ports.catalogPort
    );
    appServer = await startApplicationServer(ports);
    browser = await chromium.launch({
      args: ['--lang=en-US'],
      env: { ...process.env, LANG: 'en_US.UTF-8', LC_ALL: 'en_US.UTF-8' },
    });
    await assertDecimalPointLocale(browser);
    await shootAll({ browser, ports, catalogServer, selected, flags, checking });
  } finally {
    await browser?.close();
    await appServer?.close();
    await catalogServer?.close();
    restoreLocale();
  }
}

// ---------------------------------------------------------------------------
// Number locale
// ---------------------------------------------------------------------------

/**
 * The CFPreferences domain of the browser `chromium.launch()` runs headless:
 * Playwright's `chrome-headless-shell`, an unbundled executable, whose
 * application domain is its executable name. Its renderer is the same binary,
 * so this one domain reaches the process that draws form controls.
 */
const HEADLESS_SHELL_DEFAULTS_DOMAIN = 'chrome-headless-shell';

/** A region that writes a decimal point and English month names. */
const CAPTURE_APPLE_LOCALE = 'en_US';

/**
 * Make number fields draw a decimal point, whatever this machine's region is.
 *
 * Blink formats the text of an `<input type=number>` itself. The context's
 * `locale` (navigator.language, Intl, Accept-Language) and `--lang` do not
 * reach it on macOS: there the renderer asks `NSLocale currentLocale`
 * whenever the page language and the system region share a language, so an
 * English page on an `en_TR` Mac draws "0,05". On Linux `--lang` and
 * `LANG`/`LC_ALL` (passed at launch) decide it.
 *
 * On macOS the region is overridden for the headless shell only — an
 * `AppleLocale` default in that executable's own preferences domain, never
 * the global domain — and the previous value is put back when the run ends,
 * including on Ctrl-C. `assertDecimalPointLocale` then proves the result in
 * pixels, so a machine where this does not take still fails before anything
 * is captured.
 *
 * @returns {Promise<() => void>} restores the previous preference
 */
async function pinNumberLocale() {
  if (process.platform !== 'darwin') return () => {};
  const domain = HEADLESS_SHELL_DEFAULTS_DOMAIN;
  let previous = null;
  try {
    const { stdout } = await run('defaults', ['read', domain, 'AppleLocale']);
    previous = stdout.trim();
  } catch {
    previous = null; // not set, which is the normal case
  }
  await run('defaults', ['write', domain, 'AppleLocale', CAPTURE_APPLE_LOCALE]);
  // Synchronous, and on `exit`: that event fires after a normal end, after
  // `process.exit` (the signal handlers below call it) and after an uncaught
  // exception, and only synchronous work runs in it. A run that crashes must
  // not leave the preference behind.
  let restored = false;
  const restore = () => {
    if (restored) return;
    restored = true;
    try {
      if (previous === null) {
        execFileSync('defaults', ['delete', domain, 'AppleLocale'], { stdio: 'ignore' });
      } else {
        execFileSync('defaults', ['write', domain, 'AppleLocale', previous], { stdio: 'ignore' });
      }
    } catch {
      // `delete` of a key that is already gone fails; nothing is left to undo.
    }
  };
  process.on('exit', restore);
  for (const [signal, number] of [['SIGINT', 2], ['SIGTERM', 15], ['SIGHUP', 1]]) {
    process.once(signal, () => process.exit(128 + number));
  }
  process.stdout.write(
    `number locale: ${domain} AppleLocale=${CAPTURE_APPLE_LOCALE} for this run ` +
      `(was ${previous ?? 'unset'}; restored on exit)\n`
  );
  return restore;
}

/**
 * Refuse to capture unless number fields draw a decimal point.
 *
 * The rendered text of a number field lives in a closed user-agent shadow
 * tree, so script cannot read it; pixels can. The same value is drawn in a
 * text field (which never localises) and in a number field with identical
 * metrics. Equal bytes mean the number field wrote "1234.05".
 *
 * @param {import('@playwright/test').Browser} browser
 */
async function assertDecimalPointLocale(browser) {
  const context = await browser.newContext({
    locale: 'en-US',
    deviceScaleFactor: DEVICE_SCALE_FACTOR,
  });
  try {
    const page = await context.newPage();
    await page.setContent(
      '<style>input{font:16px monospace;width:120px;border:0;padding:0;' +
        'margin:4px;outline:0;-webkit-appearance:none;appearance:textfield}' +
        'input::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}' +
        '</style><input id="t" type="text" value="1234.05">' +
        '<input id="n" type="number" value="1234.05">'
    );
    const text = await page.locator('#t').screenshot({ caret: 'hide' });
    const number = await page.locator('#n').screenshot({ caret: 'hide' });
    if (digest(text) !== digest(number)) {
      throw new Error(
        'Number fields in this browser do not draw "1234.05" (the region ' +
          'writes a decimal comma or other digits), so every slider value in a ' +
          'capture would be wrong. On macOS the tool sets AppleLocale for ' +
          `${HEADLESS_SHELL_DEFAULTS_DOMAIN} only; if that did not take, capture ` +
          'on Linux with LANG=en_US.UTF-8.'
      );
    }
  } finally {
    await context.close();
  }
}

async function shootAll({ browser, ports, catalogServer, selected, flags, checking }) {
  const manifest = await loadManifest();
  const revision = await webRepositoryRevision();
  const capturedAt = new Date().toISOString().replace(/\.\d+Z$/, 'Z');
  const failures = [];

  {
    // The qualifier must survive the abbreviation: `slice(0, 12)` on
    // `<sha>-dirty` would print exactly the string a clean tree prints.
    const [revisionSha, ...revisionQualifier] = revision.split('-');
    const shortRevision = [revisionSha.slice(0, 12), ...revisionQualifier].join(
      '-'
    );
    process.stdout.write(
      `app ${ports.origin}  catalog ${catalogServer.origin}  ` +
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
          // The web repository's permissive-CORS mirror, which answers
          // single-range requests: the host a data link reads a file from.
          sampleOrigin: `http://${ports.host}:${ports.samplePort}`,
          catalogOrigin: catalogServer.origin,
        });
      } catch (error) {
        // One broken scenario must not cost the whole sweep. The failure is
        // reported at the end and the exit status is non-zero.
        failures.push(`${scenario.id}: ${error.message}`);
        process.stdout.write(`  ${scenario.id}  FAILED: ${error.message}\n`);
        for (const logged of error.loggedErrors ?? []) {
          process.stdout.write(`      console error: ${logged}\n`);
        }
        continue;
      }
      const { bytes, record, loggedErrors } = captured;
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
        const file = path.join(directory, `${record.name}.${record.format}`);
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
          // `_static/screenshots/` must be named by a `{figure}` block, and
          // every WebP by a page's markup. An emitted image with no reference
          // fails `python -m pytest` for the whole repository until the
          // reference lands, so the two are one change.
          process.stdout.write(
            `  NOTE: this ${record.format.toUpperCase()} now fails python -m ` +
              'pytest until the reference below is pasted into a page.\n'
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
      for (const error of loggedErrors) {
        process.stdout.write(`      console error: ${error}\n`);
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
  }

  if (failures.length > 0) {
    for (const failure of failures) process.stderr.write(`${failure}\n`);
    process.exitCode = 1;
  }
}

await main();
// Scenarios may leave servers and child processes running (a loopback host,
// `cellucid serve`, an R session); each is stopped on `exit`, and nothing
// else would end the process while they hold the event loop open.
process.exit(process.exitCode ?? 0);
