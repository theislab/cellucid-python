// Open every "Open in Cellucid" link the documentation publishes and say
// whether each one still opens a dataset.
//
// The data-loading pages link real datasets on CELLxGENE, Zenodo, Hugging
// Face, GEO, Dropbox and GitHub so a reader can click and look. Those hosts
// change without telling anyone, so a link is only as good as its last check.
// This reads every `https://www.cellucid.com/?data=…` link in `docs/`, opens
// each one, in turn, in a real browser against the application at `--app`
// (default: the published site), and prints what it saw: the dataset's name
// and cell count, how long it took and what the Transfers list says was read,
// or the message that refused it. It exits non-zero when any link fails.
//
//   node check.mjs                                     # against www.cellucid.com
//   node check.mjs --app http://127.0.0.1:4173         # against a local build
//   node check.mjs --only zenodo                       # links whose ?data= contains it
//   node check.mjs --link geo:GSE123 --link zenodo:1/a.h5ad  # candidates, before publishing
//
// Playwright is resolved out of the web repository's node_modules, as the
// screenshot tool does. One page at a time: the links are large datasets.

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { createRequire } from 'node:module';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const DOCS_ROOT = path.resolve(HERE, '..', '..');
const WORKSPACE = path.resolve(DOCS_ROOT, '..', '..');
const WEB_REPO = process.env.CELLUCID_WEB_REPO ?? path.join(WORKSPACE, 'cellucid');
const require = createRequire(pathToFileURL(path.join(WEB_REPO, 'package.json')));

const LINK = /https:\/\/www\.cellucid\.com\/\?data=[^\s)"'<>\]]+/g;
const OPEN_TIMEOUT_MS = 240_000;

function argument(name, fallback) {
  const index = process.argv.indexOf(name);
  return index === -1 ? fallback : process.argv[index + 1];
}

function markdownFiles(directory) {
  return readdirSync(directory).flatMap(entry => {
    const full = path.join(directory, entry);
    if (entry === '_build' || entry === '_tooling' || entry.startsWith('.')) return [];
    if (statSync(full).isDirectory()) return markdownFiles(full);
    return entry.endsWith('.md') ? [full] : [];
  });
}

/** Every distinct published link, with the pages that carry it. */
function publishedLinks() {
  const links = new Map();
  for (const file of markdownFiles(DOCS_ROOT)) {
    for (const match of readFileSync(file, 'utf8').matchAll(LINK)) {
      const pages = links.get(match[0]) ?? new Set();
      pages.add(path.relative(DOCS_ROOT, file));
      links.set(match[0], pages);
    }
  }
  return links;
}

async function open(browser, app, link) {
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  const started = Date.now();
  try {
    const target = new URL(link);
    await page.goto(`${app}/${target.search}`);
    const outcome = await page.waitForFunction(() => {
      const failure = document.getElementById('cellucid-startup-failure');
      if (failure) return `refused: ${failure.textContent.replace(/\s+/g, ' ').trim()}`;
      const errors = [...document.querySelectorAll('.notification-error')]
        .map(note => note.textContent.replace(/\s+/g, ' ').trim());
      if (errors.length > 0) return `refused: ${errors.join(' | ')}`;
      const name = document.getElementById('dataset-name')?.textContent.trim() ?? '';
      const ready = document.documentElement.dataset.startupPhase === 'ready';
      return ready && name !== '' && name !== '–' && name !== '—' ? 'ok' : false;
    }, null, { timeout: OPEN_TIMEOUT_MS, polling: 250 }).then(handle => handle.jsonValue());
    const read = await page.evaluate(() => document
      .querySelector('#transfers-list .transfer-row')
      ?.innerText.match(/[\d.]+ [KMG]?B read(?: of a [\d.]+ [KMG]?B file)?, in \d+ requests?/)?.[0] ?? '');
    return {
      ok: outcome === 'ok',
      detail: outcome === 'ok'
        ? `${(await page.textContent('#dataset-name')).trim()} · ` +
          `${await page.textContent('#dataset-cells')} cells · ${read}`
        : outcome,
      seconds: Math.round((Date.now() - started) / 1000),
    };
  } catch (error) {
    return { ok: false, detail: `no outcome: ${error.message.split('\n')[0]}`, seconds: Math.round((Date.now() - started) / 1000) };
  } finally {
    await page.close();
  }
}

/** The `?data=` values given with `--link`, as the links a page would publish. */
function givenLinks() {
  return process.argv.flatMap((value, index) => (process.argv[index - 1] === '--link'
    ? [[`https://www.cellucid.com/?data=${encodeURIComponent(value)}`, new Set(['(given)'])]]
    : []));
}

const app = argument('--app', 'https://www.cellucid.com').replace(/\/$/, '');
const only = argument('--only', null);
const given = givenLinks();
const links = (given.length > 0 ? given : [...publishedLinks()])
  .filter(([link]) => only === null || link.includes(only));
if (links.length === 0) {
  console.error('No published ?data= links matched.');
  process.exit(1);
}

const { chromium } = require('playwright');
const browser = await chromium.launch();
let failures = 0;
try {
  for (const [link, pages] of links) {
    const { ok, detail, seconds } = await open(browser, app, link);
    if (!ok) failures += 1;
    const data = decodeURIComponent(new URL(link).searchParams.get('data'));
    console.log(`${ok ? 'ok  ' : 'FAIL'} ${String(seconds).padStart(4)}s  ${data}\n       ${detail}\n       in ${[...pages].join(', ')}`);
  }
} finally {
  await browser.close();
}
console.log(`\n${links.length - failures} of ${links.length} links open a dataset.`);
process.exit(failures === 0 ? 0 : 1);
