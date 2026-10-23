// A read-only HTTP host for the AnnData sample catalog the screenshots load.
//
// Why this exists rather than reusing the web repository's server: the
// browser-test server (`cellucid/scripts/serve-browser-tests.mjs`) is confined
// to the web repository root by design, so it can serve the application but
// cannot serve `cellucid-datasets/anndata/`. Documentation screenshots want the
// real published datasets — a reader has to recognise the picture as the app
// they are using — so the capture run needs one more origin that hosts that
// folder. The application itself is still served by the repository's own
// script; this only supplies the catalog it points `?sampleCatalog=` at.
//
// The catalog is `datasets.json` plus the AnnData stores it lists: `.h5ad`
// files and zipped Zarr stores, which the viewer reads in place with HTTP range
// requests, and Zarr store folders, which it reads key by key. So this answers
// what a static web host (GitHub Pages) answers: a single `bytes=a-b` or
// `bytes=a-` range with 206 and that slice, `HEAD` with the size, and a Zarr
// key with no extension as bytes. The headers mirror the browser-test server's
// permissive-CORS mirror, the host the web app's own catalog specs read from.
//
// Permissive CORS is unconditional. This binds loopback only and serves a
// directory the caller names, so the exposure is the caller's own filesystem
// selection, not a network surface.

import { createReadStream } from 'node:fs';
import { lstat } from 'node:fs/promises';
import { createServer } from 'node:http';
import path from 'node:path';

const HOST = '127.0.0.1';

const CONTENT_TYPES = new Map([
  ['.h5ad', 'application/x-hdf5'],
  ['.json', 'application/json; charset=utf-8'],
  ['.txt', 'text/plain; charset=utf-8'],
  ['.zip', 'application/zip'],
]);

// Zarr keys (`.zgroup`, `X/0.0`, `zarr.json`) and session files
// (`.cellucid-session`) are bytes to the viewer.
const BYTES_CONTENT_TYPE = 'application/octet-stream';

function reply(response, statusCode, message) {
  const body = `${message}\n`;
  response.writeHead(statusCode, {
    'Access-Control-Allow-Origin': '*',
    'Cache-Control': 'no-store',
    'Content-Length': Buffer.byteLength(body),
    'Content-Type': 'text/plain; charset=utf-8',
  });
  response.end(body);
}

/**
 * The byte range one request asks for, `null` for the whole file, or `false`
 * when it cannot be satisfied.
 *
 * @param {string|undefined} header
 * @param {number} size
 * @returns {{start: number, end: number}|null|false}
 */
function requestedRange(header, size) {
  const match = /^bytes=(\d+)-(\d*)$/.exec(header ?? '');
  if (match === null) return null;
  const start = Number(match[1]);
  const end = match[2] === '' ? size - 1 : Math.min(Number(match[2]), size - 1);
  return start > end ? false : { start, end };
}

/**
 * Serve one catalog folder over loopback HTTP with permissive CORS and range
 * requests.
 *
 * @param {string} root Absolute path of the catalog folder (it holds
 *   `datasets.json`).
 * @param {number} port Loopback port to bind.
 * @returns {Promise<{origin: string, close: () => Promise<void>}>}
 */
export function startCatalogServer(root, port) {
  if (!path.isAbsolute(root)) {
    throw new TypeError(`Catalog root must be absolute; received ${root}.`);
  }

  const server = createServer(async (request, response) => {
    if (request.method !== 'GET' && request.method !== 'HEAD') {
      reply(response, 405, 'Method Not Allowed');
      return;
    }

    let pathname;
    try {
      pathname = decodeURIComponent(
        new URL(request.url, `http://${HOST}:${port}`).pathname
      );
    } catch {
      reply(response, 400, 'Invalid request URL');
      return;
    }

    // Confinement is checked on the resolved path, so `..` inside the request
    // cannot walk out of the published directory.
    const candidate = path.resolve(root, `.${pathname}`);
    const relative = path.relative(root, candidate);
    if (relative.startsWith('..') || path.isAbsolute(relative)) {
      reply(response, 403, 'Forbidden');
      return;
    }

    let metadata;
    try {
      metadata = await lstat(candidate);
    } catch {
      reply(response, 404, 'Not Found');
      return;
    }
    if (!metadata.isFile() || metadata.isSymbolicLink()) {
      reply(response, 404, 'Not Found');
      return;
    }

    const range = requestedRange(request.headers.range, metadata.size);
    if (range === false) {
      response.writeHead(416, {
        'Access-Control-Allow-Origin': '*',
        'Content-Range': `bytes */${metadata.size}`,
      });
      response.end();
      return;
    }
    const start = range === null ? 0 : range.start;
    const end = range === null ? metadata.size - 1 : range.end;
    const extension = path.extname(candidate).toLowerCase();
    response.writeHead(range === null ? 200 : 206, {
      'Access-Control-Allow-Origin': '*',
      'Accept-Ranges': 'bytes',
      'Cache-Control': 'no-store',
      'Content-Length': end - start + 1,
      'Content-Type': CONTENT_TYPES.get(extension) ?? BYTES_CONTENT_TYPE,
      ...(range === null ? {} : { 'Content-Range': `bytes ${start}-${end}/${metadata.size}` }),
    });
    if (request.method === 'HEAD' || metadata.size === 0) {
      response.end();
      return;
    }

    const stream = createReadStream(candidate, { start, end });
    stream.on('error', () => {
      if (response.headersSent) response.destroy();
      else reply(response, 500, 'Read failure');
    });
    stream.pipe(response);
  });

  return new Promise((resolve, reject) => {
    server.once('error', error => {
      reject(
        error.code === 'EADDRINUSE'
          ? new Error(
              `Screenshot catalog server cannot bind ${HOST}:${port}: the ` +
                'address is already in use. Set CELLUCID_DOCS_CATALOG_PORT to ' +
                'a free port and run again.'
            )
          : error
      );
    });
    server.listen(port, HOST, () => {
      resolve({
        origin: `http://${HOST}:${port}`,
        close: () =>
          new Promise(done => {
            server.closeAllConnections();
            server.close(() => done());
          }),
      });
    });
  });
}
