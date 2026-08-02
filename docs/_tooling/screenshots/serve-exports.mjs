// A read-only HTTP host for a prepared Cellucid export tree.
//
// Why this exists rather than reusing the web repository's server: the
// browser-test server (`cellucid/scripts/serve-browser-tests.mjs`) is confined
// to the web repository root by design, so it can serve the application but
// cannot serve `cellucid-datasets/exports/`. Documentation screenshots want the
// real published datasets — a reader has to recognise the picture as the app
// they are using — so the capture run needs one more origin that hosts that
// tree. The application itself is still served by the repository's own script;
// this only supplies the data it points `exportsBaseUrl` at.
//
// The content-type table mirrors the browser-test server exactly, including the
// decision to serve `.gz` as `application/gzip` with **no** `Content-Encoding`
// header: the application decompresses those payloads itself, and letting the
// browser transparently inflate them would hand the reader a double-decoded
// buffer.
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
  ['.bin', 'application/octet-stream'],
  ['.f32', 'application/octet-stream'],
  ['.gz', 'application/gzip'],
  ['.json', 'application/json; charset=utf-8'],
  ['.txt', 'text/plain; charset=utf-8'],
  ['.u8', 'application/octet-stream'],
  ['.u16', 'application/octet-stream'],
  ['.zip', 'application/zip'],
]);

// Session presets ship without a conventional extension.
const FALLBACK_CONTENT_TYPE = 'application/octet-stream';

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
 * Serve one directory over loopback HTTP with permissive CORS.
 *
 * @param {string} root Absolute path of the directory to publish.
 * @param {number} port Loopback port to bind.
 * @returns {Promise<{origin: string, close: () => Promise<void>}>}
 */
export function startExportsServer(root, port) {
  if (!path.isAbsolute(root)) {
    throw new TypeError(`Exports root must be absolute; received ${root}.`);
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

    const extension = path.extname(candidate).toLowerCase();
    response.writeHead(200, {
      'Access-Control-Allow-Origin': '*',
      'Cache-Control': 'no-store',
      'Content-Length': metadata.size,
      'Content-Type': CONTENT_TYPES.get(extension) ?? FALLBACK_CONTENT_TYPE,
    });
    if (request.method === 'HEAD') {
      response.end();
      return;
    }

    const stream = createReadStream(candidate);
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
              `Screenshot exports server cannot bind ${HOST}:${port}: the ` +
                'address is already in use. Set CELLUCID_DOCS_EXPORTS_PORT to ' +
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
