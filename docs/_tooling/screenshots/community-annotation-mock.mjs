// A fully local stand-in for the Cellucid GitHub Worker (wire contract 3).
//
// Community annotation is the one documented feature that signs in to a third
// party. Capturing it must not touch a real account or write to a real
// repository, so the application is pointed at a synthetic Worker origin that
// exists only inside Playwright's route table. Every request the feature makes
// is answered from this file; nothing leaves the machine.
//
// The response shapes are taken from
// `cellucid/tests/browser/community-annotation-startup-lifecycle.spec.mjs`
// (its `installGitHubApi`), which is the contract the product is tested
// against, so a scenario captured here shows the same screens a real
// connection produces. A Pull asks, in order: the repository metadata
// (`GET /api/repos/{owner}/{repo}`), the snapshot of the branch
// (`POST /graphql/snapshot`), the users folder by its tree id
// (`GET /api/repos/.../git/trees/{oid}`) and the changed files by blob id
// (`POST /graphql/blobs`). Every route reads the same files, so the blob ids
// of the snapshot, the users tree and `/graphql/blobs` always agree.
//
// The token below is a fixed placeholder. It is never rendered into any
// captured pixel: no scenario opens a storage inspector, and Playwright
// screenshots carry no browser chrome, so no URL bar or token can appear.

import { createHash } from 'node:crypto';

export const WORKER_ORIGIN = 'https://worker.example';

/** The signed-in user every scenario runs as. Deliberately fictional. */
export const USER = Object.freeze({ id: 4242, login: 'r-okafor' });
export const USER_KEY = `ghid_${USER.id}`;

/** Deliberately fictional repository names; no such repository is contacted. */
export const REPO = 'your-lab/pancreas-annotations';
export const OTHER_REPOS = Object.freeze([
  'your-lab/atlas-figures',
  'your-lab/pilot-labels',
]);

/**
 * The Git blob ids of the three schemas of `cellucid-annotation`, which the
 * app accepts by content (`ANNOTATION_SCHEMA_BLOB_OIDS` in
 * `cellucid/assets/js/app/community-annotations/wire-contract.js`). A Pull
 * reads them from the snapshot and never downloads the schemas.
 */
const SCHEMA_BLOB_OIDS = Object.freeze({
  user: '9c7263612ce91052cddca9bc5b25680bcaf09194',
  config: '967093fe91eb99234239399ab067d63305952a6b',
  merges: '102bae917d6fe7810acb5c494cee498ceadf6e27',
});

/** The template's `annotations/users/.gitkeep` (one LF byte). */
const EMPTY_USERS_SENTINEL_SHA = '8b137891791fe96927ad78e64b0aad7bded08bdc';

/** The response headers the Worker exposes on every relayed answer (`worker-contract.js`). */
const WORKER_EXPOSED_RESPONSE_HEADERS = Object.freeze([
  'ETag',
  'Retry-After',
  'X-RateLimit-Limit',
  'X-RateLimit-Remaining',
  'X-RateLimit-Reset',
  'X-RateLimit-Used',
  'X-RateLimit-Resource',
  'X-GitHub-Request-Id',
  'X-Cellucid-Operation-Id',
  'X-Cellucid-Operation-Outcome',
]);

/** One GitHub page of a Worker collection (`/auth/installations`, `/auth/installation-repos`). */
const COLLECTION_PAGE_ITEMS = 100;

const AT = '2026-07-30T09:00:00.000Z';
const RATE_LIMIT = Object.freeze({
  cost: 1, limit: 5000, remaining: 4999, resetAt: '2030-01-01T00:00:00Z',
});

function capability() {
  return {
    status: 'ok',
    service: 'Cellucid GitHub Auth',
    contractVersion: 3,
    endpoints: [
      '/auth/login',
      '/auth/callback',
      '/auth/revoke',
      '/auth/user',
      '/auth/installations',
      '/auth/installation-repos',
      '/cap/lookup-cells',
      '/graphql/snapshot',
      '/graphql/blobs',
      '/api/repos/*',
    ],
  };
}

/** The Git blob id of a text, as GitHub names the file. */
function blobOid(text) {
  const bytes = Buffer.from(text, 'utf8');
  return createHash('sha1')
    .update(`blob ${bytes.length}\0`)
    .update(bytes)
    .digest('hex');
}

/** A synthetic object id for a tree or a commit of the stand-in repository. */
function objectId(text) {
  return createHash('sha1').update(text).digest('hex');
}

function repoInfo(fullName, { author, push, allowForking }) {
  return {
    full_name: fullName,
    default_branch: 'main',
    private: false,
    allow_forking: allowForking,
    permissions: {
      pull: true, triage: false, push, maintain: author, admin: author,
    },
  };
}

/**
 * Repository configuration for the pancreas `clusters` column.
 *
 * @param {{minAnnotators?: number, threshold?: number, closed?: string[]}} [options]
 */
export function configDocument({
  minAnnotators = 3,
  threshold = 0.5,
  closed = [],
} = {}) {
  return {
    version: 1,
    supportedDatasets: [{
      datasetId: 'pancreas',
      name: 'Pancreatic endocrinogenesis (scVelo)',
      fieldsToAnnotate: ['clusters'],
      annotatableSettings: { clusters: { minAnnotators, threshold } },
      closedFields: closed,
    }],
  };
}

// ---------------------------------------------------------------------------
// Annotator files, in the exact shape `annotations/schema.json` (user file
// version 3) requires: everything a contributor wrote is kept per dataset id;
// `suggestions` is keyed `<fieldKey>:<categoryLabel>`, and `votes` maps a
// suggestion id to `up` or `down`. A suggestion id starts with its proposer's
// `ghid_<id>_`. Cellucid merges these on Pull, so the vote tallies and the
// consensus arithmetic in a captured image are computed by the product from
// these documents, not drawn by the capture tool.
// ---------------------------------------------------------------------------

export const SUGGESTION = Object.freeze({
  betaCell: 'ghid_4242_beta-cell',
  betaImmature: 'ghid_51_beta-immature',
});

function user(id, login, displayName, title, suggestions, votes) {
  return {
    version: 3,
    username: `ghid_${id}`,
    githubUserId: id,
    login,
    displayName,
    ...(title ? { title } : {}),
    updatedAt: AT,
    datasets: {
      pancreas: {
        fieldsToAnnotate: ['clusters'],
        lastAccessedAt: AT,
        suggestions,
        votes,
      },
    },
  };
}

const BETA_CELL = {
  id: SUGGESTION.betaCell,
  label: 'Beta cell',
  ontologyId: 'CL:0000169',
  evidence: 'Ins1/Ins2 are the top markers here and the cluster sits at the '
    + 'terminal end of the endocrine trajectory.',
  markers: [{ gene: 'Ins1', logFC: 4.1, pval: 0.0001 }],
  proposedBy: 'ghid_4242',
  proposedAt: AT,
};

const BETA_IMMATURE = {
  id: SUGGESTION.betaImmature,
  label: 'Beta (immature)',
  evidence: 'Ins2-positive but Ucn3-low, which reads as immature rather than '
    + 'a separate identity.',
  proposedBy: 'ghid_51',
  proposedAt: AT,
};

/** One voter only — below a three-annotator minimum, so the bucket is Pending. */
export function usersPending() {
  return {
    [USER_KEY]: user(4242, 'r-okafor', 'R. Okafor', 'Postdoc',
      { 'clusters:Beta': [BETA_CELL] },
      { [SUGGESTION.betaCell]: 'up' }),
  };
}

/** Four voters, unanimous on one label — Consensus at 100%. */
export function usersConsensus() {
  const votes = { [SUGGESTION.betaCell]: 'up' };
  return {
    [USER_KEY]: user(4242, 'r-okafor', 'R. Okafor', 'Postdoc',
      { 'clusters:Beta': [BETA_CELL] }, votes),
    ghid_51: user(51, 'l-marchetti', 'L. Marchetti', 'PhD student', {}, votes),
    ghid_77: user(77, 's-abadi', 'S. Abadi', null, {}, votes),
    ghid_88: user(88, 'm-haddad', 'M. Haddad', 'Group leader', {}, votes),
  };
}

/** Four voters split two/two on two labels — an exact tie, so Disputed. */
export function usersDisputed() {
  return {
    [USER_KEY]: user(4242, 'r-okafor', 'R. Okafor', 'Postdoc',
      { 'clusters:Beta': [BETA_CELL] },
      { [SUGGESTION.betaCell]: 'up', [SUGGESTION.betaImmature]: 'down' }),
    ghid_51: user(51, 'l-marchetti', 'L. Marchetti', 'PhD student',
      { 'clusters:Beta': [BETA_IMMATURE] },
      { [SUGGESTION.betaImmature]: 'up', [SUGGESTION.betaCell]: 'down' }),
    ghid_77: user(77, 's-abadi', 'S. Abadi', null, {},
      { [SUGGESTION.betaImmature]: 'up', [SUGGESTION.betaCell]: 'down' }),
    ghid_88: user(88, 'm-haddad', 'M. Haddad', 'Group leader', {},
      { [SUGGESTION.betaCell]: 'up', [SUGGESTION.betaImmature]: 'down' }),
  };
}

/**
 * Install the mock. Call from a scenario's `setup`, before navigation.
 *
 * @param {object} options
 * @param {import('@playwright/test').Page} options.page
 * @param {boolean} [options.signedIn]   Seed a stored session token.
 * @param {boolean} [options.bindRepo]   Pre-bind this dataset to the repository.
 * @param {boolean} [options.author]     Grant maintain/admin.
 * @param {boolean} [options.push]       Grant direct push.
 * @param {boolean} [options.allowForking]
 * @param {number|null} [options.repoStatus] Force a status for the repo lookup.
 * @param {object|null} [options.config]
 * @param {object|null} [options.users]  Map of `ghid_<id>` to user document.
 * @param {number} [options.repoCount]   How many repositories discovery returns.
 */
export async function installGitHubMock({
  page,
  signedIn = true,
  bindRepo = true,
  author = true,
  push = true,
  allowForking = true,
  repoStatus = null,
  config = null,
  users = null,
  repoCount = 3,
}) {
  await page.addInitScript(({ origin, user: seed, repoRef, userKey, bind, signIn }) => {
    window.__CELLUCID_GITHUB_WORKER_ORIGIN__ = origin;
    if (signIn) {
      // A GitHub App user token expires after 8 hours; a stored session
      // carries when it expires.
      sessionStorage.setItem(
        'cellucid:github-app-auth:session',
        JSON.stringify({
          expiresAt: Date.now() + 8 * 60 * 60 * 1000,
          token: 'local-mock-token',
          user: seed,
        })
      );
    } else {
      sessionStorage.removeItem('cellucid:github-app-auth:session');
    }
    if (bind) {
      localStorage.setItem(
        'cellucid:community-annotations:repo-map',
        JSON.stringify({
          [`pancreas::${userKey}`]: { repoRef, branchMode: 'default' },
        })
      );
    }
  }, {
    origin: WORKER_ORIGIN,
    user: USER,
    userKey: USER_KEY,
    repoRef: `${REPO}@main`,
    bind: bindRepo,
    signIn: signedIn,
  });

  // The repository's files by path: `{text, oid}`. A Publish writes into it,
  // so the next Pull reads what was published.
  const files = new Map();
  const putFile = (path, text) => files.set(path, { text, oid: blobOid(text) });
  putFile('annotations/config.json', `${JSON.stringify(config ?? configDocument(), null, 2)}\n`);
  putFile(
    'annotations/moderation/merges.json',
    `${JSON.stringify({ version: 3, updatedAt: AT, merges: [] }, null, 2)}\n`
  );
  for (const [key, document] of Object.entries(users ?? {})) {
    putFile(`annotations/users/${key}.json`, `${JSON.stringify(document, null, 2)}\n`);
  }

  const blobEntry = (name, file) => ({
    name,
    type: 'blob',
    oid: file.oid,
    object: { byteSize: Buffer.byteLength(file.text, 'utf8') },
  });
  /** The users folder: the template's `.gitkeep` and one file per contributor. */
  const usersTree = () => {
    const tree = [{
      path: '.gitkeep', mode: '100644', type: 'blob', sha: EMPTY_USERS_SENTINEL_SHA, size: 1,
    }];
    for (const [path, file] of [...files].sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))) {
      if (!path.startsWith('annotations/users/')) continue;
      tree.push({
        path: path.slice('annotations/users/'.length),
        mode: '100644',
        type: 'blob',
        sha: file.oid,
        size: Buffer.byteLength(file.text, 'utf8'),
      });
    }
    return { oid: objectId(`users\n${JSON.stringify(tree)}`), tree };
  };
  const snapshot = () => {
    const users = usersTree();
    const moderationEntries = [
      blobEntry('merges.json', files.get('annotations/moderation/merges.json')),
      { name: 'merges.schema.json', type: 'blob', oid: SCHEMA_BLOB_OIDS.merges, object: { byteSize: 100 } },
    ];
    const moderationOid = objectId(`moderation\n${JSON.stringify(moderationEntries)}`);
    const entries = [
      blobEntry('config.json', files.get('annotations/config.json')),
      { name: 'config.schema.json', type: 'blob', oid: SCHEMA_BLOB_OIDS.config, object: { byteSize: 100 } },
      { name: 'moderation', type: 'tree', oid: moderationOid, object: {} },
      { name: 'schema.json', type: 'blob', oid: SCHEMA_BLOB_OIDS.user, object: { byteSize: 100 } },
      { name: 'users', type: 'tree', oid: users.oid, object: {} },
    ];
    const annotationsOid = objectId(`annotations\n${JSON.stringify(entries)}`);
    return {
      data: {
        repository: {
          ref: {
            target: {
              oid: objectId(`commit\n${annotationsOid}`),
              annotations: { oid: annotationsOid, type: 'tree', object: { entries } },
              moderation: { oid: moderationOid, type: 'tree', object: { entries: moderationEntries } },
            },
          },
        },
        rateLimit: RATE_LIMIT,
      },
    };
  };
  const blobs = oids => {
    const repository = {};
    oids.forEach((oid, index) => {
      const file = [...files.values()].find(candidate => candidate.oid === oid) ?? null;
      repository[`b${index}`] = file === null
        ? null
        : {
            oid,
            byteSize: Buffer.byteLength(file.text, 'utf8'),
            isBinary: false,
            isTruncated: false,
            text: file.text,
          };
    });
    return { data: { repository, rateLimit: RATE_LIMIT } };
  };

  // Every full name must be unique: the chooser keys its cards by it, and a
  // duplicate makes the grid render nothing at all.
  const repositories = Array.from({ length: repoCount }, (_, index) => {
    if (index === 0) return { id: 1, full_name: REPO, private: false };
    const named = OTHER_REPOS[index - 1];
    return {
      id: index + 1,
      full_name: named ?? `your-lab/project-${String(index).padStart(3, '0')}`,
      private: index % 3 === 2,
    };
  });
  /** One page of a collection, as the Worker answers it: `{items, nextPage}`. */
  const collectionPage = (items, url) => {
    const page = Number(url.searchParams.get('page') ?? '1');
    const start = (page - 1) * COLLECTION_PAGE_ITEMS;
    return {
      items: items.slice(start, start + COLLECTION_PAGE_ITEMS),
      nextPage: start + COLLECTION_PAGE_ITEMS < items.length ? page + 1 : null,
    };
  };

  await page.route(`${WORKER_ORIGIN}/**`, async route => {
    const request = route.request();
    const url = new URL(request.url());
    const method = request.method();
    const json = (body, status = 200, headers = {}) => route.fulfill({
      status,
      contentType: 'application/json',
      headers: {
        'access-control-expose-headers': WORKER_EXPOSED_RESPONSE_HEADERS.join(', '),
        ...headers,
      },
      body: JSON.stringify(body),
    });

    if (method === 'GET' && url.pathname === '/') return json(capability());
    if (method === 'GET' && url.pathname === '/auth/user') return json(USER);
    if (method === 'POST' && url.pathname === '/auth/revoke') {
      return route.fulfill({ status: 204 });
    }
    if (method === 'GET' && url.pathname === '/auth/installations') {
      return json(collectionPage([{ id: 71, account: { login: 'your-lab' } }], url));
    }
    if (method === 'POST' && url.pathname === '/auth/installation-repos') {
      return json(collectionPage(repositories, url));
    }
    if (method === 'POST' && url.pathname === '/graphql/snapshot') {
      return json(snapshot());
    }
    if (method === 'POST' && url.pathname === '/graphql/blobs') {
      return json(blobs(JSON.parse(request.postData() ?? '{}').oids ?? []));
    }

    const contentMatch = url.pathname.match(
      /^\/api\/repos\/[^/]+\/[^/]+\/contents\/(.+)$/
    );
    if (method !== 'GET') {
      // A Publish writes the contributor's own file with the Contents API;
      // the Worker answers a mutation with its operation id and outcome.
      const headers = await request.allHeaders();
      const operation = {
        'x-cellucid-operation-id': headers['x-cellucid-operation-id'] ?? '',
      };
      if (method === 'PUT' && contentMatch) {
        const body = JSON.parse(request.postData() ?? '{}');
        const path = decodeURIComponent(contentMatch[1]);
        putFile(path, Buffer.from(body.content, 'base64').toString('utf8'));
        return json(
          { content: { sha: files.get(path).oid } },
          200,
          { ...operation, 'x-cellucid-operation-outcome': 'applied' }
        );
      }
      // No scenario publishes through a fork; one that tries is told so.
      return json(
        { error: `The screenshot mock does not answer ${method} ${url.pathname}` },
        501,
        { ...operation, 'x-cellucid-operation-outcome': 'not-applied' }
      );
    }

    const repoMatch = url.pathname.match(/^\/api\/repos\/([^/]+)\/([^/]+)$/);
    if (repoMatch) {
      if (repoStatus !== null) return json({ error: 'Not Found' }, repoStatus);
      const fullName =
        `${decodeURIComponent(repoMatch[1])}/${decodeURIComponent(repoMatch[2])}`;
      return json(repoInfo(fullName, { author, push, allowForking }));
    }
    const treeMatch = url.pathname.match(
      /^\/api\/repos\/[^/]+\/[^/]+\/git\/trees\/([0-9a-f]{40})$/
    );
    if (treeMatch) {
      const users = usersTree();
      if (treeMatch[1] !== users.oid) return json({ error: 'Not Found' }, 404);
      return json({ sha: users.oid, url: url.href, tree: users.tree, truncated: false });
    }
    if (contentMatch) {
      const file = files.get(decodeURIComponent(contentMatch[1])) ?? null;
      if (file === null) return json({ error: 'Not Found' }, 404);
      return json({
        type: 'file',
        encoding: 'base64',
        content: Buffer.from(file.text, 'utf8').toString('base64'),
        sha: file.oid,
      });
    }
    return json({ error: 'Not Found' }, 404);
  });
}
