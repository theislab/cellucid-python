// A fully local stand-in for the Cellucid GitHub Worker.
//
// Community annotation is the one documented feature that signs in to a third
// party. Capturing it must not touch a real account or write to a real
// repository, so the application is pointed at a synthetic Worker origin that
// exists only inside Playwright's route table. Every request the feature makes
// is answered from this file; nothing leaves the machine.
//
// The response shapes are taken from
// `cellucid/tests/browser/community-annotation-startup-lifecycle.spec.mjs`,
// which is the contract the product is tested against, so a scenario captured
// here shows the same screens a real connection produces.
//
// The token below is a fixed placeholder. It is never rendered into any
// captured pixel: no scenario opens a storage inspector, and Playwright
// screenshots carry no browser chrome, so no URL bar or token can appear.

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

const EMPTY_USERS_SENTINEL_SHA = '8b137891791fe96927ad78e64b0aad7bded08bdc';
const AT = '2026-07-30T09:00:00.000Z';

function capability() {
  return {
    status: 'ok',
    service: 'Cellucid GitHub Auth',
    contractVersion: 1,
    endpoints: [
      '/auth/login', '/auth/callback', '/auth/user', '/auth/installations',
      '/auth/installation-repos', '/cap/lookup-cells', '/cap/search-datasets',
      '/api/repos/*',
    ],
  };
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

function schemaDocument(kind) {
  const ids = {
    user: 'https://cellucid.com/contracts/community-annotation/user-v1.schema.json',
    config: 'https://cellucid.com/contracts/community-annotation/config-v1.schema.json',
    merges: 'https://cellucid.com/contracts/community-annotation/merges-v1.schema.json',
  };
  return { $schema: 'https://json-schema.org/draft/2020-12/schema', $id: ids[kind] };
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
// Annotator files, in the exact shape `annotations/schema.json` requires:
// `suggestions` is keyed `<fieldKey>:<categoryLabel>`, and `votes` maps a
// suggestion id to `up` or `down`. Cellucid merges these on Pull, so the vote
// tallies and the consensus arithmetic in a captured image are computed by the
// product from these documents, not drawn by the capture tool.
// ---------------------------------------------------------------------------

export const SUGGESTION = Object.freeze({
  betaCell: 'sug-beta-cell',
  betaImmature: 'sug-beta-immature',
  deltaCell: 'sug-delta-cell',
});

function user(id, login, displayName, title, suggestions, votes) {
  return {
    version: 1,
    username: `ghid_${id}`,
    githubUserId: id,
    login,
    displayName,
    ...(title ? { title } : {}),
    updatedAt: AT,
    suggestions,
    votes,
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

/** Four voters, unanimous on one label — Consensus at 100 %. */
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
      sessionStorage.setItem(
        'cellucid:github-app-auth:session',
        JSON.stringify({ token: 'local-mock-token', user: seed })
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

  const documents = new Map(Object.entries(users ?? {}));
  const blobs = [...documents].map(([key, document], index) => {
    const text = `${JSON.stringify(document, null, 2)}\n`;
    return {
      path: `annotations/users/${key}.json`,
      sha: String(index + 1).repeat(40).slice(0, 40),
      size: Buffer.byteLength(text, 'utf8'),
      text,
    };
  });

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

  await page.route(`${WORKER_ORIGIN}/**`, async route => {
    const request = route.request();
    const url = new URL(request.url());
    const method = request.method();
    const json = (body, status = 200) => route.fulfill({
      status, contentType: 'application/json', body: JSON.stringify(body),
    });

    if (url.pathname === '/') return json(capability());
    if (method === 'GET' && url.pathname === '/auth/user') return json(USER);
    if (method === 'GET' && url.pathname === '/auth/installations') {
      return json({ installations: [{ id: 71, account: { login: 'your-lab' } }] });
    }
    if (method === 'POST' && url.pathname === '/auth/installation-repos') {
      return json({ repositories });
    }
    if (method !== 'GET') {
      const headers = await request.allHeaders();
      return route.fulfill({
        status: 200,
        contentType: 'application/json',
        headers: {
          'access-control-expose-headers':
            'X-Cellucid-Operation-Id, X-Cellucid-Operation-Outcome',
          'x-cellucid-operation-id': headers['x-cellucid-operation-id'] ?? '',
          'x-cellucid-operation-outcome': 'applied',
        },
        body: JSON.stringify({ content: { sha: '4'.repeat(40) } }),
      });
    }

    const repoMatch = url.pathname.match(/^\/api\/repos\/([^/]+)\/([^/]+)$/);
    if (repoMatch) {
      if (repoStatus !== null) return json({ error: 'Not Found' }, repoStatus);
      const fullName =
        `${decodeURIComponent(repoMatch[1])}/${decodeURIComponent(repoMatch[2])}`;
      return json(repoInfo(fullName, { author, push, allowForking }));
    }
    if (/^\/api\/repos\/[^/]+\/[^/]+\/git\/trees\/[^/]+$/.test(url.pathname)) {
      const tree = [
        { type: 'tree', path: 'annotations/users', sha: 'd'.repeat(40) },
      ];
      if (blobs.length === 0) {
        tree.push({
          type: 'blob',
          path: 'annotations/users/.gitkeep',
          sha: EMPTY_USERS_SENTINEL_SHA,
          size: 1,
        });
      } else {
        for (const blob of blobs) {
          tree.push({
            type: 'blob', path: blob.path, sha: blob.sha, size: blob.size,
          });
        }
      }
      return json({ tree, truncated: false });
    }
    const blobMatch = url.pathname.match(
      /^\/api\/repos\/[^/]+\/[^/]+\/git\/blobs\/([0-9a-f]{40})$/
    );
    if (blobMatch) {
      const found = blobs.find(entry => entry.sha === blobMatch[1]);
      if (found === undefined) return json({ error: 'Not Found' }, 404);
      return json({
        sha: found.sha,
        size: found.size,
        encoding: 'base64',
        content: Buffer.from(found.text, 'utf8').toString('base64'),
      });
    }

    const contentPath = decodeURIComponent(
      url.pathname.replace(/^\/api\/repos\/[^/]+\/[^/]+\/contents\//, '')
    );
    const encoded = (document, shaCharacter) => json({
      type: 'file',
      encoding: 'base64',
      content: Buffer.from(JSON.stringify(document, null, 2), 'utf8')
        .toString('base64'),
      sha: shaCharacter.repeat(40),
    });
    if (contentPath === 'annotations/schema.json') {
      return encoded(schemaDocument('user'), 'a');
    }
    if (contentPath === 'annotations/config.schema.json') {
      return encoded(schemaDocument('config'), 'b');
    }
    if (contentPath === 'annotations/moderation/merges.schema.json') {
      return encoded(schemaDocument('merges'), 'c');
    }
    if (contentPath === 'annotations/config.json') {
      return encoded(config ?? configDocument(), 'd');
    }
    if (contentPath === 'annotations/moderation/merges.json') {
      return encoded({ version: 1, updatedAt: AT, merges: [] }, 'f');
    }
    const blob = blobs.find(entry => entry.path === contentPath);
    if (blob !== undefined) {
      return json({
        type: 'file',
        encoding: 'base64',
        content: Buffer.from(blob.text, 'utf8').toString('base64'),
        sha: blob.sha,
      });
    }
    return json({ error: 'Not Found' }, 404);
  });
}
