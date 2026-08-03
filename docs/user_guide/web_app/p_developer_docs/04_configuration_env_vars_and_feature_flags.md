# Configuration, environment variables, and feature flags

Cellucid is primarily a **static browser app**, so most configuration is done via:

- **URL parameters** (deep-linking to datasets / connections / annotation repos)
- **browser storage** (preferences like theme; safe caches; debug toggles)
- **explicit global overrides** (developer-only `window.__…` flags for local testing)

There is no general “dotenv environment variable” configuration for the frontend, because there is no build step required to run it.

:::{note}
Privacy/security expectations for storage are documented for users in:
{doc}`../o_accessibility_privacy_security/02_privacy_model`.
This page focuses on **developer-facing** “what keys exist and where they are read/written”.
:::

---

## URL parameters (deep links)

### Data loading params (dataset/source selection)

Handled by:
- `cellucid/assets/js/app/main.js`
- `cellucid/assets/js/app/url-state.js` (keeps the URL in sync with UI connections)

| Param | Example | Meaning |
|---|---|---|
| `dataset` | `?dataset=suo` | Selects a local-demo dataset by id |
| `source` | `?source=local-demo&dataset=suo` | Selects the source type for `dataset` (default: `local-demo`) |
| `remote` | `?remote=https://data.example.org&dataset=study-a` | Connects to a remote data server; `dataset` is required when its manifest contains multiple datasets |
| `anndata` | `?remote=http://127.0.0.1:8765&anndata=true` | Hint: remote server is serving live AnnData (UI shows performance warning) |
| `github` | `?github=owner/repo/path&dataset=study-a` | Connects to a GitHub-hosted exports path; `dataset` is required when its manifest contains multiple datasets |

Notes:
- A remote or GitHub manifest containing exactly one dataset selects that sole
  dataset when `dataset` is omitted. A manifest containing multiple datasets
  requires `dataset=<exact-dataset-id>`; Cellucid does not choose an arbitrary
  catalog entry.
- **Mixed content:** when the viewer origin is HTTPS (e.g. `https://www.cellucid.com`), browsers will block `remote=http://127.0.0.1:<port>` (HTTPS page fetching HTTP). For local prepared-data servers, open the exact server-backed Viewer URL (`http://127.0.0.1:<port>/?source=remote`) served by `cellucid-python`; direct AnnData uses `?anndata=true`.
- For **local-user** (browser file picker), the URL is intentionally kept clean (no local paths in the URL).
- `url-state.js` uses `history.replaceState()` so UI changes don’t spam browser history.

### Community annotation param

Handled by:
- `cellucid/assets/js/app/url-state.js`
- `cellucid/assets/js/app/community-annotations/REPO_SETUP.md` (behavior + sharing)

| Param | Example | Meaning |
|---|---|---|
| `annotations` | `?annotations=owner/repo` | Pre-selects the annotation repo for the current dataset |
| `annotations` | `?annotations=owner/repo@branch` | Same, but pins to a branch |

:::{important}
Tokens are **never** stored in the URL.
Only the repo reference may be stored (and can also be persisted in `localStorage` for convenience).
:::

### Debug / diagnostics params

| Param | Example | Meaning |
|---|---|---|
| `debug` | `?debug=1` | Enables analysis-module debug logging (see `cellucid/assets/js/app/analysis/shared/debug-utils.js`) |

:::{note}
There is no analytics debug parameter. Analytics is decided entirely by the
page's hostname in `cellucid/assets/js/app/ui/core/ga-init.js`: it initialises on
`cellucid.com`, `www.cellucid.com` and `theislab.github.io`, and nowhere else.
The resulting `window.cellucidAnalyticsEnabled` flag is defined non-writable and
non-configurable, so it cannot be flipped from the console or by a query
parameter — by design. Check it to find out whether the current origin reports.
:::

---

## Browser storage keys (preferences + safe caches)

### `localStorage` keys

| Key | Purpose | Where used |
|---|---|---|
| `CELLUCID_DEBUG` | Enables app-wide debug logger | `cellucid/assets/js/utils/debug.js` |
| `cellucid_theme` | Theme preference (`light`/`dark`) | `cellucid/assets/js/utils/theme-manager.js`, `cellucid/assets/js/app/ui/core/theme-init.js` |
| `cellucid_viewer_background` | Viewer background mode (`grid`, `grid-dark`, `white`, `black`) | `cellucid/assets/js/app/ui/core/theme-init.js`, viewer/render controls |
| `cellucid_antialias` | Antialiasing preference — `auto` (also what an absent key means), `on`, or `off`. `auto` resolves from the dataset's cell count on every dataset publication | `cellucid/assets/js/app/ui/core/antialias-preference.js`, `cellucid/assets/js/app/ui/modules/render-controls.js` |
| `cellucid_last_quote_index` | Welcome modal “quote rotation” state | `cellucid/assets/js/app/ui/onboarding/welcome-modal.js` |
| `debug` | Analysis-module debug toggle | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |

Community annotation (localStorage-backed):

| Key / prefix | Purpose | Notes |
|---|---|---|
| `cellucid:community-annotations:repo-map` | “datasetId + user → chosen annotation repo” mapping | Convenience only; no secrets |
| `cellucid:community-annotations:repo-meta` | Small per-repo preferences (e.g. branch mode) | Convenience only |
| `cellucid:community-annotations:auto-pull:v1` | Auto-pull preference | Convenience only |
| `cellucid:community-annotations:lock:*` | Cross-tab “single editor per scope” lease locks | Safety feature to prevent silent data loss |
| `cellucid:community-annotations:cache:*:files:shas` | SHA index for cached GitHub files | Content is in IndexedDB; this is only the small index |
| `cellucid:community-annotations:cache:*:session` | **The working annotation session**: profile, votes, comments, suggestions, deleted suggestions, moderation merges, last sync time | **Dataset-derived content, not a preference.** Scope is datasetId + repo@branch + GitHub user id |

:::{warning}
The `:session` entry is written to **`localStorage`**, so it survives closing the
tab. The helper that builds its key is named `toSessionStorageKey`, which refers
to the *annotation session*, not to `sessionStorage`. Do not infer a lifetime
from that name; earlier documentation did and got it wrong.

Because it holds annotation content rather than preferences, it belongs in the
privacy model: {doc}`../o_accessibility_privacy_security/02_privacy_model`.
:::

### IndexedDB databases

| Database | Store | Contents |
|---|---|---|
| `cellucid_community_annotation_file_cache` | `files` | Raw annotation JSON downloaded from GitHub |
| `cellucid_marker_cache` | `markers` | Marker-discovery results, including computed statistic arrays |

Both hold data derived from the user's dataset. The marker cache degrades
**silently** when IndexedDB is unavailable — no notification is raised.

Code: `cellucid/assets/js/app/community-annotations/file-cache.js`,
`cellucid/assets/js/app/analysis/genes-panel/marker-cache.js`.

### `sessionStorage` keys

Session storage is used for **short-lived secrets and per-tab identity**.

| Key | Purpose | Notes |
|---|---|---|
| `cellucid:github-app-auth:session` | GitHub access token **and** the signed-in user payload, in one record | Explicitly not persisted across tab close; the code refuses to write it to `localStorage` |
| `cellucid:community-annotations:tab-id:v1` | Unique tab id for scope locks | Used by scope-lock lease mechanism |

`sessionStorage` holds exactly one secret, and closing the tab destroys it.

---

## Global overrides (developer-only)

These are intended for local development and troubleshooting only.

| Global | Meaning | Where used |
|---|---|---|
| `window.CELLUCID_DEBUG = true` | Alternative way to enable the main debug logger | `cellucid/assets/js/utils/debug.js` |
| `window.__CELLUCID_DEBUG__ = true` | Enables analysis-module debug logging | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |
| `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` | Override GitHub auth worker origin (local dev only) | `cellucid/assets/js/app/community-annotations/github-auth.js` |
| `window._cellucidViewer`, `window._cellucidState` | The live viewer and state, for console debugging | `cellucid/assets/js/app/main.js` |
| `window._cellucidBenchmarkHarness` | The configuration-matrix benchmark harness, published when the Performance Benchmark panel is first opened | `cellucid/assets/js/app/main.js` |
| `window.cellucidAnalyticsEnabled` | Read-only: whether this origin reports analytics | `cellucid/assets/js/app/ui/core/ga-init.js` |

:::{important}
For security, `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` is rejected on non-local hosts unless it matches the compiled-in default origin.
This is designed to prevent accidental token exfiltration via an untrusted proxy.
:::

---

## Troubleshooting (when “configuration” causes weird behavior)

### Symptom: “The app behaves strangely after I changed something weeks ago”

Likely causes:
- Stale `localStorage` preferences or caches.

How to confirm:
- DevTools → Application → Local Storage / Session Storage.

Fix (safe reset):
- Clear `localStorage` keys that start with `cellucid` or `CELLUCID_` (theme/debug).
- For community annotation debugging, you may also clear `cellucid:community-annotations:*` keys.

### Symptom: “Community annotation sign-in works locally but fails on production”

Likely causes:
- A `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` override is being used on a non-local host (blocked by design).

Fix:
- Remove the override and use the compiled-in worker origin, or update the compiled-in default for your deployment (see `cellucid/assets/js/app/community-annotations/REPO_SETUP.md`).

---

Next: {doc}`05_app_architecture_overview`.
