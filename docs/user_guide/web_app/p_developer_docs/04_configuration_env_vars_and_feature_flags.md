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
| `remote` | `?remote=http://127.0.0.1:8765` | Connects to a remote data server and loads its first dataset |
| `anndata` | `?remote=http://127.0.0.1:8765&anndata=true` | Hint: remote server is serving live AnnData (UI shows performance warning) |
| `github` | `?github=owner/repo/path` | Connects to a GitHub-hosted exports path and loads its first dataset |

Notes:
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
| `ga_debug` | `?ga_debug=1` | Enables Google Analytics debug mode (see `cellucid/assets/js/app/ui/core/ga-init.js`) |

---

## Browser storage keys (preferences + safe caches)

### `localStorage` keys

| Key | Purpose | Where used |
|---|---|---|
| `CELLUCID_DEBUG` | Enables app-wide debug logger | `cellucid/assets/js/utils/debug.js` |
| `cellucid_theme` | Theme preference (`light`/`dark`) | `cellucid/assets/js/utils/theme-manager.js`, `cellucid/assets/js/app/ui/core/theme-init.js` |
| `cellucid_viewer_background` | Viewer background mode (`grid`, `grid-dark`, `white`, `black`) | `cellucid/assets/js/app/ui/core/theme-init.js`, viewer/render controls |
| `cellucid_last_quote_index` | Welcome modal “quote rotation” state | `cellucid/assets/js/app/ui/onboarding/welcome-modal.js` |
| `debug` | Analysis-module debug toggle | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |

Community annotation (localStorage-backed):

| Key / prefix | Purpose | Notes |
|---|---|---|
| `cellucid:community-annotations:repo-map` | “datasetId + user → chosen annotation repo” mapping | Convenience only; no secrets |
| `cellucid:community-annotations:last-repo-map` | Remembers the last used repo per dataset | Convenience only |
| `cellucid:community-annotations:repo-meta` | Small per-repo preferences (e.g. branch mode) | Convenience only |
| `cellucid:community-annotations:lock:*` | Cross-tab “single editor per scope” lease locks | Safety feature to prevent silent data loss |
| `cellucid:community-annotations:cache:*:files:shas` | SHA index for cached GitHub files | Content is in IndexedDB; this is only the small index |
| `cellucid:community-annotations:last-github-user-key` | Remembers last GitHub numeric user key (`ghid_…`) | Convenience only |

:::{note}
Community annotation also uses IndexedDB for larger cached payloads (raw user JSON files).
See `cellucid/assets/js/app/community-annotations/file-cache.js`.
:::

### `sessionStorage` keys

Session storage is used for **short-lived secrets and per-tab identity**.

| Key | Purpose | Notes |
|---|---|---|
| `cellucid:github-app-auth:token:v1` | GitHub OAuth token | Explicitly not persisted across tab close |
| `cellucid:github-app-auth:user:v1` | GitHub user identity payload | Used to show “signed in as…” |
| `cellucid:community-annotations:tab-id:v1` | Unique tab id for scope locks | Used by scope-lock lease mechanism |
| `cellucid:community-annotations:cache:*:session` | Per-scope UI/session state | Scope includes datasetId + repo@branch + GitHub user id |

---

## Global overrides (developer-only)

These are intended for local development and troubleshooting only.

| Global | Meaning | Where used |
|---|---|---|
| `window.CELLUCID_DEBUG = true` | Alternative way to enable the main debug logger | `cellucid/assets/js/utils/debug.js` |
| `window.__CELLUCID_DEBUG__ = true` | Enables analysis-module debug logging | `cellucid/assets/js/app/analysis/shared/debug-utils.js` |
| `window.__CELLUCID_DEV__` | Dev helper namespace (self-tests, etc.) | Set by UI modules when debug enabled |
| `window.__CELLUCID_GITHUB_WORKER_ORIGIN__` | Override GitHub auth worker origin (local dev only) | `cellucid/assets/js/app/community-annotations/github-auth.js` |

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
