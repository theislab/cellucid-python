# Testing, CI, and release process

This page defines how to **validate changes** to the Cellucid web app before you share them with others.

Cellucid is a high-performance interactive app; many regressions are not caught by “does it run?” checks.
The goal here is to make validation repeatable.

## At a glance

**Audience**
- Contributors: follow the manual checklist before opening PRs or sharing builds.
- Maintainers: use the release checklist to avoid “it worked locally” deployments.

---

## Current reality (dev-phase)

- The web app is shipped as static files (no required bundler build).
- Automated test coverage is limited; **manual smoke testing is mandatory** for UI/state changes.
- Some validation scripts exist (CSS token/type checks).

This page focuses on “what you can do today” to validate changes.

---

## Manual smoke test checklist (recommended)

Run these steps on at least one small dataset (fast) and optionally one larger dataset (stress).

### A) Boot + basic UI

1) Open the app (local or hosted) and confirm:
   - no console errors at startup
   - WebGL renders a point cloud after dataset load

2) Toggle sidebar open/closed; resize sidebar.

3) Switch theme (light/dark) and confirm plots and UI update.

### B) Data + fields

1) Load a dataset via at least one relevant path (local-demo / local-user / remote / GitHub).
2) Switch between:
   - a categorical field
   - a continuous field
   - a gene expression field (if present)
3) Confirm legend updates and color buffers update (no stale colors).

### C) Filtering + visibility

1) Apply at least one filter (category toggle or continuous range).
2) Confirm:
   - visible point count changes
   - points actually hide/show (alpha mask updates)
3) Toggle outlier filtering (if available) and confirm it affects visibility.

### D) Highlights + pages

1) Create a highlight group; add/remove cells.
2) Create a second highlight page; switch between pages.
3) Confirm highlight overlay and counts update correctly.

### E) Multiview + dimension

1) Click “Keep view” to create a snapshot view.
2) Switch view layout (grid/single) if available.
3) Change dimension (1D/2D/3D) on the live view and (if supported) on a snapshot view.

### F) Sessions

1) Save a session (`.cellucid-session`).
2) Reload the page and load the session.
3) Confirm:
   - camera/layout restored
   - active fields and filters restored
   - highlight pages/groups restored (or lazy-restored shortly after)

### G) Analysis

1) Open Page Analysis and run at least one plot.
2) Switch highlight pages and confirm analysis UI reacts.
3) (Optional) Open a copied analysis window and close it; confirm resources clean up.

### H) Figure export

1) Export a PNG and a small SVG.
2) For larger datasets, confirm large-export prompts work and hybrid/optimized exports complete.

### I) Community annotation (only if touched)

1) Sign in and confirm identity is shown.
2) Connect to a repo where the GitHub App is installed.
3) Pull and confirm cached files update.
4) Publish a small change (or verify the PR flow when lacking push access).

---

## Lightweight automated checks (available today)

### CSS design system validation (recommended when touching CSS)

Docs:
- `cellucid/assets/css/README.md`

Scripts:

```bash
node cellucid/scripts/validate-token-types.js
node cellucid/scripts/validate-tokens.js
```

These catch:
- token/type drift
- accidental raw hex usage outside token files
- theme mapping contract issues

### Dev-only dataset reload self-test (analysis cache reset)

When `CELLUCID_DEBUG=true`, the dataset UI exposes a dev helper:
- local-user in-place dataset switching self-test

It is also exposed as:
- `window.__CELLUCID_DEV__.runLocalUserInPlaceSwitchSelfTest(...)`

This is useful to catch:
- analysis caches not resetting on dataset switch
- leaked pending requests

Code:
- `cellucid/assets/js/app/ui/modules/dataset-controls.js`

---

## Documentation build sanity check (optional but recommended)

If you changed docs pages, you can build the `cellucid-python` docs locally:

```bash
make -C cellucid-python/docs html
```

This ensures MyST links/toctrees are valid.

---

## Release checklist (web app)

Because Cellucid is a static web app, releases often fail due to “small” deployment details.
Before deploying:

1) Run the manual smoke tests on a clean browser profile.
2) Confirm no console errors on load.
3) Confirm dataset loading works in the target deployment environment (CORS, mixed content, CSP).
4) Confirm `.gz` exports are not double-decompressed by the host (see {doc}`03_build_run_and_deployment`).
5) Confirm session save/load and export still work.
6) If community annotation is enabled, confirm worker origin + CORS allowlist are correct (see `cellucid/assets/js/app/community-annotations/REPO_SETUP.md`).

---

Next: {doc}`15_extension_points_overview`.
