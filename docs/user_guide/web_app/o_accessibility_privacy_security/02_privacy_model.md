# Privacy model

**Audience:** everyone using private data + computational users + institutional reviewers (IRB/IT/security)  
**Time:** 20–45 minutes  
**What you’ll learn:**
- A practical “what leaves my machine?” model for each loading workflow
- What Cellucid stores locally (browser memory, cache, `localStorage`, `sessionStorage`, IndexedDB)
- Where privacy surprises happen (session bundles, share links, figure metadata)
- Safe defaults for sensitive datasets (clinical / proprietary)

---

## One-paragraph summary (plain language)

Cellucid is primarily a **client-side web app**: your data is loaded into your browser and rendered locally.
Cellucid does not inherently “upload your dataset somewhere”, but:

- some workflows **do fetch data over the network** (remote servers, GitHub-hosted exports),
- some features **store state locally** in the browser (preferences, caches, Community Annotation),
- and some exported artifacts **embed metadata** (sessions, PNG/SVG provenance) that can reveal sensitive context.

If you are using sensitive data, the safe workflow is usually:

> use local prepared exports + avoid GitHub/community features + review exported artifact metadata before sharing.

:::{important}
This page is not legal advice. If you have compliance requirements (HIPAA/GDPR/IRB/corporate policy), use this page to build a concrete checklist for your organization.
:::

---

## The privacy “threat model” (what can leak)

In practice, information can leak from a visualization workflow through four routes:

1) **Network requests** (your browser talks to a server)
2) **Local browser storage** (state persists on your machine)
3) **Exported artifacts** (files you share)
4) **Screenshots/screen recordings** (often includes hidden identifiers like file paths)

This page focuses on (1)–(3). For screenshots, follow the redaction checklist in:
- {doc}`../r_screenshot_checklist/index`

---

## Workflow matrix: what leaves your machine?

This table is the “most useful” version of the privacy model.
If you only read one thing, read this.

| Workflow | Where the dataset is read from | Does dataset content traverse the network? | Typical privacy risk |
|---|---|---:|---|
| **Local demo** (`source=local-demo`) | `cellucid.com` hosted demo exports | Yes (demo files are downloaded) | Low (public demo), but still trackable traffic |
| **Local prepared export folder** (directory picker) | Your local disk | No (files read locally) | Medium: local folder name/path can leak via screenshots and figure metadata |
| **Local `.h5ad` / `.zarr`** (file picker) | Your local disk | No (files read locally) | Medium: large local file access + potential local caching; treat like opening a sensitive file in a browser |
| **Remote server** (`remote=...`) | A server you connect to (often `127.0.0.1`) | Yes (files or chunks downloaded from that server) | Depends: can be “local-only” or a true data transfer |
| **GitHub-hosted exports** (`github=owner/repo/...`) | GitHub raw content | Yes (downloaded from GitHub) | High if data is private; also GitHub access logs exist |
| **Jupyter embedded viewer** | Usually a local Python server + iframe bridge | Yes (browser ↔ local server) | Usually low if truly local, but can be risky if ports are exposed |
| **Community Annotation** (GitHub sync) | GitHub API + repo contents | Yes | You are explicitly publishing annotations/votes to GitHub; treat as public unless repo is private and access-controlled |

Key nuance:
- “No network for dataset content” does **not** mean “no network at all”.  
  The app itself is still loaded from a host (e.g., `cellucid.com`) unless you self-host or work offline.

---

## What Cellucid stores locally (and why)

There are three different “lifetimes” of data inside a browser:

### 1) In-memory only (clears on reload)

Examples:
- loaded point/embedding buffers
- loaded obs/var field arrays
- current filter/highlight state (unless you save a session)

This data disappears when you reload/close the tab.

### 2) HTTP cache (browser-managed)

When Cellucid fetches files over HTTP (demo datasets, remote server, GitHub), your browser may cache responses on disk.

Implication:
- If you load sensitive data from a remote server, that data may be present in your browser cache until you clear it.

### 3) Explicit browser storage (app-managed)

Cellucid uses browser storage for a few reasons:
- persist small preferences (theme/background) so the UI is stable across reloads
- keep lightweight “convenience mappings” (e.g., last-used Community Annotation repo)
- cache larger payloads (Community Annotation file cache) to avoid re-fetching from GitHub

What matters for privacy:
- **Storage is local to your machine and browser profile**, but it can persist for a long time.
- In managed environments, storage can be blocked or cleared automatically.

#### `localStorage` (persists until cleared)

Typical contents:
- theme preference (light/dark)
- viewer background preference
- debug toggles
- Community Annotation convenience mappings (e.g., “for dataset X, I used repo Y”)

`localStorage` should not contain secrets, but it can still reveal:
- which datasets you viewed (via dataset IDs in mapping keys),
- which repos you interacted with,
- and your preferences (which can be relevant in forensic settings).

#### `sessionStorage` (clears when the tab closes)

`sessionStorage` is used for **short-lived secrets** and per-tab state, most importantly:
- GitHub OAuth tokens for Community Annotation

This is intentional:
- closing the tab signs you out (token is gone),
- and tokens are not written into session bundles or the URL.

Related: {doc}`../j_community_annotation/index` (token lifetime and storage).

#### IndexedDB (larger local cache; persists until cleared)

IndexedDB can be used for larger cached payloads that would be unreasonable to store in `localStorage`.
In Cellucid, the most important example is:
- Community Annotation caching of GitHub files (for speed and offline-ish behavior)

Privacy implication:
- If you use Community Annotation on sensitive repos, cached copies of annotation JSON can exist locally until cleared.

---

## How to clear local traces (safe reset)

If you want to “wipe” Cellucid state on a shared machine, before a screen recording, or after handling sensitive data, do a full site-data clear.

### Option A (recommended): clear site data in the browser UI

In your browser settings for the Cellucid site:
- Clear **site data** (storage) and **cached images/files**.

This typically clears:
- `localStorage`
- `sessionStorage`
- IndexedDB
- HTTP cache entries for the site

### Option B (precise): DevTools → Application

1) Open DevTools
2) Go to **Application**
3) Clear:
   - Local Storage
   - Session Storage
   - IndexedDB (if present)
4) Hard reload the page

Advanced hint:
- If you only want to clear Cellucid keys (not everything on the origin), target keys starting with `cellucid` / `CELLUCID_` and (for Community Annotation) `cellucid:community-annotations:`.

---

## What goes over the network (and why)

Think of network activity in three buckets:

### 1) Loading the app itself

If you use `https://cellucid.com`, your browser downloads:
- HTML/CSS/JS assets
- optional font assets

If you self-host the app, the “app host” changes, but the model is the same.

### 2) Loading datasets (workflow-dependent)

- **Local folder / local `.h5ad` / local `.zarr`**: dataset content is read from disk, not fetched from the network.
- **Remote server**: dataset content (or chunks) is downloaded from that server.
- **GitHub**: dataset content is downloaded from GitHub raw endpoints.

### 3) Optional services (feature-dependent)

- **Community Annotation** uses GitHub APIs and repo contents.

---

## Exported artifacts (the most common privacy surprises)

### Session bundles (`.cellucid-session`)

Session bundles are designed to be shareable, but that means they can carry *derived* information.

They can contain:
- dataset identifiers
- field names and category labels
- highlight page/group names (your labels!)
- highlight memberships (cell index sets)

They do not contain the full dataset, but they can still be sensitive.

Read the dedicated page before sharing sessions from private datasets:
- {doc}`../l_sessions_sharing/08_security_privacy_and_trust`

---

### Figure exports (PNG/SVG)

Cellucid exports embed provenance metadata by default so you can reproduce figures later.
Depending on your loading workflow, that metadata can include:
- dataset id/name
- source type
- and sometimes local paths or source URLs

Before sharing figures publicly, skim:
- {doc}`../k_figure_export/05_metadata_and_provenance` (what is embedded)
- {doc}`../k_figure_export/07_troubleshooting_figure_export` (how to inspect/strip metadata)

---

### Share links (URLs)

Share links can encode:
- dataset selection (`source=...`, `dataset=...`)
- remote server URLs (`remote=http://...`)
- GitHub repo paths (`github=owner/repo/...`)
- annotation repo references (`annotations=owner/repo`)

This is convenient, but it means URLs can leak information in:
- chat logs
- browser history
- screen recordings

If you are privacy-constrained, prefer sharing:
- a session bundle via an approved channel, or
- an export folder via approved storage
and keep URLs “clean”.

---

## Security best practices (practical)

### 1) Treat “local server” as a real server

If you are using `cellucid-python` server mode (or any remote server), ensure it is not accidentally exposed.

Safe default:
- bind to `127.0.0.1` (localhost-only)

Risky default in shared networks:
- binding to `0.0.0.0` exposes to your LAN (anyone on the network may access the dataset if they can reach the port)

### 2) Prefer HTTPS for remote data

If you connect to a non-local remote server:
- use `https://` to avoid leaking data to passive network observers
- ensure your organization’s proxy/CSP policy is compatible with required endpoints

### 3) Only sign in (GitHub) on origins you trust

Community Annotation requires a GitHub OAuth flow.
Even though tokens are session-only, best practice is:
- use it on trusted deployments (official `cellucid.com` or your audited self-host)
- avoid signing in inside embedded/untrusted iframes

---

## Institutional / IRB checklist (copy/paste)

Use this as a starting point for approvals and internal documentation.

1) **Where does the app run?**
   - `cellucid.com` (public host) vs self-hosted internal deployment
2) **Where is dataset content stored?**
   - local disk only vs remote server vs GitHub
3) **What local persistence exists?**
   - `localStorage` preferences, annotation caches, browser HTTP cache
4) **What artifacts will be shared?**
   - session bundles, exported figures, URLs; and what redaction rules apply

---

## Troubleshooting (privacy / storage / “I thought it was local”)

### Symptom: “I loaded a local folder, but I still see network traffic”

**Likely causes (ordered):**
1) You are using the hosted app (`cellucid.com`), so the app assets load over the network.
2) You accidentally connected to a remote/GitHub source (dataset content truly loads over the network).

**How to confirm**
- DevTools → Network:
  - dataset requests look like `.bin`, `.json`, `.zarr`, etc.

**Fix**
1) If you need “no network after load”, use an offline/self-hosted deployment.
2) Ensure your source is “local-user” (directory/file picker), not `remote` or `github`.

---

### Symptom: “Community Annotation doesn’t work in private browsing / strict privacy mode”

**Likely causes:**
- The browser blocks `localStorage` and/or IndexedDB, which Community Annotation uses for caches and coordination.

**Fix**
- Use a normal browser profile (non-private window) for annotation sessions, or adjust the browser policy to allow site storage.

Related: {doc}`../j_community_annotation/03_ui_reference` (storage and failure modes).

---

### Symptom: “An exported PNG/SVG reveals a local path or private URL”

**Cause**
- Exported figures embed provenance metadata that can include `datasetUserPath` (local path) or remote source URLs.

**Fix**
1) Inspect metadata before sharing (see {doc}`../k_figure_export/05_metadata_and_provenance`).
2) Strip metadata if required (see {doc}`../k_figure_export/07_troubleshooting_figure_export`).

---

### Symptom: “A collaborator opened my session and saw sensitive labels”

**Cause**
- Session bundles preserve highlight/page/group names and field labels.

**Fix**
- Rename labels to privacy-safe identifiers (e.g., `donor_1`, `condition_A`) and re-save the session.

Read: {doc}`../l_sessions_sharing/08_security_privacy_and_trust`.

---

### Symptom: “My remote server is accessible from other machines”

**Likely cause**
- The server is bound to `0.0.0.0` or a public interface, not `127.0.0.1`.

**Fix**
- Re-run your server bound to localhost-only and/or add firewall rules.

---

## Next steps

- If you are sharing sessions: {doc}`../l_sessions_sharing/index`
- If you are exporting figures for publication: {doc}`../k_figure_export/index`
