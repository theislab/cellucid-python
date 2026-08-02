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
Cellucid does not upload your dataset’s *contents* anywhere, but:

- some workflows **do fetch data over the network** (remote servers, GitHub-hosted exports),
- the public deployment **sends usage analytics to Google**, including the
  **name, id and size of whatever dataset you open** — even a local one (see
  the next section; this is the surprise most people miss),
- some features **store state locally** in the browser (preferences, caches, Community Annotation),
- and some exported artifacts **embed metadata** (sessions, PNG/SVG provenance) that can reveal sensitive context.

If you are using sensitive data, the safe workflow is usually:

> self-host or run the app locally + use local prepared exports + avoid
> GitHub/community features + review exported artifact metadata before sharing.

:::{important}
This page is not legal advice. If you have compliance requirements (HIPAA/GDPR/IRB/corporate policy), use this page to build a concrete checklist for your organization.
:::

---

## The privacy “threat model” (what can leak)

In practice, information can leak from a visualization workflow through five routes:

1) **Network requests** (your browser talks to a server)
2) **Usage analytics** (the app tells a third party what you did)
3) **Local browser storage** (state persists on your machine)
4) **Exported artifacts** (files you share)
5) **Screenshots/screen recordings** (often includes hidden identifiers like file paths)

This page focuses on (1)–(4). Before sharing a screenshot, remove patient or
sample identifiers, private repository names, usernames, tokens, URLs that
identify private infrastructure, and local filesystem paths.

---

(analytics)=
## Usage analytics (read this before using sensitive data on the public site)

Cellucid loads **Google Analytics 4**. This is the part of the privacy model
that most surprises people, so here it is in full.

### When it is active

Analytics initialises **only** when the page is served from one of three
hostnames:

- `cellucid.com`
- `www.cellucid.com`
- `theislab.github.io`

On any other origin — a self-hosted deployment, a `cellucid serve` run, a
notebook iframe, `localhost`, a file server on your institution’s network — the
flag is set to false, is not writable afterwards, and **no analytics script is
ever fetched**. There is no configuration to get wrong: the decision is made
from the hostname before anything else runs.

The measurement is configured with IP anonymisation on, and events are sent by
the browser’s background beacon mechanism so they do not block the interface.

### What is sent

Two kinds of event:

**Interface events** — which named button you clicked, the kind of click, the
pointer type, modifier keys, and the current dataset id.

The button's name here is a fixed identifier the developers assigned, such as
`session:save` or `legend:highlight-category`. It is never read from what the
button says. That distinction matters because many controls are labelled from
your data — the legend's category buttons are named after your categories, and
the highlight page buttons after names you typed — so a tracker that identified
a control by its own text would send those. A control with no assigned
identifier is counted only as its element type.

**Dataset load events** — and these are the ones that matter for privacy:

| Field | Example content |
|---|---|
| loading method | `local-user-prepared`, `local-user-h5ad`, `remote-connect`, `github-url-param`, … |
| dataset id | the id from `dataset_identity.json` |
| **dataset name** | the human-readable name your export was given |
| cell count, gene count, obs field count, edge count | the shape of your data |
| whether connectivity exists | yes/no |
| duration, success/failure, HTTP status | how the load went |
| failure reason | the first 160 characters of the error message |

Plus page-performance measurements (paint and responsiveness timings).

:::{warning}
**Loading a local file is not analytically silent.** The values in your dataset
never leave the browser — but on the public site its *name*, its *id* and its
*shape* do, and so does the text of any error it produced.

If your export is called `AML_relapse_cohort_2026_donor7`, that string is what
gets sent. If a load fails, the error message — which can contain a file name —
is sent with it.
:::

### How to avoid it entirely

Pick any one of these; each is sufficient:

- **Self-host the app.** Serve `index.html` from your own origin. Analytics
  never initialises.
- **Use the local server or notebook viewer.** `cellucid serve` and the Jupyter
  embedding are not on the production hostnames.
- **Block it at the network or browser layer** if you must use the public site —
  but prefer one of the first two, because they remove the code path rather than
  the request.

Naming exports neutrally (`cohort_a` rather than the study name) is a sensible
habit regardless, since the same string also lands in figure metadata and
session bundles.

---

## Workflow matrix: what leaves your machine?

This table is the “most useful” version of the privacy model.
If you only read one thing, read this.

| Workflow | Where the dataset is read from | Does dataset **content** traverse the network? | Does dataset **identity and shape** reach analytics on the public site? | Typical privacy risk |
|---|---|---:|---:|---|
| **Local demo** (`source=local-demo`) | Hosted demo exports | Yes (demo files are downloaded) | Yes | Low (public demo), but still trackable traffic |
| **Local prepared export folder** (directory picker) | Your local disk | No (files read locally) | **Yes** | Medium: name, id and cell counts are reported; the local folder path can also leak via screenshots and figure metadata |
| **Local `.h5ad` / `.zarr`** (file picker) | Your local disk | No (files read locally) | **Yes** | Medium: as above, plus large local file access and potential local caching |
| **Remote server** (`remote=...`) | A server you connect to (often `127.0.0.1`) | Yes (files or chunks downloaded from that server) | Yes | Depends: can be "local-only" or a true data transfer |
| **GitHub-hosted exports** (`github=owner/repo/...`) | GitHub raw content | Yes (downloaded from GitHub) | Yes | High if data is private; also GitHub access logs exist |
| **Jupyter embedded viewer** | A local Python server + iframe bridge | Yes (browser ↔ local server) | No — the notebook viewer is never on a production hostname | Usually low if truly local, but can be risky if ports are exposed |
| **Community Annotation** (GitHub sync) | GitHub API + repo contents, through Cellucid's auth proxy | Yes | Yes | You are explicitly publishing annotations/votes to GitHub; treat as public unless repo is private and access-controlled |

Two nuances that matter:

- “No network for dataset content” does **not** mean “no network at all”. The
  app itself is still loaded from a host unless you self-host or work offline,
  and on the public host it reports what you opened — see {ref}`analytics`.
- The fourth column is empty for every row if you **self-host or run locally**.
  That single choice removes the entire analytics surface.

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

Contents:
- theme preference (light/dark)
- viewer background preference
- which welcome-screen quote was shown last
- a debug toggle
- Community Annotation: an auto-pull preference, a path→revision index, and —
  the one that matters — a **per-scope annotation session**

`localStorage` holds no secrets. It is not, however, “preferences only”: the
Community Annotation session entry is **dataset-derived content**. It holds your
profile, your votes, your comments, your suggestions, deleted suggestions,
moderation merges, and the last sync time, keyed by dataset and repository.

So `localStorage` can reveal:
- which datasets you viewed (via dataset IDs in the keys),
- which repositories you interacted with,
- **the annotations and votes you made, in readable form**,
- and your preferences (which can be relevant in forensic settings).

#### `sessionStorage` (clears when the tab closes)

`sessionStorage` is used for **short-lived secrets** and per-tab state, most importantly:
- GitHub OAuth tokens for Community Annotation

This is intentional:
- closing the tab signs you out (token is gone),
- and tokens are not written into session bundles or the URL.

Related: {doc}`../j_community_annotation/index` (token lifetime and storage).

#### IndexedDB (larger local cache; persists until cleared)

IndexedDB holds larger cached payloads that would be unreasonable to store in
`localStorage`. Cellucid uses **two** databases:

- **Community Annotation file cache** — the raw annotation JSON downloaded from
  GitHub, for speed and offline-ish behaviour.
- **Marker cache** — the results of marker-discovery runs, including the
  computed statistic arrays.

Privacy implications:
- If you use Community Annotation on sensitive repos, cached copies of
  annotation JSON can exist locally until cleared.
- The marker cache is **derived from your expression data**. It is not the
  matrix, but a per-gene statistic table for named groups is not a preference
  either. On a shared machine, clear it along with everything else.

If IndexedDB is unavailable (private browsing, strict policy), the marker cache
degrades silently — you get no warning, just slower repeat analyses.

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

If you use `https://www.cellucid.com`, your browser downloads:
- HTML/CSS/JS assets
- optional font assets

Every JavaScript library the app depends on is served from the app’s own origin;
none is pulled from a third-party CDN, and the page’s content security policy
would block one that was. The single third-party script origin permitted is
Google’s tag manager, and only on the production hostnames — see
{ref}`analytics`.

If you self-host the app, the “app host” changes, but the model is the same.

### 2) Loading datasets (workflow-dependent)

- **Local folder / local `.h5ad` / local `.zarr`**: dataset content is read from disk, not fetched from the network.
- **Remote server**: dataset content (or chunks) is downloaded from that server.
- **GitHub**: dataset content is downloaded from GitHub raw endpoints.

### 3) Optional services (feature-dependent)

- **Community Annotation** does **not** talk to GitHub directly. Sign-in and all
  repository reads and writes go through a Cellucid-operated authentication
  proxy, which holds the app credentials and relays to the GitHub API. Your
  GitHub token is sent to that proxy on every annotation request — not to
  `api.github.com`. The proxy’s address is pinned in the code and the app refuses
  to send a token to any other origin, which is what stops a tampered
  deployment from exfiltrating it.
- **ORCID lookup**: typing into the ORCID field in the annotation identity form
  queries ORCID’s public search endpoint.

### 4) What is not there

No error-reporting service, no session-replay recorder, no product-analytics
vendor beyond the Google measurement described above, and no service worker.
Outside the analytics beacon, nothing is sent that you did not initiate.

---

## Exported artifacts (the most common privacy surprises)

### Session bundles (`.cellucid-session`)

Session bundles are designed to be shareable, but that means they can carry *derived* information.

They contain:
- a dataset fingerprint: source type, dataset id, cell count, variable count, and
  a **one-way digest** of the cell ordering
- field names and category labels
- highlight page/group names (your labels!)
- highlight memberships (cell index sets)
- user-defined category codes and field overlays
- a creation timestamp

They do **not** contain expression values, coordinates, obs values, your
identity, any GitHub token, or any path from the author’s machine. The cell
ordering appears only as a digest that cannot be turned back into coordinates.

So the residual risk is real but bounded: what you are sharing is *which cells*
you selected and *what you called them*, not the data itself. On a sensitive
dataset, a group named after a patient is the leak — not the file format.

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
- bind to `127.0.0.1` (localhost-only) — this **is** the default

Risky if you change it:
- binding to `0.0.0.0` exposes the server to your LAN.

Be precise about what that means, because the server has two different security
postures:

- **Reading data is not authenticated.** Health, server info, the dataset list,
  the protocol description, the web assets and every dataset artifact are served
  to anyone who can reach the port. On `0.0.0.0`, that is anyone on the network.
  This is deliberate — it is how the browser reads your data — and it is exactly
  why the bind address is the control that matters.
- **Writing back into your Python session is authenticated.** The two endpoints
  the browser uses to reach into the notebook — posting interaction events, and
  delivering a session bundle — require a per-viewer token that is generated
  fresh in the Python process, never written to disk, and handed to the browser
  along with the viewer. The server compares it with a constant-time comparison,
  so a wrong token cannot be discovered by measuring how long the rejection
  takes, and an unknown viewer is compared against a same-length dummy so the
  failure path does not reveal the token length either. A rejected event is not
  delivered to your callbacks at all, and the token is stripped from the payload
  before your code ever sees it.

The practical reading: a stranger who reaches an exposed port can **read** your
dataset, but cannot **drive** your notebook.

Binding to `127.0.0.1` is necessary but not sufficient. Any web page you have
open can publish a short-lived DNS record that re-points its own name at
`127.0.0.1` (**DNS rebinding**); the browser then treats your local Cellucid
server as same-origin, sends no `Origin`, applies no CORS check, and hands the
page your dataset.

Cellucid closes this by validating the `Host` header on every request before any
route runs:
- it accepts exactly one well-formed `Host` naming `localhost` or a loopback IP
  literal on the port it actually bound (DNS cannot rebind an IP literal),
- everything else gets `421 Misdirected Request`, with a body that never echoes
  the attacker’s value.

Consequence for reverse proxies: a proxy such as `jupyter-server-proxy` forwards
the browser’s `Host` unchanged, so the proxied request looks exactly like a
rebound one and is refused too. That is deliberate — a rebound page is
same-origin and can forge `X-Forwarded-*` headers, so no header can distinguish
the two. Declare the proxy’s host name explicitly on the Python side instead:
`cellucid serve ... --allowed-host hub.example.org`, or
`allowed_hosts=["hub.example.org"]`. Wildcards are rejected. Details:
{doc}`../../python_package/b_concepts_mental_models/06_privacy_security_and_offline_vs_online`.

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
   - `cellucid.com` (public host) vs self-hosted internal deployment.
   - This single answer determines whether Google Analytics is active at all —
     see {ref}`analytics`. For most review boards it is the first question.
2) **Where is dataset content stored?**
   - local disk only vs remote server vs GitHub
3) **If the public host is used, what identifiers does it report?**
   - dataset name, dataset id, cell/gene/field counts, and load error text.
   - Are your export names free of study, cohort or patient identifiers?
4) **What local persistence exists?**
   - `localStorage` preferences **and Community Annotation session content**,
     IndexedDB annotation and marker caches, browser HTTP cache
5) **If a local server is used, what is its bind address?**
   - `127.0.0.1` (default) vs `0.0.0.0`. Reads are unauthenticated either way.
6) **What artifacts will be shared?**
   - session bundles, exported figures, URLs; and what redaction rules apply

---

## Troubleshooting (privacy / storage / “I thought it was local”)

### Symptom: “I loaded a local folder, but I still see network traffic”

**Likely causes (ordered):**
1) You are using the hosted app (`cellucid.com`), so the app assets load over the network.
2) You are on a production hostname, so analytics beacons are being sent — see {ref}`analytics`.
3) You accidentally connected to a remote/GitHub source (dataset content truly loads over the network).

**How to confirm**
- DevTools → Network:
  - dataset requests look like `.bin`, `.json`, `.zarr`, etc.
  - analytics traffic goes to Google’s tag manager and measurement endpoints.
  - Check `window.cellucidAnalyticsEnabled` in the console: `true` means this
    origin is a production host and is reporting.

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
