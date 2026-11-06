# Privacy model

**Audience:** everyone using private data + computational users + institutional reviewers (IRB/IT/security)  
**Time:** 20–45 minutes  
**What you’ll learn:**
- A practical “what leaves my machine?” model for each loading workflow
- What Cellucid stores locally (browser memory, cache, `localStorage`, `sessionStorage`, IndexedDB, the origin-private file system)
- What AI agents (the chat bar, Claude and ChatGPT through the connector, a local agent app) are sent
- Where privacy surprises happen (session bundles, share links, figure metadata)
- Safe defaults for sensitive datasets (clinical / proprietary)

---

## One-paragraph summary (plain language)

Cellucid is primarily a **client-side web app**: your data is loaded into your browser and rendered locally.
Cellucid does not upload your dataset’s *contents* anywhere, but:

- some workflows **do fetch data over the network** (remote servers, datasets opened from a link or a GitHub repository),
- the public deployment **sends usage analytics to Google**: which features you
  use and how loads go, with your dataset described only by coarse shape
  (cell-count range, dimensions) — never its name, file name or contents (see the
  next section),
- some features **store state locally** in the browser (preferences, caches,
  Community Annotation, the history of every step you take with each dataset,
  the list of datasets you opened, and the parts of the datasets you opened,
  kept so they reopen faster until you remove them under Session › Stored),
- an **AI agent you let work in Cellucid** sends what it reads (names,
  counts, summaries and, when it asks, values of cells or a picture of the
  plot) to the provider of its model, unless the dataset is set to **Agents
  see names and summaries only** ({ref}`privacy-agents`),
- and some exported artifacts **embed metadata** (sessions, PNG/SVG provenance) that can reveal sensitive context.

If you are using sensitive data, the safe workflow is usually:

> self-host or run the app locally + open local files (any format Cellucid
> reads) + avoid links, community features and AI agents + review exported
> artifact metadata before sharing.

:::{important}
This page is not legal advice. If you have compliance requirements (HIPAA/GDPR/IRB/corporate policy), use this page to build a concrete checklist for your organization.
:::

---

## The privacy “threat model” (what can leak)

In practice, information can leak from a visualization workflow through six routes:

1) **Network requests** (your browser talks to a server)
2) **Usage analytics** (the app tells a third party what you did)
3) **Local browser storage** (state persists on your machine)
4) **Exported artifacts** (files you share)
5) **AI agents** (what an agent reads goes to its model's provider)
6) **Screenshots/screen recordings** (often includes hidden identifiers like file paths)

This page focuses on (1)–(5). Before sharing a screenshot, remove patient or
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

and only on the app's own pages: a **viewer page** — any address with the
`viewer` parameter, which is how a Python or R notebook widget frames
`https://www.cellucid.com/?viewer=embedded` and how a Python or R program opens
a browser-tab viewer (`?viewer=true…`) — never reports, whatever host serves
it. On any other origin — a self-hosted deployment, a `cellucid serve` run,
`localhost`, a file server on your institution’s network — and on every viewer
page, the flag is set to false, is not writable afterwards, and **no analytics
script is ever fetched**. There is no configuration to get wrong: the decision
is made from the address before anything else runs
(`assets/js/app/ui/core/ga-init.js`).

Events are sent by the browser's background beacon mechanism, so they do not
block the interface. Every event Cellucid sends passes one declared schema
(`assets/js/analytics/schema.js`) before it is sent: each value is a closed
list, a bounded number, or a fixed pattern — never free text, a name, a path, a
URL, or an error message — and an event that does not fit is dropped whole.

Before Google's tag loads, the page fixes what Google may do with it, so a
change in the Google Analytics property cannot widen it: consent defaults deny
advertising storage, advertising use of user data and ad personalisation
(only analytics storage is granted), Google signals and ad-personalisation
signals are off, and every hit's page location and referrer are origin and
path only. Google's tag has two features of its own that the property turns on
or off, outside this code: **enhanced measurement** (automatic page views on
address changes, scrolls, outbound clicks, file downloads, form interactions,
video) and **user-provided data collection** (the tag looking through the page
for e-mail addresses, phone numbers and postal addresses). For Cellucid's
property, enhanced measurement is on, to show the developers how the site is
used, and user-provided data collection is off. The property's query
redaction names every parameter that says what the page opens, so no dataset
link, server address or key in one reaches Google in any address the tag
sends; `scripts/check-analytics-property.mjs` in the web app's repository
checks the live property's user-provided data collection and redaction
against this page.

### What is sent

- **One page view per visit**, with the page address and the address it came
  from reduced to origin and path (`https://www.cellucid.com/`): the query
  string — dataset ids, `remote=` addresses, `data=` links, `sampleCatalog`, a
  subsample — is never part of it.
- **Google's enhanced-measurement events** — page changes, scrolls past 90 %
  of the page, clicks on links to other sites, site searches, embedded video
  plays, file downloads from links, and form starts and submissions, each with
  the address involved (the clicked link, the downloaded file's link) and every
  data parameter redacted from it.
- **What Google Analytics adds to every hit itself:** the approximate location
  (country, region and city, from the network address, which Google Analytics
  does not keep; the property's granular location and device data collection
  is on) and the browser, operating system and device type.
- **`app_ready`** — how long startup took (milliseconds and a range such as
  `1-3s`) and whether a dataset opened.
- **`startup_failure`** — where startup stopped (a fixed phase name), which kind
  of source the page asked for, and the failure's fixed error code and class.
- **`data_load`** — once per dataset load that finishes or fails:

  | Field | Example content |
  |---|---|
  | loading method | `local-user-file`, `local-user-folder`, `remote-connect`, `data-url-param`, … |
  | source kind | `local-user`, `remote`, `link`, `local-demo` (a notebook's viewer never reports) |
  | file format it was read as (success only) | `h5ad`, `zarr`, `loom`, `tenx-h5`, `h5mu`, … |
  | outcome and duration | `success`, `3812` ms, `3-10s` |
  | cell and gene count **ranges** | `10k-100k`, `1M-10M` — never the exact counts |
  | obs field count, embedding dimensions | `14`, `2+3` |
  | how its embeddings were chosen, how many it could show and how many it opened with turned on (each a **range**) | `default` or `best_guess`; `1`, `2-3`, `4-9`, `10+` |
  | connectivity / vector field present | yes/no |
  | whether the default session its catalog names was applied | `none`, `restored`, `failed`, … |
  | dataset id | **only for the deployment's own public sample datasets** |
  | on failure | a fixed error code and class, and the HTTP status |

- **`feature_used`** — the first use, per page load, of a feature you started
  yourself, from a fixed list (for example `render_mode` / `smoke`,
  `figure_export`, `highlight_tool` / `lasso`, `data_mapping` / `embedding`,
  `embedding_columns`, `default`, `knn`, `modality`, `gene_names`, `expression`, `log1p`,
  `highly_variable`, `latent` or `vector_fields`, and `dimension` / `2` when a view is switched to a 2D
  embedding). The feature's name is a fixed
  identifier the developers assigned; nothing is read from what a control says,
  so category names, page names you typed and the `obsm` or `obsp` keys or
  embedding names you chose are never sent.
- **`web_vital`** — page-speed measurements, once each per page load: LCP, CLS,
  INP, FCP, TTFB and the long-task total and maximum, rounded, with a
  good / needs-improvement / poor rating.
- **Four device classes**, attached to the events above: touch / mouse / hybrid
  input; phone / tablet / desktop / wide window; hardware / software graphics
  (the graphics card's name itself is never sent); and whether volumetric smoke
  is supported.

:::{note}
**Loading a local file sends no name.** Its values never leave the browser, and
on the public site only the descriptors in the table do: the loading method,
outcome, duration, count ranges, dimensions, whether its embeddings were a
best guess, and whether it has a graph or a vector field. Its name, its id,
its `obsm` keys and the text of any error stay in the browser.
:::

### How to avoid it entirely

Pick any one of these; each is sufficient:

- **Self-host the app.** Serve `index.html` from your own origin. Analytics
  never initialises.
- **Use the local server or a Python/R viewer.** `cellucid serve` serves the
  app from your own machine; a notebook widget (`show_anndata()`,
  `cellucid_view()` in Jupyter or Colab) frames the hosted app as a viewer page,
  and a browser-tab viewer is a viewer page on your own machine: neither ever
  loads analytics.
- **Block it at the network or browser layer** if you must use the public site —
  an ad or tracking blocker refusing Google's script or Cellucid's own
  `assets/js/analytics/` files leaves the app fully working, because analytics
  is loaded as an optional part. Prefer one of the first two all the same,
  because they remove the code path rather than the request.

Google Analytics keeps two first-party cookies of its own on the production
hostnames, `_ga` and `_ga_K4774ZWGEQ`, to tell a returning browser from a new
one. Cellucid itself sets no cookie, and no other origin gets either.

Naming datasets neutrally (`cohort_a` rather than the study name) is a sensible
habit regardless, since the same string also lands in figure metadata and
session bundles.

---

## Workflow matrix: what leaves your machine?

This table is the “most useful” version of the privacy model.
If you only read one thing, read this.

| Workflow | Where the dataset is read from | Does dataset **content** traverse the network? | What reaches analytics on the public site? | Typical privacy risk |
|---|---|---:|---:|---|
| **Sample datasets** (`source=local-demo`) | The hosted sample catalog | Yes (the parts of the sample files a view needs are downloaded) | Load descriptors and the public sample's id | Low (public samples), but still trackable traffic |
| **Local file or folder** (file picker or drop: `.h5ad`, Zarr, `.h5mu`, `.loom`, 10x, R objects, MatrixMarket) | Your local disk | No (files read locally) | Coarse shape only — no name or id | Low–medium: only count ranges and dimensions are reported; the file name can still leak via screenshots and figure metadata |
| **Remote server** (`remote=...`) | A server you connect to (often `127.0.0.1`) | Yes (files or chunks downloaded from that server) | Coarse shape only — never the server address | Depends: can be "local-only" or a true data transfer |
| **Data links** (`data=...`: GitHub, CELLxGENE, GEO, Zenodo, Hugging Face, Figshare, Dropbox, Google Drive, OneDrive, any web address) | The host the link names | Yes (the parts a view needs are downloaded) | Coarse shape only — never the link | High if data is private; the host's access logs exist. A private Google Drive or OneDrive file would be read with your own sign-in, straight from Google or Microsoft ({doc}`../b_data_loading/16_private_drive_onedrive`); until a site's owner registers it with Google and Microsoft, that site says Drive and OneDrive are "not configured on this site yet" |
| **Notebook widget** (`show_anndata()`, `cellucid_view()` in Jupyter) | The notebook kernel, through the notebook's own connection | Only between your browser and your notebook server (the app's code comes from www.cellucid.com; the data never goes there) | No — the widget's page is a viewer page (`?viewer=embedded`), which never loads analytics | Low: no port is opened; the token never appears in an address |
| **Python or R browser-tab viewer** (`AnnDataViewer`, `serve_anndata()`, `cellucid_view()` in RStudio) | The local server the program started | Yes (browser ↔ that server, plus one WebSocket for commands and events) | No — a viewer page never loads analytics, and it is not on a production hostname | Usually low on your own computer; on a shared machine every other user of it can read the served data (see {ref}`shared-machines`) |
| **Community annotation** (GitHub sync) | GitHub API + repo contents, through Cellucid's auth proxy | Yes | That the feature was used | You are explicitly publishing annotations/votes to GitHub, in a file that names your GitHub numeric id and login and any profile fields you fill in; treat as public unless repo is private and access-controlled |
| **AI agents** (the chat bar, Claude or ChatGPT through the connector, `cellucid mcp`) | Wherever the dataset is read from (any row above) | What an agent reads goes to its model's provider: names, counts, summaries and, unless the dataset is set to names and summaries only, values of cells and pictures of the plot; never the dataset's files | Nothing of the agent's work | The provider's terms apply to what it receives; see {ref}`privacy-agents` |

Two nuances that matter:

- “No network for dataset content” does **not** mean “no network at all”. The
  app itself is still loaded from a host unless you self-host or work offline,
  and on the public host it reports what you opened — see {ref}`analytics`.
- The fourth column is empty for every row if you **self-host or run locally**.
  That single choice removes the entire analytics surface.

---

(privacy-stored-locally)=
## What Cellucid stores locally (and why)

There are three different “lifetimes” of data inside a browser:

### 1) In-memory only (clears on reload)

Examples:
- loaded point/embedding buffers
- loaded obs/var field arrays
- what was read of a file, within the tab's memory budget

This data disappears when you reload/close the tab. What you *do* with a
dataset (filters, highlights, colours, analyses: every step of its history)
is not only in memory: while **Autosave in this browser** is on (the
default, in **State › History › History settings**), each step is also
written to this browser's IndexedDB
(`cellucid-history`, below) a moment after it is made, so a reload, a closed
tab or a crash loses at most the change of its last moment
({doc}`../m_state_history/03_saved_states_and_autosave`).

### 2) HTTP cache (browser-managed)

Dataset bytes are fetched with `cache: 'no-store'` (range reads of a file, a
Zarr store's keys, a served file), so the browser does not write them to its
HTTP cache, with one exception: a Google Drive file read with your sign-in
goes on one of a few numbered request lanes (its token is a header, so each
lane's request is preflighted once) with the browser's default cache mode, so
the browser may keep those answers in its HTTP cache. What it
may also keep on disk is the app's own files, a catalog's `datasets.json` and
default sessions (revalidated with the host before each use), and small
provider API answers.

On `www.cellucid.com` (and any deployed copy, never on `localhost`), the
app's own code is also kept by its code cache (`service-worker.js`): the
scripts, stylesheets and fonts whose address names their bytes, in this
site's Cache Storage, so a later visit loads them without the network. It
keeps nothing else: no dataset byte, no page, no request to another host.
Clearing this site's data removes it.

Implication:
- The values of a file read without a sign-in are not left in the HTTP cache;
  parts of a signed-in Google Drive file, catalog files and sessions you opened
  from a link can be, until you clear it (and the parts Cellucid keeps itself
  are in the origin-private file system below).

### 3) Explicit browser storage (app-managed)

Cellucid uses browser storage for a few reasons:
- persist small preferences (theme/background) so the UI is stable across reloads
- keep lightweight “convenience mappings” (e.g., last-used Community Annotation repo)
- keep the history of what you do with each dataset, and the list of datasets
  you opened, so a closed tab keeps your work and a dataset opens again with
  one click
- cache larger payloads (Community Annotation file cache) to avoid re-fetching from GitHub

What matters for privacy:
- **Storage is local to your machine and browser profile**, but it can persist for a long time.
- In managed environments, storage can be blocked or cleared automatically.
- Nothing of it is uploaded: no Cellucid server and no account is involved.

#### `localStorage` (persists until cleared)

Contents:
- theme preference (light/dark) (`cellucid_theme`)
- viewer background preference (`cellucid_viewer_background`)
- the antialiasing choice (`cellucid_antialias`)
- whether the Session panel's **Subsampling** and **Data mapping** sections
  are open (`cellucid_session_sections`)
- the history's settings (`cellucid:history-settings`, **State › History ›
  History settings**): whether camera moves and panels are recorded, whether
  new datasets are recorded, whether going back while paused asks first,
  **Autosave in this browser**, whether the printable report draws pictures
  of the chosen steps, and **Your name in shared histories** (the one value
  here a person types: it is written into the histories you share, as the
  author of your steps); nothing about a dataset
- the chat bar's choices (`cellucid:agent-chat`): the provider, the model of
  each provider, the effort of each model, your instructions for agents, the
  token budget and whether the chat was open; never a key
- an API key you chose to **Remember on this browser**
  (`cellucid:agent-key:<provider>`, also in `sessionStorage`, below), until
  **Forget**
- the agent permissions (`cellucid:agent-permissions`): Ask, Allow or Deny
  per kind
- per dataset (by its id), whether it is set to **Agents see names and
  summaries only** and whether its notice of what agents send was answered
  (`cellucid:agent-data`)
- this browser's pairing with the Cellucid relay
  (`cellucid:agent-relay:tab-secret`: 32 random bytes, made the first time
  you connect Claude or ChatGPT; the relay keeps only a value derived from
  it), until **Agents › Manage… › Forget this browser**
- which quote the welcome screen showed last (`cellucid_last_quote_index`)
- debug switches you set yourself in the browser console (`CELLUCID_DEBUG`,
  and `debug` for the analysis code); Cellucid only reads them
- the CELLxGENE Census release last read (`cellucid:cellxgene-census-release`,
  so a `cellxgene:` link opens one round trip sooner)
- what showing a cell and running each analysis costs on this device
  (`cellucid:opening-cost:v2`, `cellucid:analysis-cost:v2`): timings and cell
  counts only, so the next opening and the next analysis are sized from
  measurements
- the speed of the hosts you read data from (`cellucid:host-link:v1`): for
  each host, its address (`https://host[:port]`, never a path, a dataset, a
  query or a token), its round trip idle and under load, its rate and when
  they were measured,
  so the next opening from that host is sized from what it really
  delivered. At most 32 hosts, the least recently used dropped; removed
  under **Session › Stored › Manage…**
- the pace of a host that once blocked this network for asking too often
  (`cellucid:host-pace:v1`): by the host's address, the least time between two
  requests it was found to bear, and when and how often it blocked
- Community Annotation: an auto-pull preference
  (`cellucid:community-annotations:auto-pull:v1`), which repository each
  dataset uses (`cellucid:community-annotations:repo-map`), the tabs' short
  write leases (`cellucid:community-annotations:lock:…`), and — the one
  that matters — a **per-scope annotation session**
  (`cellucid:community-annotations:cache:…:session`)

`localStorage` holds no sign-in token. It is not, however, “preferences only”: the
Community Annotation session entry is **dataset-derived content**. It holds your
profile, your votes, your comments, your suggestions, deleted suggestions,
moderation merges, and the last sync time, keyed by dataset and repository.

An API key remembered there is a credential for your provider account: on a
shared machine, leave **Remember on this browser** off (the key is then
kept for the tab only), or **Forget** it when you are done.

So `localStorage` can reveal:
- which hosts you read data from (the host speeds, and the pace of a host that once blocked this network),
- which datasets you viewed (via dataset IDs in the keys),
- which repositories you interacted with,
- **the annotations and votes you made, in readable form**,
- and your preferences (which can be relevant in forensic settings).

#### `sessionStorage` (per tab; clears when the tab closes)

`sessionStorage` holds the sign-ins, each a **short-lived access token**, and
per-tab state:
- the GitHub token for Community Annotation (`cellucid:github-app-auth:session`):
  a GitHub App user token that **expires after 8 hours**, kept with its expiry
  and forgotten when it expires; **Sign out** revokes it at GitHub (through the
  Cellucid service that holds the app's credentials) and forgets it. It
  arrives in the page address's fragment, which a script in the page's head
  takes out of the address before anything else runs, so even a startup that
  fails leaves it in neither the address bar nor the history;
- the Google sign-in for Google Drive (`cellucid:google-auth:session`): an
  access token that **expires after about an hour**; **Sign out** revokes it at
  Google (`oauth2.googleapis.com/revoke`);
- the Microsoft sign-in for OneDrive (`cellucid:microsoft-auth:session`): the
  access token only (**about an hour**). Microsoft's refresh token, which renews
  it without a click for up to 24 hours, is held in the page's memory and never
  stored, so after a reload the stored token serves until it expires and then
  one click signs in again. Microsoft offers a single-page application no way to
  revoke its token: **Sign out** forgets both, and the access token lapses
  within its hour;
- the chat bar's API key of each provider (`cellucid:agent-key:<provider>`),
  unless **Remember on this browser** keeps it in `localStorage` too: a reload
  keeps it, closing the tab forgets it; it is sent only to its provider
  (`api.anthropic.com`, `api.openai.com`);
- the Community Annotation tab's id (`cellucid:community-annotations:tab-id:v1`);
- whether an opening was under way in this tab (`cellucid:opening-under-way`,
  nothing but the mark: a tab the browser ended while it opened a dataset
  opens fewer cells next time);
- on a phone, whether the panel was left open or as its peek bar
  (`cellucid:phone-sheet`).

A sign-out the provider could not revoke (GitHub or Google unreachable, a
Cellucid service deployed before revocation existed, or Microsoft) is said in a
notification.

It also notes, in Chrome, Edge and other Chromium-based browsers, the id of
the dataset opened from this computer that the tab shows
(`cellucid:tab-local-dataset`), so that a reload of the tab can offer to
reopen it (below); showing any other dataset removes the note.

What this means:
- a reload of the tab keeps the sign-ins until they expire; closing the tab or
  **Sign out** ends them in this browser;
- browsers write `sessionStorage` to the profile on disk to restore tabs after
  a restart or a crash (and copy it when a tab is duplicated), which is why only
  short-lived tokens are kept there;
- tokens are never written into `localStorage`, IndexedDB, a session bundle, a
  figure, the address, the console or analytics, and an API key never into
  IndexedDB, a session bundle, the history, the address or a request to
  anything but its provider.

Related: {doc}`../j_community_annotation/index` (token lifetime and storage).

#### IndexedDB (larger local cache; persists until cleared)

IndexedDB holds larger records that would be unreasonable to store in
`localStorage`. Cellucid keeps five databases there:

- **Saved states** (`cellucid-history`) — the history of each dataset, kept
  as you work, and the states you saved under a name
  ({doc}`../m_state_history/03_saved_states_and_autosave`). Ten stores,
  written together so a crash keeps a whole write or none of it: `entries`,
  one record per saved state (whether it is the autosave or a named save, its
  name, its dataset's id and name, how to open that dataset again as a
  session file says it, how many steps it holds, its size, and when it was
  made and last written); `trackers`, `steps` and `turns`, each tracker with
  its branches, notes and bookmarks, every step (what it did, its
  parameters, when it was taken and by whom: you, an agent with its model
  and turn, or a notebook) and each agent turn; `journal`, the history's own
  edits and what is in **Recently deleted**; `chats` and `messages`, the
  conversations with agents of the chat bar, message by message (without the
  providers' reasoning or pictures); `values` and `parts`, the state after
  each step, kept once however many steps share it; and `arrays`, the typed
  arrays those hold. It is the content a session file holds: field and
  category names, names you typed, highlighted cells, filters, colours,
  views, and the results of the analyses you ran (marker statistics, a gene
  signature's scores per cell, a dot plot's means per group). It holds none
  of the dataset's own expression matrix, coordinates or obs columns. Written
  at every step while **Autosave in this browser** is on (the default;
  `cellucid:history-settings` above); off, nothing is written, and saving
  under a name still works. Listed, opened and removed under **State ›
  Session › Load…** and **Session › Stored › Manage…**.
- **My datasets** (`cellucid-library`, store `datasets`) — one record per
  dataset you opened in this browser, keyed by its id: its name; where it
  opens from (a sample's id, a link in its public form with its provider and
  whether the link carried something to grant access that was left out, a
  server's address, or a file's or folder's name); its format, layout, file
  name, host and size; its cells and genes; the cells it was last open with
  (how many, and a subsample's method, size and seed); when it was first and
  last opened, and how many times. Never the file's content, never a path,
  and never what a link carried to grant access (a key, a signature, a
  token): only that it carried one ({ref}`link-access-keys`). Listed in
  **Open data › My datasets**, forgotten one by one there (**All… ›
  Forget**) or all at once under **Session › Stored › Manage…**
  ({doc}`../b_data_loading/01_loading_options_overview`). A dataset a
  notebook or a Python or R viewer hands to its viewer is never recorded,
  nor the sample Cellucid opens by itself at startup. **All… › Export…**
  writes the datasets and the saved states you tick to a `.cellucid-library`
  file, to take them to another browser: each dataset as this list keeps it
  (never what a link carried to grant access) and each saved state whole,
  with its notes, Recently deleted and conversations with agents.
  **Import…** adds a file's datasets and saved states beside what this
  browser keeps, replacing nothing. Such a file can hold private
  conversations and the names of your datasets: treat it like a session
  file.
- **Files to reopen** (`cellucid-local-files`, store `handles`; Chrome, Edge
  and other Chromium-based browsers, when Cellucid is not inside a frame) —
  per dataset opened from this computer (by its id), the browser's *handle*
  of the file or folder you chose with **Open…** or dropped: the browser's
  own reference to that file, which only this site in this browser profile
  can use, with its name, whether it is a file or a folder, its size and
  modification time (a folder's file count and a digest of every file's
  path, size and modification time) and when it was last opened. No path and
  none of the file's content is written, and the handle never leaves this
  browser: not into a session, an export or the address. It lets **My
  datasets**, a reload (*Reopen <name>*), or a session saved on that file
  open it again with one click; reading it again needs the browser's leave,
  which it keeps while a tab of the site is open (or on every visit, if you
  allowed that) and otherwise asks you for on that click. A file that moved,
  was deleted or changed since is not opened (choose it again with
  **Open…**), and its handle is forgotten. Forgetting a dataset in **My
  datasets** forgets its handle; **Remove** on the My datasets row under
  **Session › Stored › Manage…** forgets them all. Firefox and Safari give a
  page no such handle: there nothing is kept, and a local file is chosen
  again.
- **Decoder choices** (`cellucid-decoders`) — which decoder decodes each
  compression fastest in this browser: codec families, chunk size classes
  and decoder names under a version naming the decoders, the browser and its
  cores; nothing about a dataset or a link. A later visit skips measuring
  them. Forgotten under **Session › Stored › Manage…**.
- **Community Annotation file cache**
  (`cellucid_community_annotation_file_cache`) — what each annotation file
  downloaded from GitHub holds for the dataset open (its suggestions, votes,
  comments and the contributor's profile), by Git blob id, with the reason
  for any file a Pull left out, and the commit and tree ids of the last
  completed Pull, for speed and offline-ish behaviour.

No analysis result is cached apart from the history: a result is kept with
the step that ran it, so in `cellucid-history` while **Autosave in this
browser** is on, and in a session file you share with **Results › Include**.

Privacy implications:
- The saved states are a record of your work on each dataset, in readable
  form: the names you gave groups and categories, the cells you highlighted,
  the analyses you ran and their results. On a shared machine they are
  readable by this site in this browser profile until you remove them
  (**Load…**, **Stored › Manage…**) or clear site data. Switch **Autosave in
  this browser** off to write none (named saves are still kept). They also
  hold your conversations with agents of the chat bar.
- My datasets names the datasets you opened: file and folder names, links in
  their public form, server addresses. The names of the files and folders kept
  to reopen are readable the same way; reading one still needs the browser's
  leave.
- If you use Community Annotation on sensitive repos, cached copies of
  annotation JSON can exist locally until cleared.

If IndexedDB is unavailable (private browsing in some browsers, strict
policy): no saved state is kept, so nothing is offered to continue and the
saved states list is empty; My datasets is kept for the page only and says
so (*Kept for this page only: this window keeps nothing for this site.*); no
file is kept to reopen; decoders are measured again on every visit; and
Community Annotation **Pull** fails and says so, without replacing the
annotation view on screen.

#### Private and incognito windows

A private or incognito window keeps nothing of this site after it closes:
`localStorage`, IndexedDB and the origin-private file system all go with it.
While it is open, the window keeps what its browser allows:

- Where the browser gives the page IndexedDB for the window, the saved
  states, My datasets and the files to reopen work until the window closes.
- Where it gives the page no origin-private file system (below), the
  **Stored** row says *Kept for this tab only*, and, while IndexedDB works,
  My datasets says *Kept only until this window closes.*
- Where it gives no IndexedDB at all (WebKit's private browsing, below), see
  the paragraph above: no saved state, and My datasets for the page only.

#### Origin-private file system (dataset-derived, kept across visits)

Dataset-derived bytes are written to the browser's origin-private file
system, a disk area only this site can read:

- a file that cannot be read in parts (an R object, `.rds`/`.rda`, or a
  MatrixMarket folder) is read whole once and converted into an AnnData
  store there (`cellucid-transcoded/`), so the views can read it in place;
- a file whose host serves only whole files is downloaded there as it
  arrives (`cellucid-transcoded/`, removed once a store converted from it is
  kept);
- a gene index of a large cell-major matrix spills there when it outgrows its
  memory budget (`cellucid-column-index/`);
- the gene index of a matrix stored cell by cell, once built over the open
  cells, is kept there (`cellucid-transcoded/`, one folder per run of cells:
  the open cells' file rows, row pointers and expression values or their
  positions in the file, never the file's name or address), so the same
  cells of the same file show genes at once later;
- the parts of a remote file (by link, from CELLxGENE, a provider, or a
  private Google Drive or OneDrive file; never a file served from this
  machine or a notebook) that Cellucid downloaded to show it — its structure, the
  embeddings, the obs columns, the expression chunks of the genes you read —
  are kept there as they arrive (`cellucid-transcoded/`, *Downloaded
  parts*), so opening the same version of the file again downloads none of
  them. A private Drive or OneDrive file is kept like any other, under its
  real version; no token, signature or link is part of what it is filed
  under or of its record. What a link is known by (its dataset id, so of
  sessions, exported figures and annotation caches; what is kept of it; the
  page's memory of its bytes) is its public form, with known signature, token
  and key parameters and any user name and password left out, and a short
  digest of any other parameter (so two files told apart only by such a
  parameter stay two datasets); a link signed again is the same dataset. See
  {ref}`link-access-keys`.

Everything kept is filed under the file's **version**: a name of its content
where the source has one (a CELLxGENE Census release, a checksum, a commit
or revision), else its exact size with the `Last-Modified` date and strong
`ETag` its host sends (a file a viewer host — `cellucid serve`, a notebook,
R — serves by its `ETag`; a file on this computer by its size and
modification time). Each opening asks the source for the version before
anything kept is used, and what is kept of the file's other versions is
removed then; a file that changes while it is open is refused ("open it
again to read the new version") and everything kept of it removed. What
this cannot notice: a file replaced by other content of exactly the same
size within the same second, on a host that sends no `ETag` and whose source
names no content. A file whose source names no version keeps nothing for a
later opening.

What is kept is **kept across visits**, so the same file opens again
without being downloaded again: it stays until you remove it under
**Session › Stored › Manage…**, until site data is cleared, or
until the browser needs the room (kept gene indexes and downloaded parts
first, then the least recently opened). Whatever no tab holds and is not kept for later visits is removed when
the app next starts, or before another conversion is written. Where the
browser has no such file system (a Firefox or Safari private window), the
conversion stays in memory and goes with the tab. In an ordinary Safari
window it is kept like anywhere else.

What is kept is filed under a digest of the file's identity. Beside it, so
the *Stored* dialog can say what each item is, its manifest records the
file's name and its host (for a gene index also the matrix and how its cells
were chosen); no address, path, file id, signature or token is written.

Where it is: inside your browser profile, in the storage the browser keeps
for this site (`www.cellucid.com`, or the address you opened Cellucid at),
managed by the browser on every system and in every browser — not a folder
you can open from the file manager. Clearing the site's data in the browser
removes it; a private or incognito window keeps none of it after the window
closes (the Stored row then says *Kept for this tab only*). The
browser decides how much a site may store; **Session › Stored ›
Manage…** shows, under *Where it's kept*, the place, the space as the
browser reports it, and whether the browser keeps the storage when the disk
runs low, with a button that asks it to — asked only when you press it
({ref}`stored-in-this-browser`).

Privacy implication:
- On a shared machine, an R object or MatrixMarket folder you opened (from
  this computer, a link, or a private Google Drive or OneDrive file), and the
  downloaded parts of any remote file you opened, stay on
  disk, readable by this site in this browser profile, under its file name
  and host, until you remove it under **Session › Stored ›
  Manage…**, clear site data, or the browser removes it for room. Use a
  private window to keep nothing.

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
- the origin-private file system
- HTTP cache entries for the site

### Option B (precise): DevTools → Application

1) Open DevTools
2) Go to **Application**
3) Clear:
   - Local Storage
   - Session Storage
   - IndexedDB (if present)
   - Storage › the origin-private file system, or **Clear site data**
4) Hard reload the page

Advanced hint:
- If you only want to clear Cellucid keys (not everything on the origin), target keys starting with `cellucid` / `CELLUCID_` and (for Community Annotation) `cellucid:community-annotations:`, and the IndexedDB databases `cellucid-history`, `cellucid-library`, `cellucid-local-files`, `cellucid-decoders` and `cellucid_community_annotation_file_cache`.

### Option C (inside Cellucid): Session › Stored › Manage… › Remove all

**Remove all** in the *Stored in this browser* dialog removes what Cellucid
keeps to open datasets again and to pick up your work: the kept files, every
saved state, My datasets with its file handles, the connection speeds and the
decoder choices ({ref}`stored-in-this-browser`). It leaves what an open tab
is using (until that tab closes), the preferences, the Community Annotation
data and the sign-ins; Options A and B clear those too.

---

## What goes over the network (and why)

Think of network activity in three buckets:

### 1) Loading the app itself

If you use `https://www.cellucid.com`, your browser downloads:
- HTML/CSS/JS assets
- optional font assets

Every JavaScript library the app depends on is served from the app’s own origin;
none is pulled from a third-party CDN, and the page’s content security policy
would block one that was. The third-party scripts it permits are named by
path, not by host: Google's tag (`www.googletagmanager.com/gtag/`), used only
on the production hostnames (see {ref}`analytics`), and Google's sign-in
(`accounts.google.com/gsi/`) and file picker (`apis.google.com/js/api.js` and
the modules it loads from `apis.google.com/_/scs/`), loaded only when you open
a Google Drive file. No page submits a form (`form-action 'none'`). The app
page compiles one WebAssembly module of its own, the decoders of compressed
data (Zstandard, LZ4, Blosc), from its own origin; it evaluates no other code
(`'wasm-unsafe-eval'`, never `'unsafe-eval'`).

If you self-host the app, the “app host” changes, but the model is the same.

### 2) Loading datasets (workflow-dependent)

- **A local file or folder**: dataset content is read from disk, not fetched from the network.
- **Remote server**: dataset content (or chunks) is downloaded from that server.
- **A link**: dataset content is downloaded from the host the link names
  (GitHub's raw endpoints, CELLxGENE's Census bucket, NCBI GEO, Zenodo,
  Hugging Face, Figshare, Dropbox, Google Drive, OneDrive, or the web address
  itself). A link to one file on GitHub also asks GitHub's public API
  (`api.github.com`, without any sign-in or token) which commit last wrote
  that file, so the file is read at that commit.

### 3) Optional services (feature-dependent)

- **Community annotation** does **not** talk to GitHub directly. Sign-in and all
  repository reads and writes go through a Cellucid-operated authentication
  proxy, which holds the app credentials and relays to the GitHub API. Your
  GitHub token is sent to that proxy on every annotation request — not to
  `api.github.com`. The proxy’s address is pinned in the code and the app refuses
  to send a token to any other origin, which is what stops a tampered
  deployment from exfiltrating it.
- **Cell type search** (Community Annotation's voting form): **Search CAP**,
  **Search Ontology** and **Search Markers** send the term, ontology id or
  marker genes you typed, through the same Cellucid service, to the Cell
  Annotation Platform (`celltype.info`), with no GitHub token, cookie or
  annotation content.
- **ORCID lookup**: typing into the ORCID field in the annotation identity form
  queries ORCID’s public search endpoint.
- **Google Drive and OneDrive**: on a site its owner registered with Google
  and Microsoft, opening a Drive file loads Google's sign-in and file-picker
  scripts (from `accounts.google.com` and `apis.google.com`), only then; a
  private Drive file or any OneDrive file asks you to sign in to Google or
  Microsoft in their own window. The page then reads the file straight from
  Google or Microsoft with your sign-in. No Cellucid server is involved, and
  the sign-in stays in the tab's `sessionStorage`
  ({doc}`../b_data_loading/16_private_drive_onedrive`). The site's owner
  registered `www.cellucid.com` with both (the public client IDs and Google's
  API key in `assets/js/data/sources/provider-config.js`), and with Google also
  `http://127.0.0.1:8765` and `http://localhost:8765`, the default address of
  `cellucid serve`, so Drive works there too. Whatever program holds port
  8765 on a computer can therefore show Google's consent screen under
  Cellucid's name; on a shared machine, another user can take that port first
  (see {ref}`shared-machines`). A copy of the app on any other address says
  "not configured on this site yet" for Drive and OneDrive and sends nothing.
- **The chat bar** (**Ask an agent**) sends each request from the tab
  straight to the provider you chose (`api.anthropic.com` or
  `api.openai.com`), with your own API key: your message, your instructions,
  the conversation so far and the tool results. Nothing passes through a
  Cellucid server. Only when you send a message.
- **The Cellucid relay** (`mcp.cellucid.com`), only once
  you connect Claude or ChatGPT (**Agents › Claude and ChatGPT**): this tab
  keeps one WebSocket to it, proven by this browser's pairing secret, and it
  carries each tool call to the tab and the result back to Claude or
  ChatGPT, which send it on to their models. The relay keeps no call and no
  result (they exist in its memory only while in flight, and are not
  logged), no account and no chat; it keeps the connections you allowed (the
  app, the agent's name, the app's host, times) and, for at most fifteen
  minutes, a connection request waiting for your answer. No dataset file
  passes through it ({doc}`../g_agents/02_connect_claude_and_chatgpt`).
- **A local agent app** (`cellucid mcp`, `cellucid_mcp()`) talks to a viewer
  on your computer over the viewer's own channel, with no relay; the app
  sends what it reads to the provider it uses
  ({doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`).

### 4) What is not there

No error-reporting service, no session-replay recorder, no product-analytics
vendor beyond the Google measurement described above, and no service worker
other than the code cache described above (deployed copies only: it keeps the
app's own files, asking this site again only for files the page just loaded
that the browser's HTTP cache no longer holds, and contacts no other host).
Outside the analytics beacon, one request is made that you did not start: the
hosted app reads its public sample catalog at startup (`theislab.github.io`,
and opens a connection to it early), unless the address names another catalog
(`?sampleCatalog=`, as `cellucid serve`, a notebook and the tests do) or the
page is a viewer page (`?viewer=…`), which never reads it.

---

## Exported artifacts (the most common privacy surprises)

### Session bundles (`.cellucid-session`)

Session bundles are designed to be shareable, but that means they can carry *derived* information.

**State › Session › Share…** chooses what the file holds, and **What's
inside** lists each kind of content with its count before you save:

- **Share**: the **History** (the trackers and branches you tick, with every
  value), the **Current state** only (one step, *Shared state*, with no
  notes, conversations, authors or times), or the **Steps only** (a replay
  plan for any dataset: no cells and no values);
- **Include**: notes and bookmarks (on by default), **Conversations** with
  agents (off by default), Recently deleted (off by default);
- **Results**: the analyses' results kept, or computed again on open;
- **Times**: exact, dates only, or none;
- **Dataset**: its name and where it opens from, or its name only.

A file shared as the history contains:
- a dataset fingerprint: dataset id and name, cell count, variable count, the
  data mapping (which `obsm`, `obsp`, layer and `var` keys are used), and a
  **one-way digest** of the cell ordering
- where the dataset was opened from: a sample id, a link in its **public
  form** (never the key, signature or token it carried; see
  {ref}`link-access-keys`), a server's address, or a local file's or folder's
  **name** (never its path)
- the keys of the fields and genes the views use
- field names and category labels
- highlight page/group names (your labels!)
- highlight memberships (cell index sets)
- user-defined category codes and field overlays
- the history (session format 4): every step, with its parameters, when it
  was taken and who took it (you, under **Your name in shared histories**
  when you set one; an agent, with its model and turn; a notebook), and the
  state after each step
- the results of the analyses you ran, kept with their steps (marker
  statistics, a dot plot's means per group, a gene signature's score for each
  cell), unless **Results** is set to compute them on open
- notes and bookmarks, and, only when you ticked them, the conversations with
  agents and Recently deleted
- a creation timestamp (rounded to the day or left out with **Times**)

They do **not** contain the dataset's expression matrix, coordinates or obs
columns, your identity (unless you set **Your name in shared histories**), any
sign-in token or API key, or any path from the author’s machine (a local
file's name, yes). The cell ordering appears only as a
digest that cannot be turned back into coordinates. An analysis result is
derived from the data, though: a signature score per cell or a group's mean
expression travels with the session.

So the residual risk is real but bounded: what you are sharing is *which cells*
you selected, *what you called them* and *what your analyses found*, not the
data itself. On a sensitive dataset, a group named after a patient is the
leak — not the file format. The saved states this browser keeps
(`cellucid-history`, above) hold the same content, and a `.cellucid-library`
file from **My datasets › Export…** holds the saved states you ticked whole,
conversations included.

Read the dedicated page before sharing sessions from private datasets:
- {doc}`../l_sessions_sharing/08_security_privacy_and_trust`

---

### Figure exports (PNG/SVG)

Cellucid exports embed provenance metadata by default so you can reproduce figures later.
Depending on your loading workflow, that metadata can include:
- dataset id/name
- source type
- and a source: the URL for a remote dataset, or, for a local one, just the
  **name** of the folder or file you chose. It is never an absolute path — a
  browser does not hand a page one

Before sharing figures publicly, skim:
- {doc}`../k_figure_export/05_metadata_and_provenance` (what is embedded)
- {doc}`../k_figure_export/07_troubleshooting_figure_export` (how to inspect/strip metadata)

---

### Share links (URLs)

Share links can encode:
- dataset selection (`source=...`, `dataset=...`)
- remote server URLs (`remote=http://...`)
- data links in their public form (`data=gh:owner/repo/...`, a GEO accession,
  a Zenodo record, a web address without its key, signature or token)
- annotation repo references (`annotations=owner/repo`)

(link-access-keys)=
#### A link's access key goes nowhere but the address bar

A link can carry what grants access to the file: a pre-signed URL's signature
(S3, Google Cloud Storage, Azure), a token or API key, a Dropbox share key
(`rlkey`), a Google Drive resource key, a user name and password. Cellucid
reads the file with the link exactly as you pasted it, and keeps that form in
two places only: this page's address bar (so a reload opens the same file; the
browser's history records it like any address) and the link field you typed it
into. Everywhere else the link goes in its
**public form**, which keeps only the parameters that name a file or its
version (`versionId`, `generation`, `alt`, `download`, `dl`, `raw`, `format`,
`ref`, `revision`, `version`, `v`, `path`, `file`, `filename`, `name`, `id`)
and drops every other one, the fragment, and any user name and password:
Share links, badges and embed snippets (the Share dialog says when something
was left out), session files, exported figures' metadata, notifications and
failure messages, the console, the opening progress, the notebook kernel (a
debug snapshot's address; a session applied to an AnnData or R object records
no `source`, unless the apply keeps the session file itself, `keep_session`,
which stores the file as it was saved) and Community Annotation. Whoever opens a shared link or
session of such a file needs their own access to it.

This is convenient, but it means URLs can leak information in:
- chat logs
- browser history
- screen recordings

If you are privacy-constrained, prefer sharing:
- a session bundle via an approved channel, or
- the dataset file via approved storage
and keep URLs “clean”.

---

(privacy-agents)=
## AI agents: what they are sent

An agent works in Cellucid through the same actions as the interface, and
what it reads goes to the provider of its model with the conversation. The
full account is {doc}`../g_agents/03_permissions_and_what_agents_see`; in
short:

- **What leaves the browser.** Names, counts and summaries and, when an
  agent asks, values of cells and pictures of the plot. Never the dataset's
  files, a sign-in token or an API key (a key goes only to its own
  provider). The chat bar sends from the tab to Anthropic or OpenAI; a
  connected Claude or ChatGPT receives results through the relay, which
  carries them through without keeping them; a local agent app receives them
  from the viewer on your computer.
- **Once per dataset** a notice says so, with the switch **Agents see names
  and summaries only** (also under **Permissions › Change…**). It is kept per
  dataset in this browser (`cellucid:agent-data`), applies to every agent
  (the chat, the connector, local apps) and no agent can turn it off.
- **With the switch on**, reads of single cells' values and pictures of the
  plot are refused (a picture counts as cell-level: the plot draws each
  cell); minimums, maximums and medians are left out of results; summaries
  over groups an agent made by naming cells are refused (these marks are
  kept for the page only, never saved in a session); a read Cellucid has not
  classified is refused; and results read earlier in a chat that could hold
  values are withheld when the conversation is sent again.
- **It is not anonymisation.** A summary of a very small group the data or
  you define (a narrow range, a filter, a rare category) can still reveal its
  cells' values. No minimum group size is set.
- **Other people's names in community annotation.** When an agent reads the
  community annotation of a dataset (`annotation_suggestions`), the
  suggested labels, their ontology terms and who suggested them (each person's
  GitHub login and the display name and title they publish) go to the agent's
  provider with the rest of what it reads (in the repository they are seen by
  everyone who can read it).
- **Text written by others is data.** Every result's data reaches an agent
  between `<cellucid-data>` marks (every `<` inside escaped), and agents are
  told never to take it as instructions.
- **Kept in this browser:** the conversations of the chat bar, with the
  history (`cellucid-history`); the chat's choices, permissions, the switch
  and the relay pairing in `localStorage`; API keys in `sessionStorage`
  (and `localStorage` with **Remember on this browser**). See
  {ref}`privacy-stored-locally`.

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
- **Writing back into your Python session is authenticated.** Interaction
  events and session bundles reach the notebook only over the viewer channel,
  which opens with the viewer's id and a per-viewer token generated fresh in
  the Python process and never written to disk. A notebook widget hands the
  token to its page in a message, posted to the app's exact origin; a
  browser-tab viewer carries it in its address, and its channel is one
  WebSocket (`/_cellucid/ws`). The server refuses
  a handshake from a page on another site (`Origin` must be the server's own
  host, a loopback name or a declared `allowed_hosts` name), and compares the
  token in constant time, so a wrong token cannot be discovered by measuring
  how long the rejection takes; an unknown viewer is compared against a
  same-length dummy so the failure path does not reveal the token length
  either. A socket that fails the check is closed before any event reaches
  your callbacks.

The practical reading: a stranger who reaches an exposed port can **read** your
dataset, but cannot **drive** your notebook.

(shared-machines)=
#### On a shared machine, `127.0.0.1` is every user of that machine

The loopback address keeps other *computers* out, not other *people*. On a
machine several people use at once — an HPC login or compute node, a
single-host JupyterHub (The Littlest JupyterHub), a multi-user RStudio Server
or Posit Workbench, a lab workstation with several accounts — every user of
that machine can connect to `127.0.0.1:<port>` and read everything a running
`cellucid serve`, `serve_anndata()`, `AnnDataViewer` or browser-tab
`cellucid_view()` serves, with `curl` or any program, for as long as it runs:
the dataset list, every array, the web app. Data routes carry no secret by
design (the browser reads them with plain range requests). What they cannot
do is drive your notebook or R session: the viewer channel needs the
per-viewer token.

When that matters: the data must not be readable by the other users of that
machine (unpublished patient-level data, data under an agreement that names
who may access it). Then use the notebook widget instead (`show_anndata()` or
`cellucid_view()` in Jupyter: no port is opened, the data travels over the
notebook's own authenticated connection), run the viewer on a machine only you
use, or open the file in the browser directly (**Open…**). A port that is
taken first by another user is a related risk: `cellucid serve` on its default
port 8765 is also where Google Drive sign-in is registered (see
{ref}`analytics` and the Google Drive notes above).

In the browser, only the hosted app's own pages (`https://www.cellucid.com` and
its bare apex) and pages on this machine (`localhost`, `127.0.0.1`, any port)
may read a server's data across origins (CORS); a page of any other site,
including other `cellucid.com` names, cannot.

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

Consequence for reverse proxies: a proxy in front of a browser-tab server (a
hub's proxy route) forwards
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
Even though its token expires after 8 hours and **Sign out** revokes it, best
practice is:
- use it on trusted deployments (official `cellucid.com` or your audited self-host)
- avoid signing in inside embedded/untrusted iframes

### 4) Report a security problem privately

Report a vulnerability to the contact in the
[security policy](https://github.com/theislab/cellucid/security/policy), not
in a public issue. The site publishes the same contact for automated tools at
[www.cellucid.com/.well-known/security.txt](https://www.cellucid.com/.well-known/security.txt)
(RFC 9116). The plain-language summaries for visitors are the site's
[privacy policy](https://www.cellucid.com/privacy.html) and
[terms of use](https://www.cellucid.com/terms.html).

---

## Institutional / IRB checklist (copy/paste)

Use this as a starting point for approvals and internal documentation.

1) **Where does the app run?**
   - `cellucid.com` (public host) vs self-hosted internal deployment vs a
     notebook widget or a Python/R viewer.
   - This answer determines whether Google Analytics is active at all — see
     {ref}`analytics`. Only the app's own pages on the public host report; a
     notebook widget (which frames www.cellucid.com) and a browser-tab viewer
     never do. For most review boards it is the first question.
2) **Where is dataset content stored?**
   - local disk only vs remote server vs a link's host (GitHub, a repository, a cloud drive)
3) **If the public host is used, what does it report about your data?**
   - No name, id, path or error text: only the loading method, outcome,
     duration, cell/gene count ranges, field count, dimensions and whether it
     has a graph or a vector field (a dataset id only for the public sample
     datasets). See
     {ref}`analytics`.
4) **What local persistence exists?**
   - `localStorage` preferences, the speed of each host data was read from
     (and the pace of one that blocked this network)
     **and Community Annotation session content**; in IndexedDB, the saved
     states (the history of every step on each dataset, with analysis results
     and conversations with agents, while **Autosave in this browser** is on),
     My datasets (names and
     public links of the datasets opened), handles of local files to reopen,
     decoder choices and the annotation file cache;
     converted R objects and MatrixMarket folders, whole-file downloads,
     downloaded parts of remote files (private Drive and OneDrive files
     included) and gene indexes in the origin-private file system (kept across
     visits until removed; Session › Stored › Manage…), short-lived
     sign-in tokens in `sessionStorage`, an AI provider's API key in
     `sessionStorage` (or `localStorage` with **Remember on this browser**),
     browser HTTP cache (app files, catalogs, sessions, signed-in Drive reads)
5) **If a local server is used, what is its bind address, and who else uses that machine?**
   - `127.0.0.1` (default) vs `0.0.0.0`. Reads are unauthenticated either way:
     on `127.0.0.1` every user of that machine can read the served data, on
     `0.0.0.0` anyone on the network ({ref}`shared-machines`).
6) **What artifacts will be shared?**
   - session bundles, `.cellucid-library` files, exported figures, URLs; and
     what redaction rules apply
7) **Will AI agents be used on the data?**
   - Which provider receives what agents read (Anthropic, OpenAI, or the
     provider behind a local agent app), and under which terms; whether each
     dataset is set to **Agents see names and summaries only**, knowing that
     this keeps values of single cells, extremes and pictures out but is not
     anonymisation ({ref}`privacy-agents`).

---

## Troubleshooting (privacy / storage / “I thought it was local”)

### Symptom: “I loaded a local folder, but I still see network traffic”

**Likely causes (ordered):**
1) You are using the hosted app (`cellucid.com`), so the app assets load over the network.
2) You are on a production hostname, so analytics beacons are being sent — see {ref}`analytics`.
3) You accidentally connected to a server or opened a link (dataset content truly loads over the network).

**How to confirm**
- DevTools → Network:
  - dataset requests are range requests to the file (`.h5ad`, `.zarr.zip`, …),
    a Zarr store's keys, or `_cellucid/anndata/…` on a server.
  - analytics traffic goes to Google’s tag manager and measurement endpoints.
  - Check `window.cellucidAnalyticsEnabled` in the console: `true` means this
    origin is a production host and is reporting.

**Fix**
1) If you need “no network after load”, use an offline/self-hosted deployment.
2) Ensure your source is “local-user” (**Open…** or a drop), not a server or a link.

---

### Symptom: “Community Annotation doesn’t work in private browsing / strict privacy mode”

**Likely causes:**
- The browser blocks `localStorage` and/or IndexedDB, which Community Annotation uses for caches and coordination.

**Fix**
- Use a normal browser profile (non-private window) for annotation sessions, or adjust the browser policy to allow site storage.

Related: {doc}`../j_community_annotation/03_ui_reference` (storage and failure modes).

---

### Symptom: “An exported PNG/SVG names my dataset or a private URL”

**Cause**
- Exported figures embed provenance metadata. For a local dataset that is
  `datasetUserPath`, which holds only the **name** of the folder or file you
  picked — never a path to it. For a remote one it is the source URL, which can
  be private.

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
