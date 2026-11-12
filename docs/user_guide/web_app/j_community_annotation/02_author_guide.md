# Community Annotation — Author Guide (Repo Setup + Operations)

This guide is for **dataset authors / maintainers** who want to run a community annotation round and keep it healthy at scale.

If you are an annotator (voting only), go to {doc}`01_annotator_guide`.

---

## What You’ll Do As an Author (Plain Language)

Community annotation has a simple division of labor:

- **Annotators** contribute suggestions, votes, and comments.
- **Authors** decide *what can be annotated*, *when the round is open/closed*, *how consensus is computed*, and (optionally) *how duplicate suggestions are merged*.

As an author, you will:

1) Pick a stable dataset id (critical).
2) Create an annotation GitHub repo (recommended: use the template).
3) Configure `annotations/config.json` to bind the repo to your dataset and annotatable columns.
4) Install the Cellucid GitHub App so the repo appears in the UI.
5) Connect from Cellucid, Pull, enable fields, and run the round.
6) Optionally moderate duplicates (merges).
7) Export a `cellucid-consensus.json` snapshot for downstream tooling.

---

## How To Read This Guide (Two Tracks)

::::{tab-set}

:::{tab-item} Quick Setup (Most Labs)

Follow these sections in order:

- Section 0 (checklists): avoid the common “we started too early” problems.
- Section 2 (dataset id): confirm your dataset id is stable.
- Section 3–6 (repo + app): create repo, configure `config.json`, install GitHub App.
- Section 7–9 (Cellucid UI): connect, enable fields, tune consensus settings.
- Section 12 (export): download `cellucid-consensus.json`.
- Section 13 (troubleshooting): when something breaks.
:::

:::{tab-item} Scaling / Governance (Large Groups)

Read everything, but pay special attention to:

- Section 0.1 (people/process): who merges PRs, who merges duplicates, and how decisions are communicated.
- Section 4 (what a file may hold, and what a Pull leaves out).
- Section 5 (GitHub settings): fork/PR flow, CI validation, and branch strategy.
- Section 9 (consensus parameters): how thresholds behave under disagreement.
- Section 10 (merges): when to merge vs when to leave disputed.
- Section 13 (rate limits, storage, and “it works for me” failures).
- *Limits*: what GitHub and Cloudflare bound, and what Cellucid cannot check.
:::

::::

---

## Starting point in the sidebar

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-01-disconnected.png
:alt: The Community annotation panel before any repository is connected: Repository with its (i), the line Not connected. and a Connect repo button, with the pointer on that button.
:width: 504px

Before repository setup, the Community Annotation panel reports that no
repository is connected and presents the explicit Connect repo action.
```

---

## 0) Before You Start (Author Checklist)

This section is intentionally “overkill”: most community-annotation failures come from skipping one of these.

### 0.1 Decide Your Round’s Rules (People + Process)

Before you open annotation, decide and write down:

- **Which columns are in scope** (one column at a time is easier for non-technical annotators).
- **What counts as “done”** (e.g., “≥5 voters and confidence ≥0.7”, plus manual review of disputed clusters).
- **Naming conventions** (e.g., “CD4 T cell” vs “CD4+ T”; “doublet” label; capitalization).
- **Evidence expectations** (marker genes? reference atlas? wet-lab rationale?).
- **Who moderates duplicates** (one or more maintainers/admins) and how often.
- **Timeline** (start/end date, reminders, when fields will be closed).

:::{tip}
If your annotators include wet-lab scientists, reduce cognitive load:

- Start with one coarse column (fewer categories).
- Provide a short “house style” label guide.
- Encourage comments as evidence rather than long debates in external chat.
:::

### 0.2 Dataset Checklist (Technical + Scientific)

- You have a **stable dataset id** (Cellucid uses `uns['cellucid']['id']` in the dataset file).
- You have at least one **categorical obs column** suitable for annotation (e.g. `leiden`, `cluster`, `cell_type_coarse`).
- The **category labels** for that column are stable (or you are ready to “freeze” them before annotation begins).
- You know your intended audience:
  - for broad groups, prefer coarse clusters and fewer categories;
  - for expert-only rounds, finer clusters can be appropriate.

:::{warning}
Changing any of the following after people have voted will fragment or invalidate prior work:

- dataset id
- annotatable column key (field name)
- category labels within that column (e.g. renaming clusters)

Cellucid intentionally locks category renaming/merging while annotation voting is enabled to prevent accidental breakage.
:::

### 0.3 GitHub Checklist (Repo + Permissions + CI)

- You can create a GitHub repository to store annotations (public or private).
- You (and/or your org) can install the **Cellucid GitHub App** on that repository’s owner (user/org).
- You decide how contributors will publish:
  - **Direct publish** (contributors have push access), or
  - **Fork + Pull Request** (contributors do not have push access).
- You have a plan for **branch consistency** (everyone must use the same branch).
- You enable validation CI (recommended): it rejects a broken file in a pull
  request before it is merged. CI cannot stop a direct write (a collaborator
  with push access publishes before CI runs), so Cellucid never trusts the
  repository to be valid: on Pull it leaves a broken file out on its own and
  says which one (Section 4).

:::{important}
Community annotation is designed so each contributor writes only their own user file. Authors should avoid manually editing `annotations/users/*.json` unless you are doing a targeted repair and you understand the schema.
:::

### 0.4 Dry Run (Highly Recommended)

Do a 10-minute dry run before inviting many people:

1) Use two browsers or two GitHub accounts.
2) Connect to the repo in Cellucid and Pull.
3) Make a few votes/suggestions.
4) Publish via your intended model (direct push or PR).
5) Pull again and confirm the merged view updates.

This catches:

- dataset id mismatch,
- wrong branch,
- app not installed on the repo owner,
- PR flow blocked by “forking disabled” policies,
- CI failures in `annotations/config.json`.

---

## 1) Understand the Three “Author Control Planes”

As an author you control annotation through three layers:

1) **The dataset** (what can be annotated)
   - Which categorical obs columns exist.
   - Whether categories are stable and meaningful.

2) **The annotation repo config** (`annotations/config.json`)
   - Which dataset ids are allowed.
   - Which categorical obs columns are annotatable (`fieldsToAnnotate`).
   - Per-column consensus rules (`annotatableSettings`: `minAnnotators`, `threshold`).
   - Whether annotation is temporarily locked (`closedFields`).

3) **Moderation merges** (`annotations/moderation/merges.json`, optional)
   - Used to merge duplicates so votes combine cleanly.

Everything else (suggestions, votes, comments) comes from the community (one file per user, with one section per dataset).

### What to edit (and what not to)

- ✅ Edit as author:
  - `annotations/config.json` (via UI or GitHub)
  - `annotations/moderation/merges.json` (via UI moderation; authors only)
- ❌ Avoid editing:
  - `annotations/users/*.json` (per-user data; conflict-free collaboration depends on “one user → one file”)

---

## 2) Choose a Stable Dataset ID (Critical)

Community annotation is scoped by dataset id. If the id changes, existing annotations will not appear (it becomes a different scope).

### How dataset id is determined

- For a file that states its id in `uns['cellucid']['id']`, Cellucid uses that
  id, whatever opened it: the sample menu, a data link, a server, a notebook,
  or **Open…** for a local file or folder. The same file
  keeps the same annotations everywhere.
- A file without one gets an id from how it was opened: a data link gives a
  stable one, such as `cellxgene:<dataset id>` or `geo:<GSE>:<file>`
  ({doc}`../b_data_loading/14_open_from_link`); a local file gets one from its
  format, its file name, its cell and gene counts and its gene names
  (`h5ad_atlas_1bd42f4e`). Writing the same data again (an observation column
  added) keeps that id; two files with the same name, cell and gene counts and
  gene list share annotations unless they carry their own id, which
  `cellucid.prepare()` and `cellucid_prepare()` write. See
  {doc}`../b_data_loading/06_dataset_identity_why_it_matters`.
- When writing the file with `cellucid.prepare(...)`, set `dataset_id=...`.

One repository can annotate several datasets. Each contributor's file keeps a
separate section per dataset, so two datasets that share a field and a
category label (for example `leiden` / `0`) never share suggestions, votes or
merges.

Example (recommended):

```python
from cellucid import prepare

prepare(
    "my_atlas.zarr",
    format="zarr",
    embeddings={"X_umap": adata.obsm["X_umap"]},
    obs=adata.obs,
    dataset_id="my_atlas_v1",  # keep this stable for the entire annotation round
    dataset_name="My Atlas (v1)",
)
```

How to confirm:

```python
import anndata as ad

ad.read_zarr("my_atlas.zarr").uns["cellucid"]["id"]  # 'my_atlas_v1'
```

:::{tip}
Treat dataset id like a **contract**:

- If you write the file again with small technical changes but the same clusters and meaning, keep the same id.
- If you change clustering, category labels, or biological meaning, use a new id and create a new round (or a new `supportedDatasets[]` entry).
:::

### Confirm the dataset id in the UI (recommended)

In Cellucid, the **Community annotation** panel's status block displays the dataset id (this is the id your `annotations/config.json` must match).

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-06-connected-status.png
:alt: The connected Community annotation panel: GitHub sync with the Dataset, GitHub and Repo rows and a GitHub sync… button, a Profile block naming the signed-in user with Edit and Clear, and the collapsed Annotated columns, Consensus column and Snapshot and cache groups.
:width: 504px

The status block at the top of the connected panel. **Dataset** is the id your
`annotations/config.json` must carry verbatim; **GitHub** and **Repo** complete
the scope your local annotation state is filed under.
```

---

## 3) Create the Annotation Repository

Cellucid expects a specific repository layout. The easiest path is to start from the template and customize it.

### Required layout (paths are case-sensitive)

```
annotations/
  config.json
  config.schema.json
  schema.json
  users/
    .gitkeep      (one LF byte; keeps the empty folder in Git)
    (one JSON file per contributor, ghid_<GitHub user id>.json)
  moderation/
    merges.json   (optional but recommended; authors publish it from the UI)
    merges.schema.json
.github/
  workflows/
    validate.yml  (recommended)
scripts/
  validate_user_files.py
.gitattributes
```

The browser requires both `annotations/config.schema.json` and
`annotations/moderation/merges.schema.json` even though the active
`annotations/moderation/merges.json` document itself is optional.

**The three schema files must be the template's, byte for byte.** Cellucid
identifies them by content: it reads each schema's Git blob id from the
repository and accepts the repository only when all three are exactly the
files published in
[theislab/cellucid-annotation](https://github.com/theislab/cellucid-annotation)
(user files version 3, merges version 3, config version 1). It never downloads
them. A schema that differs in any byte (an older version, an edit, other line
endings) stops Pull with a message that names the file and asks you to copy it
again. `.gitattributes` keeps every file's line endings LF on every platform,
so a checkout on Windows commits the same bytes; copy it with the schemas.

### Recommended “template-first” setup

The template repository is
[theislab/cellucid-annotation](https://github.com/theislab/cellucid-annotation).
Copy only its layout files, at the same paths:

- `annotations/` (with `annotations/users/.gitkeep`; replace the example
  `annotations/config.json` entry with your own)
- `scripts/validate_user_files.py`
- `.github/workflows/validate.yml`
- `.gitattributes`

Leave the template's own `tests/`, `scripts/generate_contract_corpus.py`,
`.github/workflows/tests.yml`, license and policy files behind: they test and
govern the template, not your annotation repository.

::::{tab-set}

:::{tab-item} GitHub UI (No Command Line)

1) Create a new GitHub repository (public or private).
2) Add the layout files listed above at the same paths, with their content
   unchanged (the schemas byte for byte).
3) Update `annotations/config.json` (Section 4).
4) Commit + push.

If you’re doing this in the GitHub web UI:

- create the folders/files exactly as shown,
- copy/paste file contents carefully (JSON is strict: no trailing commas);
  for the three schema files, upload the files themselves rather than pasting
  them, so not a byte changes.

:::

:::{tab-item} Git Users (CLI / IDE)

1) Create an empty repo on GitHub.
2) Clone it locally.
3) Copy the layout files into the repo root.
4) Edit `annotations/config.json`.
5) Commit + push.

Example on macOS/Linux (from a clone of `cellucid-annotation`):

```bash
cd /path/to/cellucid-annotation
rsync -aR annotations scripts/validate_user_files.py .github/workflows/validate.yml .gitattributes /path/to/your-annotation-repo/
cd /path/to/your-annotation-repo
git add -A
git commit -m "Initialize Cellucid annotation repo"
git push
```

`rsync -R` keeps each path relative to the template root, so the files land
at the same paths in your repository.

Windows PowerShell equivalent:

```powershell
Set-Location 'C:\path\to\cellucid-annotation'
$dest = 'C:\path\to\your-annotation-repo'
foreach ($p in 'annotations', 'scripts\validate_user_files.py', '.github\workflows\validate.yml', '.gitattributes') {
  $parent = Split-Path (Join-Path $dest $p)
  New-Item -ItemType Directory -Force -Path $parent | Out-Null
  Copy-Item -Recurse -Force $p $parent
}
Set-Location $dest
git add -A
git commit -m "Initialize Cellucid annotation repo"
git push
```

:::

::::

:::{note}
Keep this repo “boring”:

- Avoid large binaries (screenshots belong in documentation repos, not the annotation repo).
- Avoid data exports (the annotation repo should contain only JSON + scripts).
- Avoid rewriting history (force pushes) once annotation starts; it confuses caches and PR history.
:::

### Branch strategy (do not skip)

Everyone must be on the same `owner/repo@branch`, or they will appear to “disagree” because they are literally writing to different universes.

Common strategies:

- **Simple**: use `main` for the whole round.
- **Safer**: create a dedicated branch for the round (e.g. `round-2025-01`) and tell everyone to use it.

For large groups, a dedicated branch reduces accidental changes to `main` and makes it easier to archive rounds.

---

## 4) Configure `annotations/config.json`

`annotations/config.json` binds your annotation repo to one or more dataset ids and specifies which columns are annotatable.

### Minimal example (one dataset, one field)

```json
{
  "version": 1,
  "supportedDatasets": [
    {
      "datasetId": "my_atlas_v1",
      "name": "My Atlas (v1)",
      "fieldsToAnnotate": ["leiden"],
      "annotatableSettings": {
        "leiden": { "minAnnotators": 3, "threshold": 0.5 }
      },
      "closedFields": []
    }
  ]
}
```

### What each field means (and what can go wrong)

- `version` (must be `1`)
  - If you change this, validation will fail.

- `supportedDatasets` (must be a non-empty array)
  - Each entry is one dataset you want this repo to serve.
  - Dataset ids must be unique.

- `supportedDatasets[].datasetId`
  - Must match the dataset id of the dataset currently open in Cellucid: the id its file states in `uns['cellucid']['id']`, or the one its route gives it ({doc}`../b_data_loading/06_dataset_identity_why_it_matters`). When the repo does not list it, Cellucid names the id it looked for.
  - If the currently loaded dataset id is missing:
    - **annotators are blocked** (cannot Pull / view annotations)
    - **authors can still connect** (with a confirmation) and Publish an updated config to unblock everyone

- `supportedDatasets[].name`
  - Human-friendly name shown in the UI.

- `fieldsToAnnotate`
  - List of **categorical obs keys** (column names) that may be annotated.
  - Every listed key must exist as a categorical field in the loaded dataset.
    Any missing key fails the complete Pull without changing annotation state.
    Cellucid performs this check before selecting the raw-cache scope or
    downloading any user or moderation files, so the mismatch cannot mutate
    the raw-file cache.

- `annotatableSettings[fieldKey]`
  - Per-field consensus rules.
  - `minAnnotators` (whole number, 0 or more; the only upper bound is
    9,007,199,254,740,991, the largest integer a browser holds exactly):
    minimum unique voters required before a bucket can be anything other than
    “Pending”. A bucket nobody voted on is always “Pending”, even with
    `minAnnotators: 0`.
  - `threshold` (number -1..1): minimum `confidence` to reach “Consensus”.

- `closedFields`
  - Fields in this list are locked for annotators (no voting/suggestions/comments).
  - Validation rule: every closed field must also be in `fieldsToAnnotate`.

:::{important}
Validation rule (enforced by the template CI script):

- The keys in `annotatableSettings` must be exactly the same keys as
  `fieldsToAnnotate`; neither side may contain an extra or missing key.
- Every key in `closedFields` must also appear in `fieldsToAnnotate`.

If you violate this, GitHub Actions will fail and authors may be blocked from publishing updates cleanly.
:::

### How `threshold` behaves (do not guess)

Cellucid computes (per category bucket):

- `voters`: unique users who cast any vote in that bucket (across all suggestions)
- `netVotes`: for the current leading suggestion, `upvotes - downvotes`
- `confidence = netVotes / voters` (ranges from `-1` to `+1`)

Important edge cases:

- If `voters < minAnnotators` → status is **Pending** (even if there is a strong early leader).
- If nobody voted (`voters = 0`) → status is **Pending**, whatever `minAnnotators` and `threshold` say.
- Suggestions with exactly the same label count as one suggestion (Section 10); `voters` counts each person once.
- Suggestions are ranked by net vote first and upvote count second. A top tie
  exists only when suggestions have both equal net vote and equal upvote count;
  that tie is **Disputed**.

See Section 9 for worked examples and recommended defaults.

### Validate your repo inputs (recommended)

In the annotation repo (not in Cellucid), run:

```bash
python scripts/validate_user_files.py
```

This validates (Python 3.10 or newer; it reports every problem as an error,
never as a crash):

- that the three schemas are the canonical ones, byte for byte
- `annotations/config.json`
- every entry under `annotations/users/`: only `ghid_<id>.json` user files and
  the `.gitkeep` sentinel may be there (no other file, folder or link)
- `annotations/moderation/merges.json` (optional)
- that no path is a case variant of a layout path (such as
  `annotations/Users/`), which Cellucid never reads

If this fails, fix the file(s) it reports before inviting annotators.

### What a file may hold

- **No text has a length limit**: labels, evidence, ontology ids, marker
  genes, comments, merge notes, profile fields, dataset ids and names and
  field keys can be as long as the file allows. Every text is *exact*: its
  first and last characters are not whitespace, and it is never trimmed.
  Whitespace is one character class everywhere (the schemas, the validator and
  the app), the union of JavaScript's `\s` and Python's `str.isspace()`.
- **The only size bound is the file**: every active annotation document
  (`annotations/config.json`, each `annotations/users/ghid_<id>.json`, and
  `annotations/moderation/merges.json`) is at most 1,000,000 UTF-8 bytes,
  the most GitHub's repository Contents API returns as JSON
  ([GitHub: get repository content](https://docs.github.com/en/rest/repos/contents#get-repository-content)).
  One contributor's file holds every dataset they annotated in this
  repository, so this bound covers all of them together.
- **Integers** are at most 9,007,199,254,740,991.
- **JSON nesting** is at most 512 levels in both validators (Python's JSON
  reader recurses once per level and stops at the interpreter's recursion
  limit, 1000 by default); no valid annotation file nests deeper than 7.
- **Suggestion ids name their owner**: `ghid_<id>_<rest>`, where `ghid_<id>`
  is the proposer's own identity (`proposedBy`) and `<rest>` is exact text
  without `:` (Cellucid writes `ghid_<id>_<random UUID>`). Because the digits
  end at that `_`, no contributor can write an id that collides with, or
  claims to be, another contributor's suggestion. Votes, comments and merges
  may name anyone's suggestion; a deletion marker (`deletedSuggestions`)
  names only the file owner's own.

### Pull: what stops it and what is left out

A Pull has no limit on how many user files there are or on their total size
(GitHub's own listing limit is under *Limits* below). It reads the repository
at one commit, downloads only the files whose Git blob id changed since the
last completed Pull, keeps each batch in the browser's file cache as it
arrives (an interrupted Pull resumes there), and compiles the merged view by
streaming the cache one file at a time.

**One file never stops a Pull.** A user file or merges file that is not
valid (malformed JSON, over 1,000,000 bytes, failing its schema, a version
other than 3) is left out of the consensus on its own: Cellucid keeps the
reason with the file's blob id (so it is not downloaded again until it
changes), shows one warning that says how many files were left out, and lists
each file with its reason under the warning's details. For the open dataset,
Cellucid reads a user file's profile and that dataset's section only, so a
broken section of another dataset does not leave the file out here. A path
under `annotations/users/` that is not a `ghid_<id>.json` file, and a case
variant of a layout path (`annotations/Users/`), is ignored and listed the
same way. If your own published file is the one left out, the warning says so
first: Publish will not overwrite it until you fix or remove it on GitHub.

**These stop a Pull**, because nothing can be compiled without them:

- a schema file that is not the canonical one, a missing
  `annotations/config.json` or `annotations/users/` folder, or a missing branch;
- an invalid `annotations/config.json`;
- a listing of `annotations/users/` that GitHub reports as truncated;
- a configured categorical field the open dataset does not have. That check
  fails the complete Pull without changing annotation state, and it runs
  before selecting the raw-cache scope or downloading any user or moderation
  files, so the mismatch cannot change the file cache.

A Pull that stops leaves the annotation view as it was.

---

## 5) Configure GitHub Repo Settings (Highly Recommended)

### Decide how annotators will publish

You have two viable models:

1) **Direct publish** (annotators have write access)
   - Pros: simplest experience (Publish writes directly to `annotations/users/ghid_<id>.json`)
   - Cons: requires adding many people as collaborators; less review control.
     Anyone with push access can also write any file with Git or GitHub's web
     editor, including another contributor's `ghid_<id>.json`: Cellucid's
     Worker refuses such a write made through Cellucid, but Cellucid reads
     files by content and cannot tell who committed them. Give push access
     only to people you trust with every file, and review the commit history
     when in doubt.

2) **Fork + Pull Request publish** (annotators do not have write access)
   - Pros: reviewable contributions, no direct writes to your repo
   - Cons: you must merge PRs; depends on forking being allowed and not blocked by org policy

Cellucid chooses the best option per user:

- If the user can push → direct publish
- Else if the repo allows forking → fork + PR publish
- Else → user cannot publish (they can still vote locally, but nothing can be shared)

In other words, Cellucid uses a fork + Pull Request only when the repository
allows forking. It selects exactly one route from repository metadata before
mutation; a route failure is terminal and does not switch to the other route.

### Fork + PR model: what a contributor does first

Cellucid never creates a fork. A contributor without write access does two
things on GitHub before their first Publish:

1. **Fork the repository**: `https://github.com/<owner>/<repo>/fork` (any
   name; Cellucid finds a renamed fork among the repository's forks).
2. **Install the Cellucid GitHub App on their own account with access to the
   fork**: *All repositories*, or *Only select repositories* with the fork
   selected (a fork made after the installation is added under *Settings* →
   *Applications* → *Installed GitHub Apps* → *Configure*). A user token
   reaches a repository only in an account where the app is installed
   ([GitHub: authenticating on behalf of a user](https://docs.github.com/en/apps/creating-github-apps/authenticating-with-a-github-app/authenticating-with-a-github-app-on-behalf-of-a-user)),
   and the head branch and the user file are written in the fork. Opening the
   pull request uses the app's installation on your repository.

What Publish says otherwise:

- No fork: `You have no fork of <owner>/<repo> yet: fork it on GitHub, then
  Publish again.`, with a **Fork on GitHub** link to
  `https://github.com/<owner>/<repo>/fork` (also under *Details*). Nothing is written.
- A fork GitHub is still copying (GitHub creates forks asynchronously):
  `GitHub is still copying your fork <login>/<repo>; Publish again in a
  minute.` Cellucid asks once and does not wait.
- The app cannot write to the fork: `Cellucid cannot write to your fork
  <login>/<repo>: install Cellucid's GitHub App on your account with access to
  it, then Publish again.`, with an **Install the app** link.

How a fork publish works, so the pull requests you receive make sense:

- Each contributor has **one pull request per base branch**, from one head
  branch in their fork named `cellucid-annotations/ghid_<id>/…`. It carries
  their whole user file, so every dataset they annotate in this repository
  goes into the same pull request.
- On every Publish, Cellucid starts from the file in that open pull request
  (it holds sections not merged yet) or, without one, from the file on your
  branch, replaces the open dataset's section, then **resets the head branch
  to your branch's current tip** and writes the file on it. The pull request's
  diff is that one user file, and a pull request you squash-merged earlier
  never makes the next one conflict. An open pull request is reused; a new one
  is opened only when none is open.

### Enable validation CI

The template includes:

- `scripts/validate_user_files.py`
- `.github/workflows/validate.yml`

The shipped `validate.yml` runs one command, `python scripts/validate_user_files.py`,
on every push and pull request, with read-only repository access
(`permissions: contents: read`). It needs nothing but the layout files, so it
works in a repository that copied only those.

Suggested GitHub settings:

- Require the validation check to pass before merging PRs.
- For fork-based contributions, allow GitHub Actions to run on PRs (org policies may apply).

CI is a gate for pull requests and a report for everything else: a direct
publish lands before CI runs. A file CI would reject is left out on Pull with
its reason (Section 4), never silently used.

### Branch protection (advanced, but important)

Branch protection can break direct publishing:

- If direct pushes are blocked, users with “write” permissions may still see Publish fail.

There is no UI route override. Branch protection alone does not select fork +
Pull Request: when repository metadata reports push permission, Cellucid
selects direct publication and reports a terminal error if GitHub rejects it.

For direct publication, configure the connected branch to accept the exact
Cellucid commit. For an annotator to use fork + PR instead, that annotator must
have no source-repository push permission and the repository must allow
forking. Authors with maintain/admin permission remain on the direct route; use
GitHub outside Cellucid for an all-PR governance policy that blocks their direct
commits.

---

## 6) Install the GitHub App (Required for Repo Discovery)

Cellucid’s UI lists only repositories where the **Cellucid GitHub App is installed**.

1) Install the app on the user/org that owns the annotation repo.
2) If you choose “Only select repositories”, make sure the annotation repo is selected.

:::{note}
Org repos often require an org admin to approve the installation.
:::

### Optional: self-host the GitHub OAuth + API proxy (org deployments)

Cellucid's community annotation UI uses a small server component (the
checked-in deployment is a Cloudflare Worker) to run GitHub App OAuth, proxy
only the repository operations used by Cellucid, run two GraphQL reads it
writes itself, and relay fixed CAP persisted queries.

If you are using **cellucid.com**, you typically do *not* need to do anything here.

If your organization requires owning the auth infrastructure (recommended for many orgs), you can self-host.

#### Exact version-3 route contract

The Worker exposes the following ordered routes:

- `/auth/login`
- `/auth/callback`
- `/auth/revoke`
- `/auth/user`
- `/auth/installations`
- `/auth/installation-repos`
- `/cap/lookup-cells`
- `/graphql/snapshot`
- `/graphql/blobs`
- `/api/repos/*`

Opening `/` must return the exact compatible service identity and ordered route
inventory, beginning
`{"status":"ok","service":"Cellucid GitHub Auth","contractVersion":3,`
and continuing with an `endpoints` array containing the routes above. Before
OAuth or token disclosure, the browser checks that health document and stops
on a stale or incompatible Worker. An older Worker is reported as "The
community-annotation service at `<origin>` runs an older version (contract N)
than this Cellucid (contract 3); its owner must redeploy the Worker."; a newer
one asks you to reload Cellucid. A Worker that cannot be reached at all is
reported as unreachable, not as incompatible.

What the routes do, in short (the setup guide has the exact contract):

- `/graphql/snapshot` reads a branch's commit and the listings of
  `annotations/` and `annotations/moderation/` in one GraphQL query the
  Worker writes itself; `/graphql/blobs` reads many files by blob id in one
  query. The browser never sends a GraphQL document.
- `/auth/installations` and `/auth/installation-repos` answer one GitHub page
  (100 items) per request (`?page=N`), so a list of any length is read page by
  page.
- `/api/repos/*` relays the repository requests Cellucid uses, among them a
  users folder read by its tree id and the forced move of a contributor's own
  pull-request head branch (`PATCH git/refs/heads/cellucid-annotations/ghid_<id>/…`
  in the contributor's own fork).
- Every write is bound to the signed-in user, whom the Worker looks up at
  GitHub first: a user file `annotations/users/ghid_<n>.json` only by user
  `n`, a branch under `cellucid-annotations/` only as that user's own
  `cellucid-annotations/ghid_<id>/…` (created or moved only in a repository the
  user owns), a pull request only from such a branch. A refused write is
  `403` with outcome `not-applied`.
- Every answer that relays GitHub carries GitHub's `ETag`, `Retry-After` and
  `X-RateLimit-*` headers, so the browser paces its requests by GitHub's
  budget; reads accept `If-None-Match`, and a `304` costs no transfer. Every
  route caches a CORS preflight for 7200 seconds
  (`Access-Control-Max-Age: 7200`, the longest Chromium keeps one).
- A request has no total deadline: it is stopped (`504`) only when nothing
  moved for 13 seconds (GitHub's own 10-second request timeout plus a
  connection set-up), so a slow but steady upload or download is never cut.

#### Exact Worker and GitHub App configuration

The Worker requires exactly `ALLOWED_ORIGINS`, `GITHUB_CLIENT_ID`, and
`GITHUB_CLIENT_SECRET`. It does not use an App id, private key, App JWT, or
installation access token. It deploys on the Cloudflare Workers Free plan
(about 10 ms of CPU and 50 subrequests per request). The checked-in Worker
stays well inside both at any repository size: it relays file contents, trees
and GraphQL answers without parsing them, checks a write with native JSON
parsing, and needs at most two GitHub subrequests per request. A request that
did exceed the CPU budget would be answered by Cloudflare with an HTML `503`
page without CORS headers, which the browser can only report as an unreachable
Worker; the setup guide explains how to diagnose that and when to move to
Workers Paid.

On the Free plan, abuse of the Worker cannot create a bill: it can only use up
the account's daily request allowance (100,000 requests, reset at midnight
UTC), after which Cloudflare answers with an error until the next day.
Cloudflare's rate-limiting rules belong to a zone (a domain you have on
Cloudflare), so they cover a custom domain and never a `*.workers.dev`
hostname. cellucid.com's Worker is served at `https://api.cellucid.com` (its
GitHub App callback is `https://api.cellucid.com/auth/callback`); a
self-hosted deployment that wants a limit serves its Worker on its own
Cloudflare domain the same way and adds a WAF rate-limiting rule for `/cap/*`
there
([Workers limits](https://developers.cloudflare.com/workers/platform/limits/), [rate limiting rules](https://developers.cloudflare.com/waf/rate-limiting-rules/)).

Configure the GitHub App with **Expire user authorization tokens** kept on
(GitHub's default for a new app; an existing app opts in under *Optional
Features* → *User-to-server token expiration*). The Worker accepts only an
expiring user token and refuses a token that never expires, so sign-in fails
with token expiration disabled. Also select **Request user authorization
(OAuth) during installation**, disable webhooks, leave account permissions
unset, and set these repository permissions:

- Metadata: read-only
- Contents: read and write
- Pull requests: read and write

No Administration permission: Cellucid never creates a fork. Creating one
through the API needs Administration write
([GitHub: permissions required for GitHub Apps](https://docs.github.com/en/rest/authentication/permissions-required-for-github-apps)),
so a contributor without write access forks the repository on GitHub
themselves, and Cellucid finds that fork when they publish.

Two values are compiled into the Cellucid client, so a self-hosted deployment
builds its own copy of the site:

- the Worker origin, `DEFAULT_WORKER_ORIGIN` in
  `assets/js/app/community-annotations/github-auth.js`. A different runtime
  override is accepted only from a recognized local-development host. The
  exact Cellucid page origin must also be present in `ALLOWED_ORIGINS`.
- the GitHub App's slug in the **Add repo** link of the connection wizard,
  `https://github.com/apps/cellucid-community-annotations/installations/new`
  (`GITHUB_APP_INSTALL_URL` in
  `assets/js/app/ui/modules/community-annotation/github-connection-flow.js`;
  Publish's **Install the app** link uses it too).
  Set it to your own app's slug; otherwise **Add repo** opens the installation
  page of cellucid.com's app, not yours.

For commands, precise origin rules, the complete health JSON, and live
verification, follow `cellucid/docs/github-oauth-cloudflare-setup.md` in the
Cellucid web repository. A source edit does not update a running Worker; deploy
the accepted Worker source and then verify the live root and browser lifecycle.

---

## 7) Connect the Repo From Cellucid (Author Bootstrap)

1) Load your dataset in Cellucid.
2) Open the **Community annotation** accordion.
3) Click **Connect repo**. If this dataset already has a repository connection,
   use **GitHub sync…**.
4) Complete **Continue with GitHub**, enable/reload repositories if needed,
   select the repo, click **Connect**, then **Pull latest**.
5) Confirm you see author-only controls (e.g. **Annotated columns**).

### Dataset mismatch (the most common “why can’t annotators Pull?” issue)

If the dataset loaded in Cellucid is not present in `annotations/config.json`:

- Annotators are blocked (they cannot Pull).
- Authors can still connect (with a warning) and then **Publish** to write an updated config that adds/updates `supportedDatasets[]` for the current dataset id.

This “author override” exists to make first-time bootstrapping smooth.

### “Am I actually an author?” (role sanity check)

Role is derived from GitHub permissions on the annotation repo:

- author = `maintain` or `admin`
- annotator = everything else

If you are an author, you should see author-only UI blocks such as **Annotated columns**.

If you do not:

1) Confirm you are signed into the expected GitHub account in the GitHub sync modal.
2) Confirm your permission level on the repo is maintain/admin.
3) Disconnect/reconnect and Pull again.

---

## 8) Enable/Disable Annotatable Columns (Author UI)

Once connected, open **Annotated columns** inside the **Community annotation** accordion:

1) Select the categorical obs field you want to control (dropdown labeled **Categorical obs**).
2) Click **Add** to include it in annotation.
3) Optionally adjust consensus settings (Section 9).
4) Click **Publish** so others receive the settings on Pull.

To stop annotation on a field:

- **Close** locks voting/suggestions/comments for annotators (you can reopen later).
- **Remove** removes it from the annotatable list entirely.

:::{warning}
Once annotation is enabled for a categorical field:

- category renaming and category merging are disabled in the legend UI
- field renaming is disabled

Plan your cluster names and field keys before opening the annotation round.
:::

### Choosing which columns to open (practical guidance)

For mixed audiences (computational + wet-lab):

- Start with one column that has a manageable number of categories (e.g., 10–50).
- Avoid columns that are “not biological” (e.g., `batch`, `donor`) unless your project specifically wants that.
- Use a stable clustering label column rather than something that changes with filtering.

---

## 9) Tune Consensus Rules Per Column

Inside **Annotated columns**, after you select a column that is already annotatable, you will see **Consensus settings**:

- **Threshold** slider (maps to `threshold` in `annotations/config.json`)
- **Min annotators** input (maps to `minAnnotators`)
- **Apply** (apply locally)
- **Reset** (discard local edits)

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-07-manage-annotation.png
:alt: Annotated columns open in the Community annotation panel: the Categorical obs picker on the voting column clusters with Add, Remove and Close, and its consensus settings with a Threshold slider and a Min annotators field.
:width: 504px

**Consensus settings** appear only once you select a column that is
already annotatable — the block belongs to the selected column, not to the
repository. The **Threshold** slider is in whole percent; the stored value is a
signed fraction spanning −1 to +1, because confidence itself is signed.
```

After applying locally, **Publish** to write the settings to `annotations/config.json` so all annotators receive them on Pull.

### Settings created when a field is enabled

The keys in `annotatableSettings` must be exactly the same keys as
`fieldsToAnnotate`. When an author enables a field in the UI, Cellucid creates
its required entry with:

- `minAnnotators = 1`
- `threshold = 0.5`

You may then edit these values before Publish. A config file may not omit the
settings entry for an enabled field.

### Worked examples (to build intuition)

Per category bucket, Cellucid computes:

- `voters` = number of unique users who voted on any suggestion in this category
- `netVotes` = (`upvotes - downvotes`) for the current leading suggestion
- `confidence = netVotes / voters`

Examples:

| Votes in the bucket | voters | Best netVotes | confidence | Notes |
|---|---:|---:|---:|---|
| 1 user upvotes “B cell” | 1 | 1 | 1.0 | Still **Pending** if `minAnnotators > 1` |
| 3 users upvote “B cell”, nobody downvotes | 3 | 3 | 1.0 | Strong consensus |
| 3 users: 2 upvote “B cell”, 1 downvotes “B cell” | 3 | 1 | 0.33 | Often disputed unless threshold is low |
| 4 users: 3 up, 1 down | 4 | 2 | 0.5 | Exactly meets `threshold=0.5` (if not tied) |
| 4 users split: 2 upvote A, 2 upvote B | 4 | 2 | 0.5 | **Disputed** because net votes and upvotes both tie |

:::{important}
`threshold` is not “percent upvotes”.

It is a net-support share: `(upvotes - downvotes) / voters`.

Downvotes matter and reduce confidence.
:::

### Practical parameter choices

Starting points:

- Conservative rounds (high confidence): `minAnnotators=5`, `threshold=0.6–0.8`
- Fast iteration rounds: `minAnnotators=2`, `threshold=0.5`

If you expect disagreement (hard biology, rare types):

- keep `minAnnotators` relatively high (so one loud person can’t dominate),
- use comments as evidence,
- accept that some clusters will remain Disputed and require manual review.

---

## 10) Merge Duplicate Suggestions (Moderation)

Even with careful contributors, duplicate labels happen (e.g. `T cell` vs `T-cell`, synonyms, or capitalization variants).

Suggestions with **exactly the same label** need no merge. Inside a category,
Cellucid treats them as one suggestion: labels are compared after Unicode
normalization (NFC), so the same text typed with composed or decomposed
accents counts as equal, but nothing else is folded (case, spacing and
punctuation still differ). The earliest proposal is shown as the suggestion,
and each person's votes across the equal labels count once (their majority; a
tie counts as no vote). Cellucid also refuses to add, or rename a suggestion
to, a label already shown in that category. Two people proposing "T cell"
independently therefore support one label instead of splitting it into a tie.

For labels that differ but mean the same thing, an author merges them so votes are combined:

1) Open a category’s voting modal.
2) Drag one suggestion card onto another.
3) Confirm the merge and optionally write a merge note. Like comment edits, a
   merge-note edit in progress survives the modal refreshing around it.
4) **Publish** (authors publish to `annotations/moderation/merges.json`).

Merge behavior (what annotators experience):

- Votes are combined and **de-duplicated** (at most one vote per user in the merged bundle total).
- The UI keeps originals accessible via **View merged** (so evidence and history aren’t lost).

### When to merge vs when not to merge

Merge when the suggestions are truly the same meaning:

- formatting variants (`T cell` vs `T-cell`)
- spelling variants
- the same ontology id with different labels
- one label is a strict synonym of the other in your project

Do not merge when meaning differs:

- “CD4 T cell” vs “T cell” (one is more specific)
- “Doublet” vs “Cycling T cell” (different biological claims)
- unresolved debates (leave Disputed; use comments)

### What is stored in `merges.json` (advanced)

`annotations/moderation/merges.json` (version 3) is an author-only log of
merges: `{"version": 3, "updatedAt": "<UTC timestamp>", "merges": [...]}`.

Each merge entry includes:

- `datasetId`: the dataset the merge belongs to. One repository can serve
  several datasets that share a field and a category label, so a merge
  applies only to the bucket of its own dataset.
- `bucket`: which category bucket the merge applies to. Its format is
  `<fieldKey>:<categoryLabel>` when the field key contains no colon. When it
  does, Cellucid uses
  `fk~${encodeURIComponent(fieldKey)}:<categoryLabel>`; the `fk~` prefix and
  percent-encoding make the embedded colon unambiguous.
- `fromSuggestionId` → `intoSuggestionId`: the mapping (merge “from” into
  “into”). Both are suggestion ids in the owner-prefixed form
  `ghid_<proposer id>_<rest>`; they may name any contributor's suggestions.
- `by`: author identity (stored as `ghid_<githubUserId>`)
- `at`: timestamp
- `editedAt` (optional): if you later edit the merge note
- `note` (optional): a short rationale shown in the UI

Example of one entry (illustrative):

```json
{
  "datasetId": "my_atlas_v1",
  "bucket": "leiden:7",
  "fromSuggestionId": "ghid_7_6f1c2c0e-3d1b-4e37-9a83-2f0b54c7d1a9",
  "intoSuggestionId": "ghid_8_0b9e4a51-7c2d-4f8e-b3a6-5d1e9c2f7a40",
  "by": "ghid_12345",
  "at": "2025-01-01T12:34:56Z",
  "note": "Merged formatting variants: 'T cell' and 'T-cell'."
}
```

### Undoing a merge (if you merged the wrong thing)

Open **View merged** on the target suggestion. Authors can detach a merged
member there. Only the GitHub identity that created a merge record can edit or
delete its note. Publish afterward to write the changed merge log. For recovery
outside the UI, revert the commit that changed
`annotations/moderation/merges.json`; direct JSON editing is advanced and must
be validated before Pull.

---

## 11) Build a Consensus Column (Optional, Local)

Cellucid can build a local derived categorical obs column for visualization:

1) Open **Consensus column**.
2) Choose an annotatable source column (e.g. `leiden`; a column you renamed
   here is listed by its name here, and its votes are read under its dataset
   name).
3) Choose a **New column key** (it starts as `community_<column>`).
4) Threshold / min annotators start at the column's rule in the repository
   (`annotations/config.json`); a change applies to this derived column only
   and is remembered for the column.
5) Build the derived column.

Result:

- Each category becomes a label (if consensus), or `Disputed`, or `Pending`.
- This does **not** change voting rules and does **not** publish anything to GitHub.

A built column is a step of the history (`field.consensus`), like any other
change: Undo takes it back, and a session or a shared history carries it. The
step records, for each category of the source column, its outcome (the agreed
label, disputed or pending, with its voter count), the rule it was built with
(threshold and minimum annotators) and where the votes came from (repository,
branch and when it was last pulled); the cells' labels are derived from the
source column when the step is applied. So it replays by name: on another
dataset with a column and categories of the same names, each category gets
its recorded outcome. A step recorded before this version kept one code per
cell, so it is put back only on the dataset it was built on.

A column built is a snapshot of the consensus at that moment. When a Pull (or
a vote) changes the consensus of its source column, the panel says
`<column>: consensus changed since built (N categories)` with **Rebuild**,
which builds it again under its name with the rule it was built with.

A category agreed as literally `Pending` or `Disputed` stays a label of its
own; the status category is then named `Pending (status)` (or `Disputed
(status)`).

Categories that agree on the same label are one category of the derived
column, compared as suggestions are (after Unicode NFC normalization; case,
spacing and punctuation are not folded), and it is spelled as the first of
them spells it.

This is useful to:

- color the atlas by the current community consensus
- quickly spot which clusters still need attention

:::{note}
The derived column uses the threshold/minAnnotators set in the derived-column
UI (the repository's rule unless you changed it); the history step records
them, and the methods paragraph of the history states them.
:::

---

## 12) Export a Consensus Snapshot (For Downstream Tools)

Inside **Snapshot and cache**:

- **Snapshot and cache** → **Download** (`cellucid-consensus.json`)

This produces a JSON snapshot built in your browser from the locally cached raw GitHub files (it is not written back to GitHub).

Recommended author workflow:

1) Pull latest (ensure you have everyone’s newest files)
2) Download `cellucid-consensus.json`
3) Use it downstream (e.g., build an official cell-type column)

### What is inside `cellucid-consensus.json`? (Structure)

The snapshot contains:

- `suggestions`: merged suggestion cards per bucket, including `upvotes`/`downvotes` arrays
- `consensus`: per-bucket summary objects with:
  - `status`: `pending` | `disputed` | `consensus`
  - `label`: best label (or comma-joined labels in ties)
  - `confidence`: number in -1..1
  - `voters`: unique voter count
  - `netVotes`: best net vote count
  - `suggestionId`: winning suggestion id (null in ties)

Buckets use `<fieldKey>:<categoryLabel>` when the field key contains no colon.
For a field key containing `:`, Cellucid uses
`fk~${encodeURIComponent(fieldKey)}:<categoryLabel>`. The category label
follows the first delimiter verbatim. Treat `fk~` as the escape prefix only
with `%3A` in the prefixed component; a colon-free field key may itself
begin with `fk~`.

### Example downstream usage (computational)

This is one simple pattern: map consensus labels back onto an `AnnData` cluster column.

```python
import json
from urllib.parse import unquote

import pandas as pd

doc = json.load(open("cellucid-consensus.json"))

target_field = "leiden"

mapping = {}
for bucket, summary in doc["consensus"].items():
    field_key_part, category_label = bucket.split(":", 1)
    encoded_field_key = field_key_part[3:]
    field_key = (
        unquote(encoded_field_key)
        if field_key_part.startswith("fk~") and "%3a" in encoded_field_key.lower()
        else field_key_part
    )
    if field_key != target_field:
        continue
    if summary.get("status") == "consensus":
        mapping[str(category_label)] = summary.get("label")

# Example: apply to a Series (or adata.obs[target_field])
cluster = pd.Series(["0", "1", "7", "7"], name=target_field)
cluster_consensus = cluster.astype(str).map(mapping)
```

Notes:

- For disputed or pending buckets, choose an explicit display label such as
  `"Disputed"`, or leave the value missing.
- If your cluster labels are integers, cast consistently to strings.

---

## 13) Author Troubleshooting (Massive)

If you don’t find your issue here, also check {doc}`03_ui_reference` (it includes additional UI-specific guidance and error messages).

### Before you debug: capture your “scope”

Most “it doesn’t work” reports are scope mismatches. When troubleshooting, always record:

- dataset id (from the status panel)
- repo + branch (`owner/repo@branch`)
- your GitHub login (which account you’re signed into)
- the exact error message text (copy/paste if possible)

For non-technical collaborators, a screenshot of the status panel + error message is often the fastest way to debug.

### Repo setup / structure

- **“Repo missing annotations/config.json / annotations/schema.json / annotations/users/”**
  - Cause: repo not created from template or paths renamed.
  - Fix: ensure the required layout exists exactly (case-sensitive).

- **“annotations/… is not the schema this Cellucid reads (Git blob …, expected …)”**
  - Cause: one of the three schema files differs from the template's: an
    older version (for example the version-2 schemas), an edit, or changed
    line endings.
  - Fix: copy the named file again, unchanged, from
    [theislab/cellucid-annotation](https://github.com/theislab/cellucid-annotation),
    together with `scripts/validate_user_files.py` and `.gitattributes`, and
    commit it.

- **Pull warns that files were left out of the consensus**
  - Cause: those files are not valid annotation files (the warning's details
    list each one with its reason), or paths under `annotations/` are not
    annotation files.
  - Fix: run `python scripts/validate_user_files.py`, then fix or remove the
    files it names. Everything else in the repository is used meanwhile.
    A file written for the earlier contract (`"version": 2`, suggestion ids
    without the owner prefix) is one of them. There is no migration: Publish
    never overwrites a contributor's own file that is not valid, so the file
    is removed on GitHub and its contributor publishes again from the current
    Cellucid; the old file's suggestions and votes are not carried over.

- **CI validation fails immediately**
  - Cause: invalid JSON, wrong field types, or policy rules (e.g. `annotatableSettings` contains keys not in `fieldsToAnnotate`).
  - Fix: run `python scripts/validate_user_files.py` locally in the annotation repo and follow the error output.

- **Annotators report “Pull works but nothing shows up”**
  - Common causes:
    - They are on the wrong repo/branch
    - Dataset id mismatch (see next section)
    - They haven’t Published (their local work isn’t shared)

### Dataset mismatch / blocked annotators

- **Annotators cannot Pull and see a dataset mismatch error**
  - Cause: the current dataset id is not in `annotations/config.json`.
  - Fix: as author, connect anyway and Publish; this updates `supportedDatasets[]` and unblocks annotators.
  - Prevention: add the dataset id to config before inviting annotators.

### “I’m an author but the UI says I’m not”

- Cause: author role is derived from GitHub permissions (`maintain` or `admin`).
- Fix: ensure your GitHub account has maintain/admin on the annotation repo, then reconnect and Pull.
- If role remains “unknown”: this usually indicates a GitHub API reachability/auth issue; see GitHub auth troubleshooting below.

### GitHub App install / repo not appearing

- **Repo does not appear under “Choose repository”**
  - Causes:
    - the Cellucid GitHub App is not installed for the repo owner
    - the app was installed for “Only selected repositories” and the repo is not selected
    - you are signed into a different GitHub account than expected
  - Fix:
    - install/adjust the app installation and try **Reload** in the GitHub sync modal
    - verify your GitHub username in the modal matches your intended account

- **Fork + PR flow fails for annotators**
  - "You have no fork of …": they fork the repository on GitHub (the message links to it), then Publish again. Cellucid never creates a fork.
  - "Cellucid cannot write to your fork …": the GitHub App is not installed on their own account with access to the fork. Fix: install it there (*All repositories*, or *Only select repositories* with the fork selected), then Publish again.
  - "GitHub is still copying your fork …": wait a minute and Publish again.

### Publish failures (authors)

- **Publishing fails with “Sign in required.”**
  - Fix: sign in again (a token expires after 8 hours; closing the tab clears it).

- **Publishing fails but you have write access**
  - Common causes:
    - branch protection blocks direct writes
    - required status checks are configured but GitHub API rejects direct commit
  - Cellucid does not retry as a fork or PR. Adjust the connected branch rules
    to permit the direct commit, or make the change through an external GitHub
    workflow.

- **Publishing fails for annotators and you disabled forking**
  - Cause: users without push cannot publish if `allow_forking` is disabled.
  - Fix: enable forking, or grant write access to annotators.

- **“Your published file `annotations/users/ghid_<id>.json` … is not a valid annotation file”**
  - Cause: the contributor's own file on GitHub fails the contract (for
    example a version-2 file). Cellucid composes a Publish onto that file and
    never overwrites one it cannot read.
  - Fix: fix or remove the file on GitHub, then Pull and Publish again.

- **“Remote user file changed since your last Pull” / “… changed on GitHub while publishing”**
  - Cause: the contributor's file on GitHub is newer than the one their last
    Pull read (for example after publishing from another device).
  - An annotator's Publish pulls (merging the two) and publishes again by
    itself. It stops with this message only for an author whose repository
    settings or merges are not published yet, since a Pull would bring the
    repository's over them.
  - Fix: note the staged settings, Pull latest, check them, then Publish
    again. Cellucid never replaces a remote file it has not read.

### Pull is slow / rate-limited

What a Pull costs: when nothing under `annotations/` changed since your last
completed Pull, one or two requests (the repository metadata, usually a
`304` that transfers nothing, and the snapshot). A first Pull reads the
changed files by blob id in batches packed by their sizes (a batch doubles
while GitHub answers and halves when GitHub refuses one as too large or too
slow), so the number of requests follows the bytes to read, not the number of
files: 5,000 user files of about 1 KB each take about six batch requests.

- Causes:
  - GitHub's hourly request budget of the signed-in user is spent (see
    *Limits* below)
  - frequent auto-pulls across many users
  - a slow or unsteady connection
- What Cellucid does: every request is paced by the budget GitHub reports
  (`X-RateLimit-*`, `Retry-After`). When GitHub names a wait, the progress
  notification says so and the Pull continues after it; cancel it to stop
  (what was read is kept). A Pull that still meets the limit stops with a
  `GitHub rate limit:` message, and the next Pull resumes from what was read. A failed request is asked again; nothing is
  stopped by a fixed time limit, only by a connection on which nothing
  arrives at all.
- Mitigations:
  - keep the annotation repo “clean” (only JSON + scripts; avoid large binaries)
  - ask annotators to Pull on demand instead of using aggressive auto-pull intervals

### Local cache corruption / storage restrictions

- **Error that the local raw cache is unavailable**
  - Cause: browser storage policies (private mode, strict settings, embedded iframe restrictions).
  - Impact: the files a Pull read are kept in IndexedDB; annotation state and
    the cross-tab lock use local storage. The persistent raw-file cache
    requires IndexedDB. Cellucid never switches to an in-memory cache. Pull
    fails without replacing the current annotation view.
  - Fix: use a normal browser profile, allow site storage, avoid restrictive privacy modes for the annotation session.

- **Error that local annotation state was written by another Cellucid version**
  - Cause: the local annotation state of this repository was saved by an
    earlier Cellucid in a format this one does not read.
  - Fix: clear that repository's local annotation state (**Clear session**
    under **Snapshot and cache**) and Pull again. Published annotations stay
    on GitHub; local work that was never published is not carried over.

- **Error about local cache being corrupted**
  - Fix: clear site data for the Cellucid origin and Pull again.
  - Caution: clearing site data removes unsynced local changes; publish anything important first.

### CAP (Cell Annotation Platform) search issues

- The browser sends a bounded object to the configured Cellucid Worker at
  `/cap/lookup-cells`. The Worker chooses a pinned persisted query and relays it
  to CAP; it does not accept a caller-provided query or upstream.
- Search text is visible to the Worker operator, CAP, and the network path. The
  GitHub bearer token and OAuth cookies are not sent to CAP.
- A successful response is
  `{ "contractVersion": 3, "results": [...], "omittedInvalidCount": number }`.
- If your org blocks either network hop, CAP helper searches fail without
  blocking manual annotation.

### The Worker or GitHub cannot be reached

- **“The community-annotation service at … runs an older version (contract N) than this Cellucid (contract 3)”**
  - Cause: the Worker deployment is older than the Cellucid page.
  - Fix: its owner redeploys the Worker from the matching Cellucid source
    (cellucid.com users: report it). No GitHub token was sent.
- **The Worker is reported as unreachable**
  - Causes: the Worker is down, it does not allow this site's origin
    (`ALLOWED_ORIGINS`), the network blocks it, or Cloudflare answered for it
    with an error page without CORS headers (which a browser cannot read).
  - Fix: retry; a self-hosted deployment checks its origin list and its
    Cloudflare logs (see *Limits*).

### Security / privacy review questions

- **Where are GitHub tokens stored?**
  - In browser `sessionStorage` only, with their expiry: a GitHub App user
    token expires after 8 hours (the GitHub App must keep **Expire user
    authorization tokens** on; the Worker refuses a token that never expires).
    Closing the tab clears it; **Disconnect GitHub** revokes it at GitHub
    through the Worker's `/auth/revoke` (and disconnects the dataset's
    repository). The refresh token GitHub issues never leaves the
    Worker.
- **What personal data ends up in the annotation repo?**
  - Each user file contains the contributor's GitHub numeric id (in the file
    name and as `username`/`githubUserId`), their GitHub login (`login`,
    written from the account they sign in with), and the optional profile
    fields they fill in (display name, title, ORCID iD, LinkedIn handle).
    Suggestions, votes and comments name their author by the numeric identity
    `ghid_<id>`. The commits and pull requests that carry the file also show
    the GitHub account that made them.
  - The template validation disallows email fields.

---

## Limits (What Cellucid Cannot Do or Check)

These come from GitHub, Cloudflare or the browser, not from Cellucid. They are
listed so a large round can plan around them.

- **One file per contributor, at most 1,000,000 UTF-8 bytes.** GitHub's
  Contents API returns a file's content as JSON only up to 1 MB
  ([GitHub: get repository content](https://docs.github.com/en/rest/repos/contents#get-repository-content)).
  A contributor's one file holds every dataset they annotate in the repository.
- **Listing of `annotations/users/`.** GitHub documents a limit of 100,000
  entries or 7 MB for a tree listing
  ([GitHub: get a tree](https://docs.github.com/en/rest/git/trees#get-a-tree)).
  When GitHub reports the listing as truncated, Pull stops and says so rather
  than compiling the part it saw.
- **GitHub's hourly budgets.** Every request counts against the signed-in
  user's REST and GraphQL budgets
  ([REST](https://docs.github.com/en/rest/using-the-rest-api/rate-limits-for-the-rest-api),
  [GraphQL](https://docs.github.com/en/graphql/overview/rate-limits-and-node-limits-for-the-graphql-api)).
  Cellucid paces by what GitHub reports and waits when GitHub says so;
  frequent auto-pulls by many people spend each person's own budget.
- **Direct-mode trust.** Anyone with push access can write any file with Git
  or GitHub's web editor, including another contributor's
  `annotations/users/ghid_<id>.json`. The Worker refuses such a write made
  through Cellucid, but Cellucid reads files by content and cannot verify who
  committed them. Review the commit history when in doubt.
- **The Worker cannot bind a write to the connected branch.** Which branch a
  page has connected lives in the browser; the Worker binds every write to the
  signed-in user (their own file, their own `cellucid-annotations/ghid_<id>/`
  branches), not to that branch.
- **Cloudflare Workers Free plan** (cellucid.com's Worker and the checked-in
  configuration): about 10 ms of CPU and 50 subrequests per request
  ([Cloudflare: Workers limits](https://developers.cloudflare.com/workers/platform/limits/)).
  The checked-in Worker stays well inside both (it relays reads without
  parsing them, and a 1,000,000-byte publish costs a few milliseconds). A
  request over the CPU budget is answered by Cloudflare with an HTML `503`
  page without CORS headers, which the browser can only report as an
  unreachable Worker; Cloudflare's logs show it as `exceededCpu` (the Worker
  setup guide's failure diagnosis covers it).

---

## Appendix: Copy/Paste “Author Announcement” Template

If you want a ready-to-send message for annotators, adapt this.

> We’re running a Cellucid community annotation round for **<dataset name>**.
>
> - Start here: `<cellucid link>?annotations=<owner/repo>@<branch>`
> - Please annotate the column: **<fieldKey>** (look for 🗳️ in the field dropdown).
> - Please Publish your work when you’re done so others can see it (PRs must be merged).
> - Use comments to add evidence (markers, references, rationale).
> - If you can’t find the repo in the UI, the most common fix is to install the Cellucid GitHub App and then reload repos in the GitHub sync modal.
> - Deadline: <date>. We will close voting after that and export a consensus snapshot.
