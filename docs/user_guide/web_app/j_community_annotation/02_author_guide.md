# Community Annotation — Author Guide (Repo Setup + Operations)

This guide is for **dataset authors / maintainers** who want to run a community annotation round and keep it healthy at scale.

If you are an annotator (voting only), go to `01_annotator_guide`.

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

- Section 0.3 (people/process): who merges PRs, who merges duplicates, and how decisions are communicated.
- Section 5 (GitHub settings): fork/PR flow, CI validation, and branch strategy.
- Section 9 (consensus parameters): how thresholds behave under disagreement.
- Section 10 (merges): when to merge vs when to leave disputed.
- Section 13 (rate limits, storage, and “it works for me” failures).
:::

::::

---

## Starting point in the sidebar

```{figure} ../../../_static/screenshots/community_annotation/disconnected-panel.png
:alt: Community Annotation panel before an annotation repository is connected.
:width: 246px

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

- You have a **stable dataset id** (Cellucid uses `dataset_identity.json["id"]`).
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
- You enable validation CI (recommended): it prevents broken JSON from breaking Pull for everyone.

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

Everything else (suggestions, votes, comments) comes from the community (one file per user).

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

- For pre-exported datasets, Cellucid reads `dataset_identity.json["id"]`.
- When exporting via `cellucid.prepare(...)`, you can set `dataset_id=...`.

Example (recommended):

```python
from cellucid import prepare

prepare(
    # ... data args ...
    out_dir="./my_export",
    dataset_id="my_atlas_v1",  # keep this stable for the entire annotation round
    dataset_name="My Atlas (v1)",
    obs_categorical_dtype="uint16",
)
```

How to confirm:

1) Open `my_export/dataset_identity.json`
2) Verify the `id` field is what you expect.

:::{tip}
Treat dataset id like a **contract**:

- If you re-export with small technical changes but the same clusters and meaning, keep the same id.
- If you change clustering, category labels, or biological meaning, use a new id and create a new round (or a new `supportedDatasets[]` entry).
:::

### Confirm the dataset id in the UI (recommended)

In Cellucid, the Community Annotation status panel displays the dataset id (this is the id your `annotations/config.json` must match).



---

## 3) Create the Annotation Repository

Cellucid expects a specific repository layout. The easiest path is to start from the template and customize it.

### Required layout (paths are case-sensitive)

```
annotations/
  config.json
  schema.json
  users/
    (one JSON file per contributor)
  moderation/
    merges.json   (optional but recommended; authors publish it from the UI)
.github/
  workflows/
    validate.yml  (recommended)
scripts/
  validate_user_files.py
```

### Recommended “template-first” setup

The workspace contains a ready-to-copy template repo:

- `cellucid-annotation/`

::::{tab-set}

:::{tab-item} GitHub UI (No Command Line)

1) Create a new GitHub repository (public or private).
2) Copy the contents of `cellucid-annotation/` into the root of the new repo.
3) Update `annotations/config.json` (Section 4).
4) Commit + push.

If you’re doing this in the GitHub web UI:

- create the folders/files exactly as shown,
- copy/paste file contents carefully (JSON is strict: no trailing commas).

:::

:::{tab-item} Git Users (CLI / IDE)

1) Create an empty repo on GitHub.
2) Clone it locally.
3) Copy the template contents into the repo root.
4) Edit `annotations/config.json`.
5) Commit + push.

Example (conceptual):

```bash
cp -R /path/to/cellucid-annotation/* /path/to/your-annotation-repo/
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
  - Must match `dataset_identity.json["id"]` of the dataset currently open in Cellucid.
  - If the currently loaded dataset id is missing:
    - **annotators are blocked** (cannot Pull / view annotations)
    - **authors can still connect** (with a confirmation) and Publish an updated config to unblock everyone

- `supportedDatasets[].name`
  - Human-friendly name shown in the UI.

- `fieldsToAnnotate`
  - List of **categorical obs keys** (column names) that may be annotated.
  - If you list a key that does not exist in the dataset, Cellucid will warn and ignore it.

- `annotatableSettings[fieldKey]`
  - Per-field consensus rules.
  - `minAnnotators` (integer 0–50): minimum unique voters required before a bucket can be anything other than “Pending”.
  - `threshold` (number -1..1): minimum `confidence` to reach “Consensus”.

- `closedFields`
  - Fields in this list are locked for annotators (no voting/suggestions/comments).
  - Validation rule: every closed field must also be in `fieldsToAnnotate`.

:::{important}
Validation rule (enforced by the template CI script):

- Every key in `annotatableSettings` must also appear in `fieldsToAnnotate`.
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
- If the top net-vote score is tied across multiple suggestions → status is always **Disputed**.

See Section 9 for worked examples and recommended defaults.

### Validate your repo inputs (recommended)

In the annotation repo (not in Cellucid), run:

```bash
python scripts/validate_user_files.py
```

This validates:

- `annotations/config.json`
- `annotations/users/*.json`
- `annotations/moderation/merges.json` (optional)

If this fails, fix the file(s) it reports before inviting annotators.

---

## 5) Configure GitHub Repo Settings (Highly Recommended)

### Decide how annotators will publish

You have two viable models:

1) **Direct publish** (annotators have write access)
   - Pros: simplest experience (Publish writes directly to `annotations/users/ghid_<id>.json`)
   - Cons: requires adding many people as collaborators; less review control

2) **Fork + Pull Request publish** (annotators do not have write access)
   - Pros: reviewable contributions, no direct writes to your repo
   - Cons: you must merge PRs; depends on forking being allowed and not blocked by org policy

Cellucid chooses the best option per user:

- If the user can push → direct publish
- Else if the repo allows forking → fork + PR publish
- Else → user cannot publish (they can still vote locally, but nothing can be shared)

### Fork + PR model: one extra requirement most teams miss

For PR-based publishing, contributors need their fork to be accessible to the GitHub App token.

Practical recommendation you can tell annotators:

- Install the Cellucid GitHub App on your **personal GitHub account** with access to **all repositories** (so newly created forks are included automatically).

If they do not do this, the PR flow may fail in confusing ways (the UI can’t see the fork).

### Enable validation CI

The template includes:

- `scripts/validate_user_files.py`
- `.github/workflows/validate.yml`

Suggested GitHub settings:

- Require the validation check to pass before merging PRs.
- For fork-based contributions, allow GitHub Actions to run on PRs (org policies may apply).

### Branch protection (advanced, but important)

Branch protection can break direct publishing:

- If direct pushes are blocked, users with “write” permissions may still see Publish fail.

Options:

- relax protection for the annotation repo (common for JSON-only repos), or
- require PR flow for all changes (requires policy + reviewer time), or
- publish to a dedicated “staging” branch and merge via PR (more overhead, more control).

---

## 6) Install the GitHub App (Required for Repo Discovery)

Cellucid’s UI lists only repositories where the **Cellucid GitHub App is installed**.

1) Install the app on the user/org that owns the annotation repo.
2) If you choose “Only select repositories”, make sure the annotation repo is selected.

:::{note}
Org repos often require an org admin to approve the installation.
:::

### Optional: self-host the GitHub OAuth + API proxy (org deployments)

Cellucid’s community annotation UI uses a small server component (typically a Cloudflare Worker) to:

- run the GitHub OAuth flow for a GitHub App
- proxy GitHub API requests so the frontend never needs GitHub secrets

If you are using **cellucid.com**, you typically do *not* need to do anything here.

If your organization requires owning the auth infrastructure (recommended for many orgs), you can self-host.

#### What you need (conceptual)

- A **GitHub App** (not a Personal Access Token) with OAuth enabled.
- A Worker deployment (Cloudflare Worker or equivalent) that exposes:
  - `/auth/login`, `/auth/callback`
  - `/auth/user`, `/auth/installations`, `/auth/installation-repos`
  - `/api/*` proxy to `https://api.github.com/*`

#### Worker configuration (conceptual)

The Worker must be configured with:

- `ALLOWED_ORIGINS`: comma-separated allowlist (CORS), e.g. `https://your.cellucid.site,https://staging.your.cellucid.site`
- `GITHUB_APP_ID`: numeric GitHub App id
- `GITHUB_PRIVATE_KEY`: GitHub App private key (PEM)
- `GITHUB_CLIENT_ID`: GitHub App OAuth client id
- `GITHUB_CLIENT_SECRET`: GitHub App OAuth client secret

#### Pointing Cellucid at your worker (deployment-specific)

- Some builds support a runtime override (e.g. `window.__CELLUCID_GITHUB_WORKER_ORIGIN__`).
- Production deployments usually set the worker origin at build time.

See also:

- `cellucid/docs/github-oauth-cloudflare-setup.md` in the Cellucid web repo

---

## 7) Connect the Repo From Cellucid (Author Bootstrap)

1) Load your dataset in Cellucid.
2) Open the **Community Annotation** accordion.
3) Click **Connect repo**. If this dataset already has a repository connection,
   use **Connect GitHub…** or **GitHub sync…** as shown.
4) Complete **Continue with GitHub**, enable/reload repositories if needed,
   select the repo, click **Connect**, then **Pull latest**.
5) Confirm you see author-only controls (e.g. **MANAGE ANNOTATION**).

### Dataset mismatch (the most common “why can’t annotators Pull?” issue)

If the dataset loaded in Cellucid is not present in `annotations/config.json`:

- Annotators are blocked (they cannot Pull).
- Authors can still connect (with a warning) and then **Publish** to write an updated config that adds/updates `supportedDatasets[]` for the current dataset id.

This “author override” exists to make first-time bootstrapping smooth.

### “Am I actually an author?” (role sanity check)

Role is derived from GitHub permissions on the annotation repo:

- author = `maintain` or `admin`
- annotator = everything else

If you are an author, you should see author-only UI blocks such as **MANAGE ANNOTATION**.

If you do not:

1) Confirm you are signed into the expected GitHub account in the GitHub sync modal.
2) Confirm your permission level on the repo is maintain/admin.
3) Disconnect/reconnect and Pull again.

---

## 8) Enable/Disable Annotatable Columns (Author UI)

Once connected, open **MANAGE ANNOTATION** inside the Community Annotation accordion:

1) Select the categorical obs field you want to control (dropdown labeled **Categorical obs:**).
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

Inside **MANAGE ANNOTATION**, after you select a column that is already annotatable, you will see **Annotatable consensus settings**:

- **Threshold** slider (maps to `threshold` in `annotations/config.json`)
- **Min annotators** input (maps to `minAnnotators`)
- **Apply** (apply locally)
- **Reset** (discard local edits)

After applying locally, **Publish** to write the settings to `annotations/config.json` so all annotators receive them on Pull.

### Default settings (if you do not set anything)

If a field has no explicit entry in `annotatableSettings`, Cellucid uses defaults:

- `minAnnotators = 1`
- `threshold = 0.5`

In practice, for predictable behavior across devices/users, you should set `annotatableSettings` explicitly for every annotatable field.

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
| 4 users split: 2 upvote A, 2 upvote B | 4 | 2 | 0.5 | **Disputed** due to tie (two top suggestions) |

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

As an author, you can merge duplicates so votes are combined:

1) Open a category’s voting modal.
2) Drag one suggestion card onto another.
3) Confirm the merge and optionally write a merge note.
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

`annotations/moderation/merges.json` is an author-only log of merges.

Each merge entry includes:

- `bucket`: which category bucket the merge applies to (format `<fieldKey>:<categoryLabel>`)
- `fromSuggestionId` → `intoSuggestionId`: the mapping (merge “from” into “into”)
- `by`: author identity (stored as `ghid_<githubUserId>`)
- `at`: timestamp
- `editedAt` (optional): if you later edit the merge note
- `note` (optional): a short rationale shown in the UI

Example (illustrative):

```json
{
  "bucket": "leiden:7",
  "fromSuggestionId": "sug_01HXYZ...",
  "intoSuggestionId": "sug_01HABC...",
  "by": "ghid_12345",
  "at": "2025-01-01T12:34:56Z",
  "note": "Merged formatting variants: 'T cell' and 'T-cell'."
}
```

### Undoing a merge (if you merged the wrong thing)

Depending on the UI version, you may have an author-only option to detach/undo a merge from the **View merged** dialog.

If the UI does not expose an undo:

- revert the commit that changed `annotations/moderation/merges.json`, or
- edit `merges.json` manually (advanced; validate after editing).

---

## 11) Build a Derived Consensus Column (Optional, Local)

Cellucid can build a local derived categorical obs column for visualization:

1) Open **DERIVED CONSENSUS COLUMN**.
2) Choose an annotatable source column (e.g. `leiden`).
3) Choose a **New column key** (e.g. `community_cell_type`).
4) Set threshold / min annotators (used only for this derived column).
5) Build the derived column.

Result:

- Each category becomes a label (if consensus), or `Disputed`, or `Pending`.
- This does **not** change voting rules and does **not** publish anything to GitHub.

This is useful to:

- color the atlas by the current community consensus
- quickly spot which clusters still need attention

:::{note}
The derived column uses the threshold/minAnnotators you set in the derived-column UI, which may differ from the annotatable field’s official settings. If you export screenshots or figures, record which settings you used.
:::

---

## 12) Export a Consensus Snapshot (For Downstream Tools)

Inside **CONSENSUS SNAPSHOT + LOCAL CACHE**:

- **Consensus snapshot (cellucid-consensus.json)** → **Download**

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

Buckets are keyed by `<fieldKey>:<categoryLabel>`.

### Example downstream usage (computational)

This is one simple pattern: map consensus labels back onto an `AnnData` cluster column.

```python
import json

import pandas as pd

doc = json.load(open("cellucid-consensus.json"))

target_field = "leiden"

mapping = {}
for bucket, summary in doc["consensus"].items():
    field_key, category_label = bucket.split(":", 1)
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

If you don’t find your issue here, also check `03_ui_reference` (it includes additional UI-specific guidance and error messages).

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

- **Repo does not appear under “Select an annotation repository”**
  - Causes:
    - the Cellucid GitHub App is not installed for the repo owner
    - the app was installed for “Only selected repositories” and the repo is not selected
    - you are signed into a different GitHub account than expected
  - Fix:
    - install/adjust the app installation and try **Reload** in the GitHub sync modal
    - verify your GitHub username in the modal matches your intended account

- **Fork + PR flow fails for annotators**
  - Common cause: annotator did not install the GitHub App on their personal account (their fork isn’t visible to the token).
  - Fix: ask them to install the app for their personal account (ideally “all repositories”), then retry Publish.

### Publish failures (authors)

- **Publishing fails with “Sign in required.”**
  - Fix: sign in again (tokens are session-only; closing the tab clears them).

- **Publishing fails but you have write access**
  - Common causes:
    - branch protection blocks direct writes
    - required status checks are configured but GitHub API rejects direct commit
  - Fix options:
    - relax branch protection for the annotation repo, or
    - require PR flow for all changes, or
    - publish to a dedicated branch and merge via PR

- **Publishing fails for annotators and you disabled forking**
  - Cause: users without push cannot publish if `allow_forking` is disabled.
  - Fix: enable forking, or grant write access to annotators.

### Pull is slow / rate-limited

- Causes:
  - very large number of user files
  - frequent auto-pulls across many users
  - GitHub rate limits for your org
- Mitigations:
  - keep the annotation repo “clean” (only JSON + scripts; avoid large binaries)
  - avoid massive numbers of branches with many files (each branch multiplies history/tree size)
  - ask annotators to Pull on demand instead of using aggressive auto-pull intervals

### Local cache corruption / storage restrictions

- **Warning about IndexedDB unavailable (“downloads will be cached in-memory…”)**
  - Cause: browser storage policies (private mode, strict settings, embedded iframe restrictions).
  - Impact: Pull will re-download after reload; large repos will feel slower.
  - Fix: use a normal browser profile, allow site storage, avoid restrictive privacy modes for the annotation session.

- **Error about local cache being corrupted**
  - Fix: clear site data for the Cellucid origin and Pull again.
  - Caution: clearing site data removes unsynced local changes; publish anything important first.

### CAP (Cell Annotation Platform) search issues

- CAP search queries are sent to `https://celltype.info/graphql`.
- If your org blocks outbound calls, CAP helper searches will fail. This does not block manual annotation.

### Security / privacy review questions

- **Where are GitHub tokens stored?**
  - In browser `sessionStorage` only (cleared when the tab closes).
- **What personal data ends up in the annotation repo?**
  - User files contain GitHub numeric id and optional profile fields (display name/title/orcid/linkedin handle).
  - The template validation disallows email fields.

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
