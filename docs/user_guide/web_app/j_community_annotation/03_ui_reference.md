# Community Annotation — UI Reference + Troubleshooting

This page is a **button-by-button** reference for the Community Annotation UI, plus a deep troubleshooting catalog for both annotators and authors.

If you want the guided flow first:

- Annotators: {doc}`01_annotator_guide`
- Authors: {doc}`02_author_guide`
- Every screen in order, with pictures: {doc}`04_lifecycle_walkthrough`

---

## Community Annotation Accordion (Sidebar)

You’ll find **Community Annotation** in the left sidebar.
The accordion starts collapsed.

### Entry button states

Before a repository is connected, the compact panel reads
`No annotation repo connected.` and exposes one action:

- **Connect repo**: opens the four-step connection wizard.

After a repository is connected, the **GitHub sync** action reflects
authentication state:

- **Connect GitHub…**: the repository connection exists locally, but the tab is
  not signed in.
- **GitHub sync…**: the tab is signed in and the repository is connected.

When offline, GitHub actions are disabled. Existing annotation local state can
still be edited, but it cannot be Pulled or Published.

```{figure} ../../../_static/screenshots/community_annotation/disconnected-panel.png
:alt: Community Annotation panel before an annotation repository is connected.
:width: 246px

The Community Annotation entry point reports that no repository is connected and offers the explicit Connect repo action.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-06-connected-status.png
:alt: The connected panel showing Dataset, GitHub and Repo rows, a GitHub sync button, a Profile block reporting the signed-in user, and collapsed sections for Manage annotation, Derived consensus column, and Consensus snapshot and local cache.
:width: 488px

The same panel after a repository is connected. The three header rows —
**Dataset**, **Github**, **Repo** — are the annotation scope; the sections below
them are the working surface.
```

### Status panel fields (what to read when debugging)

Depending on your state, the panel can show:

- **Dataset**: the current dataset id (this must match `annotations/config.json`).
- **GitHub**: whether you’re connected, and your GitHub login.
- **Repo**: the connected repo (optionally including `@branch`).
- **Copy share link** (icon action): copies a link with `?annotations=owner/repo@branch` prefilled.

If you are reporting a problem, a screenshot of this status panel often contains the information needed to debug (dataset id + repo/branch + login).

---

## “GitHub sync” Modal (Repo Discovery + Pull/Publish)

Open it from the Community Annotation accordion.

This modal is a 4-step wizard:

### 1) Sign in with GitHub

- Button: **Continue with GitHub**

Notes:

- This uses a GitHub App OAuth flow (no token paste).
- Tokens are stored only in `sessionStorage` (closing the tab signs you out).

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-02-signin.png
:alt: Step 1 of 4 of the GitHub sync wizard, headed "Sign in with GitHub", with a Continue with GitHub button, and a header showing Dataset, GitHub not connected, and Repo not connected.
:width: 1040px

Step 1 of 4. The four step names stay visible down the left, and the header
tracks the dataset, the GitHub identity, and the bound repository throughout.
```

### 2) Install the GitHub App

Buttons:

- **Add repo** (opens GitHub App installation flow)
- **Reload** (refresh repo list after installing)

If your annotation repo does not appear later, this step is the usual cause.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-03-install-app.png
:alt: Step 2 of 4 listing three repositories as cards marked Public or Private, with Add repo and Reload buttons above and a status line reading "Showing 1-3 of 3 repositories. Page 1 of 1." below.
:width: 1040px

Step 2 of 4. The status line always states the exact count and page, because an
organisation can expose thousands of repositories.
```

### 3) Select an annotation repository

- Filter input: **Filter repositories…**
- Repos are shown as cards (public/private).
- Select a repo card, then use the wizard navigation button **Connect**.
- If replacing an existing repository, the same button reads **Switch repo**.

```{figure} ../../../_static/screenshots/community_annotation/state-no-repository-matches-filter.png
:alt: Step 3 of 4 with a filter applied that matches no repository, showing "No repositories match this filter." above a status line reading "No repositories to display."
:width: 1040px

The empty filter result on step 3. `No repositories match this filter.` is a
filter outcome, not a permission problem — clear the filter to see the full
list.
```

What “connect” means:

- the repo choice is saved locally per dataset and per GitHub account
- Cellucid validates required paths:
  - `annotations/config.json`
  - `annotations/config.schema.json`
  - `annotations/schema.json`
  - `annotations/moderation/merges.schema.json`
  - `annotations/users/`

### 4) Sync (Pull / Publish)

Buttons:

- **Pull latest**
- **Publish**

Optional:

- **Auto pull** toggle
- interval selector (10/15/60 minutes)

Navigation:

- **Back** / **Next** buttons navigate the wizard.
- A status line at the bottom shows progress, warnings, and error text.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-05-sync.png
:alt: Step 4 of 4, headed "Sync (pull / publish)", with Pull latest and Publish buttons, an Auto pull selector offering 10, 15 and 60 minute intervals, and a status line reading "Pulled latest annotations."
:width: 1040px

Step 4 of 4 after a successful **Pull latest**. The status line under the
buttons is where every outcome is reported verbatim.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-17-publish.png
:alt: The same step 4 panel with a status line reading "Publish complete."
:width: 1040px

The same panel after **Publish**. `Publish complete.` means GitHub accepted the
write.
```

### What Pull does (technical)

Pull:

- lists `annotations/users/*.json` and downloads only files whose GitHub `sha` changed since your last Pull
- downloads optional `annotations/moderation/merges.json` when present; only
  authors may publish changes to that document
- writes validated raw-file bodies to the required persistent IndexedDB cache
  and their SHA index to local storage; Cellucid never switches to an in-memory
  cache
- compiles the merged view in the browser
- applies author settings from `annotations/config.json` (annotatable columns, thresholds, closed fields)

Pull is fail-closed and atomic. A missing configured categorical field, an
invalid or oversized active document, more than 10,000 active user files, or
more than 64,000,000 aggregate decoded UTF-8 bytes across those active user
files fails the complete Pull without changing annotation state. Cellucid
checks configured fields before selecting the raw-cache scope or downloading
user/moderation files, so that mismatch cannot mutate the raw-file cache. No
partial merged annotation view is applied.

### What Publish does (technical)

Publish always writes **your** user file:

- `annotations/users/ghid_<id>.json`

If you are an author, Publish may also update:

- `annotations/config.json`
- `annotations/moderation/merges.json`

Publishing modes:

- **Direct push** if GitHub reports you can push.
- **Fork + Pull Request** if you cannot push but the repo allows forking.

Cellucid chooses one route from current permissions before publishing. If that
route fails, the failure is reported; a different route is not substituted.

Common publish pitfalls:

- branch protection blocks direct pushes,
- forking disabled blocks PR flow,
- fork exists but is not accessible to the app token (install the app on your personal account).

---

## Profile (Optional) and “Your identity” Modal

In the Community Annotation accordion, you will see **Profile (optional)**.

Buttons:

- **Edit**: opens the “Your identity” modal
- **Clear**: clears local profile fields (you still need to Publish to update GitHub)

The “Your identity” modal includes:

- **Name:** free text used to search ORCID (e.g. “Alice Smith”)
- **Affiliation / role:** free text (e.g. “Theis Lab, Postdoc”)
- **LinkedIn:** handle only (no URL; lowercase `a-z0-9-`)
- **ORCID:** accepts a checksum-valid ORCID iD and auto-suggests matches
- **Save** / **Cancel**

Notes:

- These fields are optional and are written to your GitHub user file only on Publish.
- Annotation repo validation disallows email fields (privacy).
- For eligible exact input with no surrounding whitespace and at most 240
  Unicode code points, typing at least three Unicode code points in either
  **Name** or **ORCID** starts a 250 ms-debounced direct browser request to
  `https://pub.orcid.org`, unless the exact query is served from the 60-second
  per-modal cache. An exact ORCID iD uses `/v3.0/{id}/person`; other text uses
  `/v3.0/expanded-search/` with at most eight results. The request omits
  credentials and referrer information, but the typed text still leaves
  Cellucid and is visible to ORCID and the network path.
- Save accepts an empty ORCID or an exact checksum-valid ORCID iD. The name
  text is search input, not a valid saved ORCID iD.

---

## MANAGE ANNOTATION (Author Only)

Accordion title: **MANAGE ANNOTATION**

You must be an **author** (GitHub maintain/admin on the repo) to use these controls.

Controls:

- Dropdown: **Categorical obs:** (choose a categorical column)
- Buttons:
  - **Add**: add the selected column to the annotatable list
  - **Remove**: remove the selected column from the annotatable list
  - **Close** / **Reopen**: lock/unlock annotator interaction for this column

Per-column consensus settings (shown only when the column is annotatable):

- **Threshold** slider (maps to `annotations/config.json` threshold)
- **Min annotators** numeric input
- **Apply** (apply locally)
- **Reset** (discard local edits)

Remember: config changes become shared only after the author clicks **Publish** in the GitHub sync modal.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-07-manage-annotation.png
:alt: The Manage annotation section expanded, showing a categorical obs dropdown with Add, Remove and Close buttons, and an Annotatable consensus settings block with a Threshold slider reading 50 percent, a Min annotators field, and Apply and Reset buttons.
:width: 488px

**Manage annotation** with a column already annotatable, so the per-column
consensus settings are shown. **Threshold** is edited in whole percent; it is
stored as a signed fraction whose full range is −100 % to +100 %, because
confidence is itself signed.
```

The defaults, for a column with no explicit settings, are **1** minimum
annotator and a **0.5** threshold.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-08-votable-badge.png
:alt: The Coloring and Filtering panel with the categorical obs dropdown showing an entry prefixed by a ballot-box emoji.
:width: 458px

A column that is open for annotation is prefixed with 🗳️ wherever it is
offered. A column an author has closed is prefixed with 🗳️🏁 instead, and a
column that was never opened has no prefix at all.
```

---

## DERIVED CONSENSUS COLUMN (Optional)

Accordion title: **DERIVED CONSENSUS COLUMN**

This section builds a local categorical obs column for visualization:

- Dropdown: **Annotatable column:**
- Input: **New column key:**
- Settings:
  - **Consensus threshold:**
  - **Min annotators:**
- Button: **Build derived column**

This does not publish anything; it is a local view helper.

---

## CONSENSUS SNAPSHOT + LOCAL CACHE

Accordion title: **CONSENSUS SNAPSHOT + LOCAL CACHE**

### Consensus snapshot

- Label: **Consensus snapshot (cellucid-consensus.json)**
- Button: **Download**

Downloads a locally built `cellucid-consensus.json` snapshot (not written back
to GitHub). Pull first when the export must include the newest shared files.

The snapshot contains:

- `suggestions`: merged suggestions per bucket (with `upvotes` and `downvotes` lists)
- `consensus`: per-bucket `status/label/confidence/voters/netVotes/suggestionId`

### Local cache (two clears, different meanings)

Two different clears exist (they do different things):

- **Clear session**
  - clears local votes/suggestions/comments and annotatable-column selections
    for the current dataset/repo/branch scope
  - does not touch GitHub
  - use this if your local state is confused and you want to start fresh
  - risk: you can lose unpublished work

- **Clear downloads**
  - clears locally cached raw GitHub files (`annotations/users/*`, `annotations/moderation/merges.json`)
  - use this if Pull behavior seems wrong due to a corrupted cache
  - safe: does not remove your local votes/suggestions/comments

The annotation “session” above is separate from Cellucid's portable
`.cellucid-session` files. Community Annotation is intentionally excluded from
ordinary Save State/Load State bundles. GitHub tokens remain tab-scoped in
`sessionStorage`; annotation intent and pulled-file caches retain their own
dataset/repo/branch/user scope. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

---

## Voting Modal (Per Category)

You open this by clicking a category label in the categorical legend when the selected field is annotatable.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-11-consensus.png
:alt: The full application window with the Community voting modal open over the viewer, showing a disputed status banner, two suggestion cards with vote buttons, attribution lines and comment boxes, and a New suggestion form at the bottom.
:width: 996px

The complete voting modal for one category, over the viewer it belongs to.
```

### What you see at a glance

- A **consensus status line** (Consensus / Disputed / Pending)
- A list of suggestion cards
- A **New suggestion** form

The status line is exact about which of the three states applies and why:

| State | Status line | When |
|---|---|---|
| Consensus | `Consensus: "<label>" (<confidence>% • <n> voters)` | at least `minAnnotators` voters, one clear leader, and confidence ≥ threshold |
| Disputed | `Disputed: "<label(s)>" (<confidence>% • <n> voters)` | enough voters, but the leader is tied or confidence is below threshold |
| Pending | `Pending (<n> voters)` | fewer than `minAnnotators` voters; no confidence is reported |

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-12-disputed.png
:alt: A voting modal whose status banner reads Disputed, "Alpha cell, Alpha/PP doublet", 0 percent, 4 voters, above two suggestion cards each at net 0 with two upvotes and two downvotes.
:width: 996px

**Disputed by a tie.** Both suggestions have identical net votes and identical
upvote counts, so neither leads. A tie is disputed regardless of confidence, and
the banner names every tied label alphabetically.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-09-pending.png
:alt: A voting modal whose status banner reads "Pending (2 voters)" above one suggestion card showing net 2.
:width: 996px

**Pending.** Two voters against a three-annotator minimum. No confidence figure
is reported, because the column's own rule says there is not yet enough
participation to compute one.
```

### Suggestion card anatomy

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-10-vote.png
:alt: One suggestion card headed "Alpha cell" with net 0 at the right, Ontology, Markers and Evidence lines, up and down vote buttons at 2 each, Edit, Delete and Merge buttons, an attribution line naming the proposer and their title, and a comment input showing a 0 of 500 counter.
:width: 1028px

One card, with every element named below present: label, `net` score, the three
evidence rows, the vote buttons, the owner-only actions, the attribution line,
and the comment box with its character counter.
```

Each suggestion card shows:

- label (the proposed annotation)
- `net …` score (`up - down`)
- **Ontology:** value or `—`
- **Markers:** summary or `—`
- **Evidence:** preview (with **View full evidence** if long)
- vote buttons:
  - **▲ <count>**
  - **▼ <count>**

Vote toggling:

- click ▲ or ▼ to vote
- click the same arrow again to clear your vote

Duplicates:

- cards with duplicate labels may be highlighted
- votes can split across duplicates until an author merges them

Merged bundles (author moderation):

- a merged suggestion shows a “Merged bundle …” row and a **View merged**
  button; any author can detach a member there, but only the GitHub identity
  that created the merge record can edit or delete its note
- bundle totals are de-duplicated to one effective vote per user
- a direct vote on the merged target wins; otherwise Cellucid delegates the
  majority of that user's votes across the merged members, with a tie producing
  no effective vote

Comments:

- comment box accepts up to 500 characters
- **Enter** submits a comment
- **Shift+Enter** inserts a newline without submitting



---

## Troubleshooting (Massive)

This section is organized by “what you’re trying to do”. Always start by confirming your scope:

- dataset id (status panel)
- `owner/repo@branch`
- GitHub login

### Troubleshooting — Annotators

#### “Dataset mismatch” / “Ask an author to Publish updated settings”

- Cause: the connected repo does not list the current dataset id in `annotations/config.json`.
- Fix: ask an author to connect and Publish (this updates `supportedDatasets[]` and unblocks annotators).

#### “Another browser tab/window is already connected…”

- Cause: Cellucid may enforce a single active tab per scope (dataset + repo@branch + GitHub user id) to prevent silent data loss.
- Fix:
  - close the other Cellucid tab/window for the same scope, then retry
  - if you can’t find it, close all other tabs for this origin and reopen one tab

#### I can’t publish (and I’m not allowed to fork)

- Cause: you have no push permission and the repo has `allow_forking` disabled.
- Fix: ask the author to enable forking or grant write access.

#### Fork + PR publish fails (common hidden cause)

- Cause: your fork is not accessible to the GitHub App token (app not installed on your personal account).
- Fix: install the Cellucid GitHub App on your personal account (ideally “all repositories”), then retry.

#### “Possible conflict” on Publish

- Cause: your remote user file appears newer than your last sync (common if you published from another device).
- Recommended fix: Pull latest, confirm your intended state, then Publish again.
- If you choose “Overwrite remote”: your local intent wins and remote is replaced.

#### Pull succeeds but I don’t see new votes

Checklist:

- confirm you’re on the same `owner/repo@branch` as the rest of the group
- confirm the author did not close the column (🗳️🏁)
- confirm you’re clicking the correct category in the correct column
- confirm the other person actually Published (or their PR merged)

### Troubleshooting — Authors

#### Role shows as “annotator” but I am the dataset author

- Role is based on GitHub permissions on the **annotation repo**:
  - author = maintain/admin
  - annotator = everything else
- Fix: ensure your GitHub account has maintain/admin access, then reconnect and Pull.

#### Repo does not show up under “Select an annotation repository”

- Ensure the GitHub App is installed on the repo owner and includes this repo.
- If installed with “Only select repositories”, add the repo explicitly.
- In the GitHub sync modal: use **Reload**.

#### Pull fails on an invalid or oversized active file

- Cause: one or more active annotation documents are invalid JSON, violate
  their schema, exceed 1,000,000 decoded UTF-8 bytes, or make the Pull exceed
  its active-user-file count or aggregate-user-file byte limit.
- Result: the whole Pull fails; Cellucid does not expose a partial merged view
  or replace the last valid downloaded cache.
- Fix:
  - run the template validation script in the annotation repo (`python scripts/validate_user_files.py`)
  - fix or revert the offending files
  - communicate with the contributor if needed

#### Pull/Publish fails due to storage restrictions

- Some environments block `localStorage` and/or IndexedDB (private browsing, strict privacy settings, embedded contexts).
- Symptoms:
  - cannot acquire cross-tab lock
  - `LOCAL_RAW_CACHE_UNAVAILABLE` or another local persistence error during Pull
  - local persistence/integrity errors that disconnect the repo
- Result: Pull fails without replacing the current annotation view; it does not
  continue with a transient raw-file cache.
- Fix:
  - use a standard browser profile
  - allow site storage
  - avoid running in restricted embedded iframes for annotation work

#### CAP search concerns (privacy / firewall)

- The browser calls the configured Cellucid Worker at `/cap/lookup-cells`; the
  Worker relays a pinned persisted query to CAP.
- The search text is visible to the Worker operator and CAP. The GitHub bearer
  token and OAuth cookies are not sent to CAP.
- If your environment blocks either hop, CAP helper buttons fail without
  blocking manual suggestions or voting.

#### Self-hosted worker issues (CORS / OAuth redirect)

If you self-host the GitHub auth worker:

- CORS failures usually mean your site origin is not in `ALLOWED_ORIGINS`.
- OAuth redirect failures usually mean the GitHub App callback URL does not match `https://<worker-origin>/auth/callback`.
