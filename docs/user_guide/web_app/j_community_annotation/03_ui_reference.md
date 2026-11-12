# Community Annotation — UI Reference + Troubleshooting

This page is a **button-by-button** reference for the Community Annotation UI, plus a deep troubleshooting catalog for both annotators and authors.

If you want the guided flow first:

- Annotators: {doc}`01_annotator_guide`
- Authors: {doc}`02_author_guide`
- Every screen in order, with pictures: {doc}`04_lifecycle_walkthrough`

---

## Community annotation accordion (sidebar)

You’ll find **Community annotation** in the left sidebar.
The accordion starts collapsed.

### Entry button states

Before a repository is connected, the compact panel reads **Repository**,
`Not connected.`, and exposes one action:

- **Connect repo**: opens the four-step connection wizard.

After a repository is connected (which needs this tab signed in to GitHub),
the panel shows the dataset by its name (its id in the repository on hover),
your GitHub user and the repository, and one **GitHub sync…** button that
opens the dialog with **Pull latest** and **Publish**. A tab that is not
signed in shows `Not connected.` and **Connect repo**; signing in again as the
same GitHub user brings the connection back.

When offline, GitHub actions are disabled. Existing annotation local state can
still be edited, but it cannot be Pulled or Published.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-01-disconnected.png
:alt: The Community annotation panel before any repository is connected: Repository with its (i), the line Not connected. and a Connect repo button, with the pointer on that button.
:width: 504px

The Community Annotation entry point reports that no repository is connected and offers the explicit Connect repo action.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-06-connected-status.png
:alt: The connected Community annotation panel: GitHub sync with the Dataset, GitHub and Repo rows and a GitHub sync… button, a Profile block naming the signed-in user with Edit and Clear, and the collapsed Annotated columns, Consensus column and Snapshot and cache groups.
:width: 504px

The same panel after a repository is connected. The three header rows —
**Dataset**, **GitHub**, **Repo** — are the annotation scope; the sections below
them are the working surface.
```

### Status panel fields (what to read when debugging)

Depending on your state, the panel can show:

- **Dataset**: the current dataset by its name; hover it for the dataset id
  the repository files it under (this id must match `annotations/config.json`).
- **GitHub**: whether you’re connected, and your GitHub login.
- **Repo**: the connected repo (optionally including `@branch`).
- **Copy share link** (icon action): copies a link with `?annotations=owner/repo@branch` prefilled.

Below the rows, one line at a time says what needs you:

- **`3 changes not published`** with **Publish**: your votes, suggestions,
  comments, deletions or profile edits (and, for an author, repository
  settings or merges) that this browser has and GitHub does not. Hover the
  line for the kinds ("2 votes, 1 suggestion"). It counts against your file as
  this browser last pulled or published it; after clearing the local cache it
  is not shown until the next Pull.
- **`Open in another tab — annotation paused here`** with **Use here**: the
  same dataset, repository and GitHub user are being annotated in another tab
  of this browser. One tab annotates at a time, so nothing is overwritten;
  everything else of the dataset works in both tabs. **Use here** moves
  annotation to this tab, and the other tab shows the same line (what it had
  is kept).
- a failure in one plain sentence (GitHub's own words behind **Details**).

If you are reporting a problem, a screenshot of this status panel often contains the information needed to debug (dataset + repo/branch + login; hover the dataset for its id).

---

## “GitHub sync” Modal (Repo Discovery + Pull/Publish)

Open it from the **Community annotation** accordion.

This modal is a 4-step wizard:

### 1) Sign in

- Button: **Continue with GitHub**

Notes:

- This uses a GitHub App OAuth flow (no token paste).
- The token expires after 8 hours and is stored only in this tab's
  `sessionStorage` (closing the tab signs you out; a reload does not).
  **Disconnect GitHub** (the icon beside GitHub in this dialog) revokes it at
  GitHub and disconnects this dataset's repository too; if GitHub cannot be
  reached, a notification says the token stays valid until it expires.
- A token GitHub refuses (it expired or was revoked elsewhere) signs you out
  here with one line: "Your GitHub sign-in expired or was revoked. Sign in
  again to Pull or Publish; your unpublished changes stay in this browser."
  The repository connection is kept: signing in again as the same GitHub user
  connects it again, with your unpublished changes.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-02-signin.png
:alt: Step 1 of 4 of the GitHub sync wizard: a header showing the dataset with GitHub and Repo both not connected, the four steps Sign in, Install app, Choose repository and Sync with Sign in current, the line "Cellucid gets no repository access until you install its GitHub App on repositories you choose.", and a Continue with GitHub button.
:width: 1040px

Step 1 of 4. The four step names stay visible down the left, and the header
tracks the dataset, the GitHub identity, and the bound repository throughout.
```

### 2) Install app

Buttons:

- **Add repo** (opens GitHub App installation flow)
- **Reload** (refresh repo list after installing)

If your annotation repo does not appear later, this step is the usual cause.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-03-install-app.png
:alt: Step 2 of 4 of the GitHub sync wizard, Install app: the header shows the dataset, the signed-in GitHub account and Repo not connected; the line "Add repo opens GitHub to choose repositories; then Reload.", Add repo (under the pointer) and Reload, three repository cards marked Public or Private, and "Showing 1–3 of 3 repositories. Page 1 of 1."
:width: 1040px

Step 2 of 4. The status line always states the exact count and page, because an
organization can expose thousands of repositories.
```

### 3) Choose repository

- Filter input: **Filter repositories…**
- Repos are shown as cards (public/private).
- Select a repo card, then use the wizard navigation button **Connect**.
- If replacing an existing repository, the same button reads **Switch repo**.

```{figure} ../../../_static/screenshots/community_annotation/state-no-repository-matches-filter.png
:alt: Step 3 of 4 of the GitHub sync wizard, Choose repository, with a filter typed that matches no repository and the line No repositories match this filter. under it.
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
- the three schema files must be exactly the
  [cellucid-annotation](https://github.com/theislab/cellucid-annotation)
  template's, which Cellucid checks by their Git blob ids without downloading
  them; a schema that differs stops the Pull with a message naming the file to
  copy again

### 4) Sync

Buttons:

- **Pull latest**
- **Publish**

Optional:

- **Auto pull** toggle
- interval selector (10/15/60 minutes)

Navigation:

- **Back** / **Next** buttons navigate the wizard.
- A status line at the bottom shows progress, warnings, and error text.

While any sync is running — including a background **Auto pull** — **Pull
latest** and **Publish** are disabled until it finishes. A Pull or Publish that
still reaches a running sync does not run and reports `Another GitHub sync is
still running. Wait for it to finish, then try again.`; one you cancel reports
`Canceled.`. A Pull or Publish stopped because the dataset, the repository or
the GitHub user changed, the repository was disconnected, or annotation moved
to another tab is not reported as failed: its notification is withdrawn and
what stopped it says why. A Publish stopped after GitHub had already accepted
part of it says so instead, naming the files GitHub accepted, and asks you to
Pull before publishing again. One stopped because the connection was lost is
reported as failed. Only a completed transaction reports its success text.

When GitHub cannot be reached while Cellucid checks your role on the
repository (each request already asked again while that could help), the
status line says so and the connection is kept; the role is checked again with
the next Pull (yours or an automatic one) and as soon as this computer comes
back online.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-05-sync.png
:alt: Step 4 of 4 of the GitHub sync wizard, Sync: the dataset, GitHub account and repository in the header, the line "Pull fetches the latest votes; Publish writes them for everyone.", Pull latest and Publish, an Auto pull checkbox with its interval, and the status line "Pulled latest annotations."
:width: 1040px

Step 4 of 4 after a successful **Pull latest**. The status line under the
buttons is where every outcome is reported verbatim.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-17-publish.png
:alt: The same Sync step after Publish, its status line reading "Published." under the Pull latest and Publish buttons.
:width: 1040px

The same panel after **Publish**. `Published.` means GitHub accepted the write
to the repository. When you cannot write to the repository, Publish opens (or
updates) a pull request instead and says `Opened a pull request: <link>. The
repository has your changes once its authors merge it.`; nothing new to write
says `Nothing new to publish: GitHub already has your changes.`
```

### What Pull does (technical)

Pull:

- reads the repository at one commit: the branch's commit and the listings of
  `annotations/` and `annotations/moderation/` come in one request, and
  everything after is read by Git object id from that commit, so one Pull
  never mixes two commits
- stops there when nothing under `annotations/` changed since your last
  completed Pull (one or two requests in all)
- checks the three schemas by blob id, reads `annotations/config.json` and
  lists `annotations/users/` by its tree id
- downloads only the user files (and the optional
  `annotations/moderation/merges.json`; only authors may publish changes to
  that document) whose blob id changed, many files per request
- validates each file as it arrives and keeps it in the required persistent
  IndexedDB cache, batch by batch, so an interrupted Pull resumes where it
  stopped; Cellucid never switches to an in-memory cache
- compiles the merged view in the browser
- applies author settings from `annotations/config.json` (annotatable columns, thresholds, closed fields)

One file never stops a Pull. A user file or merges file that is not valid is
left out of the consensus on its own and reported in one warning (how many
files, with each file and its reason under its details); a path under
`annotations/users/` that is not a `ghid_<id>.json` file, or a case variant of
a layout path, is ignored and listed the same way. A missing configured
categorical field, an invalid `annotations/config.json`, a schema that is not
the template's, or a users listing GitHub reports as truncated fails the
complete Pull without changing annotation state. There is no limit on how many
user files a repository holds or on their total size: changed files are cached
a batch at a time, and the merged view is compiled by streaming the cache one
file at a time. Cellucid checks configured fields before selecting the
raw-cache scope or downloading user/moderation files, so that mismatch cannot
mutate the raw-file cache.

Requests are paced by the request budget GitHub reports; when GitHub names a
wait, the progress notification says so. A request is never cut by a fixed
time limit: only a connection on which nothing arrives at all is given up and
asked again.

### What Publish does (technical)

Publish always writes **your** user file, with every change of yours not
published yet:

- `annotations/users/ghid_<id>.json`

If you are an author, Publish may also update:

- `annotations/config.json`
- `annotations/moderation/merges.json`

Publishing routes (named here; the panel does not show them, its result line
says which one ran: *Published.*, or *Opened a pull request: …*):

- **Direct publish** if GitHub reports you can push.
- **Fork + Pull Request** if you cannot push but the repo allows forking.
  You make the fork on GitHub (Cellucid never creates one) and install the
  Cellucid GitHub App on your account with access to it (*All repositories*,
  or *Only select repositories* with the fork selected).

Cellucid chooses one route from current permissions before publishing. If that
route fails, the failure is reported; a different route is not substituted.

With **Fork + Pull Request**, you have one pull request per base branch, from
the head branch `cellucid-annotations/ghid_<id>/…` in your fork. Every Publish
resets that branch to the base branch's current tip and writes your whole user
file on it (every dataset you annotate there, the open pull request's pending
sections included), so the pull request's diff is your one file and an earlier
squash-merged pull request never makes it conflict. An open pull request is
reused. An agent's vote or suggestion, published alone, is added to the file
the open pull request already proposes, so it never drops that pull request's
other changes.

If the connection drops while a write is under way, GitHub may or may not
have it: Cellucid says `GitHub may already have applied this change;
Cellucid pulls now to see.` and pulls at once, so what it knows of your file is
GitHub's again before the next Publish. **Stop** on an agent's vote while it
is written works the same way (see Agents below).

Common publish pitfalls:

- branch protection blocks direct pushes,
- forking disabled blocks PR flow,
- no fork yet: Publish says `You have no fork of <owner>/<repo> yet: fork it on
  GitHub, then Publish again.` under the sync line, with a **Fork on
  GitHub** link,
- a fork GitHub is still copying: `GitHub is still copying your fork …;
  Publish again in a minute.`,
- a fork the app cannot write to: `Cellucid cannot write to your fork …`, with
  an **Install the app** link (install it on your account with access to the
  fork).

---

## Profile (Optional) and “Your identity” Modal

In the **Community annotation** accordion, you will see the optional **Profile** block.

Buttons:

- **Edit**: opens the “Your identity” modal
- **Clear**: clears local profile fields (you still need to Publish to update GitHub)

The “Your identity” modal includes:

- **Name**: free text used to search ORCID (e.g. “Alice Smith”)
- **Affiliation / role**: free text (e.g. “Theis Lab, Postdoc”)
- **LinkedIn**: handle only (no URL; lowercase `a-z0-9-`)
- **ORCID**: accepts a checksum-valid ORCID iD and auto-suggests matches
- **Save** / **Cancel**

Notes:

- These fields are optional and are written to your GitHub user file only on Publish.
- Annotation repo validation disallows email fields (privacy).
- For eligible exact input with no surrounding whitespace (of any length),
  typing at least three Unicode code points in either
  **Name** or **ORCID** starts a 250 ms-debounced direct browser request to
  `https://pub.orcid.org`, unless the exact query is served from the 60-second
  per-modal cache. An exact ORCID iD uses `/v3.0/{id}/person`; other text uses
  `/v3.0/expanded-search/` with at most eight results. The request omits
  credentials and referrer information, but the typed text still leaves
  Cellucid and is visible to ORCID and the network path.
- Save accepts an empty ORCID or an exact checksum-valid ORCID iD. The name
  text is search input, not a valid saved ORCID iD.

---

## Annotated columns (author only)

Collapsible title: **Annotated columns**

You must be an **author** (GitHub maintain/admin on the repo) to use these controls.

Controls:

- Dropdown: **Categorical obs** (choose a categorical column)
- Buttons:
  - **Add**: add the selected column to the annotatable list
  - **Remove**: remove the selected column from the annotatable list
  - **Close** / **Reopen**: lock/unlock annotator interaction for this column

Per-column **Consensus settings** (shown only when the column is annotatable):

- **Threshold** slider (maps to `annotations/config.json` threshold)
- **Min annotators** numeric input
- **Apply** (apply locally)
- **Reset** (discard local edits)

Remember: config changes become shared only after the author clicks **Publish** in the GitHub sync modal.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-07-manage-annotation.png
:alt: Annotated columns open in the Community annotation panel: the Categorical obs picker on the voting column clusters with Add, Remove and Close, and its consensus settings with a Threshold slider and a Min annotators field.
:width: 504px

**Annotated columns** with a column already annotatable, so the per-column
consensus settings are shown. **Threshold** is edited in whole percent; it is
stored as a signed fraction whose full range is −100% to +100%, because
confidence is itself signed.
```

The defaults, for a column with no explicit settings, are **1** minimum
annotator and a **0.5** threshold.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-08-votable-badge.png
:alt: The Categorical obs select holding clusters, prefixed by a ballot-box badge marking it open for voting, and the legend beneath it: each category with its consensus label (Beta cell under Beta, a dash under the others) and its up and down vote counts.
:width: 496px

A column that is open for annotation is prefixed with 🗳️ wherever it is
offered. A column an author has closed is prefixed with 🗳️🏁 instead, and a
column that was never opened has no prefix at all.
```

---

## Consensus column (optional)

Collapsible title: **Consensus column**

This section builds a local categorical obs column for visualization:

- Dropdown: **Annotatable column**
- Input: **New column key**
- Settings:
  - **Consensus threshold**
  - **Min annotators**
- Button: **Build derived column**
- For each built column whose consensus changed since (a Pull, a vote): one
  line `<column>: consensus changed since built (N categories)` with
  **Rebuild**

The key starts as `community_<column>` and the settings at the column's rule
in the repository. This does not publish anything; it is a local view helper
and a step of the history that records its outcomes, rule and source, and
replays by name ({doc}`02_author_guide`, section 11).

---

## Consensus snapshot + local cache

Collapsible title: **Snapshot and cache**

### Consensus snapshot

- Button: **Download** (its tooltip reads `Download cellucid-consensus.json`)

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
  - clears what Pull keeps from GitHub for this scope: the open dataset's part
    of each file read from `annotations/users/*` and
    `annotations/moderation/merges.json`, the reasons files were left out, and
    the record of the last completed Pull; the next Pull reads everything again
  - use this if Pull behavior seems wrong due to a corrupted cache
  - safe: does not remove your local votes/suggestions/comments

The annotation “session” above is separate from Cellucid's portable
`.cellucid-session` files. Community Annotation is intentionally excluded from
the session files **Share…** and **Load…** handle. GitHub tokens remain tab-scoped in
`sessionStorage` and expire after 8 hours; annotation intent and pulled-file caches retain their own
dataset/repo/branch/user scope. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

---

## Voting Modal (Per Category)

You open this by clicking a category label in the categorical legend when the selected field is annotatable.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-11-consensus.png
:alt: The voting view for the Beta cluster with a banner reading "Consensus: "Beta cell" (100% • 4 voters)", the suggestion card with its ontology id, marker, evidence and four upvotes, and the New suggestion form below.
:width: 988px

The voting modal for one category: the status line, the suggestion cards, and
the **New suggestion** form beneath them.
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
| Pending | `Pending (<n> voters)` | no voter yet, or fewer than `minAnnotators` voters; no confidence is reported |

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-12-disputed.png
:alt: The voting view for the Beta cluster with a banner reading "Disputed: "Beta (immature), Beta cell" (0% • 4 voters)", two suggestion cards each at net 0 with two upvotes and two downvotes, and the New suggestion form below.
:width: 988px

**Disputed by a tie.** Both suggestions have identical net votes and identical
upvote counts, so neither leads. A tie is disputed regardless of confidence, and
the banner names every tied label alphabetically.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-09-pending.png
:alt: The voting view for the Beta cluster with a banner reading "Pending (1 voter)", one suggestion card (Beta cell, its ontology id, marker and evidence, one upvote), and the New suggestion form below.
:width: 988px

**Pending.** One voter against a three-annotator minimum. No confidence figure
is reported, because the column's own rule says there is not yet enough
participation to compute one.
```

### Suggestion card anatomy

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-10-vote.png
:alt: One suggestion card showing its label, net score, ontology, marker and evidence lines, its up and down vote buttons, the proposer attribution and the comment box, with the pointer on the upvote button.
:width: 1020px

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

- suggestions with exactly the same label (compared after Unicode NFC
  normalization; case, spacing and punctuation are not folded) are one
  suggestion: the earliest proposal is shown, and each person's votes across
  them count once (their majority; a tie is no vote)
- adding a suggestion, or renaming one, to a label already shown in the
  category is refused
- labels that differ (`T cell` and `T-cell`) can split votes until an author
  merges them

Merged bundles (author moderation):

- a merged suggestion shows a “Merged bundle …” row and a **View merged**
  button; any author can detach a member there, but only the GitHub identity
  that created the merge record can edit or delete its note
- bundle totals are de-duplicated to one effective vote per user
- a direct vote on the merged target wins; otherwise Cellucid delegates the
  majority of that user's votes across the merged members, with a tie producing
  no effective vote
- an author merges a suggestion into another by dragging its card onto the
  other, or from the keyboard with the card's **Merge…** button: choose the
  target suggestion, add an optional note, and confirm

Comments:

- the annotation contract sets no length limit on comments, labels, evidence
  or any other text; only the contributor's whole file is bounded
  (1,000,000 UTF-8 bytes)
- **Enter** submits a comment
- **Shift+Enter** inserts a newline without submitting
- unsent comments and suggestions, and open edits of your own comments,
  suggestions, and merge notes, keep their text, cursor, and focus when the
  modal refreshes (a Pull, a change from another tab, your own vote)
- if the item you were editing is deleted, merged away, or no longer on the
  page you are viewing, or voting closes, the editor closes and warns
  `Edit not saved: the … you were editing is no longer shown here (deleted,
  merged, closed, or on another page).`



---

## Renamed columns and categories

A column or category you renamed in this browser (the field toolbar, the
legend, a history step) keeps its dataset name in the repository: votes,
suggestions and the consensus are filed under the dataset's names, so other
people see them, and a rename (or a history jump, a session or the autosave
that brings one back) never moves anything in the repository. The legend and
the voting dialog show your name; the dialog's header adds the repository's
(`Alpha renamed (in the repository: cell_type • alpha)`). While a column is
open for voting, it cannot be renamed or edited, whatever name it shows.

---

(community-annotation-agents)=
## Agents and community annotation

An agent (the chat bar, Claude or ChatGPT through the connector, a local app
through `cellucid mcp`) and a Python or R notebook can read the community
annotation of the dataset shown, vote and suggest labels, through the same
session as the voting dialog. Connecting a repository and signing in to GitHub stay yours.

**Reading** (`annotation_suggestions`, no permission needed): for each
category of a column, its suggestions (label, ontology term, who suggested it
as the dialog names them, votes, your own vote, and whether it is delegated
from votes on merged suggestions) and its consensus (status, leading label,
voters) under the repository's settings, plus when this browser last pulled.
Columns and categories are named as you see them. A large column comes in
pages, each whole (`next` starts the following one). What the agent reads
goes to its provider, other people's GitHub logins and names included
({ref}`privacy-agents`).

**Voting and suggesting** (`annotation_publish`) publishes the votes and new
labels of one call together, as one commit of your file, and nothing else.
The agent is asked to put everything of a turn in one call, so a turn of
eight votes is one commit in the repository's history and one question to
you:

- you are asked first (Publishing is **Ask** by default), whoever asks: the
  question says how many votes and suggestions it publishes ("It publishes 3
  votes and 2 suggestions (each with your upvote) to team/labels@main as one
  commit"), the repository it goes to (or that it goes in a pull request when
  you cannot write to the repository), and how many of your other changes
  stay unpublished; each vote and suggestion is under **Details**;
- only your file as GitHub holds it plus these changes is written: your
  drafts, your other votes not yet published and, for an author, the
  repository settings and merges are not published (they wait for your own
  **Publish**);
- each new label comes with your upvote, as **New suggestion** adds one; a
  label already suggested for the category is refused with its id (vote on
  it instead);
- every change is checked before you are asked: if one is refused (a label
  already there, a suggestion that does not exist, the same suggestion voted
  twice differently, a field not open to voting), nothing of the call is
  published, nothing is asked, and the refusal names each one;
- votes that are already so (published and here) ask nothing and write
  nothing;
- an agent's call never opens a dialog (your profile's fields are optional;
  the file is named by your GitHub login), waits for a Pull under way, and
  pulls first when your file changed on GitHub since the last Pull;
- a suggestion is named by its id or its label; suggestions with the same
  label are one suggestion (the earliest proposal), so a label names one;
- the agent is told where the changes went: the repository, a pull request
  (with its link), or nowhere (the votes were already so), with each vote's
  counts now and each new suggestion's id.

**Undo turn cannot take back a publication.** The agent's line counts what it
published ("1 annotation commit published, not undoable"), the turn keeps it
in its record in the history, and Undo turn says it stays published. Change a
vote back by voting again.

**Stop while a publication is written**: if GitHub may already have it,
Cellucid pulls at once and the agent is told what the Pull shows: published,
or kept in this browser and not published. Your next Publish works either
way.

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

#### Fork + PR publish fails

- "You have no fork of …": fork the repository on GitHub (the **Fork on GitHub** link), then Publish again.
- "Cellucid cannot write to your fork …": the app is not installed on your account with access to the fork. Install it (**Install the app**; *All repositories*, or *Only select repositories* with the fork selected), then Publish again.
- "GitHub is still copying your fork …": Publish again in a minute.

#### Your file changed on GitHub since your last Pull

- Cause: your user file on GitHub is newer than the one your last Pull read
  (common if you published from another device or tab).
- What Publish does: it pulls (your file on GitHub and what is here are
  merged), then publishes once more, by itself; an agent's publication does
  the same. Cellucid never replaces a remote file it has not read.
- Only when you are an author with repository settings or merges not yet
  published does Publish stop with "Remote user file changed since your last
  Pull", because a Pull would bring the repository's settings over yours:
  note them, Pull latest, check them, then Publish again.

#### “Your published file … is not a valid annotation file”

- Cause: your own file on GitHub fails the annotation contract (for example a
  file written for the earlier, version-2 contract). Pull leaves it out, and
  Publish will not overwrite it.
- Fix: fix or remove the file on GitHub (or ask an author to), then Pull and
  Publish again.

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

#### Repo does not show up under “Choose repository”

- Ensure the GitHub App is installed on the repo owner and includes this repo.
- If installed with “Only select repositories”, add the repo explicitly.
- In the GitHub sync modal: use **Reload**.

#### Pull warns that files were left out

- Cause: one or more user files or the merges file are invalid JSON, violate
  their schema, exceed 1,000,000 UTF-8 bytes, or were written for the earlier
  contract (version 2); or paths under `annotations/` are not annotation files.
- Result: the Pull completes without those files; the warning says how many
  were left out and its details list each file with its reason. A file is not
  downloaded again until it changes.
- Fix:
  - run the template validation script in the annotation repo (`python scripts/validate_user_files.py`)
  - fix or remove the offending files
  - communicate with the contributor if needed

#### Pull fails: a schema “is not the schema this Cellucid reads”

- Cause: `annotations/schema.json`, `annotations/config.schema.json` or
  `annotations/moderation/merges.schema.json` differs from the template's
  file (an older version, an edit, or other line endings).
- Fix: copy the named file again, unchanged, from
  [theislab/cellucid-annotation](https://github.com/theislab/cellucid-annotation)
  (with `.gitattributes`), commit it, and Pull again.

#### Pull fails: GitHub truncated the listing of `annotations/users/`

- Cause: GitHub lists at most 100,000 entries or 7 MB of a tree
  ([GitHub: get a tree](https://docs.github.com/en/rest/git/trees#get-a-tree)).
- Result: Cellucid stops rather than compiling the part it saw.

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
- Sign-in that always fails after GitHub's authorization page usually means
  the GitHub App does not expire user tokens: keep **Expire user
  authorization tokens** on; the Worker refuses a token that never expires.
- “… runs an older version (contract N) than this Cellucid (contract 3); its
  owner must redeploy the Worker.”: deploy the Worker from the same Cellucid
  source as the page. A newer Worker asks you to reload Cellucid instead.
- **Add repo** opening the wrong GitHub App: the app's slug is compiled into
  the connection wizard (`GITHUB_APP_INSTALL_URL` in
  `assets/js/app/ui/modules/community-annotation/github-connection-flow.js`);
  set it to your app.
