# Community Annotation — Annotator Guide (UI + Voting)

This guide is for people contributing **votes, suggestions, and comments** in Cellucid.

If you are setting up the annotation repo or changing which columns are annotatable, read {doc}`02_author_guide`.

:::{tip}
If you would rather see each screen before you read about it, work through
{doc}`04_lifecycle_walkthrough` first. It is the same round in order, with a
picture of every step and of every failure state.
:::

---

## Quickstart (Pick the Path That Matches You)

::::{tab-set}

:::{tab-item} Wet-Lab / New to GitHub

You can do this without knowing Git, branches, or JSON.

1) Open Cellucid and load the dataset.
2) Open **Community annotation** in the left sidebar.
3) Click **Connect repo**, then **Continue with GitHub**.
4) At **Choose repository**, choose the repo your author gave
   you and click **Connect**.
5) Click **Pull latest** (this downloads everyone’s shared contributions).
6) Choose the 🗳️-marked column (your author will usually tell you which one).
7) Click a cluster/category name to vote.
8) When you’re done, click **Publish** so others can see your work.

If anything fails, jump to “Troubleshooting (Massive)” at the end and start with “I can’t find the repo…”.
:::

:::{tab-item} Computational / Git-Savvy

1) Verify you are on the intended scope: dataset id + `owner/repo@branch`.
2) Pull (reads one commit; downloads only files whose Git blob id changed) and confirm the merged view is current.
3) Vote/suggest; include ontology ids and evidence when possible.
4) Publish directly if permitted, or by fork + Pull Request only when the
   repository allows forking, and confirm any PR is merged.
5) Ask others to Pull to refresh their merged view.

This guide includes schema limits, caching behavior, and edge cases (tabs, devices, forks, branch mismatches).
:::

::::

---

## 0) What You Need (Annotator Checklist)

- A GitHub account you can use to sign in.
- Access to the annotation repository (public repo, or you’ve been added to a private repo).
- The **Cellucid GitHub App** installed on the repo owner (so the repo appears in the UI).

If you cannot push to the repository (fork + Pull Request), do two things on GitHub first:

- **Fork the repository** (`https://github.com/<owner>/<repo>/fork`). Cellucid never creates a fork for you; without one, Publish says so and links there.
- **Install the Cellucid GitHub App on your own GitHub account with access to that fork**: *All repositories*, or *Only select repositories* with the fork selected. Cellucid writes your user file in your fork with your sign-in, which reaches only repositories the app is installed on.

:::{important}
Your work is saved in two places:

- **Local first**: votes/suggestions/comments are saved in your browser immediately.
- **Shared**: other people only see your work after you **Publish** to GitHub (direct push or PR).

If you clear browser site data before publishing, you can lose unpublished work.
:::

---

## 1) Connect to the Annotation Repo (First Time)

1) Load the dataset in Cellucid.
2) Open the **Community annotation** accordion in the left sidebar.
3) Click **Connect repo**.
4) In the four-step **GitHub sync** wizard:
   - use **Continue with GitHub**
   - use **Add repo** and **Reload** if the required repository is not enabled
   - select the repo and click **Connect**
   - at **Sync**, click **Pull latest**

The picker resolves the repository's default branch. If the annotation round
uses a different branch, open the author's canonical
`?annotations=owner/repo@branch` link rather than guessing.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-01-disconnected.png
:alt: The Community annotation panel before any repository is connected: Repository with its (i), the line Not connected. and a Connect repo button, with the pointer on that button.
:width: 504px

The Community Annotation entry point reports that no repository is connected and offers the explicit Connect repo action.
```

### If you were given a pre-filled link

Authors can share a link like:

- `?annotations=owner/repo`
- `?annotations=owner/repo@branch`

This can pre-select the annotation repo, but you still must sign in and Pull.

---

## 2) Pull Latest (What It Does and Why It Matters)

When you click **Pull latest**, Cellucid:

- reads the repository at one commit and lists `annotations/users/*.json` (and optional `annotations/moderation/merges.json`)
- downloads only files whose Git blob id changed since your last Pull, many
  per request; when nothing changed, it stops after one or two requests
- rebuilds the merged suggestion/vote view in your browser

The persistent raw-file cache stores what each file holds for this dataset in
IndexedDB, kept as each batch arrives, so an interrupted Pull resumes where it
stopped. Annotation state and cross-tab coordination require local storage.
Cellucid never switches to an in-memory cache: if either storage boundary is
unavailable, Pull fails without replacing the current annotation view.

A contributor file that is not valid never stops your Pull: it is left out of
the consensus, and a warning says how many files were left out (each file and
its reason are under the warning's details). If it is your own file, the
warning says so first, and Publish will not overwrite it until it is fixed or
removed on GitHub; tell your author.

If you do not Pull:

- you may be voting on an outdated view (missing others’ new suggestions),
- you may accidentally create duplicates that already exist on GitHub,
- you may publish “stale” assumptions (your local view differs from the group’s current state).

:::{tip}
For everyday work, a good pattern is:

1) Pull latest
2) Annotate for 10–30 minutes
3) Pull latest again (quick sanity check)
4) Publish
:::

---

## 3) What You Can and Cannot Change

### You can

- vote **▲ up** or **▼ down** on suggestions
- propose new suggestions (label + optional ontology id/markers/evidence)
- comment on suggestions (and edit/delete your own comments)
- edit/delete your own suggestions
- Publish your own user file directly if permitted, or through a fork + Pull
  Request only when the repository allows forking

### You cannot (by design)

- decide which obs columns are annotatable
- change consensus thresholds / min annotators
- close/reopen a column
- merge duplicates

Those are author-only controls derived from GitHub repo permissions (maintain/admin).

---

## 4) Find an Annotatable Column in the UI (🗳️ Badge)

After you Pull, the author-configured annotatable columns appear as a 🗳️ badge in the categorical field dropdown.

- 🗳️ = annotation is enabled for this column
- 🗳️🏁 = the author has **closed** the column (voting disabled for annotators)

---

## 5) Open the Voting Modal for a Category (Cluster)

1) Select the annotatable categorical column (e.g. `leiden`).
2) In the categorical legend, click a category label (cluster).
3) A voting modal opens for that category.

In annotation mode:

- category labels are **locked** (no renaming/merging categories)
- each category row shows a small consensus summary and vote counts

:::{tip}
If a category is not clickable, check:

- whether you are signed in and still connected to the repo
- whether the category currently has **0 cells** (after filters) — categories with no available cells are disabled
:::

---

## 6) Read the Consensus Line (Pending / Disputed / Consensus)

At the top of the voting modal you will see one of:

- **Consensus:** a single leading label has enough support
- **Disputed:** there is no clear winner (often due to close votes or ties)
- **Pending:** not enough unique voters yet (author-configured `minAnnotators`)

### What “confidence” means (so you don’t misread it)

Cellucid computes:

- `voters`: unique users who cast any vote in this category (across all suggestions)
- `netVotes`: for the current leading suggestion, `upvotes - downvotes`
- `confidence = netVotes / voters` (range -1..1)

So:

- downvotes reduce confidence,
- one person voting multiple times still counts as **one voter** (denominator is unique users),
- suggestions are ranked first by net vote and then by upvote count; a top tie
  exists only at equal net vote and equal upvote count, and is Disputed.

If you are not sure what to do:

- add a comment with evidence,
- vote only when you feel confident the label is right/wrong.

---

## 7) Vote on Suggestions (▲ / ▼)

Each suggestion card has:

- **▲** upvote button (support)
- **▼** downvote button (oppose)
- a net score (`net up-down`)

Voting behavior:

- Clicking ▲ casts an upvote.
- Clicking ▼ casts a downvote.
- Clicking the same button again removes your vote (returns to “no vote”).
- A removed vote, or a deleted comment, stays removed when you Pull before
  publishing again: your own published file only fills in local state when it
  changed since your last Publish or Pull (for example, after publishing from
  another device).

:::{important}
Downvotes are part of the consensus math. Use them deliberately:

- Downvote when a label is clearly wrong for this cluster.
- Prefer comments for nuance (“this might fit if…”) especially early in a round.
:::

### Advanced: Can you vote on multiple suggestions?

Yes. If you want, you can upvote one label and downvote another.

Remember:

- The consensus denominator counts unique voters across the whole category, not per suggestion.
- Strategic downvoting can suppress incorrect labels, but it can also create Disputed states if the group is split.

---

## 8) Comment on Suggestions (and Edit/Delete Your Own)

Under each suggestion card:

- type a comment
- press **Enter** to submit
- use **Shift+Enter** to insert a newline without submitting

You can edit or delete your own comments (other people’s comments are read-only).
An edit in progress survives the modal refreshing around it; if the comment you
are editing disappears, the editor closes and tells you the edit was not saved.

Good comments include:

- evidence (“marker genes: …”, “matches reference atlas …”)
- reasons for disagreement (“this cluster expresses … so not …”)
- requests (“can someone check …?”)

:::{note}
Avoid personal data in comments (emails, patient info, private sample identifiers). Annotation repos are meant to be shareable within a team and may be public.
:::

---

## 9) Propose a New Suggestion (Label + Optional Evidence)

At the bottom of the modal you can add a **New suggestion**:

- **Label (required)**: short human-readable label
- **Ontology id (optional)**: e.g. `CL:0000625`
- **Marker genes (optional)**: comma-separated (e.g. `MS4A1, CD79A`)
- **Evidence (optional)**: free text

No text has a length limit of its own; your whole user file (every dataset you
annotate in the repository) is limited to 1,000,000 UTF-8 bytes, the most
GitHub serves as one JSON file. Text cannot start or end with whitespace.

Buttons you may see:

- **Add**: create the suggestion (saved locally until Publish)
- **Clear**: reset the form
- **Search CAP / Search ontology / Search markers**: helper searches (network-dependent)

### Best practices to avoid duplicates (high impact)

Before adding a new suggestion:

1) Read the existing suggestion labels.
2) If your label already exists but needs evidence, add a comment instead of duplicating.
3) If your label exists but is spelled differently, prefer commenting and ask an author to merge duplicates later.

Exactly the same label cannot be added twice in one category (Cellucid
refuses it), and if two people propose it independently, for example offline,
Cellucid counts it as one suggestion: each person's votes on it count once.
Labels that differ even slightly (`T cell`, `T-cell`, `T cells`) are separate
suggestions and split votes until an author merges them.

### About CAP search (privacy + failure modes)

Cellucid sends a bounded request to the configured Cellucid Worker at
`/cap/lookup-cells`. The Worker chooses a pinned CAP persisted query and relays
only the required variables; the browser does not send an arbitrary GraphQL
document or caller-selected upstream.

- CAP searches can help you find ontology ids, synonyms, and marker genes.
- Search text is visible to the Cellucid Worker operator, CAP, and the network
  path. The GitHub bearer token and OAuth cookies are not sent to CAP.
- A successful response has the bounded shape
  `{ "contractVersion": 3, "results": [...], "omittedInvalidCount": number }`.
- If this disclosure is not acceptable, or if you are offline or CAP is
  blocked, use manual suggestions; voting and manual annotation still work.

---

## 10) Publish Your Work (Direct Push vs Pull Request)

In the **GitHub sync** modal, click **Publish**.

What happens depends on your GitHub permissions:

- If GitHub reports push access, Cellucid writes directly to
  `annotations/users/ghid_<yourNumericId>.json`.
- Without push access, Cellucid uses a fork + Pull Request only when the
  repository allows forking. Otherwise Publish is disabled.

Cellucid selects exactly one route from the repository metadata before making
any mutation. If that route fails, the operation stops and does not switch to
the other route.

### Before you publish (recommended)

1) Pull latest (so you’re not publishing against a stale view).
2) Check whether your intended label already exists (avoid duplicates).
3) Publish.

### If a Pull Request is opened

You have one pull request per branch of the annotation repository. Publishing
again updates the same pull request (it carries your whole file, every
dataset you annotate there) instead of opening another.

1) Confirm the PR was created successfully.
2) Follow your lab/org process:
   - request review if needed
   - respond to feedback
   - wait for merge
3) Your contribution becomes visible to others after the PR is merged and they Pull.

:::{tip}
If others say “I can’t see your votes”, the top three causes are:

1) you didn’t Publish,
2) your PR is not merged,
3) they didn’t Pull after your publish/merge.
:::

---

## 11) “I’m Done” (End of an Annotation Session)

Recommended wrap-up:

- Publish your changes (or confirm your PR exists and is correct).
- Pull latest once more (sanity check that you’re up to date).
- Leave a final comment if you want to summarize your reasoning.

Remember:

- your GitHub token expires after 8 hours, and closing the tab or **Disconnect
  GitHub** ends it sooner (**Disconnect GitHub** also revokes it at GitHub and
  disconnects this dataset's repository),
- your local annotation session may persist in browser storage, but it is not a substitute for publishing.

---

## 12) Edge Cases (Read When Something Feels Weird)

### “I voted on my laptop, but my desktop doesn’t show it”

Local work is stored per browser/device. To share across devices:

1) Publish from the device where you made the votes.
2) Merge PR if applicable.
3) Pull on the other device.

### “Open in another tab — annotation paused here”

One tab of a browser annotates a dataset in a repository at a time (the same
dataset, repository@branch and GitHub user), so two tabs never overwrite each
other. A second tab says `Open in another tab — annotation paused here`;
everything else of the dataset works there.

- Fix: click **Use here** in the tab you want to annotate in. The other tab
  pauses in turn; what it had is kept.

### “Everything is disabled”

Common causes:

- the author closed the field (🗳️🏁)
- you are signed out (token expired: the panel says `Your GitHub sign-in
  expired or was revoked…`; sign in again and the repository is connected
  again, your unpublished changes still there)
- you are offline (Pull/Publish disabled)
- another tab annotates the same dataset (`Open in another tab — annotation
  paused here`, **Use here**)

### An agent voted for me

An agent you work with (the chat bar, Claude or ChatGPT, a local app) or a
Python or R notebook can vote only after you allow it (you are asked, unless
you set Publishing to Allow), and its vote is published alone: your drafts
and other unpublished changes stay yours. Undo turn cannot take a published
vote back; vote again to change it. See
{ref}`the UI reference <community-annotation-agents>`.

---

## 13) Troubleshooting (Massive)

If you don’t find your issue here, also check {doc}`03_ui_reference` (it catalogs UI controls and more error messages).

### I can’t find the repo at “Choose repository”

Most common causes:

- The Cellucid GitHub App is not installed on the repo owner, or the repo wasn’t selected during installation.
- You are signed into the wrong GitHub account.
- The repo is private and you do not have access.

Fix:

1) In the GitHub sync modal, verify your GitHub username (top of the modal).
2) Ask the author to confirm the app installation includes the repo.
3) Click **Reload** in the GitHub sync modal and retry.

### “Dataset mismatch” / “Ask an author to Publish updated settings”

- Cause: the connected repo does not list the current dataset id in `annotations/config.json`.
- Fix: ask an author to connect and Publish (this updates `supportedDatasets[]` and unblocks annotators).

### “Pull latest” fails

Common causes:

- offline / flaky network / corporate proxy
- repo structure errors (missing template files, or a schema file that is not
  the template's)
- the annotation service is older than this Cellucid ("… runs an older
  version …; its owner must redeploy the Worker")
- GitHub's hourly request budget for your account is spent. Cellucid paces
  its requests by that budget and waits when GitHub names a wait (the
  progress notification says so); a Pull that still meets the limit stops
  with a `GitHub rate limit:` message.

What you can do:

1) Retry after a minute (or after the time a rate-limit message names).
2) Verify you can access GitHub in a normal browser tab.
3) Tell the author the exact error message; authors can validate repo structure and CI.

### I can’t vote / everything is disabled

- The author may have **closed** the column (🗳️🏁).
- You may not be connected to the repo anymore (session expired, signed out).

Fix:

1) If closed: contact the author (only authors can reopen).
2) If signed out: sign in again and Pull.

### I clicked a category but nothing happens

- The category may currently have **0 cells** (due to filtering). Disabled categories cannot be opened.

Fix:

- relax filters or switch to a category with cells.

### Publish is not possible for me

- If you have no push access and the repo disables forking, Cellucid cannot create a PR.

Fix:

- ask the author to enable forking or grant write access.

### Fork + PR publish fails (but the author says forking is enabled)

What Publish says, and what to do:

- `You have no fork of <owner>/<repo> yet: fork it on GitHub, then Publish again.`: open **Fork on GitHub**, create the fork, then Publish again.
- `GitHub is still copying your fork …; Publish again in a minute.`: GitHub creates a fork in the background; wait a minute, then Publish again.
- `Cellucid cannot write to your fork …`: install the Cellucid GitHub App on your own account with access to the fork (**Install the app**; choose *All repositories*, or *Only select repositories* with the fork selected), then Publish again.

### I published, but others don’t see my votes

Checklist:

- If you published via PR: it must be merged first.
- Others must Pull after your publish/merge.
- Ensure you and others are on the same repo **and branch**.

### CAP search doesn’t work

- You might be offline, behind a restrictive firewall, or CAP may be blocked.

Fix:

- proceed with manual entry; optionally ask your author if CAP is allowed in your environment.
