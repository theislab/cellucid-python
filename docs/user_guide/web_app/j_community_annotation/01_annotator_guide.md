# Community Annotation — Annotator Guide (UI + Voting)

This guide is for people contributing **votes, suggestions, and comments** in Cellucid.

If you are setting up the annotation repo or changing which columns are annotatable, read `02_author_guide`.

---

## Quickstart (Pick the Path That Matches You)

::::{tab-set}

:::{tab-item} Wet-Lab / New to GitHub

You can do this without knowing Git, branches, or JSON.

1) Open Cellucid and load the dataset.
2) Open **Community Annotation** in the left sidebar.
3) Click **Connect GitHub…** and sign in.
4) Choose the repo your author told you to use.
5) Click **Pull latest** (this downloads everyone’s shared contributions).
6) Choose the 🗳️-marked column (your author will usually tell you which one).
7) Click a cluster/category name to vote.
8) When you’re done, click **Publish** so others can see your work.

If anything fails, jump to “Troubleshooting (Massive)” at the end and start with “I can’t find the repo…”.
:::

:::{tab-item} Computational / Git-Savvy

1) Verify you are on the intended scope: dataset id + `owner/repo@branch`.
2) Pull (SHA-based incremental download) and confirm the merged view is current.
3) Vote/suggest; include ontology ids and evidence when possible.
4) Publish (direct push if permitted; otherwise fork + PR) and confirm the PR is merged.
5) Ask others to Pull to refresh their merged view.

This guide includes schema limits, caching behavior, and edge cases (tabs, devices, forks, branch mismatches).
:::

::::

---

## 0) What You Need (Annotator Checklist)

- A GitHub account you can use to sign in.
- Access to the annotation repository (public repo, or you’ve been added to a private repo).
- The **Cellucid GitHub App** installed on the repo owner (so the repo appears in the UI).

If your author uses the fork + Pull Request model, one extra thing helps a lot:

- Install the Cellucid GitHub App on **your personal GitHub account** with access to **all repositories**, so newly created forks are included automatically.

:::{important}
Your work is saved in two places:

- **Local first**: votes/suggestions/comments are saved in your browser immediately.
- **Shared**: other people only see your work after you **Publish** to GitHub (direct push or PR).

If you clear browser site data before publishing, you can lose unpublished work.
:::

---

## 1) Connect to the Annotation Repo (First Time)

1) Load the dataset in Cellucid.
2) Open the **Community Annotation** accordion in the left sidebar.
3) Click **Connect GitHub…**
4) In the **GitHub sync** modal:
   - sign in
   - choose the repo (and branch, if your group uses a non-default branch)
   - click **Pull latest**

<!-- SCREENSHOT PLACEHOLDER
TYPE: Screenshot (Entry point)
Suggested filename: community_annotation/10_connect_github_button.png

Capture:
- Community Annotation accordion with the "Connect GitHub…" button visible (or "Choose repo…" / "GitHub sync…" depending on state).

Goal:
- Make it obvious where to start.

Crop / framing:
- Include the left sidebar and the Community Annotation section header.

Redact:
- Private dataset ids or repo names if visible.

Figure caption:
- Example: "Start from the Community Annotation accordion and click Connect GitHub…"

Alt text:
- Example: "Community Annotation section in the sidebar with the Connect GitHub button."
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing the Connect GitHub button in the Community Annotation section.
:width: 100%

Start from the Community Annotation accordion and click Connect GitHub…
```

<!-- SCREENSHOT PLACEHOLDER
TYPE: Screenshot (GitHub sync modal - repo selection and Pull)
Suggested filename: community_annotation/11_github_sync_repo_pull.png

Capture:
- GitHub sync modal after sign-in.
- Repo list visible (or a selected repo card), and the "Pull latest" button visible.

Goal:
- Show how to pick a repo + branch and Pull.

Crop / framing:
- Crop to the modal; include the branch selector if present (important for groups using non-default branches).

Redact:
- Private repo/org names if needed.

Figure caption:
- Example: "Use Pull latest to download the current community files from GitHub."

Alt text:
- Example: "GitHub sync modal showing repo selection and Pull latest."
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing repo selection and Pull latest in the GitHub sync modal.
:width: 100%

Use Pull latest to download the current community files from GitHub.
```

### If you were given a pre-filled link

Authors can share a link like:

- `?annotations=owner/repo`
- `?annotations=owner/repo@branch`

This can pre-select the annotation repo, but you still must sign in and Pull.

---

## 2) Pull Latest (What It Does and Why It Matters)

When you click **Pull latest**, Cellucid:

- lists `annotations/users/*.json` (and optional `annotations/moderation/merges.json`)
- downloads only files that changed since your last Pull (SHA-based incremental)
- rebuilds the merged suggestion/vote view in your browser

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
- Publish your own user file (direct push if permitted; otherwise PR)

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

<!-- SCREENSHOT PLACEHOLDER
TYPE: Screenshot (Field dropdown with 🗳️)
Suggested filename: community_annotation/12_field_dropdown_ballot_badge.png

Capture:
- The categorical field dropdown open, showing at least one field with a 🗳️ prefix.
- If possible, also show a 🗳️🏁 closed field to teach the difference.

Goal:
- Teach annotators how to identify annotatable columns.

Crop / framing:
- Tight crop around the dropdown and its labels.

Redact:
- Sensitive field names if needed.

Figure caption:
- Example: "Annotatable categorical columns are marked with a 🗳️ badge."

Alt text:
- Example: "Dropdown list of categorical fields with a ballot icon marking annotatable fields."
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing the 🗳️ badge in the categorical field dropdown.
:width: 100%

Annotatable categorical columns are marked with a 🗳️ badge.
```

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

<!-- SCREENSHOT PLACEHOLDER
TYPE: Screenshot (Legend in annotation mode)
Suggested filename: community_annotation/13_legend_annotation_mode.png

Capture:
- Categorical legend in annotation mode showing:
  - a consensus summary (Pending/Disputed/Consensus)
  - ▲/▼ counts (or a visible vote summary)
  - at least one clickable category label

Goal:
- Show annotators where to click to open the voting modal.

Crop / framing:
- Include the legend and the active field label, so readers know which column they are annotating.

Redact:
- Sensitive cluster labels if needed.

Figure caption:
- Example: "In annotation mode, click a category label to open the voting modal."

Alt text:
- Example: "Categorical legend showing vote counts and clickable category labels."
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing the annotating categorical legend.
:width: 100%

In annotation mode, click a category label to open the voting modal.
```

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
- ties between top net-vote suggestions are always Disputed.

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

- type a comment (max 500 characters)
- press **Enter** to submit
- use **Shift+Enter** to insert a newline without submitting

You can edit or delete your own comments (other people’s comments are read-only).

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

- **Label (required)**: short human-readable label (max 120 chars)
- **Ontology id (optional)**: e.g. `CL:0000625` (max 64 chars)
- **Marker genes (optional)**: comma-separated (e.g. `MS4A1, CD79A`)
- **Evidence (optional)**: free text (max 2000 chars)

Buttons you may see:

- **Add**: create the suggestion (saved locally until Publish)
- **Clear**: reset the form
- **Search CAP / Search Ontology / Search Markers**: helper searches (network-dependent)

<!-- SCREENSHOT PLACEHOLDER
TYPE: Screenshot (New suggestion form)
Suggested filename: community_annotation/14_new_suggestion_form.png

Capture:
- The "New suggestion" section showing:
  - Label input
  - Ontology id input (if present)
  - Marker genes input
  - Evidence input
  - Add/Clear buttons (and CAP search buttons if visible)

Goal:
- Make it easy for annotators to fill in the form correctly.

Crop / framing:
- Tight crop around the new suggestion form; include field labels.

Redact:
- None usually needed, but remove private labels if present.

Figure caption:
- Example: "Use New suggestion to propose a label with optional ontology/markers/evidence."

Alt text:
- Example: "New suggestion form with label, ontology id, markers, evidence, and Add button."
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot showing the New suggestion form.
:width: 100%

Use New suggestion to propose a label with optional ontology/markers/evidence.
```

### Best practices to avoid duplicates (high impact)

Before adding a new suggestion:

1) Read the existing suggestion labels.
2) If your label already exists but needs evidence, add a comment instead of duplicating.
3) If your label exists but is spelled differently, prefer commenting and ask an author to merge duplicates later.

Duplicates split votes until an author merges them.

### About CAP search (privacy + failure modes)

Cellucid can query the Cell Annotation Platform (CAP) at `https://celltype.info/graphql`:

- CAP searches can help you find ontology ids, synonyms, and marker genes.
- Your search terms are sent to CAP. If that is not acceptable for your environment, do not use CAP search.
- If you are offline (or CAP is blocked), searches will fail, but manual suggestions still work.

---

## 10) Publish Your Work (Direct Push vs Pull Request)

In the **GitHub sync** modal, click **Publish**.

What happens depends on your GitHub permissions:

- If you have push access: Cellucid writes directly to `annotations/users/ghid_<yourNumericId>.json`.
- If you do not have push access: Cellucid uses a fork + Pull Request flow.

### Before you publish (recommended)

1) Pull latest (so you’re not publishing against a stale view).
2) Check whether your intended label already exists (avoid duplicates).
3) Publish.

### If a Pull Request is opened

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

- closing the tab clears your GitHub token (session-only),
- your local annotation session may persist in browser storage, but it is not a substitute for publishing.

---

## 12) Edge Cases (Read When Something Feels Weird)

### “I voted on my laptop, but my desktop doesn’t show it”

Local work is stored per browser/device. To share across devices:

1) Publish from the device where you made the votes.
2) Merge PR if applicable.
3) Pull on the other device.

### “I opened two tabs and now I’m blocked”

Cellucid may enforce a single active tab per scope (dataset + repo@branch + GitHub user id) to prevent silent data loss.

- Fix: close the other tab/window for the same scope.

### “Everything is disabled”

Common causes:

- the author closed the field (🗳️🏁)
- you are signed out (token expired)
- you are offline (Pull/Publish disabled)

---

## 13) Troubleshooting (Massive)

If you don’t find your issue here, also check `03_ui_reference` (it catalogs UI controls and more error messages).

### I can’t find the repo in “Choose repo…”

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
- GitHub rate limiting
- repo structure errors (missing template files)

What you can do:

1) Retry after a minute (rate limits).
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

Common causes:

- Your fork is not accessible to the GitHub App token (you didn’t install the app on your personal account).

Fix:

1) Install the Cellucid GitHub App on your personal account.
2) Prefer “All repositories” so forks are included automatically.
3) Retry Publish.

### I published, but others don’t see my votes

Checklist:

- If you published via PR: it must be merged first.
- Others must Pull after your publish/merge.
- Ensure you and others are on the same repo **and branch**.

### CAP search doesn’t work

- You might be offline, behind a restrictive firewall, or CAP may be blocked.

Fix:

- proceed with manual entry; optionally ask your author if CAP is allowed in your environment.
