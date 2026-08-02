# Lifecycle Walkthrough (Screens, Step by Step)

**Audience:** everyone — this page assumes no GitHub experience  
**Time:** 20–30 minutes  
**What you'll learn:** every screen of one complete annotation round, in order,
with a picture of each

The other pages in this section explain *why*. This page shows *what you will
see*, from an empty panel to a published vote. Follow it top to bottom the first
time; afterwards use it as a picture index.

Every figure was captured from the running application against the published
**Pancreas** sample, with its `clusters` column opened for annotation and the
`Beta` cluster used as the worked example. Your column and category names will
differ; every control, badge, and message is the same.

The GitHub side is a local stand-in, so the repository and people named in the
figures are fictional. What the application computes from them — the vote
tallies, the confidence percentages, the consensus verdicts — is the product's
own arithmetic on those inputs, not something drawn onto the pictures.

:::{tip}
Reading this before you have a repository? Steps 1–4 are the only ones that need
anything set up. If someone has already sent you a repository name, you can skip
to step 3.
:::

---

## Step 1 — Find the panel

Community Annotation is an accordion in the left sidebar. It starts collapsed,
and before a repository is connected it holds exactly one action.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-01-disconnected.png
:alt: The Community Annotation sidebar panel reading "No annotation repo connected." above a Connect repo button.
:width: 488px

Before a repository is connected the panel states `No annotation repo
connected.` and offers only **Connect repo**. Nothing else in the feature is
reachable yet, and nothing is hidden behind a disabled control.
```

This is also what you see when your role cannot be established — see **Step 13**
below.

---

## Step 2 — Start the connection wizard

The pointer in the figure above is on **Connect repo**, which is the only thing
to press here. It opens a four-step wizard. It does not sign you in by itself,
and it does not contact GitHub until you ask it to.

---

## Step 3 — Sign in with GitHub

The wizard's four steps are always visible down the left of the dialog, so you
can see how far along you are: **Sign in with GitHub**, **Install the GitHub
App**, **Select an annotation repository**, **Sync (pull / publish)**.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-02-signin.png
:alt: The GitHub sync dialog on step 1 of 4, headed "Sign in with GitHub", with a Continue with GitHub button and a header reading GitHub not connected, Repo not connected.
:width: 1040px

Step 1 of 4. The dialog header tracks three facts throughout: the **dataset**
you have open, your **GitHub** identity, and the **repo** you are bound to. Here
the last two both read `not connected`.
```

**Continue with GitHub** sends you to GitHub to sign in and returns you here.
As the dialog says, signing in does not by itself give Cellucid access to any of
your repositories; that is the next step, and you choose which ones.

:::{important}
The sign-in token is held in `sessionStorage` only. Closing the browser tab
discards it, and a shareable `?annotations=owner/repo` link never carries it.
Every collaborator signs in as themselves.
:::

---

## Step 4 — Choose the annotation repository

Once the Cellucid GitHub App is installed on the account or organisation that
owns the repository, the repositories it can see are listed for you.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-03-install-app.png
:alt: Step 2 of 4 of the wizard listing three repositories as cards, each marked Public or Private, above a pagination status reading "Showing 1-3 of 3 repositories. Page 1 of 1."
:width: 1040px

Step 2 of 4. Each repository is a card marked **Public** or **Private**, and the
status line under the list always says exactly how many there are and which page
you are on. **Add repo** opens GitHub so you can enable more; **Reload** picks up
what you changed without restarting the wizard.
```

The list is paged and filterable because an organisation can expose thousands of
repositories. If your repository is missing, it is almost always because the
GitHub App has not been given access to it — use **Add repo**, then **Reload**.

Step 3, **Select an annotation repository**, is where you actually pick one and
press **Connect**. When the list is long enough to need it, step 3 also offers a
filter box.

### When the filter matches nothing

```{figure} ../../../_static/screenshots/community_annotation/state-no-repository-matches-filter.png
:alt: Step 3 of 4 of the wizard with a filter applied that matches no repository, showing the message "No repositories match this filter." above a status line reading "No repositories to display."
:width: 1040px

An empty result says `No repositories match this filter.` and the status line
reads `No repositories to display.` This is a filter result, not a permission
problem: clear the filter to get the full list back.
```

---

## Step 5 — Confirm you are connected

After **Connect**, the sidebar panel changes from one button to the working
surface for the whole feature.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-06-connected-status.png
:alt: The connected Community Annotation panel showing Dataset, GitHub and Repo rows, a GitHub sync button, a Profile block, and the collapsed Manage annotation, Derived consensus column, and Consensus snapshot sections.
:width: 488px

The connected panel. **Github sync** re-opens the wizard at step 4. **Profile**
holds the optional name, title, ORCID, and LinkedIn that get written into your
user file when you publish — they are saved locally first, like your votes.
```

The three rows at the top — dataset, GitHub user, repo — are the *scope*. If any
one of them changes, Cellucid opens a different annotation workspace. That is
why switching datasets appears to "lose" your votes: they are still there, filed
under the other dataset.

---

## Step 6 — Open a column for annotation (authors only)

An author decides which categorical columns the group may vote on.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-07-manage-annotation.png
:alt: The Manage annotation section expanded, showing a categorical column picker with Add, Remove and Close buttons, and an Annotatable consensus settings block with a Threshold slider reading 50% and a Min annotators field.
:width: 488px

**Manage annotation** is present for everyone but only operable by an author.
Pick a categorical column, press **Add**, and set the two numbers that decide
consensus for it: **Threshold** and **Min annotators**. **Apply** stages the
change locally; publishing writes it to `annotations/config.json`.
```

Threshold is stored as a signed fraction and edited as whole percent. Its full
range is −100 % to +100 %, because the underlying confidence score is itself
signed: a column can be configured to accept a label that more people opposed
than supported, which is almost never what you want but is representable.

The out-of-the-box defaults, when a column has no explicit settings, are
**1 minimum annotator** and a **0.5 threshold**.

---

## Step 7 — Spot a votable column

Once a column is open for annotation it is badged everywhere it is offered.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-08-votable-badge.png
:alt: The Coloring and Filtering categorical obs dropdown showing the clusters entry prefixed by a ballot-box badge, above the legend listing the cluster categories.
:width: 458px

The 🗳️ badge in the **Categorical obs** dropdown marks a column that is open for
voting. A column with no badge is an ordinary column; a column the author has
closed carries a finish-flag badge instead.
```

---

## Step 8 — Pull before you work

**Pull latest** downloads every contributor's file and rebuilds the merged view
in your browser. Do this at the start of a session, so you are voting on what
everyone else has already proposed rather than duplicating it.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-05-sync.png
:alt: Step 4 of 4 of the wizard, headed "Sync (pull / publish)", with Pull latest and Publish buttons, an Auto pull selector, and a status line reading "Pulled latest annotations."
:width: 1040px

Step 4 of 4. The status line under the buttons reports the outcome verbatim —
here, `Pulled latest annotations.` **Auto pull** can repeat this every 10, 15, or
60 minutes while the tab stays open.
```

A Pull is transactional. If any file is unreadable, oversized, or fails schema
validation, the pull is refused and the annotation view you already had is left
untouched — you never end up with half of someone's votes.

---

## Step 9 — Read the consensus

The merged result appears against each category in the legend, so you can see
the state of the whole column at a glance.

### The three states, exactly

Cellucid computes, for one column and one category:

- **voters** — how many distinct people cast any vote in that bucket, counting a
  person once no matter how many suggestions they voted on;
- **net votes** — for the leading suggestion, upvotes minus downvotes;
- **confidence** — net votes divided by voters, so it runs from −1 to +1.

Suggestions are ranked by net votes first, then by raw upvote count. The state
follows from three comparisons, in this order.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-09-pending.png
:alt: The voting modal for the Beta cluster with a status banner reading Pending and a voter count below the column minimum, above a single suggestion card.
:width: 996px

**Pending** — fewer voters than the column's minimum. Here two people have voted
and the column requires three, so no confidence figure is reported at all; the
banner names only the shortfall. The label is still shown, because the work is
real, but it is not yet a decision.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-11-consensus.png
:alt: The voting modal for category beta with a status banner reading Consensus, "Beta cell", 100 percent, 4 voters, above a suggestion card showing net 4 with ontology, marker and evidence lines.
:width: 996px

**Consensus** — enough voters, one clear leader, and confidence at or above the
threshold. Four people voted, all four supported *Beta cell* and none opposed
it, so confidence is 100 %.
```

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-12-disputed.png
:alt: The voting modal for the Beta cluster with a status banner reading Disputed and naming two tied labels at 0 percent with four voters, above two suggestion cards each showing net 0 with two upvotes and two downvotes.
:width: 996px

**Disputed** — everything else. Here two suggestions are tied on both net votes
and upvotes, so there is no single leader; a tie is always disputed, whatever
the confidence figure says. When the top is tied, the banner names every tied
label, alphabetically, separated by commas.
```

:::{note}
Disputed is a normal, healthy result. It means the group has looked at the
category and does not agree yet — which is information. Resolve it with marker
evidence, a comment thread, or an author merge, not by lowering the threshold.
:::

---

## Step 10 — Vote, discuss, and propose

Click a category name in the legend to open its voting modal.

### A suggestion card

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-10-vote.png
:alt: A single suggestion card headed "Beta cell" with its net score at the right, lines for Ontology, Markers and Evidence, up and down vote buttons, Edit, Delete and Merge buttons, an attribution line naming the proposer, and a comment input, with the pointer on the upvote button.
:width: 1028px

One suggestion. The **▲** and **▼** buttons cast your vote — pressing the one you
already chose withdraws it. **Edit** and **Delete** appear only on suggestions
you proposed. The attribution line names the proposer, with their profile title
if they set one, so a label can be discussed with the person who suggested it.
```

Every card carries a comment box. Comments are how a disputed bucket gets
resolved: state the evidence, not the conclusion.

### Proposing a label

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-13-new-suggestion.png
:alt: The New suggestion form filled in with a label, an ontology id and two marker genes, showing Search CAP, Search Ontology and Search Markers buttons beside the fields and Add and Clear buttons below.
:width: 984px

**New suggestion** takes a label plus three optional pieces of evidence: an
ontology id, marker genes, and free text. The three search buttons look the
proposal up in the Cell Annotation Platform by name, by ontology id, or by
marker set, and fill the fields in for you.
```

Only **Label** is required. Filling in an ontology id is what makes a label
comparable across studies later, and it costs one click if you use **Search
CAP**.

---

## Step 11 — Publish

Nothing you have done so far is visible to anyone else. Everything up to this
point lives in your browser.

```{figure} ../../../_static/screenshots/community_annotation/lifecycle-17-publish.png
:alt: Step 4 of 4 of the wizard with a status line reading "Publish complete."
:width: 1040px

**Publish** writes your work to GitHub and reports `Publish complete.` when the
write has been accepted.
```

Publishing writes exactly one file: your own
`annotations/users/ghid_<your-github-id>.json`. That is the whole reason the
repository is laid out one-file-per-person — two people annotating at the same
time never touch the same file, so there is nothing to merge and nothing to
conflict.

An author publishing repository settings also writes `annotations/config.json`
and, if merges changed, `annotations/moderation/merges.json`. Files that did not
change are not rewritten.

Cellucid chooses **one** publication route before it writes anything, from the
permissions GitHub reports for you: a direct write if you have push access, or a
fork plus pull request if you do not and the repository allows forking. If the
chosen route fails, that failure is final and is reported — Cellucid never
quietly tries the other route.

---

## Step 12 — Moderate duplicates (authors only)

Independent annotators reliably propose the same label twice with different
wording. An author can fold one into the other.

**Merge…** on a suggestion card opens the merge form inline: choose the
suggestion to merge into, optionally record why, and confirm. Merges are
recorded in `annotations/moderation/merges.json` and are reversible — a merged
bundle can be opened with **View merged** and detached again.

A merge changes how votes are counted, not what anyone voted. Nobody's file is
edited, and the merge is a separate, auditable record. Votes inside a merged
bundle are de-duplicated when the totals are rendered, so a person who voted for
both of two merged suggestions still counts once.

---

## Step 13 — What you see when something is wrong

### You are not signed in

You get the step 1 panel: it reports `No annotation repo connected.` and offers
**Connect repo**. There is no half-signed-in state in which some controls work.

### You can read the repository but not write to it

### Your role cannot be determined

```{figure} ../../../_static/screenshots/community_annotation/state-repository-unreadable.png
:alt: The full application window with the Community Annotation panel back in its disconnected state and a notification explaining that the user's role for the repository could not be determined and that the annotation repo was disconnected.
:width: 1440px

If GitHub cannot tell Cellucid what your permissions are — the repository was
renamed, deleted, or made private, or access was withdrawn — the repository is
disconnected rather than left in an ambiguous state, and the notification names
the repository and repeats what GitHub said.
```

Disconnecting is deliberate. A guess at your role would mean either offering an
author control that will fail on write, or hiding one you are entitled to.

### Everything looks disabled

Three causes, in order of likelihood:

1. The column is **closed** — its badge is 🗳️🏁 rather than 🗳️. Ask the author.
2. You are signed out. The panel will say so.
3. The column was never opened for annotation. It will have no badge at all.

---

## Where things are stored

| What | Where it lives | Survives a tab close? |
|---|---|---|
| Your GitHub sign-in token | `sessionStorage` | No |
| Your votes, suggestions, comments, profile | Annotation-specific local storage | Yes |
| Files downloaded by Pull | IndexedDB, with the SHA index in local storage | Yes |
| The shared record | The GitHub repository | Yes, for everyone |

Community annotation state is deliberately **not** part of a
`.cellucid-session` bundle. Saving a session does not carry your votes, and
sending someone a session file does not send them your annotations. Use
**Publish** for that. See
{doc}`../l_sessions_sharing/02_what_gets_saved_and_restored`.

---

## The round, end to end

1. Connect the repository (steps 1–5), once per dataset and browser.
2. **Pull latest** before you start.
3. Vote, comment, and propose (steps 9–10).
4. **Publish** (step 11) so the group can see your work.
5. Repeat 2–4 during the round; an author merges duplicates (step 12) and tunes
   the consensus settings (step 6) as evidence accumulates.
6. At the end of the round the author closes the columns, everyone pulls one
   last time, and the consensus is exported as a snapshot or built into a
   derived obs column.

---

## Next steps

- {doc}`01_annotator_guide` — the contributor workflow in prose, with edge cases
- {doc}`02_author_guide` — setting up the repository and running the round
- {doc}`03_ui_reference` — every control, and the troubleshooting tables
