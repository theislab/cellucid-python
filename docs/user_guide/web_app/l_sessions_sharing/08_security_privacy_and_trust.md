# Security, privacy, and trust

**Audience:** everyone sharing sessions (especially for private or clinical data)  
**Time:** 15–35 minutes  
**What you'll learn:**
- What information a session file can contain (and therefore leak)
- What session files intentionally do *not* contain (tokens, full datasets, etc.)
- How Cellucid treats sessions as untrusted input (and what that protects you from)
- Practical "safe sharing" rules for teams

---

## Treat session files like data artifacts

Even though `.cellucid-session` files do not contain the full dataset, they can still encode sensitive information.

Rule of thumb:

> If the dataset is sensitive, assume sessions derived from it are sensitive too.

---

## What a session file can contain (privacy perspective)

A session file holds the **whole history** of the dataset, not only what is on
screen: steps you undid, other branches and other trackers are in it too.
Depending on what you did, it may include:

- **The steps themselves**
  - each step's action and settings, its label, and when it was taken
  - who took it: you, an agent (by its name), or a notebook (Python or R)
  - tracker and branch names you typed
- **Dataset identifiers**
  - the dataset id, name, counts and cell-order digest of the dataset it was
    saved on
  - where the dataset was opened from: a sample, a data link (in its public
    form: never a key, signature or token the link carries) or a server (its
    address); for a file opened from your own computer, its file or folder
    name
  - which cells were open (a subsample's specification and its cells) and how
    many cells genes were read for
- **Field names and category labels**
  - including any renames you applied, and the names before each rename
  - including deleted fields
- **Highlight structure**
  - page names (e.g., `Patient A`, `Responder`; avoid private IDs)
  - group names and on/off state
- **Highlight cells**
  - which cells belong to each group (can be large)
- **Fields you made, and gene sets**
  - their definitions and exact per-cell values
  - gene sets with their names and descriptions
- **Analysis settings and results**
  - which analyses were open (sidebar panel and floating windows) and their
    settings
  - their results: differential expression, marker genes, signatures,
    correlations, dot plots and compositions, i.e. statistics derived from the
    data
- **Figure Export settings**
  - including the figure title and axis labels you typed

Practical implication:
- Do not put patient identifiers or private sample barcodes into names you type
  (highlight pages and groups, categories, trackers) if you intend to share the
  session. Renaming later is not enough: the earlier name stays in the step
  before the rename. Share the **Current state** only (it holds what is shown,
  without the steps before), or leave the branches that hold such steps
  unticked in **Share…**, and read **What's inside** before you save
  ({doc}`../m_state_history/06_sharing_and_replay`).

---

## What a session file does NOT contain (important)

By design, sessions do not include:

- the full dataset (cell coordinates, full observation table, full expression
  matrix)
- sign-in tokens: Google Drive and OneDrive sign-ins and the GitHub sign-in for
  Community Annotation are kept in the tab's `sessionStorage` and never written
  to a session file (a session saved on a signed-in file records its link, and
  whoever opens it signs in with their own account)
- a key, signature or token a data link carried: a recipient without access of
  their own opens the link you sent them first, then loads the session
  ({ref}`link-access-keys`)
- the path of a local file, or the copy of it a Chrome or Edge browser keeps to
  reopen it (that stays in the browser)
- Community Annotation votes, comments and moderation state
- Benchmarking state

This reduces risk, but does not eliminate it (sessions can still contain meaningful derived information).

The history this browser keeps (**Load…**, {doc}`../m_state_history/03_saved_states_and_autosave`) holds the same information
as a session file, but it stays in this browser and is never uploaded; see
{doc}`../m_state_history/03_saved_states_and_autosave`.

---

## Trust model: sessions are treated as untrusted input

Cellucid loads `.cellucid-session` files with safety checks:

- the signature, the manifest, every tracker and every chunk entry are checked
  exactly; unknown, missing, duplicate or reordered chunks are refused
- every chunk is decompressed to exactly its declared size, within bounds, and
  its single gzip stream is checked before it is inflated
- every array is checked against its digest and every value against its hash,
  and every step tree is checked whole
- the complete dataset fingerprint must match, and the dataset must have every
  field the shown step uses
- a failure at any point, before or after the step shown appears, puts back
  what was shown before

This is meant to protect you from:
- corrupted files,
- accidental "wrong dataset" loads,
- extreme memory allocations from maliciously crafted gzip payloads.

A load is one all-or-nothing operation. Starting a newer load or pressing
**Cancel** stops the older one and dismisses its progress without presenting
the cancellation as a failure.

:::{important}
Even with safety checks, best practice is still:
- only open sessions from sources you trust.
:::

---

## Safe sharing rules (practical)

### 1) Assume sessions are shareable only as widely as the dataset is shareable

- Public dataset → sessions can be public (if they contain no extra sensitive labels).
- Private dataset → sessions should remain private.

### 2) Redact sensitive semantics from labels

Avoid:
- patient identifiers
- clinician notes
- internal project codes

Prefer:
- `donor_1`, `donor_2`
- `condition_A`, `condition_B`
- `cluster_7` → plus separate private mapping if needed

### 3) Check what the history holds before you send it

A history file holds the steps of the branches you tick, undone ones
included. In **Share…**, tick only what should leave your computer, choose
**Current state** when the way you got there does not matter, set **Times** and
**Dataset** to what you want to say, and open **What's inside**: it lists the
branches, cell sets, results, notes, conversations, agents and models, names
of people, times and dataset the file will hold.

### 4) Don't publish sessions by accident

If you are distributing a dataset file or a dataset catalog folder publicly:
- do not include session files unless you intend them to be opened by anyone.

### 5) For high-sensitivity environments, share encrypted archives

If your organization requires it:
- share the dataset file + sessions via encrypted storage or approved systems.

---

## Next steps

- {doc}`05_share_workflows_links_sessions_data_files` (recommended sharing patterns)
- {doc}`10_troubleshooting_sessions` (when a load is refused)
- {doc}`../o_accessibility_privacy_security/02_privacy_model` (everything Cellucid keeps, and where)
