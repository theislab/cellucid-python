# Security, privacy, and trust

**Audience:** everyone sharing sessions (especially for private or clinical data)  
**Time:** 15–35 minutes  
**What you’ll learn:**
- What information a session bundle can contain (and therefore leak)
- What session bundles intentionally do *not* contain (tokens, full datasets, etc.)
- How Cellucid treats sessions as untrusted input (and what that protects you from)
- Practical “safe sharing” rules for teams

---

## Treat session bundles like data artifacts

Even though `.cellucid-session` files do not contain the full dataset, they can still encode sensitive information.

Rule of thumb:

> If the dataset is sensitive, assume sessions derived from it are sensitive too.

---

## What a session bundle can contain (privacy perspective)

Depending on what you did in the UI, a session bundle may include:

- **Dataset identifiers**
  - the dataset id, counts and cell-order digest of the dataset it was saved on
  - where the dataset was opened from, for a sample, a data link (in its
    public form: never a key, signature or token the link carries) or a
    server (its address); for a file opened from your own computer, its file
    or folder name
- **Field names and category labels**
  - including any renames you applied
  - including deleted/purged field registry state
- **Highlight structure**
  - page names (e.g., `Patient A`, `Responder`, etc. — avoid private IDs)
  - group names and enabled/disabled state
- **Highlight memberships**
  - sets of cell indices belonging to groups (can be large)
- **User-defined fields**
  - derived field definitions and exact per-cell categorical codes
- **Analysis configuration and results**
  - which analyses were open (sidebar panel and floating windows), their settings
  - their results: differential expression, marker genes, signatures,
    correlations, dot plots and compositions, i.e. statistics derived from the
    data
- **Figure Export settings**
  - including the figure title and axis labels you typed
- **Where the dataset came from**
  - the sample id, data link in its public form (never a key, signature or
    token it carries: a recipient who has no access of their own opens the
    link you sent them first, then loads the session; see
    {ref}`link-access-keys`), or server address; for a local file, its file
    or folder name (never its path, nor the handle a Chromium browser keeps
    to reopen it, which stays in that browser)

Practical implication:
- Do not put patient identifiers or private sample barcodes into highlight/page names if you intend to share sessions.

---

## What a session bundle does NOT contain (important)

By design, sessions do not include:

- the full dataset (points arrays, full obs table, full expression matrix)
- sign-in tokens: Google Drive and OneDrive sign-ins and the GitHub sign-in for
  Community Annotation are kept in the tab's `sessionStorage` and never written
  to a session file (a session saved on a signed-in file records its link, and
  whoever opens it signs in with their own account)
- Community Annotation votes/comments/moderation state
- Benchmarking state

This reduces risk, but does not eliminate it (sessions can still contain meaningful derived information).

---

## Trust model: sessions are treated as untrusted input

Cellucid loads `.cellucid-session` bundles with safety checks:

- exact-validates framing, manifest roots, chunk profiles, ordering, and byte
  lengths before success
- bounds stored and decoded sizes and preflights the single-member
  gzip/DEFLATE structure before native decompression
- rejects unknown, missing, duplicate, aliased, or dishonest contributors
- requires the complete dataset fingerprint to match
- applies state in a transaction and rolls back earlier mutations if any later
  decode, contributor, commit, or final-refresh step fails

This is meant to protect you from:
- corrupted files,
- accidental “wrong dataset” restores,
- extreme memory allocations from maliciously crafted gzip payloads.

The strict reader exposes one all-or-nothing current operation. Starting a
newer restore or pressing **Cancel** aborts the older operation and dismisses
its progress without presenting the intentional cancellation as a product
failure.

:::{important}
Even with safety checks, best practice is still:
- only open sessions from sources you trust.
:::

---

## Safe sharing rules (practical)

### 1) Assume sessions are shareable only as widely as the dataset is shareable

- Public dataset → sessions can be public (if they contain no extra sensitive labels).
- Private dataset → sessions should remain private.

### 2) Redact sensitive semantics from UI labels

Avoid:
- patient identifiers
- clinician notes
- internal project codes

Prefer:
- `donor_1`, `donor_2`
- `condition_A`, `condition_B`
- `cluster_7` → plus separate private mapping if needed

### 3) Don’t publish sessions by accident

If you are distributing a dataset file or a dataset catalog folder publicly:
- do not include session bundles unless you intend them to be opened by anyone.

### 4) For high-sensitivity environments, share encrypted archives

If your organization requires it:
- share the dataset file + sessions via encrypted storage or approved systems.

---

## Next steps

- {doc}`05_share_workflows_links_sessions_data_files` (recommended sharing patterns)
- {doc}`10_troubleshooting_sessions` (when a restore is rejected)
