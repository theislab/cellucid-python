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
  - source type and dataset id used to load the dataset
  - (for local folder loads) a display name for the folder
- **Field names and category labels**
  - including any renames you applied
  - including deleted/purged field registry state
- **Highlight structure**
  - page names (e.g., `Patient A`, `Responder`, etc. — avoid private IDs)
  - group names and enabled/disabled state
- **Highlight memberships**
  - sets of cell indices belonging to groups (can be large)
- **User-defined fields**
  - derived field definitions and (sometimes) per-cell codes
- **Analysis window configuration**
  - which analysis modes were open and their settings
- **Some caches/artifacts**
  - sessions can include “accelerator” caches that may encode partial data-derived information

Practical implication:
- Do not put patient identifiers or private sample barcodes into highlight/page names if you intend to share sessions.

---

## What a session bundle does NOT contain (important)

By design, sessions do not include:

- the full dataset (points arrays, full obs table, full expression matrix)
- dataset loading secrets/credentials (there aren’t any in normal Cellucid usage)
- GitHub OAuth tokens for Community Annotation (stored in browser session state, not in session bundles)
- Community Annotation votes/comments/moderation state
- Figure Export UI state
- Benchmarking state

This reduces risk, but does not eliminate it (sessions can still contain meaningful derived information).

---

## Trust model: sessions are treated as untrusted input

Cellucid loads `.cellucid-session` bundles with safety checks:

- validates manifest structure
- enforces size limits (including decompression caps)
- skips unknown chunk contributors
- isolates chunk restore failures so one bad chunk doesn’t brick the app
- skips dataset-dependent chunks on dataset mismatch

This is meant to protect you from:
- corrupted files,
- accidental “wrong dataset” restores,
- extreme memory allocations from maliciously crafted gzip payloads.

:::{important}
Even with safety checks, best practice is still:
- only open sessions from sources you trust, and
- avoid auto-loading sessions from untrusted datasets.
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

If you are distributing an export folder publicly:
- do not include session bundles unless you intend them to be opened by anyone,
- be careful with `state-snapshots.json` (auto-restore makes sessions “execute on open”).

### 4) For high-sensitivity environments, share encrypted archives

If your organization requires it:
- share the export folder + sessions via encrypted storage or approved systems.

---

## Next steps

- {doc}`05_share_workflows_links_bundles_exports` (recommended sharing patterns)
- {doc}`10_troubleshooting_sessions` (when restores fail or warn)
