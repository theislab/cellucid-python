# Collaboration best practices

**Audience:** teams collaborating across wet lab + computational + PI/reviewers  
**Time:** 20–45 minutes  
**What you’ll learn:**
- Practical conventions that make session files usable months later
- How to use sessions as “review packets” and “milestone checkpoints”
- How to avoid collaboration failure modes (dataset mismatch, missing context, privacy leaks)
- When to use sessions vs Community Annotation vs exported figures

**Prerequisites:**
- You can save and load sessions (see {doc}`03_save_restore_ux`)

---

## The goal: a collaborator should never have to guess

The best session sharing feels like:

> “Open this dataset, load this session, and you will see exactly the thing I’m asking about.”

The worst session sharing feels like:

> “I sent you a file. It might work. If it doesn’t, we’ll spend an hour screensharing.”

The rest of this page is how you get the first outcome reliably.

---

## Best practice #1: always ship the “context triple”

Whenever you share a session, include three pieces of context:

1) **Dataset identity**
   - “Which dataset export is this?”
   - Provide a stable dataset label (and ideally a version/commit/date).
2) **Session intent**
   - “What should the recipient look at?”
   - Provide a 1–3 sentence question or hypothesis.
3) **Expected success state**
   - “What should it look like when it worked?”
   - Example: “you should see 4 snapshots and one highlight page named ‘Treated’.”

This can be a text message, but it’s better as a small `README.md` alongside the session files.

---

## Best practice #2: use a directory structure that scales

Recommended folder layout for a project:

```
project/
  exports/
    dataset_v1/
    dataset_v2/
  sessions/
    dataset_v1/
      2025-12-30__qc__v1.cellucid-session
      2026-01-02__fig1_candidate__v3.cellucid-session
    dataset_v2/
      ...
  figures/
    fig1/
      fig1_candidate_v3.svg
  notes/
    sessions_log.md
```

Why this works:
- sessions don’t contain data, so you keep them paired by dataset version;
- session filenames stay short but meaningful;
- you can quickly answer “which session goes with which export?”

---

## Best practice #3: treat sessions like “milestones”, not “autosave”

Sessions can store large artifacts (highlight memberships, caches).
Saving constantly produces:
- huge files,
- noisy timelines,
- confusion about which one matters.

Recommended milestones:
- after QC filtering
- after first-pass labeling/highlights
- after final highlight pages
- before figure export
- “review packet” sessions for collaborators

---

## Best practice #4: encode meaning in highlight names (but avoid private IDs)

Highlights and pages often become the “story” of your analysis.

Good names:
- `Control`
- `Treated`
- `Cycling`
- `Potential doublets`
- `Suspected contaminant cluster`

Names to avoid in shared artifacts:
- patient identifiers
- sample barcodes
- internal project codes that leak sensitive information

If you need privacy-safe labels, use:
- `sample_A`, `sample_B`
- `donor_1`, `donor_2` (if allowed)
- `condition_1`, `condition_2`

See {doc}`08_security_privacy_and_trust`.

---

## Best practice #5: send “review packets” (sessions + one screenshot + one question)

For PI review, wet-lab review, or async collaboration, a good packet is:

- a `.cellucid-session` file
- one screenshot (or exported figure) that shows the “expected view”
- one question

Example:

> “Load `pbmc_demo__qc__v2.cellucid-session`. In snapshot 3, highlight group ‘Cycling’ looks split across UMAP; is this likely cell cycle vs batch? Please inspect.”

This reduces the chance that:
- the recipient loads the wrong dataset,
- the recipient doesn’t know what you want them to examine,
- the recipient assumes success when the session partially restored.

---

## Best practice #6: for multi-user labeling, use Community Annotation (not sessions)

Sessions are great for **sharing your personal state**, but they are not a multi-user synchronization mechanism.

If many people need to vote on labels and converge to consensus, use:
- {doc}`../j_community_annotation/index`

Why:
- Community Annotation is offline-first and conflict-free (each user writes their own file).
- Sessions do not store annotation repo auth or votes.

---

## Best practice #7: validate your share before you send it

Before sending a session to someone else:

1) Save the session.
2) Open a fresh tab (or a different browser profile).
3) Load the dataset again.
4) Load the session.
5) Confirm there is no dataset mismatch warning.

If you can’t reproduce your own session in a clean tab, your collaborator won’t be able to either.

---

## Best practice #8: choose the right “final artifact”

Sessions are a reproducibility tool, but for “final outputs” you often want:
- exported figures (SVG/PNG),
- exported tables (analysis outputs),
- a stable dataset export.

Rule of thumb:
- Use sessions to preserve *interactive context*.
- Use figures/tables to preserve *final results*.

Related: {doc}`../k_figure_export/index`, {doc}`../h_analysis/index`

---

## Next steps

- {doc}`05_share_workflows_links_bundles_exports` (share workflow patterns)
- {doc}`07_versioning_compatibility_and_dataset_identity` (why sessions sometimes don’t transfer across machines)
- {doc}`10_troubleshooting_sessions` (what to do when collaboration fails anyway)
