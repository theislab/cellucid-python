# Save/restore UX (manual save + restore)

**Audience:** everyone (wet lab + computational)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- Where the session controls live in the UI
- The exact steps to Save State and Load State
- How to confirm a restore succeeded (and what “partial restore” looks like)
- Best practices for naming, organizing, and sharing session files

**Prerequisites:**
- A dataset loaded in the web app

---

## Where the session controls live

In the left sidebar:

1) Open **User data** (the accordion that also contains dataset loading controls).
2) Find the block labeled **Session state**.
3) Use:
   - **Save State** (downloads a `.cellucid-session` file)
   - **Load State** (opens a file picker and restores a `.cellucid-session` file)

---

## Fast path (wet lab / non-technical): “save my work”

### Save State

1) Click **Save State**.
2) Wait for a “Session saved successfully” notification.
3) Find the downloaded file in your browser’s downloads.
   - It will end in `.cellucid-session`.
   - The default name includes a timestamp (so files don’t overwrite each other).

### Load State

1) Load the same dataset first.
2) Click **Load State**.
3) Choose the `.cellucid-session` file.
4) Wait for “Session loaded successfully”.

If nothing changes, don’t panic: see “How to confirm it worked” below.

---

## Practical path (computational): what happens under the hood

### Save State is a file download

When you click **Save State**:
- Cellucid captures state from multiple subsystems (camera, filters, highlights, etc.).
- It writes a single-file bundle and triggers a browser download.

There is no server upload by default.

### Load State is progressive restore (eager then lazy)

When you click **Load State** and pick a file:

1) **Eager restore** runs first (fast):
   - restores camera/layout/active field/filters quickly.
2) **Lazy restore** continues in the background (can be slow for large sessions):
   - restores large highlight memberships, some caches, user-defined codes, etc.

Practical implication:
- If you see the correct view but highlights “appear later”, that is expected for big sessions.

---

## How to confirm a restore worked (quick checklist)

After loading:

1) Look for a session notification:
   - “Session loaded successfully” or “Session restored (eager stage complete)”
2) Check you did **not** get a dataset mismatch warning.
   - If you did, see {doc}`07_versioning_compatibility_and_dataset_identity`.
3) Confirm the high-signal state:
   - active field (Coloring & Filtering section)
   - Active filters panel (should match what you remember)
   - snapshot layout (grid vs single)
4) If the session is large, wait for lazy restore:
   - some highlight memberships and caches can finish later.

If anything is off, jump to {doc}`10_troubleshooting_sessions`.

---

## Best practices (these save you hours later)

### 1) Treat the dataset + session as a pair

Because sessions do not contain the dataset:
- store the `.cellucid-session` file next to the dataset export folder,
- or store a clear pointer to where the dataset is hosted (URL/GitHub path).

### 2) Use a naming convention that encodes intent

Recommended naming format:

`<dataset>__<goal>__<who>__<YYYY-MM-DD>__v<N>.cellucid-session`

Examples:
- `pbmc_demo__qc-filtered__kelly__2025-12-30__v1.cellucid-session`
- `suo__figure1_candidate__team__2025-12-30__v3.cellucid-session`

Why:
- session files are binary; the filename is your metadata unless you add a separate README.

### 3) Save “milestones”, not every minute

Sessions can become large when they include:
- many highlight groups with huge memberships,
- many user-defined categorical fields,
- analysis caches.

Save milestones like:
- “after QC filters + first labeling”
- “after final highlight pages”
- “before figure export settings”

### 4) Save before you load someone else’s session

Loading is destructive to the current UI state.

Workflow:
1) Save your current session.
2) Load the other session.
3) If needed, switch back by loading your saved session.

---

## Screenshots (placeholders)

<!-- SCREENSHOT PLACEHOLDER
ID: sessions-session-controls-ui-save-load
Suggested filename: sessions_sharing/03_save-load-buttons.png
Where it appears: Sessions → Where the session controls live
Capture:
  - Show User data accordion with Session state block visible
Alt text:
  - Save State and Load State buttons in the Session state block.
Caption:
  - Session controls live under User data → Session state.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Save State / Load State buttons.
:width: 100%

Session controls live under User data → Session state.
```

<!-- SCREENSHOT PLACEHOLDER
ID: sessions-load-file-picker
Suggested filename: sessions_sharing/04_load-file-picker.png
Where it appears: Sessions → Load State
Capture:
  - Trigger the file picker (any OS is fine)
  - If your OS reveals private paths/usernames, blur or crop
Alt text:
  - File picker dialog filtering for .cellucid-session files.
Caption:
  - Load State opens a file picker restricted to `.cellucid-session` bundles.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Load State file picker.
:width: 100%

Load State opens a file picker restricted to `.cellucid-session` bundles.
```

---

## Next steps

- {doc}`04_auto_restore_latest_from_dataset_exports` (auto-loading sessions from `state-snapshots.json`)
- {doc}`05_share_workflows_links_bundles_exports` (how to send sessions to humans safely)
