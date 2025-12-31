# Troubleshooting (sessions) — fast fixes

**Audience:** everyone (especially anyone sharing sessions across machines)  
**Time:** 15–60 minutes (depending on the failure mode)  
**What you’ll learn:**
- A symptom → diagnosis → fix map for session saving/restoring
- How to distinguish “expected exclusions” from real restore failures
- How to debug auto-restore (`state-snapshots.json`) problems using DevTools

---

## Fast fix map (start here)

| Symptom | Most likely cause | First thing to try | Deep dive |
|---|---|---|---|
| “Load State did nothing” | dataset mismatch or no dataset loaded | load dataset first; watch for mismatch warning | {doc}`07_versioning_compatibility_and_dataset_identity` |
| “File picker didn’t open” | browser/iframe restriction | try in a standalone tab; try another browser | “Load State button does nothing” below |
| “Session loaded but looks different” | dataset changed, missing fields, or lazy restore not finished | wait for progress; check Active filters + active field | “Partial restore” below |
| “Auto-restore didn’t happen” | missing/invalid `state-snapshots.json` | check Network for `state-snapshots.json` | {doc}`04_auto_restore_latest_from_dataset_exports` |
| “Restore says dataset mismatch” | different dataset id or load method | load the exact dataset (same source type + id) | {doc}`07_versioning_compatibility_and_dataset_identity` |
| “Save State is extremely slow / huge files” | huge highlights/caches | save milestone sessions; reduce groups/snapshots | {doc}`09_edge_cases` |

---

## Before you debug anything (2-minute checklist)

1) **Confirm a dataset is loaded**
   - Sessions don’t contain data; restoring without a dataset can’t work as expected.
2) **Look for a dataset mismatch warning**
   - If you see it, stop: you must fix identity first.
3) **Confirm eager restore finished**
   - Look for “Session loaded successfully” or “eager stage complete”.
4) **Wait for lazy restore if the session is large**
   - Large highlight memberships and caches can finish later.

If you still have a problem, use the sections below.

---

## Save State problems

### Symptom: “Save State” does nothing

**Likely causes (ordered)**
- The app is still initializing (session serializer not ready yet).
- A browser download restriction blocked the programmatic download.
- A JavaScript error occurred during session capture.

**How to confirm**
- Look for a notification:
  - “Failed to save session” indicates capture failed.
- Open DevTools → Console and look for:
  - `Failed to save state:` errors.

**Fix**
1) Wait until the dataset is fully loaded and the UI is responsive.
2) Try **Save State** again.
3) If it still fails, try another browser (Chrome/Firefox are usually most reliable for file downloads).
4) If you are inside an embedded context (iframe/notebook), open Cellucid in a standalone tab and try again.

**Prevention**
- Save after major interactions, not during heavy loading.

---

### Symptom: “Save State” is very slow or freezes the tab

**Likely causes**
- You have extremely large highlight memberships (huge groups).
- You have many snapshots/views.
- You have large session artifacts (user-defined codes, caches).

**How to confirm**
- Does the dataset have very large `n_cells` (hundreds of thousands to millions)?
- Do you have many highlight groups/pages?

**Fix**
1) Save a smaller “milestone” session:
   - remove non-essential snapshots,
   - reduce the number of huge highlight groups,
   - close unnecessary analysis windows.
2) Try saving again.

**Prevention**
- Use the “milestone session” approach in {doc}`06_collaboration_best_practices`.

---

## Load State problems (manual file)

### Symptom: “Load State” button does nothing / file picker doesn’t appear

**Likely causes**
- Browser blocked file picker due to embedded/managed restrictions.
- The click was not considered a user gesture (rare, but can happen in some embedded contexts).
- A JavaScript error prevented the handler from running.

**How to confirm**
- Try in a fresh tab with no other modal dialogs open.
- Check DevTools → Console for errors around session loading.

**Fix**
1) Try in a standalone browser tab (not inside Jupyter/iframe).
2) Try another browser.
3) If your organization uses strict browser policies, ask for File System Access / downloads permission.

---

### Symptom: “Load State” file picker opens but you can’t select your file

**Likely causes**
- The file does not end in `.cellucid-session`.
- The file is still zipped (`.zip`) or has an extra extension (`.cellucid-session.txt`).
- Your OS file picker is filtering for “supported types” and hiding the file.

**How to confirm**
- Check the filename suffix in your file browser.

**Fix**
1) Ensure the file ends exactly in `.cellucid-session`.
2) If it’s zipped, unzip first.
3) If your OS appended `.txt`, rename the file (remove `.txt`).

---

### Symptom: “Session loaded successfully” but nothing changes

This is almost always one of two things.

**Likely causes**
1) You loaded the session into the wrong dataset (dataset mismatch → most state skipped).
2) The session contained minimal changes (you saved very early / near defaults).

**How to confirm**
- Did you see a dataset mismatch warning?
- Check whether your current view already matched the expected camera/filters.

**Fix**
1) If you saw a dataset mismatch warning:
   - load the correct dataset first, then load the session again.
2) If no mismatch warning:
   - try loading a different session file (to confirm the load path works),
   - or ask the sender what visible change to expect (snapshots? highlights? filters?).

---

### Symptom: “Session loaded, but it looks different than the sender’s”

**Likely causes (ordered)**
- Dataset mismatch (explicit warning).
- Dataset content changed even though the id/fingerprint matches (cell order differences, export differences).
- Missing fields/genes in the dataset export (session references something that doesn’t exist).
- Lazy restore not finished (large highlights/caches).
- Different screen size caused panel/layout differences (expected).

**How to confirm**
1) Check for dataset mismatch warning.
2) Check **Active filters** and the **active field** (highest-signal state).
3) Wait for lazy restore notifications to finish.
4) Compare:
   - number of snapshots,
   - highlight page/group names,
   - whether gene expression is available (if the session depends on it).

**Fix**
- If mismatch: fix dataset identity first ({doc}`07_versioning_compatibility_and_dataset_identity`).
- If missing fields: load the export version used to create the session, or re-create the session.
- If lazy restore: wait; if it never completes, see “Restore is stuck” below.

---

### Symptom: “I see my highlight groups, but no points are highlighted”

**Likely causes**
- Highlight memberships are still restoring lazily (expected for large sessions).
- Filters currently hide the highlighted cells (visibility vs highlighting interaction).
- The highlight overlay is visually subtle in the current theme/point size.

**How to confirm**
- Wait for session restore progress to finish.
- Disable filters one-by-one in Active filters.
- Zoom in or increase point size to make highlight overlay more visible.

**Fix**
1) Wait for lazy restore to finish.
2) Temporarily disable filters to verify highlight memberships exist.
3) Adjust point size/background for visibility.

Related: {doc}`../f_highlighting_selection/index`, {doc}`../e_filtering/index`

---

### Symptom: “Restore is stuck / progress never completes”

**Likely causes**
- The session file is extremely large and the browser is CPU/memory constrained.
- A specific lazy chunk restore is expensive (huge highlights, huge caches).
- The restore was canceled (you loaded another session, navigated away, or the tab lost resources).

**How to confirm**
- Does the browser tab become sluggish or memory-heavy?
- Do you see warnings/errors in DevTools console about chunk restore failures?

**Fix**
1) Give it time (large sessions can take a while).
2) If it never finishes:
   - reload the page,
   - load the dataset first,
   - load the session again,
   - avoid interacting heavily until eager restore completes.
3) If it still fails, try a smaller session (or remove huge highlight groups and resave).

---

## Auto-restore problems (`state-snapshots.json`)

### Symptom: “I opened the dataset but no session auto-loaded”

**Likely causes (ordered)**
- `state-snapshots.json` is missing.
- `state-snapshots.json` exists but has no `.cellucid-session` entries.
- `state-snapshots.json` is served as HTML (SPA fallback) and fails JSON parse.
- The `.cellucid-session` file URL 404s or is blocked by CORS.

**How to confirm**
1) DevTools → Network:
   - find `state-snapshots.json`
   - confirm status `200` and JSON content
2) DevTools → Console:
   - look for “No session bundle auto-loaded…” or state-snapshots parse warnings

**Fix**
- Follow the debugging steps in {doc}`04_auto_restore_latest_from_dataset_exports`.

---

### Symptom: “Auto-restore finds the manifest but session fetch fails”

**Likely causes**
- The session filename in `state-snapshots.json` is wrong.
- The relative path resolves incorrectly.
- Hosting requires auth (sessions must be accessible to the browser).

**How to confirm**
- In Network tab, click the `.cellucid-session` request and inspect the resolved URL.

**Fix**
1) Correct the entry in `state-snapshots.json`.
2) Prefer relative paths that are correct relative to the manifest URL.

---

## Dataset mismatch warning problems

### Symptom: “Session dataset mismatch (…) Restoring only dataset-agnostic layout.”

**What it means**
- The session was saved under a different dataset identity than what you currently loaded.
- Cellucid is refusing to apply dataset-dependent state for safety.

**Fix**
1) Load the dataset the session was created from (same source type + dataset id).
2) Then load the session again.

If you’re collaborating:
- agree on a single dataset access method (GitHub vs local folder vs hosted URL),
- or regenerate the session under the recipient’s load method.

Deep dive: {doc}`07_versioning_compatibility_and_dataset_identity`

---

## “Is this a bug?” (how to report effectively)

If you think you found a session bug, collect these before reporting:

1) What dataset source type you used (local folder / GitHub / server) and the dataset id.
2) Whether restore showed a dataset mismatch warning.
3) Which parts failed:
   - camera? filters? highlights? snapshots? analysis windows?
4) DevTools console logs around “SessionSerializer” warnings/errors.
5) The session file (if shareable) and the export folder version.

---

## Next steps

- {doc}`09_edge_cases` (expected weirdness and how to avoid it)
- {doc}`12_reference` (format notes and supported `state-snapshots.json` schemas)
