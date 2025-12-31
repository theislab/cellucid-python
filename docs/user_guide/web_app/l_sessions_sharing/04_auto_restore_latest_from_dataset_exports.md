# Auto-restore latest session (dataset exports)

**Audience:** demo builders, dataset publishers, and power users  
**Time:** 15–35 minutes  
**What you’ll learn:**
- What “auto-restore latest session” means (and when it runs)
- How `state-snapshots.json` is interpreted
- How to package session bundles inside an export folder so the app opens “already configured”
- How to debug the common “auto-restore did nothing” failures

**Prerequisites:**
- A dataset that is loaded from a directory/server/GitHub path (i.e. it has a “base URL”)
- Ability to add files into the dataset export folder

---

## What auto-restore is (and why it exists)

Cellucid supports a legacy-friendly workflow:

> “When I open this dataset, automatically restore the latest saved state.”

This is useful for:
- demos (“open the dataset and it’s already in the right view”),
- tutorials (“everyone starts from the same state”),
- sharing a default configuration (“this is the intended starting point”).

---

## When auto-restore runs

Auto-restore runs **on app startup**:

1) Cellucid loads the initial dataset.
2) Then it tries to find a session bundle list in the dataset exports directory.
3) If it finds one, it loads the **last** `.cellucid-session` entry.

:::{note}
Auto-restore is not a general “session browser”.
It always picks the last entry in the manifest.
If you want to choose among multiple sessions, use {doc}`03_save_restore_ux` and load a specific file manually.
:::

---

## The required file: `state-snapshots.json`

Auto-restore looks for a JSON file named:

`state-snapshots.json`

in the dataset exports directory (i.e., “next to” files like `obs_manifest.json` and `points_3d.bin.gz`).

That JSON file must contain one or more entries that end in:

`.cellucid-session`

Cellucid filters the entries and restores the **last** matching entry.

---

## Supported `state-snapshots.json` shapes (practical spec)

Cellucid accepts two dev-phase shapes:

### Option A (recommended): object with `states`

```json
{
  "states": [
    "cellucid-session-2025-12-30T02-00-14.cellucid-session",
    "cellucid-session-2025-12-30T03-12-01.cellucid-session"
  ]
}
```

### Option B: bare array (also accepted)

```json
[
  "cellucid-session-2025-12-30T02-00-14.cellucid-session",
  "cellucid-session-2025-12-30T03-12-01.cellucid-session"
]
```

Notes:
- Entries can be relative paths or absolute URLs.
- Query strings and hashes are ignored when checking the filename suffix.
- The match is case-insensitive (e.g., `.CELLUCID-SESSION` is still accepted).

---

## Recommended export folder layout (with sessions)

At minimum, put the session file(s) in the same folder as `state-snapshots.json`:

```
<dataset_export>/
  obs_manifest.json
  points_3d.bin.gz
  ...
  state-snapshots.json
  cellucid-session-2025-12-30T02-00-14.cellucid-session
```

If you want to group them (still works), use relative paths in `state-snapshots.json`:

```
<dataset_export>/
  state-snapshots.json
  sessions/
    v1.cellucid-session
    v2.cellucid-session
```

---

## Debugging: when auto-restore “does nothing”

This is common and usually fixable quickly.

### Step 1: check the browser console

Open DevTools → Console and look for messages like:

- `[Main] No session bundle auto-loaded (state-snapshots.json missing/empty or no .cellucid-session entries).`
- `[SessionSerializer] state-snapshots.json not available: ...`
- `[SessionSerializer] Failed to parse state-snapshots.json (expected JSON): ...`

### Step 2: check the network tab

In DevTools → Network:

1) Find the request for `state-snapshots.json`.
2) Confirm it returned:
   - HTTP `200`, and
   - real JSON (not an HTML “app shell” fallback page).

If you’re hosting exports behind an SPA router, it’s very easy for `state-snapshots.json` to return HTML.

### Step 3: verify the session file URL resolves

Once `state-snapshots.json` loads, Cellucid resolves the session file URL relative to the manifest URL (so redirects work).

Confirm:
- the session file request returns `200`,
- the response body is not empty.

### Step 4: confirm the filename really ends in `.cellucid-session`

Cellucid filters entries by filename suffix.

Common mistakes:
- `my_session.cellucid-session.json` (ends in `.json`, will be ignored)
- `my_session.cellucid-session?download=1` (fine)
- `my_session.cellucid-session#v1` (fine)
- `my_session.cellucid-session.tmp` (ignored)

### Step 5: if you see “session file truncated”, suspect hosting/proxy behavior

If you see errors like:
`Invalid chunk length ... (session file truncated?)`

Common causes:
- partial uploads,
- proxy timeouts on large files,
- servers applying unexpected content encoding.

Try:
- re-uploading the file,
- serving the file with correct static hosting settings,
- or testing with a smaller session to validate your pipeline first.

---

## Safety note: auto-loading is still “opening a file”

Auto-restore means “a dataset can cause a session file to be fetched and parsed”.

Sessions are treated as untrusted input (bounds checks + size guards), but best practice is still:
- only host session bundles you intend others to load,
- treat session bundles as potentially sensitive artifacts (see {doc}`08_security_privacy_and_trust`).

---

## Next steps

- {doc}`05_share_workflows_links_bundles_exports` (how to share a dataset+session pair)
- {doc}`10_troubleshooting_sessions` (more error patterns and fixes)
