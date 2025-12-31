# Data Flow: R → Export → Browser

**Audience:** everyone  
**Time:** 10 minutes  
**Goal:** choose the right way to load your exported data into Cellucid.

The Cellucid viewing workflow starts with **files on disk**.

R produces the files. The **browser** reads them.

## Three ways to view an export folder

### Option 1 — Load locally via the browser file picker (simplest)

This is the “no server” workflow:

1) Export a folder with `cellucid_prepare(out_dir = "...")`.
2) Open the Cellucid web app.
3) Use the dataset loader to select the export folder on disk.

Pros:
- No server to run
- Works well for small/medium datasets

Cons / limitations:
- Browser file access UX varies by OS
- Some browsers impose performance limits on very large directory trees

Start here: {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`

### Option 2 — Serve the folder locally (good for large datasets / stability)

You run a local static file server so the web app loads over `http://localhost:...`.

Pros:
- More predictable browser behavior
- Easier to share links to your own machine (same LAN)
- Often faster for large datasets because the browser can stream requests cleanly

Cons:
- You must run a server process
- You need to think about LAN exposure (security)

Start here: {doc}`../d_viewing_loading/02_host_exports_for_sharing`

### Option 3 — Host the folder remotely (for sharing and collaboration)

If you can host static files (object storage, CDN, GitHub Pages, etc.), you can share exports with collaborators without sending huge archives.

Pros:
- Best for collaboration
- Enables stable share links

Cons:
- Requires hosting and CORS setup
- Be careful with privacy (you are publishing data)

Start here: {doc}`../d_viewing_loading/02_host_exports_for_sharing`

## Where do `cellucid-python` and `cellucid-annotation` fit?

Even if you export from R, you might still use the other repos:

- **`cellucid-python`** (optional):
  - provides server mode + Jupyter embedding + hooks
  - can serve an export folder from Python if you prefer that toolchain
- **`cellucid-annotation`** (optional):
  - provides workflows for community annotation and voting
  - expects datasets to have stable identity (see {doc}`03_dataset_identity_and_reproducibility`)

## Privacy and security (do not skip)

If your dataset contains sensitive information:

- **Local file picker** is safest (data stays on your machine).
- **Local server** can expose data to your network if you bind to non-localhost.
- **Remote hosting** is publishing; treat it like making data public unless access controls exist.

See also the web app security model: {doc}`../../web_app/o_accessibility_privacy_security/index`.

## Next steps

- Export a folder: {doc}`../c_data_preparation_api/index`
- Open it in the browser: {doc}`../d_viewing_loading/01_open_exports_in_cellucid_web_app`
