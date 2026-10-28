# Troubleshooting (data loading)

**Audience:** everyone (quick checks for beginners; deeper diagnostics for computational users)  
**Time:** 5–20 minutes (depending on the symptom)  
**What you’ll learn:** how to quickly localize whether a failure is (a) browser/UI, (b) network/CORS, or (c) your data format  
**Prerequisites:** none

---

## First: identify your loading path (so you debug the right thing)

Before deep-diving, answer:

1) Are you loading:
   - a browser-selected **`.h5ad`**, `.zarr.zip`, Zarr folder or another file
     or folder Cellucid reads (**Open…**),
   - a Python-served path (the exact printed **Viewer URL**, `/?anndata=true`),
   - a notebook or R viewer (`/?viewer=true&viewerId=…`, opened by
     `show_anndata()` or `cellucid_view()`),
   - a **remote server** (`?remote=...`, only when the viewer origin can fetch that scheme),
   - a **GitHub dataset catalog** (`?data=gh:...`),
   - or another **data link** (`?data=...`: a catalog folder on a static host
     or CDN, CELLxGENE, GEO, or an `.h5ad` / Zarr address)?
2) Is it **browser-only** (no Python running) or **server-backed** (Python process running)?

If you’re unsure, start with {doc}`01_loading_options_overview`.

---

## 2-minute triage checklist (do this before anything else)

### Step 1 — Confirm the web app itself works

- Load a built-in demo dataset (if available).
- If demos fail too, this is likely an environment/browser issue:
  - see system requirements: {doc}`../a_orientation/02_system_requirements`

### Step 2 — Confirm you’re not fighting the wrong workflow

Common mismatch:
- Trying to open a file from your computer that lives on a remote machine.

A large `.h5ad` or Zarr store opens in the browser at any size: it is read in
place, as each view needs it. If the file is on a remote machine, serve it from
there instead of copying it: {doc}`04_server_tutorial`. `cellucid serve` hands
the file to the same reader in the browser.

### Step 3 — Capture the “shape” of the failure

- “Nothing happens when I click” → file picker permissions / browser issues
- “It connects but then fails to fetch something” → network/CORS/URL issue
- “It loads but is blank / missing fields” → data missing embeddings/fields

---

## Read the message before anything else

Cellucid's loading failures are written to be diagnostic, and the wording tells
you which of the three layers failed. Connection failures for the remote-server
and GitHub paths take the shape:

> *`<subject>` `<what happened>`: `<cause>` `<what to do>`*

and the middle clause is the one that identifies the layer:

| The message says | It means | Look at |
|---|---|---|
| nothing is published at that address (the server answered "not found") | The address resolved, but there is nothing there | Your path, folder, or branch |
| the server refused access | Reached, but not public | Repository visibility |
| the server rejected the request | Reached, but the request was malformed | Your address |
| the server reported a problem of its own | Reached; the far side failed | Wait, then retry |
| what is published there is larger than Cellucid opens in a browser | Reached, readable, valid — and more than this browser could allocate | The file size and the machine's memory. Serve it instead ({doc}`04_server_tutorial`), open it on a machine with more memory, or publish a smaller file. The ceilings are tabulated in {doc}`../q_troubleshooting_index/index` |
| what is published there is not in a format Cellucid can read | Reached and readable, but not an AnnData file Cellucid reads | Write it again with anndata (or `cellucid.prepare()`) |
| nothing readable came back — the connection may have failed, or that address may not hold Cellucid data | Ambiguous: no usable response arrived | Your network *and* your address |

An address the **Link or server** bar cannot use never reaches this table. A
GitHub path that is not `owner/repo/folder`, `owner/repo@branch/folder`, a
file path or a GitHub URL, and a port that is not a number from 1 to 65535
(`localhost:99999`), are refused before any request is made, with one message
about the address itself (for example `Invalid GitHub path. Use format:
owner/repo/folder (or owner/repo@branch/folder) for a folder holding
datasets.json, owner/repo/path/file.h5ad, or paste a GitHub URL.`). An empty
bar asks for an address the same way. A `?remote=` address in the page URL
must be one absolute `http(s)://` base without a trailing slash, a query or a
fragment.

:::{warning}
Read the size row and the format row as opposites. *Larger than Cellucid opens
in a browser* means the file is perfectly valid and merely too big, so
writing it again produces the same bytes and meets the same wall. That row's fix
is a different loading path, not a different file.
:::

There is one message that does **not** carry a cause clause, and its punctuation
is the tell — a full stop where the others have a colon:

> *`<subject>` `<what happened>`. Try again; if it keeps failing, check the
> address and your connection.*

That is what Cellucid says when it could not classify the failure at all, and
naming a cause anyway would be a guess. It is rare, because the ordinary
outcomes above are all classified — a connection you cancel yourself, and one
that runs past its timeout, both read as *nothing readable came back*. What is
left is a transfer the browser abandoned from underneath the request, or a
defect inside Cellucid.

Every failure notification keeps the raw error text folded under its
**Details** (and in the browser console), never in the sentence. A dataset you
chose from the list that does not open is said once, as *Could not open
"`<name>`".*, while the explanation and **Try again** stay next to the
selector. Quote the **Details** if you report the problem.

Local file failures use neither shape. They quote the specific defect instead —
for example `The selected file is not a valid HDF5 file, so it cannot be an AnnData (.h5ad) file…` or
`Cellucid found nothing to plot. It needs an obsm array with one row per cell…`.

---

## Symptom → diagnosis → fix (common issues)

### Symptom: “I chose Open… › File… or Folder…, but nothing happens”

**Likely causes (ordered)**
1) **You cancelled the file chooser.** This is by far the most common cause and
   it is deliberately silent: Cellucid shows no error, starts no spinner, and
   leaves the current dataset in place.
2) You denied the file or directory permission prompt.
3) The app is embedded in a context that blocks file access.
4) A managed-browser policy disabled local file access.
5) The selected object does not match the menu item: a folder (a `.zarr`
   store, a Cell Ranger `outs/`, a MatrixMarket folder) for **Open… ›
   Folder…**, and one file (`.h5ad`, `.zarr.zip`/`.zip`, `.rds`, `.loom`,
   `.h5`, `.h5mu`) for **Open… › File…**.

**How to confirm**
- Click the same button again and confirm the dialog. Silence after a *confirmed*
  selection is a real problem; silence after a dismissed one is expected.
- Open the standalone Cellucid page in a supported current desktop browser.
- Use a small known-current fixture with the matching control.

**Fix**
- Allow the picker permission and use the standalone page.
- If organizational policy forbids local selection, choose server mode
  explicitly: {doc}`04_server_tutorial`.

**Prevention**
- For workshops, validate **File…** and **Folder…** with a small fixture on
  the managed browser image before distributing data.

---

### Symptom: “It loads forever / spinner never ends”

:::{note}
A direct `.h5ad` has no size limit and is read from disk as each view needs it,
so a large file opens in seconds. What can take longer is the first gene of a
cell-major (CSR) matrix while it is reorganised by gene in the background; a
progress notice shows how far that is.
:::

No request waits for ever: a request the host answers as busy or failing
(HTTP 429, 503, ...), whose connection fails, or that stops receiving bytes is
asked again a few times, each repeat logged in the console as a
`[cellucid:transfer]` warning naming its URL, and then reported. When a
dataset opens at startup and nothing at all arrives for 60 s, Cellucid stops
waiting and starts without it: the panel over the empty viewer says the data
could not be opened, and its **Details** read `Startup stalled while opening
the dataset ...`, naming the reads that were open and the last request that
failed.

Any opening, at startup or later (a file, a folder, a link, a sample), that
has made no progress for 60 s says what it waits on, without stopping: its
notification reads `Nothing has changed for 1 min 0 s; waiting for` the reads
still open, with what has arrived of each (or `while` the step it is in, when
that is not a read), and the console logs `[cellucid:opening] Opening <name>
stalled while <step>: ...`. Cancel still stops it, and it goes on if what it
waits for answers. Include that console line when you report a hang.

**Likely causes (ordered)**
1) A very large cell-major (CSR) `.h5ad` opened directly in the browser while
   its matrix is still being reorganised by gene (watch the progress notice).
2) The data host is down, or a content blocker blocks some of its requests
   (the console's `[cellucid:transfer]` warnings name them).
3) Remote server URL is wrong or unreachable.
4) GitHub raw fetch blocked (corporate firewall/content blocker).

**How to confirm**
- Try loading a demo dataset.
- If using server mode, open the health endpoint:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

- If using GitHub, open the raw catalog URL in a new tab after substituting
  your exact owner, repository, branch, and path:

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json
  ```

**Fix**
- Large cell-major `.h5ad` → wait for the reorganisation notice to finish,
  or save the matrix column-major (CSC), which needs none. This is the same
  whether the file is opened in the browser or served
  ({doc}`04_server_tutorial`).
- Server mode → ensure URL/port matches the terminal output.
- GitHub → try a different network, disable blockers, or switch to server mode.

---

### Symptom: “It says Cellucid found nothing to plot”

An `.h5ad` or Zarr store without a UMAP — selected in the browser, opened from
a link, or served from Python — is not refused: it opens on a best guess (see
{doc}`13_data_mapping`). Only a file with no drawable `obsm` array at all is
refused, with this message:

```text
Cellucid found nothing to plot. It needs an obsm array with one row per cell
(<n_cells>) and at least one numeric column. obsm holds: <each entry>. Add a
2-D embedding to obsm (for example with sc.tl.umap in Scanpy), then open the
data again.
```

The list names every `obsm` entry with its shape and, for a dense array, why it
was not used: `X_df (table)`, `X_sparse (sparse 10,000 × 5)`,
`spatial (10,000 × 2, NaN at cell 17)`. A file with an empty `obsm` reads
`obsm holds: nothing.`

**Likely causes**
- AnnData has no dense numeric `obsm` array with one row per cell, or every
  such array holds a NaN, an infinity, a value too large for 32-bit floats, or
  the same coordinates for every cell.

**How to confirm**
- AnnData (Python):

  ```python
  print(adata.obsm.keys())
  ```


**Fix**
- Compute an embedding (for example `sc.tl.umap(adata)`) or fix the array the message names.
- If the file opens but shows the wrong embedding, turn the one you want on
  (and make it the default) under **Session › Data mapping**
  ({doc}`13_data_mapping`), or show it in a view with **View N shows**.

---

(no-cell-embedding)=
### Symptom: “… has no cell embedding, and Cellucid draws cells at embedding coordinates”

A count matrix with nothing to draw it with: a 10x Cell Ranger `.h5` opened
alone, or a Loom file without embedding attributes (a velocyto output). The
message names the file and what to do:

```text
filtered_feature_bc_matrix.h5 is a 10x Cell Ranger count matrix that has no
cell embedding, and Cellucid draws cells at embedding coordinates. Cell Ranger
writes UMAP and t-SNE projections beside it: open the whole outs/ folder
(Open… › Folder…), which holds this file and analysis/. Or compute one and open
the result: in Python, sc.pp.normalize_total(adata); sc.pp.log1p(adata);
sc.pp.pca(adata); sc.pp.neighbors(adata); sc.tl.umap(adata), then
adata.write_h5ad(path); in R, Seurat::NormalizeData(), FindVariableFeatures(),
ScaleData(), RunPCA() and RunUMAP(), then cellucid::cellucid_prepare(path, data = obj).
```

Opened from a link, the middle sentence says what a link allows instead: the
`outs/` folder "a link to one file cannot bring along: download the whole
outs/ folder and open it from this computer (Open… › Folder…)".

**Fix**
- Cell Ranger: open the `outs/` folder with **Open… › Folder…** (or drop it),
  so its `analysis/` projections are read with the matrix
  ({doc}`18_open_loom_and_10x`).
- A MatrixMarket folder: put a table of the cells with embedding columns
  (`UMAP_1`, `UMAP_2`) beside the matrix, or Cell Ranger's `analysis/`
  folder ({doc}`17_open_r_objects`).
- A Seurat object or SingleCellExperiment (`.rds`): run `RunPCA()` and
  `RunUMAP()` (or `scater::runUMAP()`), then `saveRDS()` again.
- Otherwise compute an embedding with the lines the message gives, and open
  the file they write.

(recognised-but-refused)=
### Symptom: “… is a Cell Ranger molecule table” or “… a MuData store whose modalities share features”

Cellucid recognised the file but does not show it:

- **`molecule_info.h5`** is a table of molecules, not a count matrix. Open
  `filtered_feature_bc_matrix.h5` from the same `outs/` folder, or the folder.
- **A MuData whose modalities share features, not cells** (`axis` 1 or −1):
  Cellucid shows cells. Open one modality as AnnData
  (`mdata.mod["rna"].write_h5ad(path)`).
- **A 10x matrix whose features of one type are stored in more than one
  place**: read it with `sc.read_10x_h5(path, gex_only=False)`, write it as
  `.h5ad` and open that.
- **A Cell Ranger 2 file whose genomes have different cells**: read one genome
  (`sc.read_10x_h5(path, genome="…")`), write it as `.h5ad` and open that.
- **A Loom file whose attributes do not match its matrix**: write it again with
  loompy.

(r-object-refused)=
### Symptom: an `.rds`, `.rda` or MatrixMarket folder is refused

The message names the object or file and the R line that fixes it
({doc}`17_open_r_objects` lists every one):

- **"… keeps the RNA assay's counts outside the file"**: a BPCells or
  HDF5Array matrix. In R, `seu[["RNA"]]$counts <- as(seu[["RNA"]]$counts,
  "dgCMatrix")` (and the same for `data`; for a SingleCellExperiment,
  `assay(sce, "counts") <- as(assay(sce, "counts"), "dgCMatrix")`), then
  `saveRDS()`.
- **"… has no layers in the file (SaveSeuratRds() moved them …)"**: load it
  with `LoadSeuratRds()` and save the object with its layers in it.
- **"This is a Seurat v2 object"**: `UpdateSeuratObject()`, then `saveRDS()`.
- **"This is an h5Seurat file"**: `SeuratDisk::LoadH5Seurat()`, then `saveRDS()`.
- **"… holds 2 single-cell objects"** (opened where no one can be asked, such
  as a server's dataset list): open the file in the viewer and choose one, or
  save the object to open on its own with `saveRDS()`.
- **"… compressed with zstd"**, **"saved as text"**, **"in R's native binary
  form"**: save it again with `saveRDS(obj, file)` and its defaults.
- **"This browser's storage has room for …"** or **"This browser's storage
  for this site is full …"**: a compressed file is read into
  the browser's storage first, and the files kept from earlier openings that
  no tab has open were already removed to make room (least recently opened
  first). In a private window the storage is small: open it in a normal
  window. Otherwise close other tabs that have a file read whole open, free
  space on the disk, or write the object once as a Cellucid-ready AnnData with
  `cellucid_prepare()` (R) or `cellucid.prepare()` (Python), which is read in
  place.
- **"This MatrixMarket matrix is symmetric"** (or complex, or a dense array):
  it is not a count matrix.
- **"… has no `…_barcodes.tsv.gz` beside it"** (or no features or genes
  list): the matrix's lists are missing or named with another prefix. Put
  the sample's three files side by side with one shared prefix
  (`GSM123_matrix.mtx.gz`, `GSM123_features.tsv.gz`,
  `GSM123_barcodes.tsv.gz`) and open the folder again.
- **"… has both `…features.tsv.gz` and `…genes.tsv.gz` beside it"**: keep
  only the list that names the matrix's rows in the folder.
- **"… are two matrices of one sample"** (`.mtx` and `.mtx.gz` with one
  prefix): keep one of them.

---

### Symptom: “Dataset loads, but I see no points (blank canvas)”

**Likely causes (ordered)**
1) WebGL/GPU issue (context lost, unsupported browser/GPU).
2) Embedding contains NaN/Inf or extreme values.
3) All points are filtered/hidden (rare on first load, but possible with a saved session).

**How to confirm**
- Try a demo dataset.
- Open browser console and look for WebGL/context errors.
- For data sanity, check embeddings for NaN/Inf (Python):

  ```python
  import numpy as np
  X = adata.obsm["X_umap_2d"]
  print(np.isfinite(X).all(), X.min(), X.max())
  ```

**Fix**
- WebGL issues → see {doc}`../a_orientation/02_system_requirements`
- NaN/Inf → clean/recompute embeddings and write the file again.

---

### Symptom: “Fields list is empty / missing expected metadata”

**Likely causes**
- `obs` is empty, or the file's launch defaults name only some columns
  (`uns['cellucid']['defaults']['obsKeys']`).
- Field names changed (session/expectations mismatch).

**How to confirm**
- AnnData: `print(adata.obs.head())`

**Fix**
- Ensure `obs` is populated before writing the file.

:::{note}
An `obs` column or gene name that starts or ends with whitespace, or that
contains a control or zero-width character (`U+200B`, `U+2060`, `U+FEFF`), is
not a name Cellucid can show exactly. Both writers refuse it, and an `.h5ad` or
Zarr opened directly is refused when it loads, with a message such as
`Invalid obs field list: field name "cell_type " cannot be displayed exactly:
Field name cannot have leading/trailing whitespace`. Rename the column
(for example `adata.obs.columns = adata.obs.columns.str.strip()`) and open the
file again.
:::

---

### Symptom: “Gene search returns nothing / genes are missing”

**Likely causes (ordered)**
1) No gene expression provided (no `X`, layer or `raw.X`), so there is nothing to query.
2) The genes were named from the wrong `var` column, so the names you search for are not the ones published.
3) The matrix in use has another gene list than the one you expect (`raw.X` often holds more genes than `X`).
4) Data is extremely large and gene fetch is slow (looks like “nothing”).

**How to confirm**
- AnnData: `print(adata.X is None)`, `print(adata.n_vars)`, `print(adata.var.head())`

**Fix**
- Choose the matrix and the gene names under **Session › Data mapping › Genes** (see [Genes and expression](13_data_mapping.md#genes-and-expression)).
- Provide a valid gene expression matrix, and name the genes from the `var` column your readers will search — `gene_id_column=...` in server/Jupyter mode, `var_gene_id_column=...` with `prepare()`.
- For a large row-major (CSR) matrix, the first gene waits while the browser
  builds a gene index, once; later genes read only their own entries. A
  column-major (CSC) matrix needs no index.

---

### Symptom: “Vector field overlay toggle is disabled / no fields appear”

**Likely causes (ordered)**
1) The dataset truly has no vector fields (common).
2) The field was refused when the overlay first showed it (a value that cannot
   be drawn); it is then listed disabled under **Session › Data mapping ›
   Vector fields** with the reason.
3) Naming mismatch: vectors are present but not discoverable (not `<field>_<basis>` after an `X_<basis>` embedding in the file with the same number of columns; e.g. `velocity_umap_2d` beside a bare `X_umap` with no `X_umap_2d`, which must be named `velocity_umap`).
4) Dimension mismatch: you only have 2D vectors but you’re looking at 3D (or vice versa).

**How to confirm**
- Open **Session ›
  Data mapping › Vector fields**, which lists every `obsm` array that fits an
  embedding and what is bound; in a notebook, `viewer.data_mapping["vectorFields"]`
  lists the bindings. See {doc}`13_data_mapping` and {doc}`05_jupyter_tutorial`.
- In Python, print `adata.obsm.keys()` and look for `<field>_<basis>` entries such as `velocity_umap` beside `X_umap`.

**Fix**
- Bind the array to its embedding under **Session › Data mapping › Vector fields**.
- Name the vector arrays by the supported convention and open the dataset again.
- Switch to the dimension that actually has vectors (2D vs 3D).

For overlay behavior and deeper diagnostics, see:
- {doc}`../i_vector_field_velocity/index`
- {doc}`../i_vector_field_velocity/07_troubleshooting_velocity_overlay`

---

### Symptom: “A link says the folder holds no `datasets.json`”

**Likely causes (ordered)**
1) Wrong repository path (pointing at another folder than the catalog folder).
2) Missing `datasets.json` (you didn’t write/commit it).
3) The files are not on the repository's default branch, but the shorthand omitted `@branch`.
4) Corporate network blocks `raw.githubusercontent.com`.

**How to confirm**
- The notification names the exact folder Cellucid read:

  ```text
  https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/ holds no datasets.json, so it is not a dataset catalog. …
  ```

- Open that folder's `datasets.json` in a browser tab. If it 404s there, the path or branch is wrong.

```{figure} ../../../_static/screenshots/data_loading/fail-github-catalog-not-found.png
:alt: A Cellucid notification saying the linked data could not be opened because the GitHub folder holds no datasets.json, after a repository path rather than its catalog folder was opened.
:width: 760px

One wrong path: one notification, *The linked data could not be opened: …
holds no datasets.json, so it is not a dataset catalog. …*, with the raw error
under **Details**. The folder it names is the exact request that failed.
```

**Fix**
- Write and commit `datasets.json` in the catalog folder
  ({doc}`07_file_format_expectations`).
- Or link to one file or Zarr store folder directly; it needs no catalog.
- Shorthand without `@branch` reads the repository's default branch; use
  `owner/repo@branch/path` for another branch.
- If blocked, use server mode or a different network.

Compare your repository with the known-good reference repository in
{doc}`11_custom_dataset_repository`, then see {doc}`02_local_demo_tutorial`.

---

### Symptom: “Google or Microsoft sign-in does not open, or does not finish”

The sign-in step sits under the **Link or server** bar and in a notification
({doc}`16_private_drive_onedrive`).

- **"The browser blocked … sign-in window."** Allow pop-ups for the site
  (the address bar's pop-up icon; in Safari also **Settings → Websites →
  Pop-up Windows**), then click the button again.
- **"… sign-in window was closed before it finished."** Click the button
  again; the window opens again. A sign-in you finish in a window Cellucid
  took for closed still counts. If you did sign in and Cellucid runs inside
  another page (a notebook), open Cellucid in its own tab.
- **"Google Drive works on a local cellucid serve only on its default port,
  8765."** Google answers only the addresses Cellucid registered, port
  included, for public and private files. Start `cellucid serve` without
  `--port`, or open the file on www.cellucid.com.
- **"Google Drive works in a notebook viewer only on port 8765."** A
  browser-tab viewer from Python (`AnnDataViewer`) serves on 8765 unless
  another server holds it, and then on the next free port (it prints a note). Stop the server holding 8765 (an earlier
  viewer: `viewer.stop()`; `cellucid serve`) and show the data again, or open
  the file on www.cellucid.com.
- **"Google Drive works on www.cellucid.com, not on …"** The page is on
  another computer's address or another site. Open the file on
  www.cellucid.com.
- **"Microsoft sign-in needs https:// or this computer's own address."** The
  page is plain `http://` on another computer's address, which Microsoft and
  the browser refuse for a sign-in. Open www.cellucid.com.
- **"Microsoft sign-in cannot reach Cellucid inside another site's page."**
  The browser keeps a page embedded in another site (a notebook at
  `localhost` showing a viewer at `127.0.0.1`) apart from Microsoft's sign-in
  window. Open Cellucid in its own tab, such as www.cellucid.com.
- **"Google Drive is not configured on this site yet"** or **"Signing in to
  Google / Microsoft is not configured on this site yet."** The site you use
  has no Google or Microsoft registration: a copy of the app hosted elsewhere
  never has one ({doc}`16_private_drive_onedrive`). Open the link on
  www.cellucid.com, or download the file and open it with **Open…**.
- **"Your organisation (…) lets only an administrator approve Cellucid"**
  (Microsoft said "Need admin approval", `AADSTS90094`, `AADSTS90095`,
  `AADSTS65001`). Your work or school account's organisation approves apps
  itself. Ask your IT administrator, with the steps in
  [Work and school accounts](16_private_drive_onedrive.md#work-and-school-accounts-administrator-approval)
  (the line under the bar links them), or click **Use another account** and
  sign in with a personal Microsoft account.
- **"Microsoft sign-in ended without Cellucid's permission"** (`AADSTS65004`).
  You declined Cellucid's request, or left Microsoft's "Need admin approval"
  page with "Return to the application without granting consent", which
  Microsoft reports the same way. Click **Sign in again** and accept, or see
  the previous item.
- **"Your organisation lets only the people it assigned use Cellucid"**
  (`AADSTS50105`), or **"Your organisation's sign-in policy does not let this
  browser or device sign in to Cellucid"** (Conditional Access, `AADSTS53003`).
  Ask your IT administrator, or use **Use another account**.
- **"Microsoft does not accept this site's sign-in (AADSTS…)"** or **"Google
  does not accept this site's sign-in (…)"**. The site's registration does not
  match the page (a redirect address, the client): open Cellucid at
  www.cellucid.com.
- **"Your Google Workspace administrator does not allow Cellucid for this
  account"**, or Google's window said **"Access blocked"**. Your organisation
  restricts third-party apps. Click **Use another account** (a personal Google
  account works), or ask your administrator, with the steps in
  [Google Workspace accounts](16_private_drive_onedrive.md#google-workspace-accounts-administrator-approval).
- **"You declined Cellucid's request to read the Google Drive files you choose
  for it."** Click **Sign in again** and allow it.
- **"Google sign-in expired" or "Microsoft sign-in expired".** Click **Sign in
  again**; the reads that waited go on.
- **"Google Drive sends this file only with a Google sign-in for now"**, or a
  read failing with **"Google Drive is not sending this file to this browser
  without a Google sign-in for now"**. The file is public, but Google is
  holding back downloads without a sign-in from your network, after many of
  them (it answers them with an "unusual traffic" page that the browser does
  not let Cellucid read), for up to about an hour. Signing in avoids it:
  click **Sign in to Google** and **Select** the file in Google's picker, and
  it is read with your sign-in (signed in already, opening the file again
  goes straight to the picker). Or wait and open it again
  ({doc}`16_private_drive_onedrive`). On a shared lab or university network,
  other people's downloads count too, and an address released from a hold is
  held again sooner (after some 30 MB). To publish a dataset for everyone,
  use Zenodo, Figshare, Hugging Face or GitHub, which need no sign-in.

---

### Symptom: “Remote server connect fails / CORS blocked / mixed content”

```{figure} ../../../_static/screenshots/data_loading/fail-remote-server-unreachable.png
:alt: A Cellucid notification saying the Cellucid server could not be reached at a loopback address, and what to check.
:width: 760px

Connecting to a port with nothing listening on it. Cellucid deliberately names
*both* possibilities, because a failed connection and a wrong address are
indistinguishable from the browser's side.
```

**Likely causes**
- Wrong URL (host/port mismatch) — this and a dead server produce the same
  message.
- You pasted the **Viewer URL** into the **Link or server** bar. Its `?…`
  query makes the bar open it as a link to data, not as a server; type the
  **Local URL** the banner prints (`http://127.0.0.1:8765`), `127.0.0.1:8765`
  or just the port.
- Server is bound to `127.0.0.1` on a remote machine, but you’re trying to access it directly.
- The Cellucid page is HTTPS and the address is `http://`, which Cellucid
  refuses before the browser does:
  `An HTTPS Cellucid page requires an explicit HTTPS remote server URL`.

**How to confirm**
- Check server banner output and copy the exact viewer URL printed by the server.
- Visit the health endpoint — it separates the two causes the message cannot:

  ```text
  http://127.0.0.1:8765/_cellucid/health
  ```

  JSON back means the address was wrong; nothing back means the server is not
  running or not reachable.

**Fix**
- For remote machines: use an SSH tunnel (recommended): {doc}`04_server_tutorial`
- Keep the server bound to `127.0.0.1` and access via `http://localhost:<port>` through the tunnel.

**Reassurance**

A failed connection does not throw away what you already had. The dataset that
was open stays open, and the control returns to its unconnected state:

```{figure} ../../../_static/screenshots/data_loading/keep-dataset-after-failed-connect.png
:alt: The whole Cellucid window after a server address with nothing behind it was opened: the Pancreas sample is still drawn and still named in the Dataset box, and one error notification sits in the lower right.
:width: 1440px

After a refused connection the previous dataset is untouched — same points, same
summary, same fields. Only the notifications changed. You can correct the
address and try again without reloading the page.
```

---

## When to stop debugging and switch workflows (high leverage)

If you hit one of these, switching workflows saves time:

- **Large `.h5ad` on a remote machine** → serve it from there over an SSH tunnel.
- **GitHub too large / blocked** → write the dataset as a Zarr store folder, switch to server mode, or host the files on an accessible HTTP server.
- **File picker blocked** → switch to server mode.

---

## What to include when asking for help (so others can debug quickly)

Copy/paste:
- your loading path (H5AD / Zarr ZIP / Zarr folder / link / server / GitHub)
- approximate dataset size: `n_cells`, `n_genes`
- browser + OS
- any console error text
- for server mode: the exact `cellucid serve ...` command and the printed URL

---

## Related pages

- Loading options map: {doc}`01_loading_options_overview`
- Local file picker: {doc}`03_browser_file_picker_tutorial`
- Server mode: {doc}`04_server_tutorial`
- Publishing a dataset catalog: {doc}`02_local_demo_tutorial`
- Complete custom repository reference: {doc}`11_custom_dataset_repository`
- Dataset identity: {doc}`06_dataset_identity_why_it_matters`
- Vector fields overlay: {doc}`../i_vector_field_velocity/index`
