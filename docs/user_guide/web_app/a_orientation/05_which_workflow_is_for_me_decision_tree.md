# “Which workflow is for me?” decision tree

**Audience:** everyone (especially mixed teams)  
**Time:** 5–10 minutes  
**Goal:** pick the simplest workflow that matches your data + sharing needs

---

## The short answer (most common choices)

- If you are a **wet-lab scientist** and someone sent you a `.h5ad`, a `.zarr`
  folder or a `.zarr.zip` → open it with **Open…** in the web app (**File…**
  for an `.h5ad` or `.zarr.zip`, **Zarr folder…** for a folder), or drop it
  on the page; nothing is installed and nothing is uploaded
  ({doc}`../b_data_loading/03_browser_file_picker_tutorial`).
- If someone sent you a **link** → paste it into the **Link or server** bar
  and click **Open** ({doc}`../b_data_loading/14_open_from_link`).
- If you are a **computational user with AnnData** → use `cellucid-python`
  (`show_anndata(...)` in a notebook, `cellucid serve` in a terminal, and
  `prepare(...)` to write a Cellucid-ready AnnData that opens fastest).
- If you work in **R** → `cellucid_prepare()` writes the same Cellucid-ready
  AnnData ({doc}`../../r_package/index`).
- If you need **multi-user labeling/voting** → use **Community annotation** (GitHub-backed).

---

## Decision tree (plain language first)

::::{tab-set}

:::{tab-item} Wet‑Lab / New to Coding

1) What did your colleague send you?
   - A `.h5ad` file → use **Open… → File…**. You do not need Python for this.
   - A folder ending in `.zarr` → use **Open… → Zarr folder…**, or drop the
     folder on the page.
   - A `.zarr.zip` → use **Open… → File…**.
   - A link → paste it into the **Link or server** bar and click **Open**.
   - Nothing yet → ask a computational colleague for an AnnData file. One
     written by `cellucid.prepare()` in Python or `cellucid_prepare()` in R
     opens fastest; both write the same thing, so it does not matter which
     they use.

2) Do you need to *share* what you did?
   - Yes → use **Save** under **Session file** to download a `.cellucid-session` file, then send it.
   - If you need many people to vote on labels → use **Community annotation** instead (it’s designed for that).

3) Does the dataset feel too big / slow?
   - Ask for the file written by `prepare()` (a gene-major layout made for
     reading one gene at a time), or for a “server mode” workflow when the data
     lives on another machine (your colleague runs a small server there and the
     browser reads only what each view needs).

:::

:::{tab-item} Computational / Power User

Start here:

```
Is the AnnData file (.h5ad, .zarr, .zarr.zip) on this computer?
  ├─ Yes → open it in the web app (Open… → File… / Zarr folder…, or drop it)
  └─ No →
       Is it on the web (GEO, Zenodo, Hugging Face, GitHub, a URL)?
         ├─ Yes → paste the link into the Link or server bar
         └─ No → is it in Python (AnnData in memory or on a remote disk)?
                  ├─ in a notebook → show_anndata(adata, dataset_name="My dataset", dataset_id="my-dataset")
                  └─ on a remote machine → cellucid serve <file> + SSH tunnel
```

For the fastest web experience, write the dataset once with
`cellucid.prepare("my-dataset.zarr", format="zarr", ...)` in Python or `cellucid_prepare()` in
R; the viewer reads either one identically ({doc}`../../r_package/index`).

Then decide on sharing:

- One-to-one sharing / reproducibility → `.cellucid-session` bundles
- A public or intranet catalog on a static host → publish a folder of AnnData
  files with a `datasets.json` catalog and open it with
  `?data=https://host/datasets/`
  ({doc}`../b_data_loading/02_local_demo_tutorial`)
- Many-to-many annotation / voting → Community Annotation repo

:::

::::

---

## Workflow table (when to use what)

| You want to… | Recommended workflow | Why | Common gotcha |
|---|---|---|---|
| Try Cellucid quickly | Demo dataset | Zero setup | Not your data |
| Explore your own data locally | Web app: Open… → File… (`.h5ad`, `.zarr.zip`) or Zarr folder… | No server needed | Files not written gene-major read more bytes per gene |
| Share a dataset with collaborators | Host the AnnData file (HTTP) and share its `?data=` link ({doc}`../b_data_loading/14_open_from_link`) | Everyone opens the same dataset URL | CORS / hosting config |
| Work in notebooks | `cellucid-python` (embedded viewer + hooks) | Tight Python ↔ UI loop | Notebook iframe restrictions (fullscreen/pointer lock) |
| Prepare a dataset without Python | `cellucid-r` (`cellucid_prepare()`) | Same Cellucid-ready AnnData, read by the same viewer | Arguments are named differently from `prepare()` |
| Compare multiple hypotheses | Multiview snapshots (“Keep view”) | Side-by-side views | Smoke mode is disabled when snapshots exist |
| Run a multi-user labeling round | Community Annotation | Conflict-free collaboration | Requires GitHub repo + app installation |

---

## Edge cases and “which path wins?”

### “I can load my H5AD directly, why run `prepare()`?”

Every AnnData file opens as it is. Writing it again with `prepare()` helps
when you need:

- fast gene reads (a gene-major matrix, so one gene is one small read),
- a stable dataset id, name and citation stored in the file itself
  (`uns['cellucid']`), so sessions and links keep matching it,
- defaults (which embedding, which `obs` columns) chosen by the author.

### “My dataset is huge (millions of cells). What should I do?”

Prefer:

- an AnnData written by `prepare()` (gene-major, chunked for range reads),
  and/or
- server mode when the data lives on a remote machine (the browser reads only what each view needs, over the tunnel).

Also plan for performance from day one:

- keep category counts reasonable,
- avoid exploding the number of views,
- be cautious with smoke mode.

### “I have sensitive data”

Choose workflows based on your constraints:

- a local file opened in the browser (no upload) is simplest,
- remote hosting requires careful access control,
- community annotation uses GitHub (be sure your repo is private if needed).

---

## Troubleshooting (choosing a workflow)

### Symptom: “I’m not sure what I have”

**How to confirm**

- If you have a `.h5ad` file → use **Open… → File…** in the browser, or
  `cellucid serve` / `show_anndata()` when you work in Python
  ({doc}`../b_data_loading/03_browser_file_picker_tutorial`,
  {doc}`../b_data_loading/04_server_tutorial`).
- If you have a `.zarr` directory → use **Open… → Zarr folder…** in the
  browser, drop the folder on the page, or point `cellucid serve` / `show_anndata()`
  at it ({doc}`../b_data_loading/03_browser_file_picker_tutorial`,
  {doc}`../b_data_loading/04_server_tutorial`).
- If you have a folder holding `datasets.json` → it is a dataset catalog;
  publish it on a static host and open its address in the **Link or server**
  bar
  ({doc}`../b_data_loading/02_local_demo_tutorial`).

### Symptom: “My colleague sent me a folder but the app can’t load it”

Likely you need the data-loading troubleshooting section ({doc}`../b_data_loading/index`) and/or the system requirements page ({doc}`02_system_requirements`).

---

## Next steps

- If you picked a web workflow: continue to {doc}`../b_data_loading/index`
- If you picked the notebook/Python workflow: go to {doc}`../../python_package/index`
- If your data lives in R: go to {doc}`../../r_package/index`
