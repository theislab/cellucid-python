# Open from a link

**Audience:** everyone who wants to look at data that already lives somewhere
public
**Time:** 3 minutes
**What you'll learn:**
- How the **Link or server** address bar decides what it opens.
- What you can paste into it, and what each one opens.
- How data is read where it lives, without downloading it first.
- What the `?data=` parameter of a share link means.
- Why a link can be refused, and what to do then.

---

## At a glance

The **Link or server** address bar opens data where it already lives: a dataset on CELLxGENE,
a GEO series, a Zenodo record, a Figshare item, a Hugging Face dataset, a
Dropbox, Google Drive or OneDrive file, a GitHub repository, or any web address
of an `.h5ad` file, a Zarr store folder, a Zarr store zipped into a `.zip`, a
Loom (`.loom`) or MuData (`.h5mu`) file, a Seurat or SingleCellExperiment
object (`.rds`, `.rda`, `.RData`), or a dataset catalog folder (one holding
`datasets.json`). Every file is read in place, except an R object or a GEO
MatrixMarket sample, which is read once as it downloads
({doc}`07_file_format_expectations`). Nothing needs a server
of your own, and only a private Google Drive file or a OneDrive file needs a
sign-in ({doc}`16_private_drive_onedrive`).

It is in the **Session** panel under **Open data**. Paste a link or an id,
then click **Open** (or press Enter). The page address gains a `?data=`
parameter naming what you opened in its provider's form (the `?data=` column
of [What you can paste](#what-you-can-paste): a pasted Google Drive link
becomes `gdrive:<id>`), also for a link the page was opened with, and never a
sign-in token or a download address. So copying the address (or using **Share**)
gives others the same view. The same bar connects to a Cellucid server
({doc}`04_server_tutorial`): see [The address bar](#the-address-bar).

Cellucid reads the data **in place**: only the bytes a view needs are
requested from the host, with HTTP range requests. A 1.6 GB CELLxGENE file
shows its embedding after reading about 26 MB.

The fastest way to see it is to click one of the datasets in
[Try it](#open-from-link-try-it) below: each opens in Cellucid straight from
where it is published.

## The address bar

One bar opens both linked data and a `cellucid serve` server. What you type
decides which, by the same rules every time:

| You type | It opens |
|---|---|
| A port number: `8765` | The server on this computer, `http://localhost:8765` |
| `localhost:8765`, `127.0.0.1:8765`, or `host:port` without `http://` | The server at that address, over `http://` |
| A provider's form: `hf:user/repo`, `GSE123456`, a DOI, `zenodo:<record>`, `owner/repo/path`, a CELLxGENE, Zenodo, Figshare, Google Drive, OneDrive, Dropbox, Hugging Face or GitHub link | That provider's data ([What you can paste](#what-you-can-paste)) |
| A web address ending in a file or store name Cellucid reads (`.h5ad`, `.h5mu`, `.loom`, `.h5`, `.rds`, `.rda`, `.zip`, `.zarr`) or `datasets.json`, or one with a `?query` or `#fragment` | That file, store or catalog |
| `http://localhost:8765` (this computer, no path) | The server |
| Any other web address | A server if one answers there (Cellucid asks its `/_cellucid/health` endpoint), otherwise the folder of datasets or the file at that address |

While something is connected, a line under the bar names it (**Server** or
**Link**, then its host and the end of its address; the whole address is its
tooltip) with **Disconnect**. The bar stays free for the next address. The
dataset on screen stays while another opens, with its progress (and Cancel)
in a notification; the new one replaces it only once it is ready, and on a
failure or Cancel the dataset on screen stays. A web address that is first
checked for a server shows that check as the opening's first step. A catalog
that holds several datasets lists them in the dataset list above, whose label
then names its host ("Dataset from theislab.github.io"); when the dataset on
screen came from another link, the catalog's default dataset replaces it.
Otherwise the dataset list offers the samples by name; while a link or a
server is connected, their files (each sample's title and cell count) are read
only once you open the list. **Disconnect** brings the samples back.

```{figure} ../../../_static/screenshots/data_loading/open-from-link-control.png
:alt: The Link or server bar with "zenodo:23288289" typed into it, the pointer on its Open button.
:width: 496px

A Zenodo record typed into the **Link or server** bar. **Open** (or Enter)
opens it, and the page's address then carries `?data=zenodo:23288289`.
```

```{figure} ../../../_static/screenshots/data_loading/connect-data-link.png
:alt: The Open data group connected to a dataset catalog folder: the dropdown labelled "Dataset from" the catalog's host and holding its datasets, Open…, and the Link or server bar holding the address, with one line under it reading Link, the address and a compact Disconnect button.
:width: 496px

A dataset catalog folder opened from a link: the dataset list is labelled with
the catalog's host and lists its datasets, and the line under the bar reads
**Link**, the host and **Disconnect**, which brings the samples back.
```

## What you can paste

| You paste | It opens | `?data=` form |
|---|---|---|
| A CELLxGENE dataset id, or its **Explore** link (`cellxgene.cziscience.com/e/<dataset id>.cxg/`) | That dataset's `.h5ad` | `cellxgene:<dataset id>` |
| A GEO series (`GSE123456`) or sample (`GSM1234567`), its accession page, or a supplementary file link | The one dataset among its supplementary files: an `.h5ad`, a gzipped R object (`.rds.gz`, `.rda.gz`, `.RData.gz`) or another file Cellucid opens, or a MatrixMarket sample (`GSM…_matrix.mtx.gz` with its `_features`/`_genes` and `_barcodes` files), named by its matrix | `geo:GSE123456`, `geo:GSE123456/<file>` or `geo:GSE123456/<GSM…_matrix.mtx.gz>` |
| A Zenodo record link, file link, DOI (`10.5281/zenodo.<id>`) or `zenodo:<record>` | The record's one file Cellucid opens (`.h5ad`, `.zip`, `.rds`, `.h5mu`, …) or the file the link names; a concept DOI opens the latest version | `zenodo:<record>[/<file>]` |
| A Figshare item page (`figshare.com/articles/<type>/<name>/<id>[/<version>]`, also on an institution's `*.figshare.com` portal), a Figshare DOI (`10.6084/m9.figshare.<id>[.v<n>]`), any DOI that resolves to a Figshare item, or `figshare:<id>` | The item's one file Cellucid opens (or the file the link names, by name or by `?file=<file id>`); a link without a version opens the latest | `figshare:<id>[v<version>][/<file>]` |
| A Google Drive file link (`drive.google.com/file/d/<id>/…`, `…/open?id=<id>`, `…/uc?id=<id>`) | That file, read through Google's Drive API; a private file after **Sign in to Google** ({doc}`16_private_drive_onedrive`); folder links are refused | `gdrive:<id>` |
| A OneDrive or SharePoint file share link (`1drv.ms/…`, `onedrive.live.com/…`, `<tenant>-my.sharepoint.com/:u:/…`) | That file, after **Sign in with Microsoft** ({doc}`16_private_drive_onedrive`); folder links are refused | `onedrive:<drive>/<item>` once opened |
| A Hugging Face dataset link (`huggingface.co/datasets/<user>/<repo>/…`) or `hf:<user>/<repo>[@<revision>][/<path>]` | A dataset catalog folder, or one file Cellucid opens; pinned to the commit the branch points at | `hf:<user>/<repo>[@<revision>][/<path>]` |
| A Dropbox file link | That file, read from Dropbox's content host (folder links need a Dropbox sign-in and are refused) | `dropbox:<share link>` (`https://www.dropbox.com/…`, without `dl=`) |
| A GitHub repository path (`owner/repo/folder`, read from the default branch, or `owner/repo@branch/folder`), a `github.com/.../tree/...`, `.../blob/...` or `.../raw/...` link, or a `raw.githubusercontent.com` link (the **Raw** button's `refs/heads/<branch>` form included) | The dataset catalog in that folder ({doc}`11_custom_dataset_repository`), the Zarr store folder (`….zarr`) or the one file the path names ({doc}`15_host_on_github`); a file stored with Git LFS is read from GitHub's media host | `gh:owner/repo/folder`, `gh:owner/repo/path/atlas.zarr`, `gh:owner/repo/path/file.h5ad` |
| A web address ending in `/`, or its `datasets.json` | The dataset catalog in that folder ({doc}`07_file_format_expectations`) | the address itself |
| A web address of a Zarr store folder (`…/name.zarr/`, Zarr v2 or v3, with consolidated metadata as anndata writes it): a bucket, a lab server, GitHub raw or a Hugging Face `resolve` folder | That store, its chunks read as views need them | the address itself |
| A web address of a file Cellucid opens (an `.h5ad`, a Zarr `.zip`, an `.rds`, …): a bucket, a lab server, a Hugging Face file link (`huggingface.co/datasets/<user>/<repo>/resolve/main/<file>`), or a Zenodo API file link (`zenodo.org/api/records/<id>/files/<file>/content`; Zenodo's ordinary download link cannot be read by other sites) | That file | the address itself |

A link that holds several datasets (a catalog) opens on the dataset its
`datasets.json` names as default and lists all of them in the dataset list
(labelled "Dataset from <host>" while the catalog is connected); a share link
adds `&dataset=<id>` to pick one. The list
reads each dataset's name and cell count from its file's metadata; a dataset's
file is opened when you choose it. A catalog that names a default session for a
dataset applies it, exactly as the built-in samples do.

(open-from-link-try-it)=
## Try it

Every dataset below is published by someone else, on the host its section
names, and opens in Cellucid read in place. Click **Open**. The page address
carries the link (`?data=…`), so it is also what you would share.

"Opening reads" is how much of the file Cellucid fetched to show its first
view, measured on 2026-10-22 (the R object and MatrixMarket rows on
2026-10-24, the Google Drive and OneDrive rows on 2026-10-25); the rest is read only when a view needs it (showing a gene,
switching an embedding). An R object or a MatrixMarket matrix is read whole,
once. The tag names the ecosystem whose tools write the file: **Python** for
AnnData `.h5ad`, **R** for a Seurat object; a Zarr store, written by either,
has none. Hosts change: if a link stops
working, the [link checker](#open-from-link-checking) says which.

### CELLxGENE

Census copies of the 2025-11-10 release. An **Explore** link works too:
[this one](https://www.cellucid.com/?data=https%3A%2F%2Fcellxgene.cziscience.com%2Fe%2Fa49d9109-1d0c-4b36-8139-19aa9a83428c.cxg%2F)
opens the first row.

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Healthy human liver: lymphocytes | Python | 2,346 | 15 MB | 4 MB | [Open](https://www.cellucid.com/?data=cellxgene:a49d9109-1d0c-4b36-8139-19aa9a83428c) |
| Tabula Sapiens: neural | Python | 2,685 | 122 MB | 11 MB | [Open](https://www.cellucid.com/?data=cellxgene:6e4acdae-59ef-452b-a6d4-260cba368586) |
| Human Brain Cell Atlas: hippocampus, rostral CA3 | Python | 8,077 | 106 MB | 7 MB | [Open](https://www.cellucid.com/?data=cellxgene:d5452b83-7c3d-4d7c-ab7a-c7fece7196c5) |
| Human Immune Health Atlas: other cells | Python | 10,937 | 56 MB | 10 MB | [Open](https://www.cellucid.com/?data=cellxgene:81ffdcd6-3f04-40bc-8bb0-9748f4ffe576) |
| Tabula Sapiens: mammary | Python | 30,936 | 1.2 GB | 18 MB | [Open](https://www.cellucid.com/?data=cellxgene:2ba40233-8576-4dec-a5f1-2adfa115e2dc) |
| Foveal and parafoveal retina | Python | 34,323 | 252 MB | 8 MB | [Open](https://www.cellucid.com/?data=cellxgene:f5be9ed2-5d41-4a52-91e4-4ff24ff84900) |
| Tabula Sapiens: thymus | Python | 42,729 | 1.6 GB | 26 MB | [Open](https://www.cellucid.com/?data=cellxgene:0ced5e76-6040-47ff-8a72-93847965afc0) |
| Krasnow Lab Human Lung Cell Atlas, 10x | Python | 65,662 | 569 MB | 12 MB | [Open](https://www.cellucid.com/?data=cellxgene:8c42cfd0-0b0a-46d5-910c-fc833d83c45e) |
| Pancreas, harmonized across studies | Python | 121,916 | 915 MB | 14 MB | [Open](https://www.cellucid.com/?data=cellxgene:31f657dc-1875-4c4b-a5ca-ce63b3ef3a82) |

### Zenodo

A record with one openable file opens it; a record with several names them,
and the link names one with `zenodo:<record>/<file>`.

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Day 3 Session 4: processed Scanpy AnnData (record 13255338) | Python | 624 | 8 MB | 2 MB | [Open](https://www.cellucid.com/?data=zenodo:13255338/processedScanpyAnndata.h5ad) |
| Uninfected fibroblasts, RNA velocity (record 7608772) | Python | 1,545 | 193 MB | 4 MB | [Open](https://www.cellucid.com/?data=zenodo:7608772/uninfected_fibroblast.h5ad) |
| AML atlas (record 15713160) | Python | 1,732 | 1.2 GB | 4 MB | [Open](https://www.cellucid.com/?data=zenodo:15713160/AML_Atlas.h5ad) |
| Unwounded skin fibroblasts, integrated (record 10198892) | Python | 7,559 | 105 MB | 8 MB | [Open](https://www.cellucid.com/?data=zenodo:10198892/unwounded_fibroblast_merged.h5ad) |
| Bradyzoite heterogeneity, with a 3-D embedding (record 23288289) | Python | 8,057 | 15 MB | 2 MB | [Open](https://www.cellucid.com/?data=zenodo:23288289) |
| Ankylosing spondylitis T cells, TCR-DeepInsight (record 12741480) | Python | 15,065 | 22 MB | 7 MB | [Open](https://www.cellucid.com/?data=zenodo:12741480/Yi_2023_Ankylosing_Spondylitis.h5ad) |
| Kang et al. 2018 PBMC, IFN-β (record 10850398) | Python | 24,562 | 264 MB | 7 MB | [Open](https://www.cellucid.com/?data=zenodo:10850398/kang18_counts_25k.h5ad) |
| Liao et al. 2020 COVID-19 BALF (record 10850398) | Python | 65,792 | 828 MB | 25 MB | [Open](https://www.cellucid.com/?data=zenodo:10850398/liao20_sub.h5ad) |
| Human embryo at Carnegie stage 9, annotated; written by anndata 0.7 (record 14730926) | Python | 2,150 | 17 MB | 3 MB | [Open](https://www.cellucid.com/?data=zenodo:14730926/cs9_scRNA_annotation.h5ad) |
| Mast cells across mouse organs, a Seurat object with three assays (record 8071138) | R | 371 | 20 MB | 20 MB, once | [Open](https://www.cellucid.com/?data=zenodo:8071138/Mice_Mast_cells.Rds) |

### Figshare

Items published on figshare.com, each holding one `.h5ad`. The file is read
from Amazon S3 through the signed link Figshare's file host gives the page,
which works for ten seconds; within those seconds every read goes straight to
S3.

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| scverse tutorial data: PBMC 3k, processed (item 22577536) | Python | 2,638 | 100 MB | 3 MB | [Open](https://www.cellucid.com/?data=figshare:22577536) |
| Scanpy tutorial: PBMC 3k (item 16447278) | Python | 2,638 | 32 MB | 4 MB | [Open](https://www.cellucid.com/?data=figshare:16447278) |
| Scanpy tutorial: trajectory, version 2, by its DOI (item 16540047) | Python | 2,730 | 12 MB | 2 MB | [Open](https://www.cellucid.com/?data=10.6084%2Fm9.figshare.16540047.v2) |
| Squidpy tutorial: Visium mouse brain (item 16566057) | Python | 704 | 28 MB | 11 MB | [Open](https://www.cellucid.com/?data=figshare:16566057) |

### Google Drive

Files shared on Google Drive with "Anyone with the link" open with no
sign-in: each is read through Google's Drive API with Cellucid's public key,
in place, one range at a time. On 2026-10-25 Drive answered each read in
about 0.65 s (0.56–0.98 s; a few kilobytes and a few megabytes alike), at about
4 MB/s a request, and Cellucid opened each file with the same requests and
bytes as from a host answering ranges after the same 0.65 s. Read with a
Google sign-in, a file takes the same requests and the same time
({doc}`16_private_drive_onedrive`): the 1.7 GB thymus below opened signed in
on 9.3 MB in 84 requests.

Google may refuse downloads without a sign-in from a network for a while
after heavy use (a few hundred megabytes within minutes; after a first such
hold, some 30 MB), and on a shared lab or university network other people's
downloads count too. Cellucid then asks you to sign in to Google, which
always works. So Drive and OneDrive suit sharing with collaborators, above all
privately; to publish a dataset for everyone, Zenodo, Figshare, Hugging Face
or GitHub serve it with no sign-in and no such hold.

The kit of test files Cellucid keeps on Google Drive, Python's and R's
formats side by side:

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Healthy human liver: lymphocytes (CELLxGENE) | Python | 2,346 | 15 MB | 4 MB | [Open](https://www.cellucid.com/?data=gdrive:1lylTBkM-sOjKMsA2gve3gKzcrI3XIyVQ) |
| Tabula Sapiens: thymus (CELLxGENE) | Python | 42,729 | 1.7 GB | 9 MB | [Open](https://www.cellucid.com/?data=gdrive:1tWxN8Ls8Ars8yQ7YYlQ1t-C6FmHXCivF) |
| Pancreas endocrinogenesis, Cellucid-ready Zarr ZIP |  | 3,696 | 3.8 MB | 0.3 MB | [Open](https://www.cellucid.com/?data=gdrive:1YHfT0PYXbcOqmBnT8xcdSLOYgrl-aqnR) |
| Pancreas endocrinogenesis, Loom | Python | 3,696 | 5.9 MB | 1.3 MB | [Open](https://www.cellucid.com/?data=gdrive:1yVD9HcN6dFIlhWhMmDiv2U8P_e1PF4El) |
| Synthetic CITE-seq, MuData (`rna`, `prot`) | Python | 60 | 0.1 MB | 0.1 MB | [Open](https://www.cellucid.com/?data=gdrive:1LvugT8VSn2ZAeqHFy2fEyWEySp6gOVLT) |
| Pancreas endocrinogenesis, Seurat v5 `.rds` | R | 3,696 | 9.8 MB | 9.8 MB, once | [Open](https://www.cellucid.com/?data=gdrive:1Mat7COHtWTvoiIlI6xeR3NxSq53UqvnP) |
| Pancreas endocrinogenesis, Seurat v5 `.RData` | R | 3,696 | 9.8 MB | 9.8 MB, once | [Open](https://www.cellucid.com/?data=gdrive:1yxS0L0LfqJXyzzhMI13ewdQWfIzzIgWp) |
| Pancreas endocrinogenesis, SingleCellExperiment `.rds` | R | 3,696 | 5.4 MB | 5.4 MB, once | [Open](https://www.cellucid.com/?data=gdrive:1Nw1PYY4WE4dUnfdhxpV9nEy2nboQMT3X) |
| Mast cells across mouse organs, Seurat with three assays (Zenodo 8071138) | R | 371 | 20 MB | 20 MB, once | [Open](https://www.cellucid.com/?data=gdrive:1A3UVhxjBjUJfD5hWRp5K1yAQo4Vk0RVc) |

A 10x Cell Ranger `.h5` on Drive ([this one](https://www.cellucid.com/?data=gdrive:1SLdms6shM74UV1leT2dsD--Qqf6X_rW5)) is refused, as it is from any host:
it holds counts and no embedding, and the message names the `outs/` folder
to open instead.

Files their authors share on Google Drive, from the tutorials and
repositories that link them:

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| RegVelo tutorial: zebrafish neural crest, preprocessed (theislab/regvelo) | Python | 697 | 42 MB | 2 MB | [Open](https://www.cellucid.com/?data=gdrive:1Nzq1F6dGw-nR9lhRLfZdHOG7dcYq7P0i) |
| ClusterMap: STARmap mouse placenta (wanglab-broad/ClusterMap) | Python | 7,224 | 104 MB | 2 MB | [Open](https://www.cellucid.com/?data=gdrive:1DDCowUuZ7PPFUSZsjvSqntWkYJMjf1Na) |
| RegVelo: Perturb-seq, with PHATE and t-SNE (theislab/regvelo) | Python | 12,393 | 379 MB | 3 MB | [Open](https://www.cellucid.com/?data=gdrive:1CKADDykadX0CLjYGPgc3YXA-AGRSaA1k) |
| COSG tutorial: SEA-AD middle temporal gyrus, 20,000-cell subsample (genecell/COSG) | Python | 20,000 | 1.7 GB | 8 MB | [Open](https://www.cellucid.com/?data=gdrive:1nH-CRaTQFxJ5pAVpy8_hUQn1nrIcakq2) |
| CyteType example: Human Lung Cell Atlas core, downsampled (NygenAnalytics/CyteTypeR) | Python | 27,887 | 245 MB | 5 MB | [Open](https://www.cellucid.com/?data=gdrive:13O0dyUnwJKLPm8fncRt597S5hs2COsxx) |
| Developing human lung (He et al.), its scVI/scANVI integration (theislab/sctram) | Python | 71,650 | 285 MB | 7 MB | [Open](https://www.cellucid.com/?data=gdrive:1Yhd-Nt5kGphZz1bYItJSRSAm-JVxm34J) |
| scHolography: human skin Visium at low resolution, a Seurat object (YiLab-SC/scHolography) | R | 659 | 12 MB | 12 MB, once | [Open](https://www.cellucid.com/?data=gdrive:1aF0zv5uEHu-42BJQWXtkLy7qrf9ELfq9) |

Any of the link forms Drive gives works the same way: the **Share** link
(`drive.google.com/file/d/<id>/view?usp=sharing`), `open?id=`, the
`uc?id=…&export=download` links that `gdown` uses (also under `/u/0/`),
`docs.google.com/uc?id=` and `drive.usercontent.google.com/download?id=`.

### OneDrive

OneDrive and SharePoint send no file without a Microsoft sign-in, not even one
shared with "Anyone with the link": Microsoft Graph and the OneDrive API
answer such a request for a share link with HTTP 401 (checked 2026-10-25). A
OneDrive link therefore opens after **Sign in with Microsoft**
({doc}`16_private_drive_onedrive`), with any personal Microsoft account, and
is then read in place from Microsoft's download address, with plain range
requests that need no sign-in of their own. The same kit, shared from a
personal OneDrive with "Anyone with the link" (read with a personal account on
2026-10-25):

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Healthy human liver: lymphocytes (CELLxGENE) | Python | 2,346 | 15 MB | 3 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQABK1r8aMTnTa67w5AFZOXNAfiIHx3QvSPmnViA_x9vfhk) (Microsoft sign-in needed) |
| Foveal and parafoveal retina (CELLxGENE) | Python | 34,323 | 264 MB | 5 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQBJiL1OCpqBRpC9Ug346BYlAT_1bXWFDRKXWMnwsqaoOU0) (Microsoft sign-in needed) |
| Pancreas endocrinogenesis, Cellucid-ready Zarr ZIP |  | 3,696 | 3.8 MB | 0.3 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQCvOqEz5Nm8SJ5aMbRRicVuATbYzopTeG1Q8WRuZBC0Sxw) (Microsoft sign-in needed) |
| Pancreas endocrinogenesis, Loom | Python | 3,696 | 5.9 MB | 1.5 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQB83DYmO8VfR53u9on1M6KrAeNbQF9uHDun_0IJZJB8sys) (Microsoft sign-in needed) |
| Synthetic CITE-seq, MuData (`rna`, `prot`) | Python | 60 | 0.1 MB | 0.1 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQCuSSGLF9jVRKAXzNh7YO0wAbgx7nWyYrdVI585JwBmrAs) (Microsoft sign-in needed) |
| Pancreas endocrinogenesis, Seurat v5 `.rds` | R | 3,696 | 9.8 MB | 9.8 MB, once | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQBifyr0Gm5eTZjx_Zyfh34CAUSaMW68VGPj-o_YTzmsMhg) (Microsoft sign-in needed) |
| Pancreas endocrinogenesis, Seurat v5 `.RData` | R | 3,696 | 9.8 MB | 9.8 MB, once | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQAizMaHYAdxRojZpIlItUejAf5N3Qcv9rY6GTqPbK1ZWJ0) (Microsoft sign-in needed) |
| Pancreas endocrinogenesis, SingleCellExperiment `.rds` | R | 3,696 | 5.4 MB | 5.4 MB, once | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQDYps4oi-SqQp3FXeAodLdjAcdFfTlmEuHRFwbHdtbb4kw) (Microsoft sign-in needed) |
| Mast cells across mouse organs, Seurat with three assays (Zenodo 8071138) | R | 371 | 20 MB | 20 MB, once | [Open](https://www.cellucid.com/?data=https%3A%2F%2F1drv.ms%2Fu%2Fc%2F479799233e1614b2%2FIQDSAf3ha2nmT7bM51LKKpXDAeRit-YJFLvG8_36ytt_x9I) (Microsoft sign-in needed) |

### Hugging Face

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Kang et al. IFN-β PBMC counts | Python | 8,000 | 38 MB | 2 MB | [Open](https://www.cellucid.com/?data=hf:Nikita1113/tb-science-scrna-kang-ifnb) |
| PBMC, preprocessed | Python | 17,041 | 431 MB | 10 MB | [Open](https://www.cellucid.com/?data=hf:datavil/pbmc3k) |

### GEO

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| GSE350831: Flex pilot, annotated | Python | 21,751 | 868 MB | 9 MB | [Open](https://www.cellucid.com/?data=geo:GSE350831/GSE350831_Flex_pilot_paper_annotated_bbknn.h5ad) |
| GSE327505: Saber et al. 2026, processed | Python | 25,651 | 2.1 GB | 12 MB | [Open](https://www.cellucid.com/?data=geo:GSE327505) |
| GSE328433: annotated scRNA-seq, 15 h APF; converted from Seurat by SeuratDisk | Python | 12,187 | 537 MB | 7 MB | [Open](https://www.cellucid.com/?data=geo:GSE328433) |
| GSE336826: zebrafish embryo at 12 hpf, a MatrixMarket sample with its UMAP and cluster tables |  | 4,278 | 161 MB | 170 MB, once | [Open](https://www.cellucid.com/?data=geo:GSE336826/GSM9843410_12hpf_scRNA_matrix.mtx.gz) |

### Dropbox

Share links to files published with the [Spateo](https://spateo-release.readthedocs.io/)
package's sample data, and an older tutorial file.

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Spateo sample: mouse embryo seqFISH | Python | 17,694 | 29 MB | 2 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fcm3uw8czhz5hu30%2FseqFISH.h5ad%3Fdl%3D1) |
| Spateo sample: mouse brain, bin 60 | Python | 7,765 | 70 MB | 3 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fc5tu4drxda01m0u%2Fmousebrain_bin60.h5ad%3Fdl%3D1) |
| Human pancreas, a tutorial file written by anndata 0.6 (compound `obs` and `obsm`) | Python | 14,693 | 176 MB | 6 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fqj1jlm9w10wmt0u%2Fpancreas.h5ad%3Fdl%3D1) |

### GitHub

From [theislab/cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets),
one small dataset in every form a link opens, a multimodal one as a MuData and
as a Seurat object, a dataset catalog of five, and Cellucid's own sample
catalog.

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Synthetic cell-type islands, 2-D, as a Zarr folder (read through its consolidated metadata) |  | 72 | folder | 13 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr) |
| The same dataset as an `.h5ad` | Python | 72 | 73 KB | 73 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad) |
| The same dataset as a Zarr ZIP |  | 72 | 38 KB | 52 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.zarr.zip) |
| Synthetic CITE-seq, a MuData Zarr folder of two modalities (`rna`, `prot`) |  | 60 | folder | — | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-citeseq.zarr) |
| The same cells as a Seurat object (assays `RNA`, `ADT`) | R | 60 | 8 KB | 8 KB, once | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-citeseq-seurat.rds) |
| A catalog of the five (1-D, 2-D, 3-D, the MuData and the Seurat object; the 1-D and 3-D ones with velocity on their own embedding), opens on the 2-D one |  | 72 | folder | — | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog) |
| Cellucid's own sample catalog, five atlases (opens on Suo et al.) |  | 561,947 | folder | — | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-datasets/anndata) |

### Any web address

| Dataset | Tag | Cells | File | Opening reads |  |
|---|---|---|---|---|---|
| Cellucid's sample catalog, from its GitHub Pages folder |  | 561,947 | folder | — | [Open](https://www.cellucid.com/?data=https%3A%2F%2Ftheislab.github.io%2Fcellucid-datasets%2Fanndata%2F) |
| A CELLxGENE Census file, by its bucket address | Python | 2,346 | 15 MB | 4 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fcellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com%2Fcell-census%2F2025-11-10%2Fh5ads%2Fa49d9109-1d0c-4b36-8139-19aa9a83428c.h5ad) |

## Where the bytes come from

- **CELLxGENE.** CELLxGENE's own site and download links can only be read by
  CELLxGENE itself, so Cellucid reads the copy of each dataset that
  [CZ CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/)
  publishes in its public bucket, from the latest Census release. A dataset
  shows at the version in that release. The browser remembers the release it
  read last and looks for the file there while it reads the current release
  list, so a link usually opens one round trip sooner; when a newer release
  is out, the file is looked up again in it. A dataset published or updated after
  the release cannot be opened from a link: download its `.h5ad` from
  CELLxGENE and open the file with **Open… › File…** instead. CELLxGENE data
  is CC BY 4.0; the Dataset box names the source and links to it.
- **GEO.** The supplementary files of a series or a sample, read from NCBI's
  file host. A series' files are read from two listings fetched together: its
  supplementary folder and its `filelist.txt`, which names every file of the
  series' `_RAW.tar` with its exact size. Each of those files is read from
  its sample's own folder, so a sample opens without the archive and without
  a request to learn its size. A dataset is one file Cellucid opens (an
  `.h5ad`, a `.loom`, an R object, which GEO stores gzipped as `.rds.gz` and
  read once as it downloads) or a **MatrixMarket sample**: a matrix and the
  feature and barcode lists that share its prefix, such as
  `GSM5008737_RNA_3P-matrix.mtx.gz`, `GSM5008737_RNA_3P-features.tsv.gz` and
  `GSM5008737_RNA_3P-barcodes.tsv.gz` (`_genes.tsv.gz` in Cell Ranger 2's
  naming). Lists under a shorter shared prefix serve every matrix below it
  (`GSE1_barcodes.tsv.gz` for `GSE1_RNA_matrix.mtx.gz` and
  `GSE1_ADT_matrix.mtx.gz`), and tables of the sample's cells with the same
  prefix are read with it for the embedding. A series with several datasets
  names each one's link, so you can open one: a sample by its matrix,
  `geo:GSE164378/GSM5008742_HTO_5P-matrix.mtx.gz`. Most GEO count matrices
  have no embedding and are then refused with the lines that compute one
  ({doc}`17_open_r_objects`).
- **Zenodo.** The record's file API, which allows other sites to read files;
  the website's own download links do not. Zenodo limits anonymous requests
  per minute, so the record listing supplies each file's size and checksum
  and no request is spent probing.
- **Figshare.** The item from Figshare's API, which names each file's size,
  checksum and licence; each file from Amazon S3, through the signed link
  Figshare's file host (`ndownloader.figshare.com`) redirects to, which works
  for ten seconds. The page asks the file host for that link, because the file
  host's bot check refuses requests from a browser's background reader in some
  browsers, and the reader reads S3 with it; when its ten seconds are over, one
  one-byte request gets the next. A published version never changes.
- **Google Drive.** Google's Drive API: the only Google address that lets
  another site read a file's bytes (its download pages refuse cross-site
  reads). A public file is read with Cellucid's public API key, a private one,
  or any file you granted to Cellucid in Google's picker, with your sign-in.
  Its first bytes are asked for together with its details, so it opens one
  round trip after the paste. When Google holds back downloads without a
  sign-in from your network address, Cellucid tells it in about eight seconds
  (the time Google takes to answer) and asks you to sign in; a read held back
  after the file opened is told in about nine.
- **OneDrive and SharePoint.** Microsoft Graph names the file and a short-lived
  download address, which is read with plain range requests; when it expires,
  Graph is asked for a new one, once for all the reads waiting on it.
- **Hugging Face.** The Hub's file addresses, which redirect to a CDN that
  allows range reads. The link is pinned to one commit when it opens, so the
  files cannot change under an open dataset.
- **GitHub and any web address.** The host itself. It must allow other sites
  to read its files (CORS). A single file is read in place when the host
  answers range requests; a host that only sends whole files has its file
  downloaded once into memory, with its progress shown, and then opened.

## What to expect

- **While it opens** you see what is happening from the first moment: over
  the empty viewer when nothing is shown yet, or as a notification while
  another dataset stays on screen. It names the link, the step under way
  (*Finding the file on NCBI GEO…*, *Reading the file's structure…*,
  *Reading embedding X_umap…*), how much has arrived, how fast, and how long
  is left when that can be told. A request that failed for a moment and is
  made again shows as *Retrying…* with why. **Cancel** stops it where it is;
  the dataset on screen stays. How it ended (opened, and how long it took; or
  what failed) is kept in the notification center.
- **Opening it again** in the same page (after **Disconnect**, or when a
  session reopens its dataset) asks the provider nothing: what it said the
  first time, and the file's first bytes, are kept until the page closes.
- **Opening** reads the file's structure, its cell metadata, its gene names
  and the default embedding, nothing else. On a typical connection this takes
  a few seconds for a small file and up to about half a minute for a large
  atlas whose metadata is spread out. Reads are 16 KiB pages that grow into
  longer requests when a pass goes forward.
- **A host that compresses a file for transfer** (`Content-Encoding`: NCBI
  does so for GEO's `.rds` files) leaves no part of it readable on its own.
  Cellucid notices, downloads the whole file once, and says so in the
  notification.
- **The first gene** of an `.h5ad` whose matrix is stored by cell (CSR, the
  usual layout) needs the matrix's index read once, because one gene's values
  are spread over all of it. Over a typical internet connection that happens
  when you first show a gene, never when the file opens; from a host that
  delivers bytes as fast as a disk (the browser measures the link while
  opening), it starts in the background as soon as the file opens. Later
  genes are immediate. A file stored by gene (CSC) reads one gene at a time.
- **A file with nothing to plot** (`obsm` empty, as many GEO files are) is
  refused, with what `obsm` holds. A file without a UMAP but with another
  embedding opens on a best guess ({doc}`13_data_mapping`). Add an embedding
  (for example `sc.tl.umap(adata)` in Scanpy), save the file, and open it
  again.
- **A session** saved while linked data is open records the link, so loading
  the session opens the same data again. A link that carries its own access
  key (a Dropbox link) is then in the session file too: share such a session
  only with people who may open the data.

```{figure} ../../../_static/screenshots/data_loading/opening-failed.png
:alt: The panel over the empty viewer after a linked file could not be opened: "The data could not be opened", one sentence saying why (the file was not found at its host, HTTP 404) and what to do, a folded Details, and the buttons Open without this link, focused, and Try again.
:width: 856px

A link given in the page's address (`?data=`) that cannot be opened, here a
file the host does not have: with nothing else to show, the reason stays over
the empty viewer, the error under **Details**, with **Open without this link**
(the app without that address) and **Try again**.
```

## How much has been read

While linked data is open, the **Dataset** box says where it is read from
(its format and host) and how much has been read so far, for example
"5.5 MB read in 43 requests"; the line reads "Reading…" while bytes are
arriving. **Details** in the same box lists the file, its size when the host
states it, and the rest of what the dataset holds. Each read is also written to
the browser console as `[cellucid:transfer]`. A request that fails for a moment
and is made again is said once, in a notification that closes when bytes
arrive again.

Linked data is read in place, so nothing is saved to disk: what was read stays
in memory and is released when the dataset closes. Only a file on a host that
sends whole files is read whole, into memory; nothing is saved to your
computer unless you export a file.

```{figure} ../../../_static/screenshots/data_loading/dataset-box-details.png
:alt: The Dataset box after an .h5ad file was opened from a web address: its name, the cells and genes line, "h5ad file" and the host, how much has been read in how many requests, and Details open, listing File, Format, Size, Read, Cells, Genes, Obs fields, Embeddings and the rest.
:width: 488px

An `.h5ad` read in place from a web address: the Dataset box says how much has
been read in how many requests, and **Details** says how it is read and what
the dataset holds.
```

## Why a link can be refused

```{figure} ../../../_static/screenshots/data_loading/fail-data-link-unrecognised.png
:alt: A Cellucid notification reading "Cellucid does not recognise "not a link"", followed by the kinds of link it opens.
:width: 760px

Text Cellucid does not recognise is refused before anything is requested,
with the kinds of link it opens.
```

| Message | Cause | What to do |
|---|---|---|
| Cellucid does not recognise … | The text is not a link or id it knows | Paste one of the forms in the table above |
| `<host>` did not let this site read … (CORS) … | The host does not let other sites read its files, or could not be reached | Use a host that allows it, or download the file and open it with **Open… › File…** |
| … is not in the Census release … | The CELLxGENE dataset is newer than the Census copy | Download it from CELLxGENE and open the file |
| A CELLxGENE collection holds several datasets … | You pasted a collection link | Open the collection, choose a dataset, paste its **Explore** link |
| GEO … ships its data only gzip-compressed | A gzipped `.h5ad` or other range-read file cannot be read in place (a gzipped R object can: `.rds.gz` is read once as it downloads) | Download and decompress it, then open the decompressed file |
| GEO … has `GSM…_matrix.mtx.gz`, but no `GSM…_barcodes.tsv.gz` for it … | The sample's feature or barcode list is missing from GEO, or under another name | Download the files from the GEO page, put the three side by side, and open the folder with **Open… › Folder…** |
| GEO … has `GSM…_features.tsv.gz` and `GSM…_genes.tsv.gz` for … | Two feature lists share the matrix's prefix; which names its rows cannot be told from the names | Download the matrix, its barcodes and the list you mean into a folder of their own and open it with **Open… › Folder…** |
| … is not a dataset Cellucid opens on its own | The link names a sample's feature or barcode list | Open the sample by its matrix (the message names the link) |
| … is a 10x Cell Ranger count matrix that has no cell embedding … | A Cell Ranger `.h5` holds counts only; its projections are in the `outs/` folder beside it | Download the `outs/` folder and open it with **Open… › Folder…**, or compute an embedding ({doc}`18_open_loom_and_10x`) |
| Cellucid found nothing to plot … | The file has no embedding | Add one, as above |
| This is a private Figshare link (figshare.com/s/…) … | Figshare offers no way for a browser to read a private link | Publish the item, or download the file and open it with **Open… › File…** |
| Figshare's file host refused this browser with a bot check … | Figshare's file host answered with an AWS WAF challenge instead of the file; a page cannot answer it for another site (seen with automated, headless browsers) | Try again in a minute or in another browser, or download the file from Figshare and open it with **Open… › File…** |
| This is a Figshare file download link … | A bare `ndownloader` link does not name the file's size and name | Paste the item's page or DOI |
| The DOI … names a record on …, which is not a Figshare item | A DOI of another repository | Paste a Zenodo or Figshare DOI, or a link to the file |
| This link is a Google Drive folder / a OneDrive folder … | Cellucid opens single files from Drive and OneDrive | Share the file itself (an `.h5ad`, a `.zarr.zip`, or another file Cellucid reads; {doc}`16_private_drive_onedrive`) |
| This is the address of a OneDrive page, not a link to a file … | You pasted the address of OneDrive's file list, or of a file open in OneDrive | In OneDrive, select the file, **Share**, **Copy link**, and paste that link |
| Google Drive sends this file only with a Google sign-in for now / Google Drive is not sending this file to this browser without a Google sign-in for now … | Google holds back downloads without a sign-in from a network that made many of them, for a while | **Sign in to Google** and **Select** the file in the picker, or open it again later ({doc}`16_private_drive_onedrive`) |
| This file is private on Google Drive, or the link is wrong … not configured on this site yet | A private Drive file on a site where Google sign-in is not set up (a copy of the app hosted elsewhere) | Open the link on www.cellucid.com, or download the file and open it with **Open… › File…** |
| Google Drive is not configured on this site yet / Signing in to Google (Microsoft) is not configured on this site yet | The site has no Google or Microsoft registration: a copy of the app hosted elsewhere | Open the link on www.cellucid.com, or download the file and open it with **Open… › File…** |
| Google Drive works on a local cellucid serve only on its default port, 8765 … / in a notebook viewer only on port 8765 … / on www.cellucid.com, not on … | Google answers only registered addresses, port included, for public and private files | Start `cellucid serve` without `--port`; for a notebook viewer, stop the server holding 8765 (an earlier viewer, a `cellucid serve`) and show the data again; or use www.cellucid.com |
| Microsoft sign-in needs https:// … / cannot reach Cellucid inside another site's page … | A plain `http://` page on another computer, or a viewer embedded in another site | Open Cellucid in its own tab, such as www.cellucid.com |
| Your Google account cannot open this file / Your Microsoft account cannot open this file … | The file is not shared with the account you signed in with | Ask its owner to share it with that account |
| Your organisation must approve Cellucid … | A work or school account opening a colleague's file, in an organisation that lets only administrators grant `Files.Read.All` | Ask your administrator ({doc}`16_private_drive_onedrive`) |
| Google has limited downloads of this file for now … | Google's per-file download quota | Try later, or ask the owner to share a copy |
| … has 2 files Cellucid opens. Open one of them: … | A record, series or repository holds several openable files | Paste one of the links the message lists |
| … holds no datasets.json, so it is not a dataset catalog | A folder link that points at a web page or the wrong folder | Link to the folder that holds `datasets.json`, to a Zarr store folder, or to one file |
| … answers, but it is not a Cellucid server … | A server address (a port, or `localhost:8000/…` without `http://`) where a web host serves a folder, not `cellucid serve` | Type the address of the folder's `datasets.json`, or start the server with `cellucid serve` |
| Invalid datasets.json: the catalog must contain exactly: datasets, default, schema (missing: schema) | A `datasets.json` without `"schema"` | Its publisher writes `"schema": 2`, `"default"` and `"datasets"` |
| This .h5ad file uses something the browser cannot read: … filter … | A compression filter the browser has no decoder for (Blosc, Zstandard, LZ4, bitshuffle, bzip2, SZIP) | Save it again with gzip: `adata.write_h5ad(new_path, compression="gzip")` |

Each refusal can be seen on real published data; paste one of these into the
**Link or server** bar:

- several openable files: `zenodo:15713160`
- nothing to plot: `geo:GSE307135`
- LZF compression: `zenodo:7838976/GSM3768152_Brain_8_dge_asStruCluster.h5ad`
- newer than the Census release: `cellxgene:8f3cb9c4-38d5-4b9c-b309-bcccd76d0075`

(open-from-link-checking)=
## Checking the links on this page

`docs/_tooling/data-links/check.mjs` in the Python repository opens every
`https://www.cellucid.com/?data=…` link the documentation publishes, one at a
time in a real browser, and prints the dataset each one opened (name, cells,
what was read) or the message that refused it:

```sh
node docs/_tooling/data-links/check.mjs                              # against www.cellucid.com
node docs/_tooling/data-links/check.mjs --app http://127.0.0.1:4173  # against a local build
node docs/_tooling/data-links/check.mjs --link geo:GSE328433          # a candidate, before publishing it
```

## Related pages

- {doc}`16_private_drive_onedrive`: Google Drive, OneDrive, and sharing private data with a collaborator.
- {doc}`11_custom_dataset_repository`: publish your own datasets on GitHub.
- {doc}`03_browser_file_picker_tutorial`: open a file from your computer.
- {doc}`13_data_mapping`: choose which embeddings, genes and graph a file shows.
