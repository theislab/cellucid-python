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
of an `.h5ad` file, a Zarr store folder, a Zarr store zipped into a `.zip`, or
a dataset catalog folder (one holding `datasets.json`). Nothing needs a server
of your own, and only a private Google Drive file or a OneDrive file needs a
sign-in ({doc}`16_private_drive_onedrive`).

It is in the **Session** panel under **Open data**. Paste a link or an id,
then click **Open** (or press Enter). The page address gains a `?data=`
parameter naming what you opened, so copying the address (or using **Share**)
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
| A provider's form: `hf:user/repo`, `GSE123456`, a DOI, `zenodo:<record>`, `owner/repo/path`, a CELLxGENE, Zenodo, Hugging Face, Dropbox or GitHub link | That provider's data ([What you can paste](#what-you-can-paste)) |
| A web address ending in `.h5ad`, `.zip`, `.zarr` or `datasets.json`, or one with a `?query` | That file, store or catalog |
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
**Disconnect** brings the samples back.

## What you can paste

| You paste | It opens | `?data=` form |
|---|---|---|
| A CELLxGENE dataset id, or its **Explore** link (`cellxgene.cziscience.com/e/<dataset id>.cxg/`) | That dataset's `.h5ad` | `cellxgene:<dataset id>` |
| A GEO series (`GSE123456`), its accession page, or a supplementary file link | The series' `.h5ad` supplementary file | `geo:GSE123456` or `geo:GSE123456/<file>.h5ad` |
| A Zenodo record link, file link, DOI (`10.5281/zenodo.<id>`) or `zenodo:<record>` | The record's one `.h5ad` or `.zip` (or the file the link names); a concept DOI opens the latest version | `zenodo:<record>[/<file>]` |
| A Figshare item page (`figshare.com/articles/<type>/<name>/<id>[/<version>]`, also on an institution's `*.figshare.com` portal), a Figshare DOI (`10.6084/m9.figshare.<id>[.v<n>]`), any DOI that resolves to a Figshare item, or `figshare:<id>` | The item's one `.h5ad` or `.zip` (or the file the link names, by name or by `?file=<file id>`); a link without a version opens the latest | `figshare:<id>[v<version>][/<file>]` |
| A Google Drive file link (`drive.google.com/file/d/<id>/…`, `…/open?id=<id>`, `…/uc?id=<id>`) | That file, read through Google's Drive API; a private file after **Sign in to Google** ({doc}`16_private_drive_onedrive`); folder links are refused | `gdrive:<id>` |
| A OneDrive or SharePoint file share link (`1drv.ms/…`, `onedrive.live.com/…`, `<tenant>-my.sharepoint.com/:u:/…`) | That file, after **Sign in with Microsoft** ({doc}`16_private_drive_onedrive`); folder links are refused | `onedrive:<drive>/<item>` once opened |
| A Hugging Face dataset link (`huggingface.co/datasets/<user>/<repo>/…`) or `hf:<user>/<repo>[@<revision>][/<path>]` | A dataset catalog folder, or one `.h5ad` / `.zip`; pinned to the commit the branch points at | `hf:<user>/<repo>[@<revision>][/<path>]` |
| A Dropbox file link | That file, read from Dropbox's content host (folder links need a Dropbox sign-in and are refused) | `dropbox:<link>` |
| A GitHub repository path (`owner/repo/folder`, read from the default branch, or `owner/repo@branch/folder`), a `github.com/.../tree/...`, `.../blob/...` or `.../raw/...` link, or a `raw.githubusercontent.com` link (the **Raw** button's `refs/heads/<branch>` form included) | The dataset catalog in that folder ({doc}`11_custom_dataset_repository`), the Zarr store folder (`….zarr`) or the one `.h5ad` or Zarr ZIP the path names ({doc}`15_host_on_github`); a file stored with Git LFS is read from GitHub's media host | `gh:owner/repo/folder`, `gh:owner/repo/path/atlas.zarr`, `gh:owner/repo/path/file.h5ad` |
| A web address ending in `/`, or its `datasets.json` | The dataset catalog in that folder ({doc}`07_file_format_expectations`) | the address itself |
| A web address of a Zarr store folder (`…/name.zarr/`, Zarr v2 or v3, with consolidated metadata as anndata writes it): a bucket, a lab server, GitHub raw or a Hugging Face `resolve` folder | That store, its chunks read as views need them | the address itself |
| A web address of an `.h5ad` file or a Zarr `.zip`: a bucket, a lab server, a Hugging Face file link (`huggingface.co/datasets/<user>/<repo>/resolve/main/<file>`), or a Zenodo API file link (`zenodo.org/api/records/<id>/files/<file>/content`; Zenodo's ordinary download link cannot be read by other sites) | That file | the address itself |

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
view, measured on 2026-10-22; the rest is read only when a view needs it
(showing a gene, switching an embedding). Hosts change: if a link stops
working, the [link checker](#open-from-link-checking) says which.

### CELLxGENE

Census copies of the 2025-11-10 release. An **Explore** link works too:
[this one](https://www.cellucid.com/?data=https%3A%2F%2Fcellxgene.cziscience.com%2Fe%2Fa49d9109-1d0c-4b36-8139-19aa9a83428c.cxg%2F)
opens the first row.

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Healthy human liver: lymphocytes | 2,346 | 15 MB | 4 MB | [Open](https://www.cellucid.com/?data=cellxgene:a49d9109-1d0c-4b36-8139-19aa9a83428c) |
| Tabula Sapiens: neural | 2,685 | 122 MB | 11 MB | [Open](https://www.cellucid.com/?data=cellxgene:6e4acdae-59ef-452b-a6d4-260cba368586) |
| Human Brain Cell Atlas: hippocampus, rostral CA3 | 8,077 | 106 MB | 7 MB | [Open](https://www.cellucid.com/?data=cellxgene:d5452b83-7c3d-4d7c-ab7a-c7fece7196c5) |
| Human Immune Health Atlas: other cells | 10,937 | 56 MB | 10 MB | [Open](https://www.cellucid.com/?data=cellxgene:81ffdcd6-3f04-40bc-8bb0-9748f4ffe576) |
| Tabula Sapiens: mammary | 30,936 | 1.2 GB | 18 MB | [Open](https://www.cellucid.com/?data=cellxgene:2ba40233-8576-4dec-a5f1-2adfa115e2dc) |
| Foveal and parafoveal retina | 34,323 | 252 MB | 8 MB | [Open](https://www.cellucid.com/?data=cellxgene:f5be9ed2-5d41-4a52-91e4-4ff24ff84900) |
| Tabula Sapiens: thymus | 42,729 | 1.6 GB | 26 MB | [Open](https://www.cellucid.com/?data=cellxgene:0ced5e76-6040-47ff-8a72-93847965afc0) |
| Krasnow Lab Human Lung Cell Atlas, 10x | 65,662 | 569 MB | 12 MB | [Open](https://www.cellucid.com/?data=cellxgene:8c42cfd0-0b0a-46d5-910c-fc833d83c45e) |
| Pancreas, harmonized across studies | 121,916 | 915 MB | 14 MB | [Open](https://www.cellucid.com/?data=cellxgene:31f657dc-1875-4c4b-a5ca-ce63b3ef3a82) |

### Zenodo

A record with one openable file opens it; a record with several names them,
and the link names one with `zenodo:<record>/<file>`.

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Day 3 Session 4: processed Scanpy AnnData (record 13255338) | 624 | 8 MB | 2 MB | [Open](https://www.cellucid.com/?data=zenodo:13255338/processedScanpyAnndata.h5ad) |
| Uninfected fibroblasts, RNA velocity (record 7608772) | 1,545 | 193 MB | 4 MB | [Open](https://www.cellucid.com/?data=zenodo:7608772/uninfected_fibroblast.h5ad) |
| AML atlas (record 15713160) | 1,732 | 1.2 GB | 4 MB | [Open](https://www.cellucid.com/?data=zenodo:15713160/AML_Atlas.h5ad) |
| Unwounded skin fibroblasts, integrated (record 10198892) | 7,559 | 105 MB | 8 MB | [Open](https://www.cellucid.com/?data=zenodo:10198892/unwounded_fibroblast_merged.h5ad) |
| Bradyzoite heterogeneity, with a 3-D embedding (record 23288289) | 8,057 | 15 MB | 2 MB | [Open](https://www.cellucid.com/?data=zenodo:23288289) |
| Ankylosing spondylitis T cells, TCR-DeepInsight (record 12741480) | 15,065 | 22 MB | 7 MB | [Open](https://www.cellucid.com/?data=zenodo:12741480/Yi_2023_Ankylosing_Spondylitis.h5ad) |
| Kang et al. 2018 PBMC, IFN-β (record 10850398) | 24,562 | 264 MB | 7 MB | [Open](https://www.cellucid.com/?data=zenodo:10850398/kang18_counts_25k.h5ad) |
| Liao et al. 2020 COVID-19 BALF (record 10850398) | 65,792 | 828 MB | 25 MB | [Open](https://www.cellucid.com/?data=zenodo:10850398/liao20_sub.h5ad) |
| Human embryo at Carnegie stage 9, annotated; written by anndata 0.7 (record 14730926) | 2,150 | 17 MB | 3 MB | [Open](https://www.cellucid.com/?data=zenodo:14730926/cs9_scRNA_annotation.h5ad) |

### Figshare

Items published on figshare.com, each holding one `.h5ad`. The file is read
from Amazon S3 through the signed link Figshare's file host gives the page,
which works for ten seconds; within those seconds every read goes straight to
S3.

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| scverse tutorial data: PBMC 3k, processed (item 22577536) | 2,638 | 100 MB | 3 MB | [Open](https://www.cellucid.com/?data=figshare:22577536) |
| Scanpy tutorial: PBMC 3k (item 16447278) | 2,638 | 32 MB | 4 MB | [Open](https://www.cellucid.com/?data=figshare:16447278) |
| Scanpy tutorial: trajectory, version 2, by its DOI (item 16540047) | 2,730 | 12 MB | 2 MB | [Open](https://www.cellucid.com/?data=10.6084%2Fm9.figshare.16540047.v2) |
| Squidpy tutorial: Visium mouse brain (item 16566057) | 704 | 28 MB | 11 MB | [Open](https://www.cellucid.com/?data=figshare:16566057) |

### Google Drive and OneDrive

A Google Drive file shared with "Anyone with the link" opens with no sign-in;
a private one, and every OneDrive file, after one sign-in. See
{doc}`16_private_drive_onedrive`.

### Hugging Face

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Kang et al. IFN-β PBMC counts | 8,000 | 38 MB | 2 MB | [Open](https://www.cellucid.com/?data=hf:Nikita1113/tb-science-scrna-kang-ifnb) |
| PBMC, preprocessed | 17,041 | 431 MB | 10 MB | [Open](https://www.cellucid.com/?data=hf:datavil/pbmc3k) |

### GEO

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| GSE350831: Flex pilot, annotated | 21,751 | 868 MB | 9 MB | [Open](https://www.cellucid.com/?data=geo:GSE350831/GSE350831_Flex_pilot_paper_annotated_bbknn.h5ad) |
| GSE327505: Saber et al. 2026, processed | 25,651 | 2.1 GB | 12 MB | [Open](https://www.cellucid.com/?data=geo:GSE327505) |
| GSE328433: annotated scRNA-seq, 15 h APF; converted from Seurat by SeuratDisk | 12,187 | 537 MB | 7 MB | [Open](https://www.cellucid.com/?data=geo:GSE328433) |

### Dropbox

Share links to files published with the [Spateo](https://spateo-release.readthedocs.io/)
package's sample data, and an older tutorial file.

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Spateo sample: mouse embryo seqFISH | 17,694 | 29 MB | 2 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fcm3uw8czhz5hu30%2FseqFISH.h5ad%3Fdl%3D1) |
| Spateo sample: mouse brain, bin 60 | 7,765 | 70 MB | 3 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fc5tu4drxda01m0u%2Fmousebrain_bin60.h5ad%3Fdl%3D1) |
| Human pancreas, a tutorial file written by anndata 0.6 (compound `obs` and `obsm`) | 14,693 | 176 MB | 6 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fwww.dropbox.com%2Fs%2Fqj1jlm9w10wmt0u%2Fpancreas.h5ad%3Fdl%3D1) |

### GitHub

From [theislab/cellucid-demo-custom-datasets](https://github.com/theislab/cellucid-demo-custom-datasets),
one small dataset in every form a link opens, a dataset catalog of three, and
Cellucid's own sample catalog.

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Synthetic cell-type islands, 2-D, as a Zarr folder (read through its consolidated metadata) | 72 | folder | 13 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr) |
| The same dataset as an `.h5ad` | 72 | 73 KB | 73 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad) |
| The same dataset as a Zarr ZIP | 72 | 38 KB | 52 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.zarr.zip) |
| A catalog of three Zarr folders (1-D, 2-D, 3-D; the 1-D and 3-D ones with velocity on their own embedding), opens on the 2-D one | 72 | folder | 13 KB | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog) |
| Cellucid's own sample catalog, five atlases (opens on Suo et al.) | 561,947 | folder | — | [Open](https://www.cellucid.com/?data=gh:theislab/cellucid-datasets/anndata) |

### Any web address

| Dataset | Cells | File | Opening reads |  |
|---|---|---|---|---|
| Cellucid's sample catalog, from its GitHub Pages folder | 561,947 | folder | — | [Open](https://www.cellucid.com/?data=https%3A%2F%2Ftheislab.github.io%2Fcellucid-datasets%2Fanndata%2F) |
| A CELLxGENE Census file, by its bucket address | 2,346 | 15 MB | 4 MB | [Open](https://www.cellucid.com/?data=https%3A%2F%2Fcellxgene-census-public-us-west-2.s3.us-west-2.amazonaws.com%2Fcell-census%2F2025-11-10%2Fh5ads%2Fa49d9109-1d0c-4b36-8139-19aa9a83428c.h5ad) |

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
- **GEO.** The `.h5ad` supplementary files of a series, read from NCBI's file
  host. A series with several `.h5ad` files names them, so you can open one.
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
  reads). A public file is read with Cellucid's public API key, a private one
  with your sign-in. Its first bytes are asked for together with its details,
  so it opens one round trip after the paste.
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
  atlas whose metadata is spread out. Reads are 256 KiB pages that grow into
  longer requests when a pass goes forward.
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

## Why a link can be refused

| Message | Cause | What to do |
|---|---|---|
| Cellucid does not recognise … | The text is not a link or id it knows | Paste one of the forms in the table above |
| `<host>` did not let this site read … (CORS) … | The host does not let other sites read its files, or could not be reached | Use a host that allows it, or download the file and open it with **Open… › File…** |
| … is not in the Census release … | The CELLxGENE dataset is newer than the Census copy | Download it from CELLxGENE and open the file |
| A CELLxGENE collection holds several datasets … | You pasted a collection link | Open the collection, choose a dataset, paste its **Explore** link |
| GEO series … ships its AnnData only gzip-compressed | A compressed file cannot be read in place | Download and decompress it, then open the `.h5ad` file |
| Cellucid found nothing to plot … | The file has no embedding | Add one, as above |
| This is a private Figshare link (figshare.com/s/…) … | Figshare offers no way for a browser to read a private link | Publish the item, or download the file and open it with **Open… › File…** |
| Figshare's file host refused this browser with a bot check … | Figshare's file host answered with an AWS WAF challenge instead of the file; a page cannot answer it for another site (seen with automated, headless browsers) | Try again in a minute or in another browser, or download the file from Figshare and open it with **Open… › File…** |
| This is a Figshare file download link … | A bare `ndownloader` link does not name the file's size and name | Paste the item's page or DOI |
| The DOI … names a record on …, which is not a Figshare item | A DOI of another repository | Paste a Zenodo or Figshare DOI, or a link to the file |
| This link is a Google Drive folder / a OneDrive folder … | Cellucid opens single files from Drive and OneDrive | Share the `.h5ad` or `.zarr.zip` file itself ({doc}`16_private_drive_onedrive`) |
| This file is private on Google Drive, or the link is wrong … not configured on this site yet | A private Drive file on a site where Google sign-in is not set up | Open it on www.cellucid.com |
| Google Drive works on a local cellucid serve only on its default port, 8765 … / in a notebook viewer only on port 8765 … / on www.cellucid.com, not on … | Google answers only registered addresses, port included, for public and private files | Start `cellucid serve` without `--port`, create the notebook viewer with `AnnDataViewer(data, port=8765)`, or use www.cellucid.com |
| Microsoft sign-in needs https:// … / cannot reach Cellucid inside another site's page … | A plain `http://` page on another computer, or a viewer embedded in another site | Open Cellucid in its own tab, such as www.cellucid.com |
| Your Google account cannot open this file / Your Microsoft account cannot open this file … | The file is not shared with the account you signed in with | Ask its owner to share it with that account |
| Your organisation must approve Cellucid … | A work or school account opening a colleague's file, in an organisation that lets only administrators grant `Files.Read.All` | Ask your administrator ({doc}`16_private_drive_onedrive`) |
| Google has limited downloads of this file for now … | Google's per-file download quota | Try later, or ask the owner to share a copy |
| … has 2 files Cellucid opens. Open one of them: … | A record, series or repository holds several openable files | Paste one of the links the message lists |
| … holds no datasets.json, so it is not a dataset catalog | A folder link that points at a web page or the wrong folder | Link to the folder that holds `datasets.json`, to a Zarr store folder, or to one file |
| Invalid datasets.json: the catalog is missing 'schema' … | A `datasets.json` without `"schema"` | Its publisher writes `"schema": 2`, `"default"` and `"datasets"` |
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
