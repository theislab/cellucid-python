# Host a Store for Sharing

**Audience:** anyone sharing datasets (lab leads, computational users, educators)
**Time:** 15–60 minutes (depends on upload size)
**Goal:** publish a store written by `cellucid_prepare()` so others open it with
one link.

A Cellucid-ready AnnData store needs no server of its own: any host that serves
files over HTTPS will do, because the web app reads the store in place with
range requests. You share one link; the person who opens it needs nothing but a
browser.

```{warning}
Publishing a store on a public host is publishing your data. Treat it like
making the dataset public. For data that must stay private, keep it on your own
machine and open the `.zarr.zip` with the **Zarr ZIP** button instead.
```

## Folder or single file?

| Write | When |
| --- | --- |
| `out_path = "my_dataset.zarr", format = "zarr"` (a folder) | GitHub (a repository or its Pages site): no file over 100 MB, no LFS, and GitHub answers many small files at once over HTTP/2. |
| `out_path = "my_dataset.zarr.zip", format = "zarr.zip"` (one file) | Everything else: a Zenodo or Figshare record, a Hugging Face file, a Dropbox file, a bucket, a lab web server, a local file, `cellucid serve` and Jupyter. The ZIP's entries are stored uncompressed, so any part of it is one range request, and ranges that lie close together are fetched as one; most such hosts answer six requests at a time over HTTP/1.1, where a folder's one request per chunk is much slower ({doc}`../c_data_preparation_api/index`). |

Both hold the same files, byte for byte.

## GitHub

1. Write the store into your repository, for example
   `data/my_dataset.zarr` (folder) or `data/my_dataset.zarr.zip` (file).
2. Commit and push it. A folder store keeps every file small (one chunk holds
   at most 262,144 values, a few MiB before compression), which suits GitHub's
   file size limit; a `.zarr.zip` larger than 100 MB needs Git LFS, which the web app
   also reads.
3. In the web app, paste into **Open from link**
   `owner/repo/data/my_dataset.zarr` (the folder) or
   `owner/repo/data/my_dataset.zarr.zip` (the file). Without `@branch` this
   reads the repository's default branch; `owner/repo@branch/data/...` reads
   another.
4. To serve the folder from GitHub Pages as well, commit an empty `.nojekyll`
   file at the root of the published branch: without it Pages drops the
   store's `.zgroup`, `.zattrs` and `.zmetadata` files and its `_index`
   folders. {doc}`../../web_app/b_data_loading/15_host_on_github` has the
   details and measured speeds.

To publish several datasets from one repository, with a default and a list the
viewer shows, follow {doc}`../../web_app/b_data_loading/11_custom_dataset_repository`.
A complete example repository, a catalog of three Zarr folders, is opened with
this connection value:

```text
theislab/cellucid-demo-custom-datasets/catalog
```

## Hugging Face, Zenodo, Dropbox, a bucket, a lab server

- **Hugging Face:** upload the `.zarr.zip` to a dataset repository and paste
  its `huggingface.co/datasets/<user>/<repo>/…` link.
- **Zenodo:** upload the `.zarr.zip` to a record and paste the record link or
  its DOI.
- **Dropbox:** share the `.zarr.zip` file and paste the file link.
- **A bucket or a web server:** serve the `.zarr.zip` file over HTTPS, with
  CORS headers that let other sites read its files and with range requests
  answered, then paste the address.

Every accepted form, and what each opens, is listed in
{doc}`../../web_app/b_data_loading/14_open_from_link`.

## Share the link

After the store opens, the page address carries it in `?data=…`. Copy that
address (or use **Share**) and send it: whoever opens it sees the same dataset,
with the launch defaults you stored in `uns['cellucid']`.

## Checklist before you share

- The store opens in the web app from its public address (try it in a private
  window).
- `dataset_name`, `dataset_description`, `citation_text`/`citation_url` and
  `source_name`/`source_url` say what the dataset is and how to cite it.
- `dataset_id` is stable: sessions and community annotations are keyed on it, so
  keep it the same when you replace the store with a new version.

## Troubleshooting

- “Nothing is published at that address” → the link points at the wrong path or
  branch; open it in a browser tab first.
- “CORS” or “the server does not allow range requests” → the host is not
  configured as above; use GitHub, Hugging Face or Zenodo, or fix the server's
  headers.
- The store opens but looks wrong → {doc}`03_validate_stores_and_debug_loading`.
