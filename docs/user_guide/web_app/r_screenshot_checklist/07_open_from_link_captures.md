# Open from link, Transfers, and the Session panel: captures to take

Working material, not published. Prepared 2026-10-22 together with the
scenarios in `docs/_tooling/screenshots/scenarios.mjs` (`buildDataLinkScenarios`).
Every scenario below was dry-run with `--out` into a scratch directory and
settles; none has been emitted. Read
{doc}`00_capture_tooling_and_conventions` first.

Every link in these scenarios points at a loopback host of the capture run (the
docs catalog host serving `cellucid-datasets/anndata/`, or the web repository's permissive-CORS mirror), never the
internet, so the images are reproducible offline.

## Run

```sh
cd cellucid-python/docs/_tooling/screenshots
node capture.mjs shoot loading-session-panel loading-link-control-typed \
  loading-link-catalog-connected loading-link-file-transfers \
  loading-link-refusal-unrecognised loading-link-refusal-not-found --emit
```

`shoot` prints each `{figure}` block with its `:width:`. Paste it where the
table below says, with the caption given here, **in the same change** — an
orphan PNG fails `python -m pytest` (the `referenced == screenshots` contract).
Then delete the two images this set replaces, with their `captures.json`
records, and run `python -m pytest tests/test_docs_current_api_contract.py` and
`sphinx-build -W`.

## The images

| Scenario | File (`_static/screenshots/data_loading/`) | Page and place | Caption |
| --- | --- | --- | --- |
| `loading-session-panel` | `data-loading-session-panel.png` (**re-captured in place**, same name) | the seven pages that already show it; nothing to paste | The figure blocks stay; check each page's alt text still matches (**Open from link** below **Server**, **Session file** with **Save** and **Load**). |
| `loading-link-control-typed` | `open-from-link-control.png` (new) | `b_data_loading/14_open_from_link.md`, end of **At a glance** | The **Open from link** control with a Zenodo record typed in. **Connect** opens it; the page address then carries `?data=zenodo:23288289`. |
| `loading-link-catalog-connected` | `connect-data-link.png` (new; **replaces** `connect-github-catalog.png`) | `14_open_from_link.md`, after the "What you can paste" table; and in place of `connect-github-catalog.png` in `02_local_demo_tutorial.md`, `09_screenshots.md` and `11_custom_dataset_repository.md` | A dataset catalog folder (`datasets.json` and the AnnData files it lists) opened from a link: **Disconnect** closes it, and **Transfers** names where its files are read from, how much has come in, and that nothing is saved to disk. |
| `loading-link-file-transfers` | `transfers-file-in-place.png` (new) | `14_open_from_link.md`, section **Transfers** | One `.h5ad` read in place from a web address: the bytes read of the file's size, in how many requests, why, and where they go. |
| `loading-link-refusal-unrecognised` | `fail-data-link-unrecognised.png` (new) | `14_open_from_link.md`, section **Why a link can be refused**, before the table | Text Cellucid does not recognise is refused before anything is requested, with the kinds of link it opens. |
| `loading-link-refusal-not-found` | `fail-data-link-not-found.png` (new; **replaces** `fail-github-catalog-not-found.png`) | `14_open_from_link.md`, same section, after the table; and in place of `fail-github-catalog-not-found.png` in `02_local_demo_tutorial.md`, `08_troubleshooting_data_loading.md`, `09_screenshots.md` and `11_custom_dataset_repository.md` | A web address that is not a dataset catalog folder: the notification says it holds no `datasets.json` and what a folder link must point at. The dataset that was open stays open. |

## What changed in the tool for this set

- `capture.mjs` passes `urls.sampleOrigin` (the CORS mirror that answers range
  requests) to scenarios, and a scenario's `query` may be a function of `urls`
  (a `?data=` link to a loopback host).
- A startup that names its data exactly (`?data=`) keeps the welcome overlay
  hidden; the tool now waits for the overlay only until startup shows it or is
  ready without it (`keepWelcome` still throws if it never opens).
- `PANCREAS_READY` and the other count waits read "Showing all 3,696 **cells**";
  the app has said "cells" since the panel was renamed, so every Pancreas
  scenario failed at its first step before this fix.
