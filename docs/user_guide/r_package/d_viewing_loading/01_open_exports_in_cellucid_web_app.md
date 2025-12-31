# Open Exports in the Cellucid Web App

**Audience:** everyone  
**Time:** 5–15 minutes  
**Goal:** load a local exported dataset folder into the Cellucid web app.

This is the simplest workflow:
- export a dataset folder with `cellucid_prepare(out_dir=...)`
- load it locally via the browser file picker

```{note}
The web app has a full file-picker tutorial with screenshots here:
{doc}`../../web_app/b_data_loading/03_browser_file_picker_tutorial`.
This page is the short R-focused checklist.
```

## Prerequisites

You have an export folder that contains at least:
- `dataset_identity.json`
- `obs_manifest.json`
- `points_2d.bin` or `points_3d.bin` (optionally `.gz`)

Quick check:

```r
list.files(out_dir)
```

## Step-by-step checklist

### 1) Open the web app

Open the Cellucid web app in a desktop browser.

```{tip}
If you have trouble loading large datasets, try Chrome or Firefox first.
```

### 2) Find the “load dataset” controls

In the left sidebar, look for the data loading panel (often called something like “Dataset Connections”).

### 3) Choose the local folder picker

Select the option that lets you browse local data / pick a folder.

### 4) Select your export folder (`out_dir`)

Select the folder you exported (the folder that contains `dataset_identity.json`).

```{warning}
Some browsers show a file picker that looks like “select a file”, but the correct choice is the **folder** that contains the dataset files.
If you select a file inside the folder, Cellucid may not be able to load the rest of the dataset.
```

### 5) Wait for the dataset to load

For small/medium exports this should be quick.
For large exports (especially with many genes), initial directory scanning may take longer.

## What success looks like

After loading, you should be able to:

- see points in the canvas
- color by an `obs` field
- (if you exported expression) search for a gene from your exported gene set and color by it

## If loading is slow or fails

Start with:
- {doc}`03_validate_exports_and_debug_loading`
- {doc}`../i_troubleshooting_index/04_web_app_loading_issues`

And the web app’s broader data-loading troubleshooting:
- {doc}`../../web_app/b_data_loading/08_troubleshooting_data_loading`
