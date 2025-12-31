# Metadata and provenance

**Audience:** computational users + anyone preparing manuscripts (reproducibility)  
**Time:** 10–20 minutes  
**What you’ll learn:**
- What provenance/metadata is attached to exported figures (if supported)
- What information you should record even if it is not embedded automatically
- A template “methods” sentence/block for papers

**Prerequisites:**
- A dataset loaded
- You can export at least one figure (see {doc}`02_export_ui_walkthrough`)

---

## What “provenance” means for a Cellucid figure

Provenance is the answer to:

> “What exactly am I looking at, and how could I recreate it?”

For a Cellucid figure, that usually means recording:

**Data identity**
- dataset name and dataset id (or an exported folder version tag)
- how the dataset was loaded (local folder, server, URL, notebook)
- *ideally* a version/commit hash for the dataset export, if you maintain one

**View identity**
- which view was exported (live view vs a snapshot view, or multiview grid)
- the dimension (1D/2D/3D)
- camera pose/projection (implicitly captured by the export, but you still want the saved state)

**Scientific state**
- active color-by field (and whether it is categorical vs continuous)
- legend mapping (category colors or colormap)
- filters/visibility (what was hidden at export time)
- highlights/selection emphasis (if used)

**Export settings**
- plot size (W×H)
- format (SVG vs PNG)
- DPI (for PNG)
- large dataset strategy (full/optimized/hybrid/raster)
- include axes/legend/title/background
- frame export crop rectangle (if used)

:::{tip}
If you saved a session bundle right before exporting, you already captured most of the view/scientific state. The remaining missing piece is export settings (which are not stored in sessions).
:::

---

## What metadata is embedded into exports (must be explicit)

Cellucid exports embed metadata by default (there is no “metadata off” toggle in the current UI).

### PNG metadata (tEXt chunks)

PNG exports include standard `tEXt` fields such as:
- `Software`: `Cellucid (cellucid.com)`
- `Website`: `https://cellucid.com`
- `Creation Time`: ISO timestamp
- `Dataset` and `Dataset ID` (when known)
- `Color Field` (field key)
- `Source File` (may be a URL or a local path, depending on how the dataset was loaded)
- `Description` (human-readable summary)
- `Comment`: a compact JSON blob that includes structured provenance:
  - dataset (name/id/source URLs/user path/citation if present)
  - view (id/label)
  - field (key/kind)
  - filters summary
  - export settings (format/size/DPI/strategy/legend/axes/background/crop)

### SVG metadata (RDF / Dublin Core + Cellucid JSON)

SVG exports include a `<metadata>` block with:
- Dublin Core fields like creator/publisher/date/source/description
- Cellucid-specific tags including:
  - the active color field
  - the source file (URL or local path)
  - a JSON provenance blob similar to the PNG `Comment`

### Filenames (also part of provenance)

Exports are named conservatively:
- `<dataset>_<color-field>_<view>_<timestamp>.svg`
- `<dataset>_<color-field>_<view>_dpi300_<timestamp>.png` (when exporting multiple DPIs)

This makes it easier to tell “what is this file?” even without opening it.

:::{important}
The embedded metadata is useful, but it does not replace version control. For long-term reproducibility, also record:
- the Cellucid app version/commit (if you have it), and
- the dataset export version you used.
:::

---

## Recommended “methods text” template (for papers)

Use this as a starting point and adapt to your norms/journal requirements:

> “Figures were generated in Cellucid (cellucid.com; version/commit: ___) from dataset ___ (dataset id/version: ___). The exported view used color-by field ___ with filters ___ applied. The interactive state (camera, fields, and filters) was saved as a session bundle (___). Figures were exported at plot size ___×___ using ___ format (___ DPI for PNG) and ___ strategy (full/optimized/hybrid).”

If you used point reduction (optimized vector), add:

> “For SVG export, points were reduced using density-preserving sampling to keep approximately ___ points for visualization.”

If you used a crop (Frame export), add:

> “Figure framing used an export crop without changing the interactive camera view.”

---

## How this relates to sessions and sharing

### Is a session bundle sufficient?

A session bundle is the best way to preserve the interactive state (views, cameras, filters, fields, highlights).
However:
- figure export UI settings are **not** stored in sessions,
- so you still need to record export settings (or recover them from embedded metadata in the exported file).

### What if the dataset export changes?

If the dataset export changes between exports (common in iterative pipelines):
- field keys can change,
- category ordering can change,
- legends can change,
- filters can target different distributions,
- and a “re-export” may not match.

For manuscript figures, treat the dataset export as immutable once you start final export production.

### How to package for collaboration

To make a figure reproducible for a collaborator, share:

1) the dataset export folder (or a stable URL), plus
2) the session bundle, plus
3) the exported artifact(s) (SVG/PNG), plus
4) (optional but recommended) a short README describing export settings and any post-processing.

:::{important}
If you need to share figures publicly, inspect embedded metadata first. It may contain local file paths or private URLs depending on your loading workflow.
:::

Related pages:
- {doc}`../l_sessions_sharing/index`
- {doc}`01_figure_export_goals_wysiwyg_and_reproducibility`

---

## Screenshot placeholder (you will replace later)

<!-- SCREENSHOT PLACEHOLDER
ID: figure-export-metadata-example
Suggested filename: figure_export/11_metadata-example.png
Where it appears: User Guide → Web App → Figure Export → 05_metadata_and_provenance.md
Capture:
  - Export a PNG (or SVG) and show a view of its embedded metadata (or a sidecar/provenance panel if the UI provides it)
  - Acceptable: a terminal `exiftool`-style view, macOS Preview “Inspector”, or a screenshot of the app’s provenance display
Redact:
  - Remove: private dataset ids/names, local file paths
Alt text:
  - An exported figure’s embedded metadata or provenance information.
Caption:
  - Provenance metadata makes it possible to recreate the same exported figure later; record dataset/version and export settings when preparing manuscripts.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for figure export metadata/provenance.
:width: 100%

Record (or embed) enough provenance to recreate the same export later, especially for manuscripts and collaboration.
```

---

## Next steps

- {doc}`07_troubleshooting_figure_export` (missing metadata, unexpected outputs)
- {doc}`06_edge_cases` (datasets/legends that stress the export pipeline)
