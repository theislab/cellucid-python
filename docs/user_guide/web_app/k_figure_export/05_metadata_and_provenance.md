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
- explicit SVG point strategy (Full vector, Optimized vector, or Hybrid)
- include axes/legend/title/background
- frame export crop rectangle (if used)

:::{tip}
If you saved a session bundle right before exporting, you already captured most of the view/scientific state. The remaining missing piece is export settings (which are not stored in sessions).
:::

---

## What metadata is embedded into exports (must be explicit)

Cellucid exports embed metadata by default (there is no “metadata off” toggle in the current UI).

Every panel’s record is read out of the same atomic snapshot the panel is drawn
from — the camera and render state, the point positions, colours and
transparency the viewer is drawing with, that view’s own filters, and its legend
— so an export cannot stamp a figure with a view it did not draw. That snapshot
carries the view itself into the figure rather than an approximation of it: the
reference grid the viewer’s background draws, the atmospheric fog that carries
depth, and, under **Background: Match viewer**, the viewer’s own paper — whose
luminance also chooses the ink for the title, axes, legend, and selection badge.
The bounded differences are listed in {doc}`06_edge_cases`.

:::{important}
**Provenance is per panel.** A multi-panel export makes one provenance claim for
a file that shows several different views, so the record carries **one entry per
rendered panel** — its own view, color field, and filters. There is no top-level
`view`, `field`, or `filters` value that could be read as applying to the whole
figure. A single-view export is simply the one-entry case.
:::

### PNG metadata (UTF-8 iTXt chunks)

PNG exports include standard uncompressed `iTXt` fields:
- `Software`: `Cellucid (cellucid.com)`
- `Website`: `https://cellucid.com`
- `Creation Time`: ISO timestamp
- `Dataset` and `Dataset ID` (when known)
- `Color Field`: **only when every exported panel uses the same color field.**
  A single-view export always carries it (when a field is active); a multi-panel
  export whose panels use different fields omits the key entirely rather than
  claiming one panel's field for the figure. The per-panel record below is then
  the only truthful answer.
- `Source File` (may be a URL or a local path, depending on how the dataset was loaded)
- `Description`: a human-readable summary, one clause per panel
- `Comment`: a compact JSON blob with the structured record

### What `Description` looks like

Every panel appears with the letter it is drawn with, its field, and its own
filters. A single-view export drops the letter (there is only one panel):

```text
Views: Live (Field: cell_type; Filters: none) • Source: https://example.org/pancreas
```

A grid names each panel separately:

```text
Views: A. Live (Field: cell_type; Filters: total_counts: 500.00 – 12000.00) | B. Endoderm (Field: leiden; Filters: none) • Source: /Users/you/exports/pancreas
```

Read it as: `Views:` then one `Panel. Label (Field: …; Filters: …)` clause per
panel joined by ` | `, then the source. `none` means exactly that — no color
field, or no active filter. The SVG `dc:description` is built by the same code,
so it is identical text for the same export.

### What `Comment` looks like (structured JSON)

```json
{
  "generator": "https://cellucid.com",
  "exporter": { "name": "Cellucid", "website": "https://cellucid.com" },
  "exportedAt": "2026-05-14T09:12:44.108Z",
  "dataset": {
    "name": "Pancreas", "id": "pancreas", "sourceType": "github-repo",
    "baseUrl": "https://example.org/pancreas", "userPath": null,
    "source": { "name": null, "url": null, "citation": null }
  },
  "views": [
    {
      "panel": "A", "id": "live", "label": "Live",
      "field": { "key": "cell_type", "kind": "category" },
      "filters": ["total_counts: 500.00 – 12000.00"]
    },
    {
      "panel": "B", "id": "view_2", "label": "Endoderm",
      "field": { "key": "leiden", "kind": "category" },
      "filters": []
    }
  ],
  "export": {
    "format": "png", "width": 1200, "height": 900, "dpi": 300,
    "strategy": null, "includeAxes": true, "includeLegend": true,
    "legendPosition": "right", "background": "white",
    "backgroundColor": "#ffffff", "crop": null
  }
}
```

Notes on reading it:
- `views` always has one entry per rendered panel, in drawing order.
- `panel` is the drawn panel letter (`A`, `B`, … `AA`) in a multi-panel export
  and `null` in a single-view export, where no letter is drawn.
- `field.key` / `field.kind` are `null` when the panel has no active color
  field; `kind` is `"category"` or `"continuous"`.
- `filters` lists only real restrictions; the app's “No filters active”
  placeholder never enters the record, so an unfiltered panel gets `[]`.
- `export.dpi` is `null` for SVG jobs, and `export.strategy` is the explicit SVG
  point strategy (`null` for PNG-only requests).

### SVG metadata (RDF / Dublin Core + Cellucid JSON)

SVG exports include a `<metadata>` block with:
- Dublin Core fields: creator, publisher, relation, date, source, identifier,
  and a `dc:description` that is the same text as the PNG `Description`
- `cellucid:colorField` — present under the same rule as PNG `Color Field`:
  only when every panel shares one field
- `cellucid:sourceFile` — the source URL or local path
- `cellucid:json` — the same structured blob as the PNG `Comment`, differing
  only in `export.format` (`"svg"`) and `export.dpi` (`null`)

Because both formats are built from one source, the SVG and the PNG of a single
export cannot disagree about what they show.

### Filenames (also part of provenance)

Exports are named conservatively:
- `<dataset>_<color-field>_<view>_<timestamp>.svg`
- `<dataset>_<color-field>_<view>_dpi300_<timestamp>.png` (when exporting multiple DPIs)
- `<dataset>_<color-field>_<view>_batch_<timestamp>.zip` (multi-file export)

Two segments behave differently for a multi-panel export:

- `<view>` is the exported view's label for a single-view export, and the
  literal `multiview` for a grid.
- `<color-field>` is **omitted** when the panels do not all share one color
  field, for the same reason `Color Field` is: a filename may not name a field
  that only some panels use. A grid where every panel shares the field keeps the
  segment.

This makes it easier to tell “what is this file?” even without opening it.

A batch ZIP adds no alternate scientific representation: each entry is the
same native artifact, with its own embedded metadata, that the corresponding
single-file renderer produces.

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

For a multi-panel figure, state each panel separately — the embedded
`Description` already lists them in that form, so you can copy it:

> “Panel A shows ___ colored by ___ (filters: ___); panel B shows ___ colored by ___ (filters: ___).”

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

## Interface reference

```{figure} ../../../_static/screenshots/figure_export/labels-annotations.png
:alt: Figure Export Labels and Annotations controls for title, axes, legend, and annotation visibility.
:width: 224px

Labels and Annotations makes title, axis, legend, and annotation choices explicit.
```

---

## Next steps

- {doc}`07_troubleshooting_figure_export` (missing metadata, unexpected outputs)
- {doc}`06_edge_cases` (datasets/legends that stress the export pipeline)
