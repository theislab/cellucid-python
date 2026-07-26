# Troubleshooting (figure export)

**Audience:** everyone  
**Time:** 10–30 minutes (depending on symptom)  
**What you’ll learn:**
- How to diagnose common export failures and mismatches
- Which checks to do first (state, format/mode, size)
- Practical fixes and prevention steps

**Prerequisites:**
- A dataset loaded (or at least a reproducible failing state)

---

## Quick triage checklist (do this first)

This checklist resolves a large fraction of figure export reports:

1) **Confirm the active view**
   - In multiview grid compare, click inside the panel you intend to export.
   - If you meant to export all panels, ensure you are in **grid compare** and have **Export all views (split-view)** enabled.

2) **Confirm visibility**
   - Check filters/visibility: exporting “0 visible points” produces an empty plot (see {doc}`06_edge_cases`).

3) **Try a minimal “sanity export”**
   - Export a small **PNG** at a small plot size (e.g., 600×450, 150 DPI).
   - If this fails, the issue is likely download permissions, environment, or an internal error (not file size).

4) **Reduce size before debugging**
   - For big exports, first cut plot size in half and export at 150 DPI.
   - If the smaller export works, scale up cautiously.

5) **Check for Export Fidelity Warnings**
   - **Export blocked** identifies the exact missing requirement. It never
     changes the requested format or strategy.

6) **Check the browser console**
   - Look for `[FigureExport]` errors/warnings.
   - If you’re filing a bug, copy the message and include it (see {doc}`09_reference_implementation_notes`).

7) **Check the expected delivery**
   - One requested file downloads natively.
   - Any multi-file choice downloads one ZIP containing every requested file.

---

## If you see “Export blocked”

Cellucid blocks an export when it cannot exactly match the active view.

Common blockers and what they mean:

- **“Exact point export unavailable”**
  WebGL2 or the required camera matrices are unavailable. Restore that capability
  before exporting.

- **“Connectivity overlay not exported”**  
  Connectivity edges are visible, but figure export covers the point layer only.
  Disable the overlay before exporting.

There is no “export anyway” action: resolve every listed blocker first.

---

## Symptom: “Export button does nothing”

### Likely causes (ordered)

- Downloads are blocked for the site (popup/download permissions).
- A browser extension blocks blob downloads (privacy tools, ad blockers).
- Export threw an error (check console; an error toast may also appear).
- Export is running but you missed the browser’s download UI indicator.

### How to confirm

- Try exporting a tiny PNG (e.g., 600×450 at 150 DPI).
- Check your browser downloads UI (downloads tray/arrow).
- Open the developer console and look for `[FigureExport]` errors.
- Try an incognito/private window (extensions disabled).

### Fix

1) Allow downloads for the site (browser settings).
2) Disable extensions temporarily (or use an incognito window).
3) Export a small PNG first; then increase plot size/DPI.
4) If you see a console error, capture it and file a bug with context (see {doc}`09_reference_implementation_notes`).

### Prevention

- Confirm WebGL2 is available before requesting PNG or Hybrid.
- For very large datasets, avoid full vector SVG; use Hybrid or PNG.

---

## Symptom: “Export failed” (error toast/notification)

### Likely causes (ordered)

- Export size/DPI is too large for available memory.
- Browser/GPU limitations (especially for very large PNG or 3D views).
- Missing required view buffers or render state.

### How to confirm

- Try a smaller export (smaller plot size; 150 DPI).
- Try exporting **SVG** (Hybrid) instead of a huge PNG, or vice versa.
- Check the browser console for the exact error message.

### Fix

1) Reduce plot size and/or DPI.
2) Switch strategy:
   - Large dataset → Hybrid SVG or PNG.
   - Medium dataset → Optimized vector.
3) If the error message mentions a missing viewer method (for example, `getViewPositions`), file a bug with the exact error text and the steps that produced it.

---

## Symptom: “SVG is enormous / crashes Illustrator”

### Likely causes (ordered)

- You exported **Full vector** with too many visible points.
- Category explosion: a huge categorical legend dominates the layout and file size.
- Multi-panel export with many views and a large plot size.

### How to confirm

- Check the SVG file size (hundreds of MB is a red flag).
- Open the SVG in a web browser first (fast validity check).
- Export again using **Hybrid** or **Optimized vector** and compare file sizes.

### Fix

1) Export **Hybrid SVG** (recommended for large datasets and most 3D exports).
2) If you need fully vector points, use **Optimized vector** and tune the “Keep” point count.
3) Reduce legend complexity:
   - move legend to bottom,
   - disable legend,
   - reduce/merge categories upstream if appropriate.
4) Export PNG instead when editability is not required.

### Prevention

- Choose the SVG strategy explicitly; Full vector has no automatic size guard.
- Use {doc}`04_quality_knobs_and_best_practices` as your sizing/strategy playbook.

---

## Symptom: “PNG looks different than the viewer”

### Likely causes (ordered)

- Export background differs (e.g., White BG vs Match viewer).
- You changed the active view/field/filters right before export and are looking at a different state than you think.
- Very small plot size + later scaling up makes the export look blurry compared to the interactive view.

### How to confirm

- Export a small PNG first (to verify state capture works).
- Export again with **Match viewer** background and compare.
- Confirm that no **Export blocked** dialog remains unresolved.

### Fix

1) Make the exported plot larger (don’t export tiny and scale up).
2) Use **PNG (300 DPI)** or **Hybrid SVG** for 3D/shader-heavy views.
3) Restore WebGL2 in the current environment by resolving the named policy,
   driver, or extension restriction.
4) If text/legend quality is the issue, export SVG (Hybrid) and edit text in a vector editor.

---

## Symptom: “PNG or Hybrid export is blocked”

### Likely causes (ordered)

- WebGL2 is unavailable.
- The active view has not published complete camera matrices.

### How to confirm

- Read the **Exact point export unavailable** details.
- Confirm the active view is fully loaded and the canvas has a WebGL2 context.

### Fix

1) Restore WebGL2 or the missing camera matrices in the current environment.
2) Submit the same PNG or Hybrid request again. Cellucid does not route the
   request to a different renderer.

---

## Symptom: “Legends/axes missing or misaligned”

### Likely causes (ordered)

- Legend/axes are disabled in export settings.
- Multi-panel export: a shared legend is only shown when all panels share the same active field (otherwise legend behavior is limited).
- Plot size is too small for the requested legend/category count, causing layout crowding.
- Font substitution in your downstream editor shifts text widths.

### How to confirm

- Ensure **Legend** and **Axes** are enabled in the export panel.
- Export a single-view figure (disable “Export all views”) to isolate multiview legend behavior.
- Increase plot size and export again.
- Open the SVG in a browser to see whether it’s a renderer issue or an editor-specific issue.

### Fix

1) Increase plot size (and for PNG, choose appropriate DPI).
2) Move legend to **Bottom** or disable legend.
3) For large datasets, use **Hybrid** (reduces SVG complexity).
4) If the issue is editor-specific, try Inkscape vs Illustrator vs browser rendering.

---

## Symptom: “Export all views didn’t export my whole grid”

### Likely causes (ordered)

- You are not in **grid compare** mode (multi-panel export is tied to grid layout mode).
- You have no snapshot views (only the live view exists).
- The “Export all views (split-view)” option was not enabled at export time.

### How to confirm

- Switch to multiview grid compare and confirm multiple panels are visible.
- Export again with “Export all views (split-view)” enabled.

### Fix

1) Create snapshot views (Keep view).
2) Switch to **Grid compare** layout.
3) Enable **Export all views (split-view)** and export.

---

## Symptom: “Frame export crop didn’t work / crop seems wrong”

### Likely causes (ordered)

- Frame export is not enabled (checkbox off).
- You expected the preview to show the cropped output, but you are still in “edit framing” preview mode.
- Your crop frame is effectively full-frame (so the export looks uncropped).

### How to confirm

- Enable **Show preview** and adjust the crop frame.
- Click **Confirm** to switch the preview into a framed preview (so you can see what you will export).

### Fix

1) Enable **Frame export**.
2) Adjust the frame and click **Confirm** to validate.
3) If you want the crop aspect to match your plot, click **Lock aspect**.
4) Use **Reset** (or double-click) to return to full frame and try again.

---

## Symptom: “Connectivity/edges are missing in the export”

This is expected in the current export implementation:
- exports include the point layer (and point-based overlays),
- connectivity lines/edges are not exported.

Workarounds:
- Export a point-only figure and add edges via a separate plotting workflow (and document it).
- If edges are critical, file a feature request with a concrete example and target use case.

---

## Symptom: “The export contains sensitive metadata (dataset path/URL)”

### What’s happening

Exports embed provenance metadata to support reproducibility. This can include:
- dataset names/ids in filenames,
- source URLs,
- local dataset paths (depending on how data was loaded).

### Fix / mitigation

See {doc}`05_metadata_and_provenance` for:
- how to inspect metadata,
- how to strip metadata from PNG/SVG if needed,
- and how to share “data + session + figures” without leaking local paths.
---

## Next steps

- {doc}`06_edge_cases` (pathological states that stress export)
- {doc}`03_export_formats_and_renderers` (choose modes to avoid crashes)
- {doc}`09_reference_implementation_notes` (developer references)
