# Field selector UX

**Audience:** everyone (users selecting fields; computational users debugging “where did my column go?”)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- Where field selection lives in the UI (and how it behaves in multiview)
- How obs field dropdowns differ from gene expression search
- What **Duplicate / Rename / Delete / Clear** actually do
- How **Deleted Fields** restore + “Confirm deletion” work (soft delete vs purge)
- What the `*` marker means and why renames don’t break loading

---

## Where the field selector lives

All field selection for coloring happens in the left sidebar:

- **Coloring & Filtering** → the three selectors:
  - **Categorical obs**
  - **Continuous obs**
  - **Gene Expression**

The legend and related controls appear underneath in the dashed **Display options** box.

:::{important}
In multiview (live + snapshots), the sidebar controls apply to the **currently active view**.

If you’re in Grid compare, click a panel first to make it the active view, then pick a field.
:::

<!-- SCREENSHOT PLACEHOLDER
ID: field-selector-three-rows-and-actions
Suggested filename: web_app/fields_legends/10_field-selector-rows-actions.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 02_field_selector_ux.md
Capture:
  - Dataset loaded
  - Sidebar open to “Coloring & Filtering”
  - Make sure all three rows are visible:
    - Categorical obs dropdown populated
    - Continuous obs dropdown populated
    - Gene Expression box visible (i.e., dataset has genes)
  - Hover a field action icon so the tooltip appears (“Duplicate field”, “Rename field”, etc.)
Crop:
  - Include: the three selector rows and the row of action icons on the right of each selector
  - Exclude: unrelated sidebar sections if they add clutter
Annotations:
  - Call out: #1 Categorical obs dropdown, #2 Continuous obs dropdown, #3 Gene Expression search, #4 action buttons (copy/rename/delete/clear)
Alt text:
  - The three field selectors in the Coloring & Filtering panel with action buttons.
Caption:
  - Explain that you can select exactly one active “color-by” source at a time, and the action buttons modify the field list (duplicate/rename/delete).
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the field selector rows and action buttons.
:width: 100%

Pick the active “color-by” field from Categorical obs, Continuous obs, or Gene Expression, then use the action buttons to duplicate/rename/delete/clear.
```

---

## The “one active field” rule (obs vs gene expression)

At any moment, **only one source drives coloring**:

- an **obs** field (categorical *or* continuous), or
- a **gene expression** (var) field.

When you select something:
- picking a gene clears the obs dropdown selections,
- picking an obs field clears gene selection,
- choosing **None** clears the active field (legend hides).

---

## Obs dropdowns (Categorical obs vs Continuous obs)

### What you’ll see

Each obs dropdown always includes:
- **None** (clears that selector)
- a list of available fields of that kind

If there are no fields of a kind, the dropdown shows an empty-state label:
- `(no categorical obs fields)`
- `(no continuous obs fields)`

### The `*` marker and “Original:” tooltips

If you rename a field, Cellucid shows:
- a trailing `*` in the dropdown label (example: `clusters_cleaned *`)
- a tooltip with the original name (`Original: clusters`)

This is not cosmetic: the app keeps the original key internally so it can still load data correctly after a rename.

---

## Gene expression selector (searchable dropdown)

Gene expression is usually a long list (thousands of genes), so the UI is a search box:

- Click in **Gene Expression** → results dropdown opens
- Type to filter genes by substring match
- Click a result to load and activate it

Keyboard shortcuts:
- `Enter` selects the **top** result in the dropdown
- `Esc` closes the dropdown (and removes focus)

Result limit:
- the dropdown shows up to **100** matches
- if there are more, you’ll see a message like `...and 2,341 more. Type to narrow results.`

:::{note}
Loading a gene can take time (especially on remote datasets). While a gene is loading, Cellucid temporarily disables other selectors so the active field state can’t race.
:::

<!-- SCREENSHOT PLACEHOLDER
ID: gene-expression-search-dropdown-open
Suggested filename: web_app/fields_legends/11_gene-search-dropdown.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 02_field_selector_ux.md
Capture:
  - Dataset with many genes loaded
  - Click Gene Expression search box to open the dropdown
  - Type a query that yields >10 results (so the reader sees scrolling/selection)
  - Optional: pick a query that yields >100 results so the “...and N more” message appears
Crop:
  - Include: the Gene Expression input and the open dropdown list
Annotations:
  - Call out: query text, one result row, “...and N more” message (if present)
Alt text:
  - Gene expression search box with a dropdown list of matching genes.
Caption:
  - Explain that genes are selected via search and loaded on demand; refine the query to narrow results.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the gene expression searchable dropdown.
:width: 100%

Gene expression uses a searchable dropdown: type to narrow matches, then click (or press Enter) to load a gene and color by it.
```

---

## Field action buttons (Copy / Rename / Delete / Clear)

Each selector row has the same four actions:

### Copy (Duplicate)

Creates a new field like:
- `clusters (copy)`
- `n_counts (copy)`
- `MS4A1 (copy)` (genes)

Why you would duplicate:
- you want to experiment with **category merges / deletions** without losing the original column,
- you want a “working copy” to rename for clarity (paper figures, screenshots),
- you want to preserve a specific set of hidden categories and colors for one workflow.

Implementation detail (useful to know for power users):
- duplicating a **categorical obs** field makes a true copy of category codes (so edits won’t affect the original),
- duplicating a **continuous** field (obs or gene) creates a user-defined “alias” field that can be reconstructed from the original key during session restore.

### Rename

Renames the field key shown in the UI.

Naming rules (enforced):
- name must be non-empty
- no leading/trailing whitespace
- cannot contain `:`
- must be unique among visible fields

To “undo” a rename:
- rename the field back to its original key (the `*` marker disappears).

:::{note}
If Community Annotation voting is enabled for a categorical field, rename may be disabled for that field (to keep votes stable). You’ll see an error message if you try.
:::

### Delete (soft delete)

Removes the field from dropdowns *without destroying it*:
- the field moves to **Deleted Fields**
- you can restore it later

This is a “soft delete” and is designed to be safe.

### Clear (clear selection)

Clears the active selection for that selector.

If all selectors are cleared, there is no active field:
- the legend hides,
- coloring returns to the default viewer state for “no field selected”.

---

## Deleted Fields panel (Restore vs Confirm deletion)

When you delete fields (or when Cellucid creates derived fields and hides the originals), a **Deleted Fields** section appears at the bottom of **Coloring & Filtering**.

It contains:
- an **Obs** group (deleted obs fields)
- a **Genes** group (deleted genes / var fields)

Each deleted item has two actions:

### Restore (undo delete)

Brings the field back into the dropdowns.

If restoring would create a name collision, Cellucid auto-renames the restored field (for example: `clusters (restored)`), and shows a notification.

### Confirm (purge restore capability)

This is intentionally destructive:
- the field becomes **non-restorable** in the current session *and* in saved session bundles.

Use Confirm when you are sure you won’t want the column back and want to reduce clutter.

<!-- SCREENSHOT PLACEHOLDER
ID: deleted-fields-restore-confirm
Suggested filename: web_app/fields_legends/12_deleted-fields-panel.png
Where it appears: User Guide → Web App → Fields, Coloring, and Legends → 02_field_selector_ux.md
Capture:
  - Delete at least one obs field and (optionally) one gene so both groups appear
  - Expand the Deleted Fields accordion
  - Ensure Restore and Confirm buttons are visible
Crop:
  - Include: “Deleted Fields” header + at least one row with its buttons
Annotations:
  - Call out: Restore vs Confirm and explain “soft delete vs permanent purge”
Alt text:
  - Deleted Fields panel showing restore and confirm deletion buttons.
Caption:
  - Emphasize safety: delete is reversible; confirm deletion is permanent.
-->
```{figure} ../../../_static/screenshots/placeholder-screenshot.svg
:alt: Placeholder screenshot for the Deleted Fields restore/confirm panel.
:width: 100%

Deleted Fields lets you restore soft-deleted fields, or permanently confirm deletion to remove restore capability.
```

---

## Advanced: Community Annotation “voting mode” indicators (optional)

If you are connected to a Community Annotation repo, categorical fields can be marked as “annotatable”:
- the field name gains a `🗳️` badge in the dropdown
- category label clicks can open a voting modal (labels are locked while voting is enabled)

Power-user gesture:
- **Right-click** the *Categorical obs* dropdown to enable/disable annotation voting for the selected field (authors only).

If you’re not using Community Annotation, you can ignore these indicators.

---

## Next steps

- For exact color semantics (missing values, log scale, “None” gray), read `d_fields_coloring_legends/03_color_by_behavior`.
- For every legend interaction (colormaps, sliders, category merges), read `d_fields_coloring_legends/04_legend_behavior`.
