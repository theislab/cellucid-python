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

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: Three stacked field controls labelled CATEGORICAL OBS with a dropdown reading clusters, CONTINUOUS OBS with a dropdown reading None, and GENE EXPRESSION with a text box reading Search genes; each is followed by four small square icon buttons showing two overlapping rectangles, a pencil, a bin and a cross, greyed out on the two unused rows, with the mouse pointer pressing the clusters dropdown.
:width: 488px

The three selectors, in the order they appear. The four small buttons to the
right of each are **Duplicate**, **Rename**, **Delete** and **Clear**; they are
greyed out until that selector holds a field, which is why they look inert on
the two empty rows here.
```

:::{tip}
The buttons are icon-only, and their meaning lives in the tooltip rather than a
label. Hovering the first one on the categorical row shows
`Duplicate the selected categorical obs field` — note the app says *duplicate*,
not *copy*.
:::

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

### The `*` marker for renamed fields

If you rename a field, Cellucid shows a trailing `*` in the dropdown label
(example: `clusters_cleaned *`).

This is not cosmetic: the app keeps the original key internally so it can still load data correctly after a rename.

:::{note}
**The obs dropdowns carry no `Original:` tooltip.** The option text is the only
thing the app sets there, so hovering a renamed obs option shows nothing. An
`Original: <old name>` tooltip does exist in two other places — on **gene
search results**, and on rows in the **Deleted Fields** panel — so if you need
to see what a renamed field used to be called, look there.
:::

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

```{figure} ../../../_static/screenshots/web_app/gene-search-dropdown-open.png
:alt: A field row labelled GENE EXPRESSION with a text box containing the letters Rp, and an open dropdown list below it showing gene names that contain Rp, one per row.
:width: 480px

Typing into the box opens the list immediately and filters it on every
keystroke. Matching is a case-insensitive substring test, so `Rp` matches
anywhere in the name, not only at the start.
```

If nothing matches, the dropdown says so rather than going blank:

```{figure} ../../../_static/screenshots/web_app/gene-search-no-match.png
:alt: The same field row with CD19 typed into the box and a dropdown below reading No gene matches, followed by a sentence stating how many gene names this dataset publishes and warning that it may be a subset of the source data and to check for typos, and a link reading Why a gene may be missing.
:width: 480px

The empty state names the number of genes **this** dataset publishes, which is
the number you need in order to tell “I mistyped it” apart from “this export
does not contain it”. See {doc}`07_genes_in_the_built_in_samples`.
```

:::{note}
Loading a gene can take time (especially on remote datasets). While a gene is loading, Cellucid temporarily disables other selectors so the active field state can’t race.
:::

### What the list contains

The dropdown holds exactly the gene names the loaded dataset published, plus any
gene fields you created inside Cellucid. There is no second identifier to search
instead and no alias table, so a name the export did not publish matches nothing
even when the gene exists in the underlying study.

That is worth knowing before you conclude a gene is missing from your data: the
four human built-in samples publish only part of their prepared gene axis, for a
documented reason. See {doc}`07_genes_in_the_built_in_samples`.



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

Brings the field back into the dropdowns, and reports `Restored "<name>"`.

:::{important}
**A restore that would collide with a visible name is refused, not renamed.**
Cellucid stops with
`Cannot restore "<name>" while that visible field name exists` and the field
stays in Deleted Fields. There is no automatic `(restored)` suffix.

To get the field back, rename the field that is currently holding the name,
then restore again.
:::

### Confirm (purge restore capability)

This is intentionally destructive:
- the field becomes **non-restorable** in the current session *and* in saved session bundles.

Use Confirm when you are sure you won’t want the column back and want to reduce clutter.

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

- For exact color semantics (missing values, log scale, “None” gray), read {doc}`03_color_by_behavior`.
- For every legend interaction (colormaps, sliders, category merges), read {doc}`04_legend_behavior`.
