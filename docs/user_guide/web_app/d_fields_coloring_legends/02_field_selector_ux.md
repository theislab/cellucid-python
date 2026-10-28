# Field selector UX

**Audience:** everyone (users selecting fields; computational users debugging “where did my column go?”)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- Where field selection lives in the UI (and how it behaves in multiview)
- How obs field dropdowns differ from gene expression search
- What **Duplicate / Rename / Delete / Clear** actually do
- How **Deleted fields** restore + “Confirm deletion” work (soft delete vs purge)
- What the `*` marker means and why renames don’t break loading

---

## Where the field selector lives

All field selection for coloring happens in the left sidebar:

- **Coloring & filtering** → the three selectors:
  - **Categorical obs**
  - **Continuous obs**
  - **Gene expression**

The legend and related controls appear underneath in the **Legend** section.

:::{important}
In multiview (live + snapshots), the sidebar controls apply to the **currently active view**.

If you’re in Grid compare, click a panel first to make it the active view, then pick a field.
:::

```{figure} ../../../_static/screenshots/web_app/field-selectors-three-routes.png
:alt: The Color by heading with its (i), above the three field selectors stacked in order — Categorical obs holding clusters, Continuous obs holding None, and an empty Gene expression search box — each followed by four small icon buttons, with the pointer pressing the first.
:width: 504px

The three selectors, in the order they appear. The small buttons to the
right of each are **Duplicate**, **Rename**, **Delete** and **Clear**, and on
the two obs rows also **Export** (after **Duplicate**: the field as a CSV of
obs names and values, {doc}`../f_highlighting_selection/03_highlight_ui`);
they are greyed out until that selector holds a field, which is why they look
inert on the two empty rows here.
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
search results**, and on rows in the **Deleted fields** panel — so if you need
to see what a renamed field used to be called, look there.
:::

---

## Gene expression selector (searchable dropdown)

Gene expression is usually a long list (thousands of genes), so the UI is a search box:

- Click in **Gene expression** → results dropdown opens
- Type to filter genes by substring match
- Click a result to load and activate it

Keyboard shortcuts:
- `↓` / `↑` move the highlighted result through the dropdown (opening it if
  it is closed); they wrap from the last result to the first and back, and
  `↑` with nothing highlighted starts at the last result
- `Enter` selects the highlighted result, or the **top** result when none is
  highlighted
- `Esc` closes the dropdown (and removes focus)

A gene you are about to choose starts loading early. When a result stays
under the pointer, or highlighted with the arrow keys, for 150 ms, its values
are fetched and decoded in the background, so choosing it usually colours the
view without a wait. Moving on to another result cancels that load, and
choosing the result uses it instead of starting a second one. A background
load that fails stays silent; choosing the gene loads it again and reports the
error then.

The search covers the gene list in use. A dataset's gene names are read
right after its first view is drawn: until then the search shows **Reading
gene names…** and cannot be typed in (and, if they could not be read, says
so, with the reason as its tooltip). Changing the gene names or the
expression matrix under **Session › Data mapping** replaces that list; a gene
a view is coloured by stays selected when the new list has it (see
{doc}`../b_data_loading/13_data_mapping`).

Results:
- every match is listed; the dropdown opens with the first **100** and builds
  the next 100 each time you scroll or arrow to its end
- **Home** and **End** jump to the first and the last match

```{figure} ../../../_static/screenshots/web_app/gene-search-dropdown-open.png
:alt: The Gene expression search box holding the letters Rp with its list open below, naming the genes that contain them.
:width: 496px

Typing into the box opens the list immediately and filters it on every
keystroke. Matching is a case-insensitive substring test, so `Rp` matches
anywhere in the name, not only at the start.
```

If nothing matches, the dropdown says so rather than going blank:

```{figure} ../../../_static/screenshots/web_app/gene-search-no-match.png
:alt: The Gene expression list after typing a symbol this sample does not publish: No gene matches "CD19", a sentence giving the 3,753 gene names the dataset publishes, and a link Why a gene may be missing.
:width: 496px

The empty state names the number of genes **this** dataset publishes, which is
the number you need in order to tell “I mistyped it” apart from “this dataset
does not contain it”. See {doc}`07_genes_in_the_built_in_samples`.
```

:::{note}
Loading a gene can take time (especially on remote datasets). While a gene is loading, Cellucid temporarily disables other selectors so the active field state can’t race.
:::

### What the list contains

The dropdown holds exactly the gene names the loaded dataset published, plus any
gene fields you created inside Cellucid. There is no second identifier to search
instead and no alias table, so a name the file does not hold matches nothing
even when the gene exists in the underlying study.

That is worth knowing before you conclude a gene is missing from your data: the
four human built-in samples publish only part of their upstream gene axis, for a
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
- no control or zero-width character (`U+200B`, `U+2060`, `U+FEFF`)
- must be unique among visible fields

A rename is dataset-wide: it applies to the live view and to every kept view at
once, including their panel labels and filter summaries, and a saved session
records it that way.

To “undo” a rename:
- rename the field back to its original key (the `*` marker disappears).

:::{note}
If Community Annotation voting is enabled for a categorical field, rename may be disabled for that field (to keep votes stable). You’ll see an error message if you try.
:::

### Delete (soft delete)

Removes the field from dropdowns *without destroying it*:
- the field moves to **Deleted fields**
- you can restore it later

This is a “soft delete” and is designed to be safe.

Like a rename, a delete is dataset-wide: when you next select a kept view,
the field is in its **Deleted fields** too, and is no longer that view's
active field. A restore and a Confirm are dataset-wide too.

### Clear (clear selection)

Clears the active selection for that selector.

If all selectors are cleared, there is no active field:
- the legend hides,
- coloring returns to the default viewer state for “no field selected”.

---

## Deleted fields panel (Restore vs Confirm deletion)

When you delete fields (or when Cellucid creates derived fields and hides the originals), a **Deleted fields** section appears at the bottom of **Coloring & filtering**.

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
stays in Deleted fields. There is no automatic `(restored)` suffix.

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
