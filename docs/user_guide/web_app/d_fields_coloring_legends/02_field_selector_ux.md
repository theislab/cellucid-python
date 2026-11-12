# Field selector UX

**Audience:** everyone (users selecting fields; computational users debugging “where did my column go?”)  
**Time:** 15–25 minutes  
**What you’ll learn:**
- Where field selection lives in the UI (and how it behaves in multiview)
- How obs field dropdowns differ from gene expression search
- What **Duplicate / Group / Rename / Delete / Clear** actually do
- How **Deleted fields** restore + **Delete for good…** work (soft delete vs purge)
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
:alt: The Color by heading with its (i), above the three field selectors stacked in order — Categorical obs holding clusters, Continuous obs holding None, and an empty Gene expression search box — the Categorical obs row followed by six small icon buttons, the Continuous obs row by five and the gene row by four, with the pointer pressing the first selector.
:width: 504px

The three selectors, in the order they appear. The small buttons to the
right of each are **Duplicate**, **Rename**, **Delete** and **Clear**, and on
the two obs rows also **Export** (after **Duplicate**: the field as a CSV of
obs names and values, {doc}`../f_highlighting_selection/03_highlight_ui`),
and on the categorical row **Group** before it (below); they are grayed out until that selector holds a field, which is why they look
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

### What a dataset opens colored by

A built-in sample opens as its starting view colors it. Any other dataset
(your own file or folder, a link, a served object) opens colored by its best
categorical field, shown in **Categorical obs** like any choice, and its
history starts from it (it is not a step). The field is chosen by its name,
compared in one spelling (lower case, any punctuation as `_`, so `Cell Type`,
`cell.type` and `CellType` are alike), in this order:

1. a cell type: `cell_type`, `celltype`, then any name holding *cell type*
   (`author_cell_type`, `predicted.celltype.l2`) and `cell_ontology_class`;
2. another annotation: `*annotation*`, `ann_level_*`, `majority_voting`,
   `predicted_labels`, `labels`, `cell_state`, `*subtype*`;
3. a lineage or compartment: `*lineage*`, `*compartment*`, `cell_class`,
   `broad_cell_type`, …;
4. clusters: `ident`, `leiden*`, `louvain*`, `seurat_clusters`, `cluster*`,
   `*_snn_res.*`, …;
5. the tissue or organ, then a condition (`condition`, `disease`, `treatment`,
   `development_stage`, …), then the cell-cycle phase.

Within one rank, the second level of a hierarchy comes first, then the first,
then deeper (`predicted.celltype.l2` before `.l1` and `.l3`), then the file's
order. A field with one category is skipped. Never chosen: a field named like
a sample or where it came from (`donor`, `patient`, `sample`, `orig.ident`,
`library`, `batch`, `plate`, `assay`, …), an id or barcode, or a quality flag
(`doublet`, `is_*`, `pass_qc`, …).

With none of those names, the first other categorical field with 2 to 72
categories is used (72 is how many colors the category palette designs;
beyond them the colors are generated and harder to tell apart). With none
at all, the dataset opens with **None**. The composition and dot plot analyses
take their **Cell types** / **Group by** and **Samples** defaults from the same
ranks.

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

A field you made (a copy, a group, a column from highlight pages) has no
original name: its name is whatever you call it now, so a rename shows no `*`
and offers no **Put name back**, before and after the session is saved and
opened again.

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
- Type (or paste) to filter genes by substring match, ignoring letter case
- Click a result to load and activate it

A gene name is matched as it is written: a dot is a character like any other
(`AC058791.1`, `RP11-34P13.7`). An Ensembl accession typed with a version
(`ENSG00000141510.17`) also finds the same accession listed without one or
with another version, since the version numbers a release of the annotation,
not another gene; gene sets pasted into Gene signature or Dot plot are matched
the same way.

Keyboard shortcuts:
- `↓` / `↑` move the highlighted result through the dropdown (opening it if
  it is closed); they wrap from the last result to the first and back, and
  `↑` with nothing highlighted starts at the last result
- `Enter` selects the highlighted result, or the **top** result for what the
  box says when none is highlighted (the gene named exactly first, then names
  that start with it). When the name matches no gene, the dropdown says **No
  gene matches** and screen readers hear "*name* is not a gene in this
  dataset."
- `Esc` closes the dropdown (and removes focus)

A gene you are about to choose starts loading early. When a result stays
under the pointer, or highlighted with the arrow keys, for 150 ms, its values
are fetched and decoded in the background, so choosing it usually colors the
view without a wait. Moving on to another result cancels that load, and
choosing the result uses it instead of starting a second one. A background
load that fails stays silent; choosing the gene loads it again and reports the
error then.

The search covers the gene list in use. A dataset's gene names are read
right after its first view is drawn: until then the search shows **Reading
gene names…** and cannot be typed in (and, if they could not be read, says
so, with the reason as its tooltip). A gene chosen from elsewhere meanwhile (an
analysis result's gene) is not shown, and a notification says that the gene
names are still being read.

A gene chosen is always shown or a notification says why not:

- chosen while more cells are being opened (a progressive opening's step, or
  **Apply** under **Cells to open**), it is read on the cells opened after
  that step and shown then;
- chosen while the gene list is replaced (**Session › Data mapping**, whether
  you applied it or the app did because the values contradict the file's
  `log1p` record), it is read from the new list and shown, under the name the
  new list gives it; if the new list does not have it, a notification says so;
- chosen while a session or a sample's saved view is being restored or saved,
  it is not shown, and a notification says which and asks you to choose it
  again when that has ended.

Choosing another gene, pressing **Cancel** on its loading notification, or
opening another dataset replaces a choice without a message. Changing the gene names or the
expression matrix under **Session › Data mapping** replaces that list; a gene
a view is colored by stays selected when the new list has it (see
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
:alt: The Gene expression list after typing a symbol this sample does not publish: No gene matches "CD19", a sentence saying it searched the 3,753 genes of this dataset named by its var index and that another var column (symbols) can name them instead, a Choose gene names… button, and a link Why a gene may be missing.
:width: 496px

The empty state names the number of genes **this** dataset publishes, which is
the number you need in order to tell “I mistyped it” apart from “this dataset
does not contain it”. See {doc}`07_genes_in_the_built_in_samples`.
```

:::{note}
Loading a field or a gene can take time (especially on remote datasets). While one is loading, the menus show the field or gene you chose and are disabled until it is shown, so the active field state can’t race; a refresh of the menus meanwhile (the gene list arriving, an annotation connection changing) keeps showing your choice.
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

Each selector row has the same four actions (the categorical row also
**Group**, below):

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

(group-categories-into-a-new-field)=
### Group categories into a new field

The categorical row has one more button after **Duplicate**: **Group**
(`Group the categories of the selected categorical obs field into a new
field`). It opens **Group categories**, one row per category with a box for
the group it joins, to make a coarser annotation in one go: on Pancreas, a
lineage from `clusters` with *Ductal* and *Ngn3 low EP* typed as
*Progenitor* and *Alpha*, *Beta*, *Delta* and *Epsilon* as *Endocrine*.

- Categories given the same group become one category of the new field, in
  the order the groups are first typed; the group labels typed so far are
  offered in every box.
- A category left empty keeps its name and comes after the groups. With
  **Others into** ticked, every category left empty joins that one instead
  (*Other* by default, or a group's label).
- **New field** names the field made (by default `clusters (grouped)`,
  numbered when taken). The field it is made from stays as it is: a grouping
  is another annotation beside it, which is what you compare it with.
- A group takes the color of the first category it holds when that color can
  be told apart from the kept categories' and the groups' before it, else the
  next palette color that can (two groups never read as one: grouping
  Kanemaru's atrial cardiomyocytes against the rest gave amber and bright
  orange). It is shown while any of its categories was; a kept category keeps
  its own color. Every color the app gives a new category or page follows the
  same rule (the [legend](04_legend_behavior.md) and
  [highlight pages](../f_highlighting_selection/03_highlight_ui.md) say where).
- The new field is shown at once. It is one step of the history, so one
  **Undo** removes it.
- One line under the rows says what **Group** makes ("2 groups from 6
  categories; 2 others keep their names"); a label two groups share, or a
  kept category's name a group also takes, is refused there in a sentence
  and the dialog stays open.

Like a merge, the grouping is a function of the file's column: cells opened
later (a larger subsample) are labeled by the same rule. An agent, a
notebook or a replay makes the same field with the action `field.group`
({doc}`../p_developer_docs/30_action_reference`), as one step.

### Rename

Renames the field key shown in the UI. The name box takes the selector's
place: press {kbd}`Enter` or click elsewhere to rename, {kbd}`Escape` to leave
it as it was. A name refused (one a field shown has: *A field with this name
already exists*) is said on the box, which stays open until you change it.

Naming rules (enforced):
- name must be non-empty
- no leading/trailing whitespace
- no control or zero-width character (`U+200B`, `U+2060`, `U+FEFF`)
- must be unique among visible fields

A rename is dataset-wide: it applies to the live view and to every kept view at
once, including their panel labels and filter summaries, and a saved session
records it that way.

A rename is one step of the history, so **Undo** (`Ctrl+Z` / `⌘Z`) takes it
back; renaming the field back to its original key does too (the `*` marker
disappears).

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

## Deleted fields panel (Restore vs Delete for good)

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
`A field named "<name>" is shown already; rename it before restoring the
deleted one.` and the field stays in Deleted fields. There is no automatic `(restored)` suffix.

To get the field back, rename the field that is currently holding the name,
then restore again.
:::

### Delete for good…

This is intentionally destructive, so it asks first (**Delete for good**, or
Cancel):
- the field leaves **Deleted fields** and can no longer be restored there,
  nor from a session saved afterwards;
- it is one step of the history, so **Undo** (`Ctrl+Z` / `⌘Z`) still brings it
  back to **Deleted fields** (the step reads `Deleted <name> for good`).

Use **Delete for good…** when you are sure you won’t want the column back and want to reduce clutter.

### Two deleted fields of one name

A name is unique among the fields shown, not among the deleted ones: delete
`score`, rename another field to `score` and delete it too, and **Deleted
fields** lists both, as `score (1 of 2)` and `score (2 of 2)`. Each row's
**Restore** and **Delete for good…** act on that row's field. An agent or a
notebook names one by its place with `which`
(`viewer.field_restore(field={"source": "obs", "key": "score"}, which=2)` in
Python, `cellucid_field_restore(viewer, list(source = "obs", key = "score"),
which = 2)` in R); named by its name alone, the action is refused with a
sentence saying so.

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
