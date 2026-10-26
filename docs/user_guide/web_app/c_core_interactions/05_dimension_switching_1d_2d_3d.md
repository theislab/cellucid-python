# Switching embeddings (1D / 2D / 3D)

**Audience:** everyone (switching embeddings is foundational)  
**Time:** 8–15 minutes  
**What you’ll learn:**
- What a view shows: one embedding, whose dimension is 1D, 2D or 3D
- How to switch it (the picker, keys `1`, `2` and `3`)
- What changes (positions, navigation mode, overlays) vs what stays

---

## What a view shows

Every view shows one **embedding**: a set of coordinates for the cells. Its
dimension follows from the embedding:

- **1D**: one coordinate per cell (often a timeline or pseudotime-like axis)
- **2D**: a 2D embedding (UMAP, t-SNE, a PCA plane, spatial coordinates)
- **3D**: a 3D embedding (a 3D UMAP, three principal components)

A dataset can carry any number of embeddings of each dimension. The ones turned
on under **Session › Data mapping** ({doc}`../b_data_loading/13_data_mapping`)
are the ones views can show; one embedding per dimension is that dimension's
**default**. Different views can show different embeddings, of the same or of
different dimensions.

---

## The picker

**Compare views › View N shows** (with a single view, **Shows**) lists every
embedding that is turned on, grouped by dimension, as `UMAP · 3D`,
`PCA (unintegrated) · 2D`. Each dimension's default comes first in its group;
hover an entry to see whether it is the default and the label the dataset
stores it under (for a file opened in the browser, the `obsm` key and its
columns, such as `X_pca[0,1,2]`). Choosing one
changes the selected view only. The picker is disabled, and says why, when the
dataset shows a single embedding.

The selected view is the one highlighted in the **Views** list above the
picker; its label names it. In **Grid compare**, clicking inside a panel
selects it too.

---

## Keys 1, 2 and 3

Pressing `1`, `2` or `3` (outside a text field) changes the selected view:

- If the view shows another dimension, it switches to that dimension's
  default embedding.
- If it already shows that dimension, it moves to the next embedding of that
  dimension that is turned on, in the picker's order, and wraps around after
  the last.

A short notice says what is shown now and its place in the cycle, for example
`2D · t-SNE (2 of 3)`. When no embedding of that dimension is turned on, the
notice says so — `No 1D embedding is on. Turn one on in Session › Data
mapping.` — and nothing changes. The shortcut sheet (`?`) lists the keys as
`1D, 2D, 3D embedding · again: next`.

---

## Per view

Each view keeps its own embedding:

- the live view can show a 3D UMAP,
- a kept view a 2D t-SNE,
- another kept view three principal components in 3D.

Cells, highlights and selections are shared; colouring, filters, camera,
navigation mode, centroid labels and vector fields stay with each view.
Sessions keep each view's embedding. See
{doc}`04_view_layout_live_snapshots_small_multiples`.

---

## What changes when a view switches embedding

### Changes

- **Point positions** change to the coordinates of the new embedding.
- **Centroid labels** move to the new embedding's category centres.
- **The navigation mode**, when the dimension changes: 3D selects **Orbit**,
  1D and 2D select **Planar** — unless you chose a mode for this view in that
  dimension before; see below.
- **Vector fields** follow: a field is drawn only in views that show the
  embedding it is bound to.
- Level of detail, picking and culling use the new embedding's spatial tree.
  The first time an embedding is shown it is loaded (with a progress notice)
  and its tree is built; switching back to an embedding shown before is
  immediate, because its positions and GPU buffers are kept for a while.

### Does not change

- The dataset, the cells, highlights and selections.
- The view's field, filters and colours.
- The camera position: every embedding is scaled to the same box.
- **A navigation mode you chose for this view.** A view you set to Free-fly in
  3D flies again whenever it shows a 3D embedding.

:::{tip}
If switching leaves you “lost in space” (empty-looking view), click **Reset camera**.
:::

---

## What switching actually looks like

The three frames below are one session on the Pancreas sample (which publishes
1D, 2D and 3D embeddings), switched with the **Shows** picker under
**Compare views**.

```{figure} ../../../_static/screenshots/web_app/window-dimension-3d.png
:alt: The whole window showing the Pancreas 3D embedding on the Studio backdrop, with the Navigation Mode select reading Orbit and, under Compare views, the Shows picker reading UMAP · 3D.
:width: 1440px

**3D**, the Pancreas default: **Shows** reads *UMAP · 3D* and **Navigation**
reads `Orbit`.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-2d.png
:alt: The same window after choosing UMAP · 2D in the Shows picker, with the pointer on it: the cells now lie on a plane seen face-on over a dotted backdrop, and the Navigation Mode select has followed the view to Planar.
:width: 1440px

**2D**. The coordinates changed — the cells are now
coplanar — and the camera followed. `Mode` moved to `Planar` on its own, so the embedding is seen
face-on rather than at an angle.
```

```{figure} ../../../_static/screenshots/web_app/window-dimension-1d.png
:alt: The same window after choosing UMAP · 1D in the Shows picker: the cells are spread along one axis, and the Navigation Mode select has followed the view to Planar.
:width: 1440px

**1D**. Every cell keeps its colour but has only one coordinate, so the whole
dataset becomes a coloured line — an ordering, not a map. `Mode` is `Planar`
here too: only a 3D embedding gets Orbit.
```

:::{important}
**Switching dimension switches navigation mode with it.** 3D uses **Orbit**
and 1D/2D use **Planar**, view by view, every time a view switches to another
dimension, as the frames above show.

A mode you choose for a view — **Free-fly**, or **Orbit** while it shows 2D —
is remembered for that view in that dimension, and it comes back whenever the
view shows that dimension again. Choose the dimension's own mode again and the
rule takes over.
:::

---

## Edge cases

- **The picker is disabled**: the dataset shows one embedding. Turn more on
  under **Session › Data mapping**, if it has more.
- **A key does nothing**: no embedding of that dimension is turned on (the
  notice says so), or a text field has focus.
- **Switch fails with an error**: the embedding's values were refused (for
  example a NaN); the notice names the cell, and the view keeps what it showed.
  See {doc}`../b_data_loading/13_data_mapping`.
- **The picture changed but the camera feels wrong**: you chose a navigation
  mode for this view in this dimension before. Set **Mode** to the dimension's
  default — **Planar** for 1D/2D, **Orbit** for 3D — and it follows the
  dimension again.

---

## Troubleshooting (switching embeddings)

For a full catalog, see {doc}`06_troubleshooting_core_interactions`.

### Symptom: “I can’t find 2D / the embedding is not in the picker”

**Likely cause:** it is not turned on, or the dataset does not have it.

**Fix:** turn it on under **Session › Data mapping › Embeddings**; if it is not
listed there, the dataset does not carry it (a data preparation issue).

### Symptom: “After switching, the plot looks empty”

**Likely causes:**

- camera is framed poorly for the new embedding,
- filters removed all visible points.

**Fix:**

1) Click **Reset camera**.
2) Clear filters/outlier thresholds.
