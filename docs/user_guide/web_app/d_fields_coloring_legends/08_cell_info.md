# Cell info

**Audience:** everyone  
**Time:** 2 minutes  
**What you'll learn:**
- How to see the values of the cell under the pointer
- What the card shows, in one view and in a grid of views
- Why it stays fast on millions of cells

---

## Turning it on

**Show cell info** is in **Coloring & filtering**, under the field selectors
and above **Legend**. It is on by default; untick it to turn the card off. A
saved session remembers it.

Rest the pointer on a cell and a card appears beside it. Move to
another cell and the card follows. Move off the cells, or start moving the
camera, and it goes away.

On a touch screen, press and hold a cell for about half a second. The card
stays until your next touch, so a quick tap still selects or moves the camera
as before.

## What the card shows

From top to bottom:

- **The value the view is coloured by**, beside a dot in the cell's own colour.
  For a category it is the category label; for a gene or a continuous field it
  is the number, to four significant figures. The field's name is under it.
- **Where that value sits on the colour scale**, for a continuous field or a
  gene: the active colour map with a marker at the cell's value and the scale's
  ends below it. A log colour scale is placed on the log axis.
- **How many cells of that category are visible in the view**, for a
  categorical field: exactly the cells that selecting the category there would
  pick, so filters, hidden categories and the outlier filter are all taken into
  account. In a grid of views, the name of the view the pointer is over.
- **The cell's other categories**, only when **Show recent fields** (under
  **Show cell info**) is ticked: every other categorical field that is
  already loaded, each with its colour dot. A cell with no value for a field
  shows *missing*. When there are more fields than the window is tall, the card
  adds columns; only when even the window's full width is used up does it name
  how many more fields there are instead of showing them.

Each view of a grid is read on its own terms: the card describes the field and
colours of the view under the pointer, which may differ from the view you are
editing. Hidden cells (filtered out, or outside the current level of detail) are
never picked.

## Cost on large datasets

The card is designed to cost nothing you can see, even on millions of cells:

- **It never loads anything.** It reads only fields that are already in
  memory. A field that is not loaded is not shown, rather than fetched.
- **It never redraws the plot.** The card is a page element moved over the
  canvas; the scene is not rendered again for it.
- **A lookup only when the pointer rests.** Moving the pointer only records
  where it is; the cell is looked up once the pointer has been still for about
  60 ms, like a tooltip. Sweeping across the plot does no work at all, so frame
  times while moving are the same with the card on or off.
- **Nothing while the camera moves.** During a drag, a touch gesture, a held
  navigation key, inertia or a zoom that is still easing, there is no lookup
  and no card. It returns when the camera settles.
- **The index is built in the background.** On a view of 250,000 cells or
  more the lookup walks a spatial index. It is built in a worker, off the main
  thread, so the page never pauses for it (on ten million cells the first hover
  answers in about a millisecond and the card appears once the index is ready).
  Level of detail, frustum culling and the selection tools use the same index.
  The worker receives the positions in slices, so the page never holds a second
  full copy of them, at any size. On a smaller view the lookup scans its cells
  directly and nothing is built.
- **Counts are made in idle time.** The per-category count is one pass over the
  view's cells (about 25 ms at ten million), run when the page is idle rather
  than on a hover, and kept until that view's colouring or visibility changes.
  The number appears on the card a moment after the card itself.

See {doc}`../n_benchmarking_performance/02_performance_considerations_what_gets_slow_and_why`
for what does make a large view slow.
