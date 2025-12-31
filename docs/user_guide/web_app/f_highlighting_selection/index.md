# Highlighting and Selection (Groups, Pages, Tools)

Selection is how you **pick cells right now** (a temporary candidate set). Highlighting is how you **persist** those picks as highlight groups inside highlight pages (so you can compare alternatives, drive analysis, and save/share via sessions).

This section is the “source of truth” for:
- what a highlight *is* (vs filters, coloring, and visibility),
- how highlight groups and highlight pages work (pages are named + colored; groups are persistent selections),
- how each selection tool behaves (especially 2D vs 3D pitfalls),
- what syncs across views/snapshots and with Python/Jupyter,
- edge cases + troubleshooting + screenshot checklist.

## Fast path

If you just want to “select a cluster and save it”:
1) Read `01_highlight_mental_model` (5–10 min, avoids most confusion)
2) Use `02_selection_tools_document_each_tool` to pick the right tool (lasso vs proximity vs KNN vs annotation-based)
3) Use `03_highlight_ui` to understand pages/groups and what **Confirm**/**Clear** actually do
4) If something feels wrong, jump to `06_troubleshooting_highlighting`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
