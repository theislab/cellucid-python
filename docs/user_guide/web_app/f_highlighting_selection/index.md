# Highlighting and Selection (Groups, Pages, Tools)

Selection is how you **pick cells right now** (a temporary candidate set). Highlighting is how you **persist** those picks as highlight groups inside highlight pages (so you can compare alternatives, drive analysis, and save/share via sessions).

This section is the “source of truth” for:
- what a highlight *is* (vs filters, coloring, and visibility),
- how highlight groups and highlight pages work (pages are named + colored; groups are persistent selections),
- how each selection tool behaves (especially 2D vs 3D pitfalls),
- what syncs across views/snapshots and with Python/Jupyter,
- edge cases, troubleshooting, and a verified highlighting capture.

## Fast path

If you just want to “select a cluster and save it”:
1) Read {doc}`01_highlight_mental_model` (5–10 min, avoids most confusion)
2) Use {doc}`02_selection_tools_document_each_tool` to pick the right tool (lasso vs proximity vs KNN vs annotation-based)
3) Use {doc}`03_highlight_ui` to understand pages/groups and what **Confirm**/**Clear** actually do
4) If something feels wrong, jump to {doc}`06_troubleshooting_highlighting`

## Every page in this section

| Page | Read it when |
|---|---|
| {doc}`01_highlight_mental_model` | You want to know how highlighting differs from filtering and colouring. |
| {doc}`02_selection_tools_document_each_tool` | You are choosing between lasso, proximity, KNN and annotation-based selection. |
| {doc}`03_highlight_ui` | You need the button-by-button map, including `Create Categorical`. |
| {doc}`04_selection_synchronization` | You work across snapshots, or drive selection from Python/Jupyter. |
| {doc}`05_edge_cases_highlighting` | Something behaved surprisingly but not wrongly (3D lasso, overlapping groups, huge pages). |
| {doc}`06_troubleshooting_highlighting` | Something is broken and you want symptom → fix. |
| {doc}`07_screenshots` | You want to see the panel as it actually looks, with verified captures. |

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
