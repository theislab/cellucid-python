# API Reference Coverage

`cellucid-r` exports the writer `cellucid_prepare()`, which writes a
Cellucid-ready AnnData (or MuData) store, `cellucid_transition_drift()`, which
computes a vector field for it from a cell-cell transition matrix, the viewer
(`cellucid_view()`, `cellucid_serve()` and the functions that talk to an open
viewer), `cellucid_apply_session()` and `cellucid_clear_web_cache()`. This
section is the flat reference for them: every exported function and the
complete argument lists on one page, and the error messages on another. The
narrative explanation of each argument lives in
{doc}`../c_data_preparation_api/index` and {doc}`../d_viewing_loading/index`.

```{note}
The R help pages (`?cellucid_prepare`, `?cellucid_transition_drift`,
`?cellucid_view`, `?cellucid_highlight`, `?cellucid_on`,
`?cellucid_apply_session`, `?cellucid_debug_connection`,
`?cellucid_clear_web_cache`) are hand-written and authoritative, and
`cellucid-r`'s own test suite rebuilds each `\usage` block into a function and
compares it to the real signature. The two pages here restate them for the
web, with links into the guide.
```

**Recommended reading order**

1) {doc}`01_public_functions_and_classes`
2) {doc}`02_error_messages_and_exceptions_document_patterns`

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Public Functions
:link: 01_public_functions_and_classes
:link-type: doc

Every exported function; every argument of `cellucid_prepare()`, grouped by
what it is for, with its default, where it lands in the store, and the page
that explains it; `cellucid_transition_drift()` and the viewer's signatures.
:::

:::{grid-item-card} {octicon}`alert;1.5em;sd-mr-1` Error Messages
:link: 02_error_messages_and_exceptions_document_patterns
:link-type: doc

How the package reports problems, and every message it raises, grouped by the
input it is about.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
