# API Reference Coverage

`cellucid-r` exposes exactly one public function, `cellucid_prepare()`. This
section is the flat reference for it: the complete argument list on one page,
and the error vocabulary on another. The narrative explanation of each argument
lives in {doc}`../c_data_preparation_api/index`.

```{note}
The R help page (`?cellucid_prepare`) is hand-written and authoritative, and
`cellucid-r`'s own test suite rebuilds its `\usage` block into a function and
compares it to the real signature. The two pages here restate it for the web,
with links into the guide.
```

**Recommended reading order**

1) {doc}`01_public_functions_and_classes`
2) {doc}`02_error_messages_and_exceptions_document_patterns`

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`list-unordered;1.5em;sd-mr-1` Public Functions and Objects
:link: 01_public_functions_and_classes
:link-type: doc

Every argument of `cellucid_prepare()`, grouped by what it is for, with its
default and the page that explains it.
:::

:::{grid-item-card} {octicon}`alert;1.5em;sd-mr-1` Error Messages and Exceptions
:link: 02_error_messages_and_exceptions_document_patterns
:link-type: doc

How the package reports problems, the common error families, and where to debug
each one.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
