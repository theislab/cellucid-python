# Landing Pages

These pages answer the “what is this?” and “how do I get started?” questions for the **R export workflow**.

This documentation is intentionally written for mixed audiences:
- **Wet lab / non-technical users**: step-by-step “do this, then that” instructions.
- **Computational users**: explicit shapes/dtypes, performance considerations, and edge cases.
- **Developers/maintainers**: file format details and reproducibility notes.

**Recommended reading order**

1) {doc}`01_what_is_cellucid_r`
2) {doc}`04_quick_start_3_levels`
3) {doc}`03_supported_inputs_and_workflows`
4) {doc}`02_installation` (or earlier if you can’t install)

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} {octicon}`info;1.5em;sd-mr-1` What is cellucid-r?
:link: 01_what_is_cellucid_r
:link-type: doc

What the R package does today, what it doesn’t do yet, and how it fits into the Cellucid ecosystem.
:::

:::{grid-item-card} {octicon}`zap;1.5em;sd-mr-1` Quick Start (3 Levels)
:link: 04_quick_start_3_levels
:link-type: doc

Copy/paste: a tiny toy export, a “real dataset” export, and an advanced export with performance options.
:::

:::{grid-item-card} {octicon}`checklist;1.5em;sd-mr-1` Supported Inputs & Workflows
:link: 03_supported_inputs_and_workflows
:link-type: doc

Exactly what inputs `cellucid_prepare()` expects (shapes, dtypes, conventions), plus common gotchas.
:::

:::{grid-item-card} {octicon}`download;1.5em;sd-mr-1` Installation
:link: 02_installation
:link-type: doc

How to install from GitHub, optional dependencies, and installation troubleshooting.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
