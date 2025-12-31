# Accessibility, Privacy, and Security

These pages make explicit what Cellucid **does and does not do** in two areas that often get skipped in visualization tooling:

- **Accessibility**: how to produce figures and workflows that work for diverse readers (color vision deficiency, low vision, keyboard-only use), plus current limitations of a WebGL canvas-based UI.
- **Privacy model**: what stays local, what can touch the network depending on your loading workflow, what is stored in the browser, and what exported artifacts can contain.

If you are in a regulated environment (clinical, IRB-controlled, corporate IP), the key idea is:

> Cellucid is “just a web app”, but your *workflow choices* determine what leaves your machine and what gets saved into shareable artifacts.

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```

---

## Fast path (pick your situation)

| You need to… | Do this first (safe default) | Then read |
|---|---|---|
| Use sensitive/clinical data without uploading it | Use **local prepared exports** + avoid Community Annotation + review figure metadata before sharing | {doc}`02_privacy_model` |
| Share sessions safely with a team | Treat `.cellucid-session` files like sensitive derived artifacts; avoid patient IDs in labels | {doc}`../l_sessions_sharing/08_security_privacy_and_trust` |
| Make figures readable for colorblind audiences | Use perceptually-uniform colormaps; check the export preview’s colorblind simulation | {doc}`01_accessibility` + {doc}`../k_figure_export/04_quality_knobs_and_best_practices` |

---

## Recommended reading order

1) {doc}`02_privacy_model` (know what leaves your machine and what gets saved)
2) {doc}`01_accessibility` (make what you publish readable to more people)

Related (often relevant):
- {doc}`../k_figure_export/05_metadata_and_provenance` (export metadata can include dataset paths/URLs)
- {doc}`../l_sessions_sharing/08_security_privacy_and_trust` (session bundle privacy/trust)
- {doc}`../a_orientation/02_system_requirements` (managed environments, blocked storage, iframe restrictions)

---

## Pages in this section

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`accessibility;1.5em;sd-mr-1` Accessibility
:link: 01_accessibility
:link-type: doc

Color/contrast guidance, keyboard shortcuts, motion sensitivity tips, and practical limits of a WebGL viewer.
:::

:::{grid-item-card} {octicon}`lock;1.5em;sd-mr-1` Privacy model
:link: 02_privacy_model
:link-type: doc

What stays local vs hits the network, what is stored in the browser, and what exported artifacts can contain.
:::

::::

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
