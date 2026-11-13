# Citing Cellucid

<!-- BEGIN cellucid-citation (generated from CITATION.cff in theislab/cellucid; do not edit) source-sha256=763df407bd03460165ac0110d8abcf61a45268dd79420297ae5ad9263dc0ab13 body-sha256=17bc575ddc1919c06fb46e6919867c737c9d11902c2d597e90222fd31a15b152 -->
```{important}
**Placeholder citation.** The Cellucid preprint is in preparation for bioRxiv and has no DOI yet; the reference below is a placeholder and will be replaced when the preprint is published.
```

If you use Cellucid — the web app, the Python package, or the R package — please cite:

Inecik, K. PLACEHOLDER: Cellucid preprint, title to be announced (forthcoming on bioRxiv). *bioRxiv* (in preparation).

```bibtex
@unpublished{cellucid,
  author = {Inecik, Kemal},
  title  = {{PLACEHOLDER: Cellucid preprint, title to be announced (forthcoming on bioRxiv)}},
  note   = {bioRxiv preprint in preparation. Placeholder citation without a DOI; it will be replaced when the preprint is published.},
}
```
<!-- END cellucid-citation -->

In the web app, **Cite** in the sidebar footer opens this reference with
**Copy reference** and **Copy BibTeX**
({doc}`user_guide/web_app/a_orientation/08_cite_cellucid`).

## Other things you may need to cite

- **A built-in sample dataset.** Cite the study it came from, not Cellucid.
  {doc}`user_guide/web_app/b_data_loading/12_sample_dataset_provenance` names
  each sample's study and the citation to use.
- **The exact package version you used.** Report it next to the reference above:
  `cellucid.__version__` in Python, and `citation("cellucid")` in R, which lists
  the R package with its version after the reference above.

## Where this citation comes from

Every Cellucid repository's `CITATION.cff` carries the same `preferred-citation`,
so GitHub's "Cite this repository" button shows this reference too. It is written
in exactly one place, `CITATION.cff` in
[theislab/cellucid](https://github.com/theislab/cellucid), and copied from there,
together with the README citations, this page, and the R package's citation
entry, by the citation sync script. Each repository's test
suite fails if any of those copies drifts from its own `CITATION.cff`.
