# Jupyter Hooks System (Python ↔ Frontend)

This chapter documents the **Jupyter hooks system**: how your notebook Python code can **drive** the Cellucid viewer (commands) and **react** to what a user does in the viewer (events).

Context (so the ecosystem is clear):
- **Cellucid** is the **web app** (the interactive viewer UI).
- **cellucid-python** is the **helper package** you use from **Python/CLI** to prepare data, run servers, embed the viewer, and integrate with notebooks.
- **cellucid-annotation** is a helper repo for community annotation workflows.
- **cellucid-r** will exist in the future (not ready yet).

## Where to start (pick one)

- “I want a minimal working loop (select → highlight)” → {doc}`02_quickstart_minimal_roundtrip`
- “Does this work in Colab / VSCode / JupyterHub?” → {doc}`03_supported_environments_matrix`
- “I need the mental model and vocabulary” → {doc}`04_glossary_and_mental_model`
- “I need message + event schemas” → {doc}`16_reference`
- “Something is broken” → {doc}`14_troubleshooting_hooks`

```{toctree}
:maxdepth: 1
:hidden:
:glob:

[0-9]*
```
