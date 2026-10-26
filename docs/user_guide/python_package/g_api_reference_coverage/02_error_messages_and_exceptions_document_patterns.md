# Error messages and exceptions

This page defines **how we document errors** across `cellucid-python`.

Goal: when an error happens, the user should be able to answer:
- What failed?
- Why did it fail?
- What should I do next?

This matters for mixed audiences:
- Wet lab / non-technical users need plain-language fixes.
- Computational users need exact shapes/keys/types to debug quickly.
- Developers need consistent exception types and message structure.

---

## Exception taxonomy (what to raise, and how to describe it)

Use these consistently and document them the same way:

| Exception | When it should happen | What the message must include |
|---|---|---|
| `TypeError` | wrong type (e.g. `vector_fields` not a dict) | expected type + received type |
| `ValueError` | wrong value/shape/range | expected shape/range + actual shape/value |
| `KeyError` | missing key in mapping (`adata.obsm`, manifest, chunk id) | missing key + where it was searched |
| `FileNotFoundError` | missing file/dir path | the resolved path |
| `TimeoutError` | waiting for UI/event too long | which event + how long waited |
| `RuntimeError` | unexpected runtime failure (e.g., session capture returned error) | short reason + next action |

For CLI:
- return non-zero exit codes (typically `1` for failure, `130` for Ctrl+C)
- print a one-line actionable message to stderr, beginning `Error:`, as the last
  line the command writes
- print **no traceback** for the condition types (`ValueError`, `TypeError`,
  `KeyError`, `RuntimeError`, `ImportError`, and every `OSError`, which
  includes `FileNotFoundError`, `TimeoutError`, a taken port and an unbindable
  host): they are the conditions a user can correct, and a traceback tells a wet-lab reader nothing
  they can act on
- print the **full traceback** for every other exception type — `AttributeError`,
  `IndexError`, `AssertionError` and the rest — because those mean Cellucid
  failed, and the traceback is the bug report; the `Error:` line then names the
  issue tracker

`--verbose` adds the traceback to the first case as well, so a condition that
turns out to have been misclassified is still reportable. See
{doc}`../h_developer_docs/05_cli_architecture_and_commands` for the exact
classification.

---

## Message structure template (recommended)

Use messages that fit this shape:

1) **What**: the object/parameter that is invalid
2) **Expected**: shape/type/range
3) **Got**: actual shape/type/value
4) **Fix hint** (optional but recommended): what to change

Example (shape mismatch):
- “Vector field 'velocity' on embedding 'X_umap' must have exactly 2 columns, got 3.”

Example (missing input):
- “prepare() needs embeddings=: at least one {obsm key: (n_cells, k) array},
  for example embeddings={'X_umap': adata.obsm['X_umap']}, or data= an AnnData
  or MuData that holds them.”

---

## Common error categories to document (by feature)

### Data preparation (`prepare`)

Document these explicitly:
- missing embeddings (an empty or absent `embeddings=` mapping)
- mismatched `n_cells` across inputs (`obs`, embeddings, `gene_expression`,
  `connectivities`, vector fields)
- NaN/Inf in embeddings, expression, the graph and vector fields, counted with
  the first affected positions
- keys and column names that cannot be a path in a store on every platform
- text the viewer prints verbatim -- `dataset_name`, `dataset_description`,
  `citation_text`, `source_name` -- carrying a character with no glyph
- an existing output (`force=True`), and a folder that is not a Zarr store

Display text is a family of its own because the defect is invisible. A name
with a trailing space renders exactly like the clean name beside it, so the
message makes the character visible and says why the writer will not remove
it, opening with the argument name: `dataset_name is displayed verbatim, so it
must not carry characters that have no glyph: 'Suo\u200b' contains U+200B ZERO
WIDTH SPACE. ...`

Where documented:
- {doc}`api/prepare`
- {doc}`../c_data_preparation_api/11_troubleshooting_prepare`

### Jupyter viewers and hooks

Document:
- “viewer does not appear” (notebook context / iframe restrictions)
- “browser cannot reach localhost” (explicit proxy or tunnel required)
- “events not arriving” (viewer not ready / its channel cannot open: a proxy
  that does not pass WebSockets, or a page origin the server refuses)
- session capture timeouts

Where documented:
- {doc}`api/jupyter`, {doc}`api/viewers`, {doc}`api/sessions`

### Servers (CLI + Python)

Document:
- port collisions
- remote server access confusion (localhost vs remote)
- configured viewer-source generation failures

Where documented:
- {doc}`api/server`, {doc}`api/cli`

### Session bundles

Document:
- invalid magic header (wrong file)
- invalid/oversized bundle manifest
- decompression limits (zip-bomb protection)
- exact dataset fingerprint mismatch rejection

Where documented:
- {doc}`api/sessions`

---

## Troubleshooting entry template (copy/paste)

Use this exact structure for each entry:

### Symptom
What the user sees (including exact error text when possible).

### Likely causes (ordered)
3–7 causes; each should be testable.

### How to confirm
Concrete checks (what to print / what URL to visit / what key to inspect).

### Fix
Step-by-step; “safe fixes first”.

### Prevention
How to avoid it next time (data conventions, naming, recommended workflow).

---

## Edge-case messaging guidelines (what to include)

When the failure is a common user mistake, the message should also mention:
- the most common fix path (e.g., “write the file with `prepare(...)` for faster reads”)
- any relevant constraints (for example, “Cellucid supports 1D, 2D, and 3D
  embeddings”)

When the failure is a security boundary (sessions, server exposure):
- clearly state the risk and safe alternative (e.g., SSH tunnel instead of binding `0.0.0.0`)
