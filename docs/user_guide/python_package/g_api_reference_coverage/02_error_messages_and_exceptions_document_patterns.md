# Error messages and exceptions (documentation patterns)

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
- print a one-line actionable message to stderr

---

## Message structure template (recommended)

Use messages that fit this shape:

1) **What**: the object/parameter that is invalid
2) **Expected**: shape/type/range
3) **Got**: actual shape/type/value
4) **Fix hint** (optional but recommended): what to change

Example (shape mismatch):
- “`X_umap_2d` must have exactly 2 columns, got 3. Shape is (10000, 3).”

Example (missing key):
- “Declare one or more exact embedding keys in `obsm`: `X_umap_1d`,
  `X_umap_2d`, or `X_umap_3d`.”

---

## Common error categories to document (by feature)

### Export (`prepare`)

Document these explicitly:
- missing embeddings (`X_umap_1d/2d/3d`)
- mismatched `n_cells` across inputs (`obs`, embeddings, `gene_expression`, `latent_space`)
- embedding keys whose suffix and column count disagree
- NaN/Inf handling in continuous quantization
- text the viewer prints verbatim -- category labels, `dataset_name`,
  `dataset_description`, `source_name`, `source_url`, `source_citation` --
  carrying a character with no glyph
- “export folder huge / slow export” performance guidance

Display text is a family of its own because the defect is invisible. A label
with a trailing space renders exactly like the clean label beside it, so the
message has to make the character visible and say why the writer will not
remove it:

```text
Categorical field 'organ' has 2 labels the viewer cannot show as written:
  - 'Gut ' ends with U+0020 SPACE
  - 'Liver ' ends with U+0020 SPACE
Labels are drawn verbatim in the legend, the field selector, and exported
figures, so these characters are invisible on screen while still making the
label a different value from the one it looks like. Cellucid does not clean
them for you: trimming a label rewrites your annotation and can merge two
categories into one, moving cells between them. Clean the column and export
again, for example:
    obs['organ'] = obs['organ'].astype('string').str.strip()
    obs['organ'] = obs['organ'].astype('category')
```

Every offending label in the field is listed in one message (up to ten, then a
count of the rest), so one edit fixes the column. Two labels that a
whitespace-collapsing renderer draws identically get their own message naming
both. Identity arguments use the same vocabulary, opening with the argument
name: `source_name is displayed verbatim, so it must not carry characters that
have no glyph: 'Suo\u200b' contains U+200B ZERO WIDTH SPACE. ...`

Where documented:
- {doc}`api/export`

### Jupyter viewers and hooks

Document:
- “viewer does not appear” (notebook context / iframe restrictions)
- “browser cannot reach localhost” (explicit proxy or tunnel required)
- “events not arriving” (viewer not ready / wrong viewer id / blocked POST)
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
- invalid/oversized manifest
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
- the most common fix path (e.g., “export with `prepare(...)` for better performance”)
- any relevant constraints (for example, “Cellucid supports 1D, 2D, and 3D
  embeddings”)

When the failure is a security boundary (sessions, server exposure):
- clearly state the risk and safe alternative (e.g., SSH tunnel instead of binding `0.0.0.0`)
