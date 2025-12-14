# Data Export

```{eval-rst}
.. currentmodule:: cellucid
```

Functions for exporting AnnData to the Cellucid binary format for static web hosting or sharing.

---

## prepare

Export an AnnData object to optimized binary files that can be served statically from any web server, CDN, or cloud storage (S3, GCS, etc.).

```python
from cellucid import prepare

# Export to local directory
prepare(adata, "./my_export")

# The export directory can then be:
# - Served locally with cellucid.serve("./my_export")
# - Uploaded to a static host (GitHub Pages, Netlify, Vercel)
# - Stored in cloud storage (S3, GCS) with public access
```

```{eval-rst}
.. autofunction:: prepare
```

---

## Output Format

The exported directory contains optimized binary files:

```
my_export/
├── metadata.json      # Dataset metadata and schema
├── positions.bin      # Cell coordinates (float32)
├── categories/        # Categorical observation data
│   └── *.bin
├── continuous/        # Continuous observation data
│   └── *.bin
└── genes/             # Gene expression data (sparse)
    └── *.bin
```

---

## See Also

- {func}`~cellucid.show` - Display exported data in Jupyter
- {func}`~cellucid.serve` - Serve exported data via HTTP
