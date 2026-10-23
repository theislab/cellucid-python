"""
Writing a Cellucid-ready AnnData from arrays and data frames.

:func:`prepare` is the package's public surface. A Cellucid-ready AnnData is a
standard AnnData (``.h5ad``, a ``.zarr`` folder, or a ``.zarr.zip``) laid out
for the web app's range reads, with ``uns['cellucid']`` describing the dataset;
cellucid-r's ``cellucid_prepare()`` writes the same. The private submodules:

``_inputs``
    Checks every input and converts it into the elements of the AnnData.
``_store``
    Writes those elements with anndata's encodings, the chunks and codecs the
    browser reads fastest, and moves the finished output into place.
``_prepare``
    :func:`prepare` itself.
"""

from ._prepare import prepare

__all__ = ["prepare"]
