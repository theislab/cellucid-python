# Makefile for Cellucid documentation

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument runs help
help:
	@echo "Cellucid Documentation Build Targets"
	@echo ""
	@echo "  html       Build HTML documentation"
	@echo "  clean      Remove build artifacts"
	@echo "  serve      Build and serve locally with auto-reload"
	@echo "  linkcheck  Check external links"
	@echo "  doctest    Run doctests in documentation"
	@echo ""

.PHONY: help html clean serve linkcheck doctest

html:
	@echo "Building HTML documentation..."
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS)
	@echo ""
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

clean:
	rm -rf $(BUILDDIR)
	rm -rf user_guide/python_package/api/_autosummary
	@echo "Cleaned build artifacts."

serve: html
	@echo "Serving documentation at http://localhost:8000"
	cd $(BUILDDIR)/html && python -m http.server 8000

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(SOURCEDIR) $(BUILDDIR)/linkcheck $(SPHINXOPTS)

doctest:
	$(SPHINXBUILD) -b doctest $(SOURCEDIR) $(BUILDDIR)/doctest $(SPHINXOPTS)

# Convenience target for ReadTheDocs
.PHONY: livehtml
livehtml:
	sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) --port 8000 --open-browser
