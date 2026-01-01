from __future__ import annotations


def test_web_cache_resolves_relative_js_imports() -> None:
    from cellucid.web_cache import _extract_candidate_paths_for_file

    text = """
    /* JSDoc-style type import that should be ignored:
       import('./plugin-contract.js')
    */
    import { foo } from '../data/data-source.js';
    async function load() {
      const module = await import('../../external/hdf5_hl.js');
      const u = new URL('./worker.js', import.meta.url);
      return module && u;
    }
    """

    paths = _extract_candidate_paths_for_file(text, base_url_path="/assets/js/data/h5ad.js")

    assert "/assets/js/data/data-source.js" in paths
    assert "/assets/external/hdf5_hl.js" in paths
    assert "/assets/js/data/worker.js" in paths
    # Should not treat comment-only type imports as real requests.
    assert "/assets/js/data/plugin-contract.js" not in paths


def test_web_cache_resolves_css_urls() -> None:
    from cellucid.web_cache import _extract_candidate_paths_for_file

    css = "body { font-family: url('../fonts/inter-latin.woff2'); }"
    paths = _extract_candidate_paths_for_file(css, base_url_path="/assets/css/main.css")
    assert "/assets/fonts/inter-latin.woff2" in paths


def test_web_cache_extracts_html_src_href() -> None:
    from cellucid.web_cache import _extract_candidate_paths_for_file

    html = '<script src="assets/js/app/main.js"></script><link href="/assets/css/main.css" />'
    paths = _extract_candidate_paths_for_file(html, base_url_path="/index.html")
    assert "/assets/js/app/main.js" in paths
    assert "/assets/css/main.css" in paths

