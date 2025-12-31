# Benchmarking and performance screenshots

This folder is the recommended destination for screenshots used in:
- `cellucid-python/docs/user_guide/web_app/n_benchmarking_performance/`

Suggested workflow:

1) Capture screenshots at ~1200–1600px width (or higher if you will crop heavily).
2) Prefer PNG for UI clarity (avoid JPG artifacts on text).
3) Redact anything sensitive (dataset IDs/names, local file paths, private repo/org names).
4) Keep aspect ratios consistent across a page when possible.

The docs currently render a generic placeholder SVG at:
- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

In the guide pages, each screenshot block includes an HTML comment describing:
- what to capture,
- what to highlight/crop,
- and suggested alt text + caption content.

To replace placeholders:

1) Add your real screenshot file(s) into this folder, and
2) Update the `{figure}` path in the corresponding doc page.

