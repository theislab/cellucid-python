# Community annotation screenshots

This folder is intentionally checked in as the recommended destination for user-guide screenshots.

Suggested workflow:

1) Capture screenshots at ~1200–1600px width.
2) Prefer PNG for UI clarity (avoid JPG artifacts on text).
3) Redact anything sensitive (dataset names, private repo names, GitHub usernames if needed).
4) Keep aspect ratios consistent across the guide when possible.

The docs currently render a generic placeholder SVG at:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

In the guide pages, each screenshot block includes an HTML comment describing:

- what to capture,
- what to highlight/crop,
- what to use as the figure caption.

Replace the placeholder by:

- adding your real screenshot file into this folder, and
- updating the `{figure}` path in the corresponding doc page.
