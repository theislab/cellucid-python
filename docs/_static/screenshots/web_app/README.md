# Web app screenshots

This folder is the recommended destination for screenshots used in the **Web App Guide**.

Suggested workflow:

1) Capture screenshots at ~1200–1600px width.
2) Prefer PNG for UI clarity (avoid JPG artifacts on text).
3) Redact anything sensitive (dataset ids/names, private repo/org names, GitHub usernames if needed).
4) Keep aspect ratios consistent across a page when possible.

The docs currently render a generic placeholder SVG at:

- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

In the guide pages, each screenshot block includes an HTML comment describing:

- what to capture,
- what to highlight/crop,
- what to use as the figure caption.

Replace the placeholder by:

- adding your real screenshot file into this folder, and
- updating the `{figure}` path in the corresponding doc page.

