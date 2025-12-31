# Vector field / velocity overlay screenshots

This folder is the recommended destination for screenshots used in:
- `cellucid-python/docs/user_guide/web_app/i_vector_field_velocity/`

Suggested workflow:

1) Capture screenshots at ~1200–1600px width (or higher if you will crop heavily).
2) Prefer PNG for UI clarity (avoid JPG artifacts on text).
3) Redact anything sensitive (dataset ids/names, private repo/org names, local file paths).
4) Keep aspect ratios consistent within a page when possible (step-by-step visuals feel cleaner).

The docs currently render a generic placeholder SVG at:
- `cellucid-python/docs/_static/screenshots/placeholder-screenshot.svg`

Each screenshot placeholder block in the docs includes an HTML comment describing:
- what to capture,
- what state is required,
- how to crop/annotate,
- suggested alt text and caption.

To replace placeholders:

1) Add your real screenshot file(s) into this folder, and
2) Update the `{figure}` path in the corresponding doc page.
