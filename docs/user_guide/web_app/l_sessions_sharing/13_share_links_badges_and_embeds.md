# Share links, badges and embeds

**Audience:** anyone publishing a dataset for others to open  
**Time:** 3 minutes  
**What you'll learn:**
- Which datasets can be shared with a link
- How to put an "Open in Cellucid" badge in a README
- How to embed the viewer in a web page

---

## Which datasets can be shared

A shared link reopens the dataset for someone else, so it only exists when the
dataset lives somewhere they can reach:

| Where the dataset comes from | Shareable |
| --- | --- |
| A built-in sample dataset | Yes |
| A dataset opened from a link (GitHub, Zenodo, CELLxGENE, GEO, Figshare, Dropbox, a web address, …) | Yes; a Google Drive or OneDrive file opens for others only as far as its sharing allows (a private file after they sign in) |
| A link to a web address on this computer (`localhost`) | No |
| A Cellucid server on another machine (`cellucid serve` on a reachable host) | Yes |
| A Cellucid server on this computer (`localhost`) | No: only this computer can reach it |
| A file or folder opened from disk | No: it never leaves your browser |
| A Jupyter notebook session | No |

For a shareable dataset, a **Share** button appears beside the **Dataset**
heading in the **Session** panel. It opens **Share this dataset** with three
snippets, each with a **Copy** button.

The snippets reopen the *dataset*, not your view of it. To share a camera,
colouring, filters and highlights as well, save a session file and send it
with the link; see {doc}`05_share_workflows_links_sessions_data_files`.

## Link

A plain link to the public site, for example
`https://www.cellucid.com/?data=gh:owner/repo/exports&dataset=atlas`. It points
at www.cellucid.com even when you are running your own copy of the app, and it
carries only the parameters that pick the dataset.

## Open in Cellucid badge

Markdown for a README, a wiki page or an issue:

```markdown
[![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:owner/repo/exports&dataset=atlas)
```

It renders as a small badge that opens the dataset in one click — the same
idea as an "Open in Colab" badge. Use it in the repository that hosts your
exports (see {doc}`../b_data_loading/11_custom_dataset_repository`).

## Embed in a web page

HTML for a lab website, a project page or a paper's supplementary site:

```html
<iframe src="https://www.cellucid.com/?data=gh:owner/repo/exports&dataset=atlas"
        title="Cellucid: atlas" width="100%" height="640"
        style="border: 0;" allow="fullscreen" loading="lazy"></iframe>
```

- **Size.** Full width and 640 px tall suits most pages. In a frame narrower
  than 720 px the viewer uses its phone layout: the menu covers the whole frame
  while it is open, and the plot has the frame to itself when it is closed.
- **Requirements.** The visitor's browser needs WebGL2, as for the app itself.
  `loading="lazy"` keeps the viewer from loading until the frame scrolls into
  view, so a long page stays fast.
- **Welcome dialog.** A dataset opened from a link or a server opens straight
  into the viewer. A built-in sample shows the short welcome dialog first, as it does
  on www.cellucid.com.
