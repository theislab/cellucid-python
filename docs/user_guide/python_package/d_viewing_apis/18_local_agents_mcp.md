# Local agents (MCP): Claude Code, Claude Desktop, Cursor

`cellucid mcp` lets an agent app on your computer work in a Cellucid viewer:
Claude Code, Claude Desktop, Cursor, or any app that takes a local MCP server.
It is a local [MCP](https://modelcontextprotocol.io) server over stdio. The app
starts it; it drives a viewer the `cellucid` Python or R package serves on the
same computer. No relay, no deployed site and no account are involved: the
app talks to `cellucid mcp`, and `cellucid mcp` talks to the viewer page over
the viewer's own channel.

The agent gets **the same tools and the same permissions** as Claude and
ChatGPT get through the Cellucid connector: every action of the viewer, each
described the same way, and every change runs in the viewer as one step of
your history, labeled with the app's name ("Claude Code"). You can stop the
agent, and undo its whole turn, from the viewer.

## At a glance

**Audience**
- Wet lab / non-technical: set it up once (below), then ask your agent app to
  look at your data.
- Computational: drive a notebook's viewer from Claude Code while you work in
  the notebook.

**Prerequisites**
- `pip install 'cellucid[mcp]'` (the `mcp` extra brings the official MCP SDK)
- An agent app that runs local MCP servers

## Set it up once

The app starts `cellucid mcp` itself whenever it needs it; you only tell it
the command once. Use the `cellucid` of the environment that has your data's
packages (an absolute path, such as `/path/to/env/bin/cellucid`, when the app
does not see that environment's `PATH`).

**Claude Code**

```bash
claude mcp add cellucid -- cellucid mcp
```

**Claude Desktop**: Settings > Developer > Edit Config, then in
`claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "cellucid": { "command": "cellucid", "args": ["mcp"] }
  }
}
```

**Cursor**: in `.cursor/mcp.json` (this project) or `~/.cursor/mcp.json`
(every project), the same `mcpServers` entry as Claude Desktop's.

To have the app open a dataset itself, add `--open` and a path:
`"args": ["mcp", "--open", "/path/to/data.h5ad"]`.

## Two ways to a viewer

**Drive a viewer that is already open** (the default). Open your data as you
do now, then ask the agent:

```bash
cellucid serve data.h5ad                     # a terminal: open the Viewer URL it prints
```

```python
import cellucid
viewer = cellucid.show_anndata(adata)        # a notebook cell
# or, in a script: cellucid.AnnDataViewer(adata)
```

```r
cellucid::cellucid_view(seurat_object)       # R (a terminal, RStudio, Positron)
```

Every such viewer is published to agents on this computer: a small record in
your user cache (`~/Library/Caches/cellucid/agents` on macOS,
`%LOCALAPPDATA%\cellucid\agents` on Windows, `~/.cache/cellucid/agents`
elsewhere; `CELLUCID_AGENTS_DIR` sets another), readable by you only, naming
the computer it is on and a secret the viewer's process checks. On macOS and
Linux the folder must be yours and writable by no one else, or nothing is
published there or read from it (one sentence in the log, and in the agent's
answer, says so). `cellucid mcp` drives the newest viewer whose page is open,
and moves to another when the page of the one it drives closes while another
viewer's is open; `cellucid mcp --viewer pbmc3k` names one by its dataset id,
dataset name or viewer id. Viewers opened after the app started are found
when the agent next calls a tool. Set `CELLUCID_AGENTS=off` before opening a
viewer to keep it from agents altogether.

The agent is never given a viewer's page address (it carries the viewer's
secret): it is told where the viewer is served, and asks you to open its page
again when none shows it. If a page is reloaded or closed while the agent
waits on a call, the agent is told at once that whether the call was done is
not known, and reads the history before trying again.

**Open a dataset for the agent**: `cellucid mcp --open data.h5ad` serves the
file (any format `cellucid serve` reads) and opens its viewer in your
browser; the viewer stops when the app ends the server. `--no-browser`
prints the address on stderr instead, `--port` chooses the port.

Not reachable: a page opened at a server's bare `/?anndata=true` address (no
viewer channel: open the printed Viewer URL instead), and an R viewer shown
in a Jupyter notebook's output (the R kernel runs no server between cells).

## What you see, and what you decide

- **Who is acting.** While the agent works, the viewer shows one line,
  "Claude Code (local) is working", with **Stop**; once done, its steps and
  **Undo turn**, which takes the agent's whole turn back.
- **Your permissions decide.** The agent permissions you set in the viewer
  (Agents panel > Permissions > Change…) apply to it as to every agent: changes run,
  heavy work (longer than about ten seconds here), writing a file, deleting
  what was saved, a large download, publishing, pausing the history and
  opening a link you have not seen in that tab are asked first, or whatever
  you chose instead. "Allow for this
  connection" holds for this app's session until the page closes.
- **Asked in the viewer, and in the app.** A question appears as the
  viewer's confirm card. When the app can ask too (an MCP elicitation, as
  Claude Code does), it asks the same question there; answer in either
  place. The model never answers it: the app shows it to you.
- **A refusal is one sentence** the agent reads ("The person's settings in
  Cellucid do not let agents delete what was saved…"). The agent's tool list
  is the same whatever your settings and never changes: the viewer tells the
  agent your settings after its first tool result, and again after a result
  once you changed them. (A changed tool list would make the app rewrite its
  cached prompt, costing several times as much and dropping the agent's
  earlier reasoning.)
- **What agents read goes to their provider.** The first time an agent acts
  on a dataset whose notice you have not answered, the viewer says so once,
  with the dataset's switch "Agents see names and summaries only"
  ({doc}`../../web_app/g_agents/03_permissions_and_what_agents_see`).
- **Results are marked as data.** The data in each result the viewer
  answers (names, labels, notes, a dataset's own text) comes between
  `<cellucid-data>` marks, and the agent is told it is data, never
  instructions.
- **Choices.** When the agent needs you to choose, it asks one short
  question with the choices as a numbered list, and waits for your answer.

## What the agent can ask

- `cellucid_status`: which viewer it drives, every viewer found, and what to
  tell you when none is open.
- Every action of the viewer: read the dataset (fields, categories, genes,
  the history), color, filter, highlight, make annotation fields, run
  analyses (markers, differential expression), move the camera, take a
  picture of the view, undo and redo. The list is the web app's action
  catalog; the reference is {doc}`../../web_app/p_developer_docs/30_action_reference`.

## Privacy

The viewer runs every call in your browser, on your computer, and answers
`cellucid mcp` there; nothing goes through a Cellucid server. What the agent
reads (names, counts, summaries, values, a picture of the view) goes to the
agent app, and from there to the model provider the app uses. A dataset's
"Agents see names and summaries only" switch, set in the viewer, keeps values
of single cells, extremes (a minimum, a maximum, a median) and pictures out of
every agent's results, this one's included, and refuses summaries over groups
an agent made by naming cells: it is enforced by the page, which runs every
agent's call the same way. It is not anonymization: a summary of a very small
group (a narrow range, a filter, a rare category) can still reveal its cells'
values.

## When something is wrong

| What the agent says, or you see | What it means, what to do |
|---|---|
| "No Cellucid viewer is open on this computer…" | Open your data in a viewer (above), then ask again. |
| "No page shows the Cellucid viewer of … now…" | The viewer runs, but no browser page shows it: open the address it names, or show the notebook cell again. |
| "The Cellucid viewer page did not take the call…" | The page is older than the package: reload it. |
| "The Cellucid viewer's page went…" | The page was reloaded or closed during the call; the agent reads the history before trying again. |
| "The agents folder … is not this user's alone…" | Another user owns the folder or others may write in it: `chmod 700` it, or set `CELLUCID_AGENTS_DIR`. |
| "cellucid mcp needs the MCP SDK…" (in the app's log) | `pip install 'cellucid[mcp]'` in the environment the app runs `cellucid` from. |
| The app lists no Cellucid tools | Check the app's MCP log (stderr carries everything `cellucid mcp` says; stdout is the protocol). |

## Related

- {doc}`../g_api_reference_coverage/api/cli` (`cellucid mcp` and its options)
  and {doc}`../g_api_reference_coverage/api/mcp` (`cellucid.mcp.run_mcp`)
- The R side: {doc}`../../r_package/d_viewing_loading/09_local_agents_mcp`
- The web app's agents (the chat bar, Claude and ChatGPT through the
  connector, permissions): {doc}`../../web_app/g_agents/index`
- How it works: {doc}`../../web_app/p_developer_docs/31_local_agent_connection`
