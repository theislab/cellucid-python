# Local Agents (MCP): Claude Code, Claude Desktop, Cursor

**Audience:** anyone who wants an agent app on their computer to work in
the Cellucid viewer of an R object
**Time:** 5 minutes, once
**What you'll do:** let Claude Code, Claude Desktop or Cursor drive a
viewer `cellucid_view()` opened, or one it opens itself.

`cellucid_mcp()` is a local [MCP](https://modelcontextprotocol.io) server
over stdin and stdout, the R peer of Python's `cellucid mcp`
({doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`, which
explains what you see and decide in the viewer, the permissions and
privacy; all of it holds here). The app starts it; the agent gets the web
app's tools and works under the permissions set in the viewer, as
"Claude Code (local)", one undoable turn at a time.

## Set it up once

```bash
claude mcp add cellucid -- Rscript -e "cellucid::cellucid_mcp()"
```

Claude Desktop (`claude_desktop_config.json`) and Cursor (`.cursor/mcp.json`):

```json
{
  "mcpServers": {
    "cellucid": { "command": "Rscript", "args": ["-e", "cellucid::cellucid_mcp()"] }
  }
}
```

Use the full path of `Rscript` when the app does not see your `PATH`.

## Drive a viewer you opened

```r
library(cellucid)
viewer <- cellucid_view(seu)       # a terminal, RStudio or Positron
```

Every R viewer with a server (`cellucid_view()` outside a notebook,
`cellucid_serve()`) is published to agents on this computer, as Python's
viewers are, in the same folder; `cellucid_mcp()` drives the newest one
whose page is open, R's or Python's, and Python's `cellucid mcp` drives R's.
`cellucid_mcp(viewer = "pbmc3k")` names one by its dataset id, dataset name
or viewer id. Not reachable: a viewer shown in a Jupyter notebook's output
(IRkernel; the R kernel runs no server between cells), and a viewer bound
to one non-loopback interface (`host`), whose server this computer's
loopback does not reach.

## Let the app open a dataset

```bash
Rscript -e 'cellucid::cellucid_mcp(open = "/path/to/object.rds")'
```

`open` serves the file in this R process and opens its viewer in the
browser (`browser = FALSE`: its address is written to stderr); the viewer
stops when the app ends the server.

## Good to know

- What you see and decide in the viewer is the same for every agent, the
  chat bar and the connector's included: your permissions, the confirm card,
  **Stop** and **Undo turn**, and a dataset's "Agents see names and summaries
  only" switch ({doc}`../../web_app/g_agents/03_permissions_and_what_agents_see`).
- The agent's tool list is the same whatever your settings and never
  changes while the server runs: the viewer tells the agent your settings
  after its first tool result, and again after a result once you changed
  them. (A changed list would make the app rewrite its cached prompt, costing
  several times as much and dropping the agent's earlier reasoning.)

- stdout carries only the MCP messages; everything else goes to stderr.
- R has no official MCP SDK; `cellucid_mcp()` implements the stdio protocol
  itself (the `initialize` handshake of protocol versions 2024-11-05 to
  2025-11-25, and the 2026-07-28 per-request protocol, whose apps are asked
  in the viewer's card only). It reads stdin without blocking while R's
  event loop runs, so the viewer `open` serves keeps answering its page.
- See {doc}`../g_api_reference_coverage/index` and `?cellucid_mcp` for
  every argument.
