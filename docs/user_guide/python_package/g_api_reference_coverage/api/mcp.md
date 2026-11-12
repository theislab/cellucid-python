# Local agents (`cellucid mcp`)

```{eval-rst}
.. currentmodule:: cellucid.mcp
```

`cellucid mcp` is a local MCP server over stdio for an agent app on this
computer (Claude Code, Claude Desktop, Cursor): the app works in a Cellucid
viewer with the web app's tools, under the permissions set in the viewer.
The app starts it; set-up, what you see and decide, and privacy are in
{doc}`../../d_viewing_apis/18_local_agents_mcp`. It needs the `mcp` extra:

```bash
pip install 'cellucid[mcp]'
```

The command (`cellucid mcp`, {doc}`cli`) and {func}`run_mcp` are the same
server; `run_mcp` is what the command calls, for an app that starts Python
itself:

```json
{"mcpServers": {"cellucid": {"command": "python", "args": ["-c", "from cellucid.mcp import run_mcp; run_mcp()"]}}}
```

- **Tools**: every action of the viewer, generated from the web app's action
  catalog exactly as the Cellucid connector lists them (names, input
  schemas, descriptions, annotations), plus `cellucid_status`. While a
  viewer page is open, its own list is used (the same tools). The list is
  the same whatever your settings and never changes while the server runs
  (`tools.listChanged` is false): the viewer tells the agent your settings
  after its tool results.
- **Calls** run in the viewer page as the app's agent ("Claude Code"),
  through a connection of kind `local` ("Claude Code (local)" in the viewer).
  A question the person is asked is the viewer's confirm card and, when the
  app declares the `elicitation` capability, an elicitation in the app.
  A call the app cancels (`notifications/cancelled`) is stopped in the
  viewer.
- **Results** of the viewer's tools are the connector's: the data between
  `<cellucid-data>` and `</cellucid-data>` (every `<` inside escaped),
  Cellucid's own sentences after it, a refusal in one sentence, and a
  dataset's "Agents see names and summaries only" switch applied by the page
  ({doc}`../../../web_app/g_agents/03_permissions_and_what_agents_see`).
- **Output**: stdout carries only MCP messages; everything else (the
  `--open` viewer's address, why no viewer is reachable) is said on stderr.

```{eval-rst}
.. autofunction:: run_mcp
```
