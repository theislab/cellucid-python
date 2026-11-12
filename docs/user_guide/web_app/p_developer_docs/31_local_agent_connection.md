# Local agent connection (`cellucid mcp`)

An agent app on the person's computer (Claude Code, Claude Desktop, Cursor)
drives a Cellucid viewer through `cellucid mcp` (cellucid-python) or
`cellucid_mcp()` (cellucid-r): a local MCP server over stdio. There is no
relay and no deployed site in the path; the page that runs the calls is the
viewer page a package serves. User documentation:
{doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`,
{doc}`../../r_package/d_viewing_loading/09_local_agents_mcp`.

## One engine for agents from outside the tab

The page runs a local agent's call exactly as it runs a connector agent's:
`app/agent/bridge/agent-calls.js` (`createAgentCalls`) is the one engine
behind both the relay client (`relay-client.js`, connection kind
`connector`) and the viewer channel (connection kind `local`). It

- dispatches the call's action through `runToolCall` as the agent
  (`origin: 'agent'`, the app's name, a turn, the connection
  `{id, kind: 'local', client}`), so the person's permissions, the confirm
  card, the activity line ("Claude Code (local)", `presence/activity.js`
  `agentLabel`), Stop and Undo turn are the connector's;
- groups a connection's calls into one turn until someone else makes a step
  or `AGENT_TURN_IDLE_MS` passes;
- answers the tool list (`mcp-tools.js` `mcpToolsOf`: the same for every
  connection, whatever the person's settings) and tells the agent those
  settings after its connection's first tool result and after any result
  once they changed, as a note appended to the result's text
  (`agent-calls.js` `withSettingsTold`, words in `settings-told.js`);
- tells the caller when the person is asked about a call, with the card's
  own sentence (`presence.questionFor`), and when they answered.

A local call is therefore guarded like any other
({ref}`agent-data-guard`): the dataset's switch "Agents see names and
summaries only" applies to it, a refusal is one sentence, and the page's
result carries its data between the data tags (`RESULT_DATA_TAGS`). What a
package answers itself, with no page in the path (`cellucid_status`, and the
sentences of `local_agents.json`), is the package's own text: the status
result's dataset name is not yet marked as data the way the relay's status
is.

The presence adds two things a local app needs: `answer(signal, allowed)`
(the person answered the same question in the agent's app, an MCP
elicitation: the card goes as if clicked, Allow once or Deny) and
`onPermissionsChanged` (the settings, or a grant for one connection,
changed: the host tells the app to list its tools again).

## The tools: one source

The tools, the instructions (`instructionsFor('local')`), the status tool
(`bridge/local-agents.js` `LOCAL_STATUS_TOOL`) and the app names
(`LOCAL_AGENT_NAMES`: `claude-code` is "Claude Code", `claude-ai` "Claude
Desktop", `cursor-vscode` "Cursor"; else the app's `clientInfo.title` or
`name`) are generated into `mcp-tools.json` (`mcpToolsDocument`) by
`scripts/generate-action-catalog.mjs`, in the same run as the action
catalog: `tests/fixtures/mcp-tools.json` here, `src/cellucid/mcp_tools.json`
and `inst/mcp_tools.json` in the packages, each made from that repository's
own `action_catalog.json` (its `catalogDigest`). Its `tools` are
`mcpToolsOf(catalogue)` with the default permissions: exactly the relay's
bundled list (`tests/mcp-tools-document-contract.test.mjs` checks them equal,
field for field). The packages use it only while no page is open; an open
page's own list (the same function, with the person's settings for this
connection) is used otherwise. Nothing is hand-copied: the packages read the
document, and three tests pin its digest.

## The viewer channel's agent messages

`viewer-protocol.json` (schema 5) adds, host to page: `agent_list
{requestId, connection}`, `agent_call {requestId, tool, args, agent,
connection}`, `agent_cancel {requestId}`, `agent_answer {requestId,
allowed}`; page to host: `agent_tools` and `agent_result` (binary: the MCP
tools, or the MCP tool result, as UTF-8 JSON; every call gets exactly one
result, an error included), `agent_asking {requestId, message, detail}`,
`agent_answered {requestId, allowed}`. No event says the tool list
changed: it never does.
`connection` is `{id, client}` (`local-` and 12 hex characters per MCP
session; the app's `clientInfo.name`). The page reads them in
`viewer-host-source.js` (`onAgent`, `sendAgentResult`, …) and `main.js`
wires them to the engine. A host lists them in `/_cellucid/protocol` like
every message, so a page sends agent events only to a host that takes them,
and an older page answers an agent command with a `command_error` that
reaches the agent ("reload the page"). A package delivers agent events to the
viewer's agent port, never to the notebook's hooks.

## Reaching a viewer in another process

A notebook's viewer runs in the kernel; `cellucid mcp` is another process.
Every package viewer with a channel is **published** (`cellucid serve` and
`serve_anndata()` included: their page holds a viewer channel of its own,
whose events reach no notebook; `CELLUCID_AGENTS=off` publishes nothing):

- **a record** per viewer in the agents folder (`$CELLUCID_AGENTS_DIR`, else
  `<user cache>/cellucid/agents`: `%LOCALAPPDATA%` on Windows,
  `~/Library/Caches` on macOS, `$XDG_CACHE_HOME` or `~/.cache` elsewhere;
  both packages compute the same folder), written atomically, readable by
  its user only (schema 2): `schema, viewer, runtime, host, pid, started,
  dataset {id, name}, url, endpoint, token`. `url` is where the viewer is
  served, never the page address (that carries the viewer token, which no
  agent is given). A reader checks every field's type and reads only
  records of its own host whose endpoint is on loopback; a record of this
  host whose process ended is removed, another host's (a home folder shared
  between machines) left alone. On POSIX the folder must be the user's and
  writable by no one else, or it is refused with one sentence.
- **an endpoint** on loopback, `<origin>/_cellucid/agent/<viewer id>`:
  Python runs one small HTTP server per process (a notebook widget has no
  server of its own); R mounts the route on the viewer's httpuv server.
  `GET` describes the viewer, `POST …/command` takes one agent command
  (`202`; `409` when no page holds the channel, with the sentence the agent
  reads), `GET …/events?consumer=<connection id>&wait=<s>` returns that
  consumer's events (`agent_tools`, `agent_result`, `agent_asking`,
  `agent_answered`, `agent_failed`, `viewer_stopped`).
  Every request carries `Authorization: Bearer <token>`; one with an
  `Origin` header (a browser), from another host, or naming another `Host`
  is refused. A consumer silent for `CONSUMER_IDLE_SECONDS` (its `cellucid
  mcp` was killed) is swept on a clock: its calls are stopped in the page,
  their cards withdrawn, and its queued events let go. Anyone on the
  computer can connect without the token, so Python's endpoint closes a
  connection idle for `CONNECTION_IDLE_SECONDS` (120 s; each holds a thread
  of the viewer's process) and answers a command nested deeper than Python
  follows `400`; R's httpuv server holds no thread per connection.

`cellucid mcp` drives the viewer it opened itself (`--open`, in-process:
the viewer's port directly), else the newest published viewer whose page is
open (`--viewer` names one), found again when the one it drives stops, or,
when none was named, when its page is not open and another viewer's is. The
endpoint is plain HTTP, not a WebSocket, so that each package's client is
its standard library (Python's `urllib`, R's `curl`) and either package's
`cellucid mcp` drives either package's viewer; both endpoints hold an
`events` request up to `EVENT_WAIT_LIMIT_SECONDS` for an event (R through a
promise httpuv answers later), and the client asks again at once. The shared fixture
`local_agents.json` (cellucid-python `tests/fixtures/`, cellucid-r
`tests/testthat/fixtures/`, byte for byte) pins the folder rule, the record,
the endpoint, the app-name cases and every sentence an agent reads.

**The page going.** When the page holding a viewer's channel closes, is
reloaded or is replaced by a newer one, every agent call and list it was
running fails at once with the page-gone sentence (whether it was done is
not known; read the history first), never a wait. A page taking the
channel again changes nothing in the tool list.

Not published: a page opened at a server's bare `/?anndata=true` (no
channel) and an IRkernel notebook's viewer (no server between cells).

## The MCP server

Python uses the official SDK (`mcp`, the `cellucid[mcp]` extra; the 2.x line
speaks the `initialize` handshake and the 2026-07-28 per-request protocol);
R has no official SDK and implements the stdio subset itself
(`R/mcp.R`). Both are thin over one core (`cellucid._mcp_core`, R's
equivalent): `initialize` gives the local instructions and
`tools.listChanged: false`; `tools/list` the tools and `cellucid_status`;
`tools/call` sends `agent_call` and returns the page's result as it is.
The local instructions tell the agent to ask the person in plain words with
numbered choices (the chat bar's `ask_person` tool is not listed).
While the person is asked, the server sends progress (when the app gave a
progress token) and, when the app declared `elicitation`, asks the card's
question as a form elicitation with no fields (an app that declared URL
elicitation only is not asked): Accept sends `agent_answer true`, Decline
`false`, dismissing it leaves the card to answer; once the card is answered
the app's question is withdrawn (`notifications/cancelled` for the
elicitation request), and an answer the app sends later changes nothing.
Progress counts 1, 2, … within a call. `notifications/cancelled` sends `agent_cancel`
(the page stops the call: a change rolled back, but a playing camera path
stopped where it is; a question withdrawn).

**The tool list never changes.** It is the same whatever the person's
settings, and neither server sends `notifications/tools/list_changed`: for
Claude 5.5 models a changed list rewrites the app's cached prompt (about
55k tokens) and drops the conversation's earlier reasoning. Measured on
Opus 5.5, the connector's next turn after a permission change cost $0.41
against about $0.10 cached, about 4×. The person's settings reach the agent
after its tool results instead ("[Cellucid: the person's settings for you
now. …]"), and the local instructions say so.

stdout carries only MCP messages: the `--open` viewer starts once the
transport holds stdout (Python's SDK sends what else is printed there to
stderr; R starts it after the handshake, its prints sunk to stderr).

## Tests

- web: `tests/agent-calls.test.mjs` (the engine as a local connection),
  `tests/viewer-host-source-contract.test.mjs` (the agent commands and
  events), `tests/mcp-tools-document-contract.test.mjs` (the generated
  document is the connector's list), `tests/viewer-protocol-contract.test.mjs`,
  and `tests/browser/mcp-local.spec.mjs`: a real Python viewer, `cellucid
  mcp` in another process and this checkout's app in Chromium: a step of
  "Claude Code", the "(local)" line, Undo turn, the list changing with the
  settings, a card asked in the app and declined there, a call the app
  cancels withdrawing its card, "This connection"; and `cellucid mcp --open`
  serving a file and stopping its viewer with it.
- cellucid-python: `tests/test_mcp.py` (pins, the wire, the port, the
  endpoint's refusals, a remote link, the session, the server over stdio
  with a fake page: initialize, lists, calls, errors, cancellation,
  elicitation, list changes, a viewer stopping).
- cellucid-r: `tests/testthat/test-mcp.R`, the same, and Python's
  `cellucid mcp` driving an R viewer.
