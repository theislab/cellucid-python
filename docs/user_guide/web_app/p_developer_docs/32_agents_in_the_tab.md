# Agents in the tab: chat bar, permissions, data guard and relay

This page documents how an AI agent acts in a Cellucid tab: the one tool
call every way in goes through, the person's permissions and the confirm
card, the chat bar and its providers, where conversations meet the history,
what agents are sent of the data, the prompt-injection guard, and the relay
that lets Claude and ChatGPT reach the tab. The local way in (`cellucid mcp`)
is {doc}`31_local_agent_connection`; the actions an agent calls are
{doc}`29_action_layer_domains_and_history` and {doc}`30_action_reference`.
The user's view of the same machinery is {doc}`../g_agents/index`.

**Audience:** contributors changing `app/agent/`, `relay/`, an action an
agent calls, or what a tool result carries

**Time:** 30–45 minutes

**Prerequisites:** {doc}`29_action_layer_domains_and_history` (the
dispatcher, origins, turns, the history)

The code is the reference: `cellucid/assets/js/app/agent/README.md` and the
READMEs of `agent/bridge/`, `agent/chat/`, `agent/presence/`,
`agent/providers/`, and `cellucid/relay/README.md` state the same contracts
beside the modules.

---

## Three ways in, one tool call

```
 the chat bar             Claude / ChatGPT            Claude Code, Claude Desktop,
 (agent/chat/,            (connector: relay/ →        Cursor (cellucid mcp →
  the person's API key)    agent/bridge/relay-client)  viewer channel, 31_…)
   │ connection 'chat'       │ connection 'connector'    │ connection 'local'
   └───────────────┬─────────┴──────────────┬────────────┘
                   ▼                        ▼
        agent/bridge/tool-call.js  runToolCall({tools, name, args, agent, turn, connection, …})
          1. the dataset's switch: refuse a cell-level read or a shaped summary (data-guard.js)
          2. getActions().dispatch(type, params, {origin: 'agent', agent, connection, turn, …})
               └─ the dispatcher passes every change to the confirmer (agent/presence/)
          3. the switch again (turned on meanwhile?), single-cell statistics and pictures left out
          4. outcomeAsJson → MCP content (mcp-content.js): the JSON between the data tags
```

Every way in runs a call with `runToolCall`: one tool list
(`bridge/mcp-tools.js` `mcpToolsOf`, from the catalogue), one result
encoding, one data guard. An agent acts only through the action dispatcher:
its changes are history steps with `origin: 'agent'`, its name, its
connection (`meta.connection: {id, kind, client, host}`) and its turn, so the
history groups them, Undo turn takes them back, and a permission granted
"for this chat" or "for this connection" keys on the connection, never on a
way in.

Either way in works alone and both work at once: the chat bar needs no
connector, the presence needs no key or chat bar, and a failure of one (a
refused key, a provider down, the relay down, a tab never paired) leaves the
other working and says in one line what to do.

| Folder or module | What it is |
|---|---|
| `agent/index.js` | `initAgents` (from `ui/core/ui-coordinator.js`): the dock at the foot of the plot, the presence, the chat bar and the share choice |
| `agent/instructions.js` | One instruction text, `instructionsFor('connector' \| 'chat' \| 'local')`; a tool it names in backticks must exist (`tests/agent-bridge.test.mjs`) |
| `agent/bridge/` | One tool call (`tool-call.js`), the catalogue as MCP tools (`mcp-tools.js`), outcomes as MCP content (`mcp-content.js`), what a result may carry (`data-guard.js`), the calls from outside the tab (`agent-calls.js`, for the relay and local connections), the relay's WebSocket and wire (`relay-client.js`, `relay-protocol.js`), the pairing (`tab-secret.js`), connect codes and the approval link (`connect-code.js`, `connect-fragment.js`) |
| `agent/presence/` | The dispatcher's confirmer: the person's permissions, the confirm card, the activity line, the notice of what agents send |
| `agent/chat/` | The chat bar: conversation, agent turn, context rule, keys, transcript, export, suggestions, the attached selection, the question with choices |
| `agent/providers/` | One module per model provider behind one contract (`contract.js`): `anthropic.js`, `openai.js`, `scripted.js` (browser specs only) |
| `agent/connect/` | The Agents panel's Claude and ChatGPT part: the connector's status line, address, approval dialog, Manage… |
| `agent/history-link.js` | The one module that knows both a chat and the history |
| `agent/data-scope.js` | The per-dataset switch "Agents see names and summaries only", the notice, and which groups an agent shaped |

## The person's permissions (`agent/presence/`)

The presence hangs on the dispatcher only (its events, `stop`,
`setConfirmer`); nothing in it knows the chat or the relay, so a new way in
gets the permissions, the confirm card and the activity line without code.

**The confirmer contract.** The dispatcher passes every change of an agent
(never a read, never another origin) to the confirmer as `{definition,
params, meta, signal, estimate, confirm, detail, stage, asking}`:

- `stage` `'start'`: before the action runs, with `estimate` (the owner's
  `provideEstimate`: `{ms, bytes, what}` from what this page measured, or
  null) and `confirm` when the definition asks before it runs;
- `stage` `'handler'`: from the implementation of a `confirmAt: 'handler'`
  action (`ctx.confirm(detail)`), once it knows the cost, with `detail` one
  sentence (a download's size);
- `asking()`: called when the person is shown the card; subscribers hear
  `asking`, then `answered` (`allowed`), so the chat says it waits and the
  relay tells a waiting agent;
- the confirmer resolves `true` (the action runs) or anything else (the
  dispatcher refuses with code `declined`: "The person declined: …"). Without
  a confirmer, an action with `confirm` is refused (`needs-confirmation`).

**Kinds and modes** (`presence/permissions.js`, pure; this browser's
`localStorage` `cellucid:agent-permissions`). `PERMISSION_KINDS`: `change`
(anything in the tab: an undoable step, another dataset opened), `heavy`
(estimated at `HEAVY_MS`, ten seconds, or longer on this device and link),
and the five `confirm` kinds (`writes-file`, `deletes-saved`,
`large-download`, `publishes`, `pauses-recording`), each Ask, Allow or
Deny; the strictest of a question's kinds decides (`decide`). The modes are
presets: Ask before heavy or irreversible (the default: changes allowed,
heavy and irreversible asked), Allow everything undoable without asking, and
Read only. Deny refuses with one sentence the agent reads (`deniedSentence`).
No threshold is made up: heavy work is judged from the page's own
measurements against `HEAVY_MS`, the attention limit the gene read's time
target is set from.

**The confirm card.** Who asks (`activity.js` `agentLabel`: "Claude (chat)",
"Claude (Claude Desktop)", "Claude Code (local)"), what (the action's title;
its parameters under Details) and why (each kind asked about). Allow
**Once**, **This chat** / **This connection** (that kind, for that
`meta.connection`, until the page closes) or **Always** (the kind set to
Allow). A card never stays once its question is moot: the agent gave up,
the person pressed Stop on that turn, another dataset is shown, a saved
state was loaded over the history, the tab closed. The first card takes the
keyboard focus unless the person is typing elsewhere.

**The activity line.** One line per agent turn (connection + agent + turn):
while it acts, what it runs and Stop (the dispatcher's `stop` of that turn,
also while one of its questions waits: Stop withdraws the card); once done
with steps, Undo turn and dismiss. Two agents at once are two lines.

**Tool lists.** A connector's or local connection's list says per tool what
the person's settings make of it ("asks", "runs without asking",
"refused"), asked anew for each list (`mcp-tools.js` `permissionOf`, from
the presence's `permissionFor({kind, connection})`); `onPermissionsChanged`
lets a host announce a changed list. The chat's list says none and never
changes within a chat: a model's reasoning is bound to the exact tools it
was given (Anthropic refuses reasoning sent back with other tools, measured
2026-11-04), so its model is told the permissions, the dataset's switch and
the person's instructions with the person's message after any of them
changed (`chat/index.js` `contextFor`).

## The chat bar (`agent/chat/`)

`initAgentChat({dom, actions, presence})`. One message of the person is one
turn (`turn-<id>`): every step the agent takes in answer carries the turn,
the agent's name and the chat as its connection (`{id: <conversation>,
kind: 'chat'}`). `agent-loop.js` `runAgentTurn` runs model call → tool
calls (one after the other, in the model's order, through `runToolCall`) →
results → …, with the effort, the token budget (checked before the next
model call and before a retry), Stop, waiting out a busy provider, and a
conversation too long.

- **Conversation** (`conversation.js`, pure, append only): `{id, createdAt,
  updatedAt, epoch, summary, messages}`; messages are `user` (`text, turn,
  at, context?, attached?`), `assistant` (`turn, provider, agent, model,
  parts, at, epoch, stop, detail?, usage?, otherCalls?`), `tool` (`turn,
  results: [{id, name, title, content, isError, step}]`) and `note`
  (`stopped`, `error`, `refusal`, `limit`, `summary`; never sent to a model).
  Nothing sent to a model is edited afterwards (reasoning and the prompt
  cache are bound to the prefix).
- **Context rule** (`context.js`): the budget is the model's own window. A
  conversation that would not fit is summarised once, at the start of a
  turn, by the same model; a new epoch starts. Every call counts: a failed
  or stopped try and the summariser's calls are carried by the next message
  (`otherCalls`), so the tokens shown and the budget count them.
- **Keys** (`keys.js`): this tab (`sessionStorage`) unless Remember
  (`localStorage`), `cellucid:agent-key:<provider>`; a key goes only to its
  provider's API origin, never into the DOM after it is entered, the
  history, a session file or the URL. **Preferences** (`preferences.js`,
  `cellucid:agent-chat`): provider, model per provider, effort per model, the
  person's instructions, the token budget, open.
- **Providers** (`providers/`): Anthropic Messages API from the browser and
  OpenAI's Responses API (`store: false`), models listed from each API
  (no model list in the code), failures mapped to the contract's kinds
  (`auth`, `permission`, `quota`, `rate`, `busy`, `network`, `unreadable`,
  `model`, `too-long`, `request`); `rate`, `busy` and `network` are waited
  out, an unreadable answer asked for once more.
- **Above the message field** (`#agent-chat-tray`): a question's choices as
  buttons while it waits, and the selection chip. **Choices**
  (`ask.js`): the chat's own tool `ask_person` (2–5 choices, strict for both
  providers) ends the turn; a click or Enter sends the choice as the
  person's message, keeping a half-typed draft; a connector does not list
  this tool, and its instructions say to ask in plain words with numbered
  choices. **Selection chip** (`selection.js`): "Selection: N cells ×" while
  cells are highlighted attaches a reference to each message
  (`message.attached.selection`: the active page's name, its groups and
  counts, never the cells, the names marked as data); × leaves it out until
  the selection changes. The agent reads the cells through the tools, as the
  dataset's switch allows.
- **Sent line** (`sentOf`): each turn's footer says what left the browser;
  with the switch on, that values and pictures were kept out, from the marks
  the tool calls set (not from the switch at send time).
- **Export** (`export.js`): the chat as Markdown, raw HTML escaped and script
  links left as text (`markdownSafe`); the data marks are stripped for people
  to read.

## Conversations and the history (`agent/history-link.js`)

The chat and the history know each other only here. After every message the
chat keeps its conversation with the tracker shown (`history.setChat(tracker,
id, transcript)`; `conversation.js` `transcriptOf`: no providers' reasoning,
each picture replaced by a line saying it was there). Another tracker,
dataset or loaded session shows its latest conversation (or a new one). A
turn running when the dataset changes stops and is kept where it began,
unless its own agent opened that dataset (an action of effect `dataset`):
then it goes on there and the conversation is kept with both histories. A
conversation kept mid-turn (the page reloaded or closed while a call ran)
holds calls with no result; continued, each is answered with what became of
it (`closeInterrupted`: its step, when the history kept one, else nothing
applied), so a model accepts it again.

- `turnSteps`, `stepsOfTurn`: a turn's steps from the history's turn records.
- `revealStep` (a message to its step; the State panel registers how) and
  `revealTurn` (a step to its message; the chat registers how).
- `noteTurn`: the chat's one-line account of a turn (`chat/export.js`
  `turnNoteOf`) kept as the turn's note through `history.note.add`
  (`from_reply`) when the agent wrote none.
- `undoAgentTurn`: the one Undo turn of the chat, the activity line and the
  State panel, through `history.turn.undo`: back to before the turn when it
  is alone at the end of the branch; when someone else acted after it began,
  the person chooses only its own steps (to Recently deleted, dependants
  applied again) or everything since.

A session file holds the conversations only when the share or save includes
them ({doc}`10_sessions_persistence_and_serialization`).

(agent-data-guard)=
## What agents are sent of the data (`data-scope.js`, `bridge/data-guard.js`)

**Free reading is the default**: names, counts, summaries and, when an agent
asks, values of cells and pictures of the plot go to its provider with the
conversation; the dataset's files never leave the browser. The presence says
so once per dataset in a compact notice (when the chat opens or sends, or
when a connector's or local agent first acts), with the switch and Got it.

**The switch** "Agents see names and summaries only" is the person's, never
an action: no agent can turn it off. It is kept per dataset identity in
`localStorage` `cellucid:agent-data` (`createAgentDataStore`; each write
merges into what is stored and other tabs follow it through `storage`
events); a dataset with no identity, or a browser that keeps no site data,
keeps it for the page only. The permission settings show it for the dataset
shown. With it on, `runToolCall` (so the chat, the connector and local
agents alike):

- refuses, before it runs, a read `data-guard.js` classifies as `cells`
  (`CELL_LEVEL_READS`: `cells.values`, and `view.picture`, since the plot
  draws each cell) or does not classify at all, with one sentence;
- refuses a value summary over groups an agent shaped by naming cells
  (`SHAPED_INPUTS`; `data-scope.js` `createShapeTracker` follows, from the
  dispatcher's events, groups an agent made with `highlight.cells`, columns it
  pushed, and fields and pages made from them, by id and every name they
  had);
- checks again once the action ran (the person may have turned the switch on
  meanwhile);
- sends what runs without pictures and without single-cell statistics
  (`SINGLE_CELL_STATISTICS`: a minimum, a maximum, a median are each one
  cell's own value; means, counts, sums and shares stay), each said after the
  data;
- withholds, when the chat sends the conversation again, earlier results that
  could hold cell values (`earlierResultContent`: cell-level, unclassified, or
  summaries with extremes or shaped groups) and their pictures.

**Every read is classified.** `resultDataOf(definition)` answers `cells`,
`summaries` (`SUMMARY_READS`, an analysis's settings, every move or edit of
the history, a change that declares no result) or null. An action nobody
classified is refused while the switch is on, so an action added later is
guarded before anyone classifies it, and `tests/agent-data-guard.test.mjs`
fails until it is classified. Adding a read means adding it to
`CELL_LEVEL_READS` or `SUMMARY_READS` (and to `SHAPED_INPUTS` when it
summarises values over groups).

**Not anonymisation.** A summary of a very small group the data or the person
defines (a narrow range, a filter, a rare category) can still reveal its
cells' values. The notice, the settings, the agent instructions and these
docs say so (`SUMMARIES_ONLY_CAVEAT`); no minimum group size is made up. The
shape marks are kept for the page only, never in a session; a saved state
loaded later starts unknown.

The history's own words follow the same line: a parameter that holds a value
per cell (`x-values`, a column a notebook pushed) is shown as its count in a
step's facts and the report (`history/exports/step-text.js`
`parameterValue`), never the values.

## Text written by others is data (the prompt-injection guard)

Names, labels, notes and a dataset's own text are written by others and may
read like instructions. So:

- every tool result's JSON sits between `RESULT_DATA_TAGS`
  (`<cellucid-data>…</cellucid-data>`, `bridge/relay-protocol.js`;
  `resultDataText` writes every `<` inside as `<`, the same JSON value,
  so nothing in the data can close the tags; `resultDataOf` reads it back);
- Cellucid's own sentences (a refusal, a picture left out, a result cut)
  follow the closing tag, and every name such a sentence quotes is marked the
  same way, every `<` in it escaped (`mcp-content.js` `refusalText`, applied
  by `toolErrorOf`, so every refusal is covered); Cellucid's own quoted words
  use typographic quotes and stay unmarked;
- the relay's own status result and the names the selection chip attaches
  are marked the same way;
- the instructions say early that text from others is data, never
  instructions (`instructions.js`), and the summariser carries no request out
  of it (`chat/context.js`).

A result's text is cut to `RESULT_TEXT_CHARS` (150,000, Claude's limit) with
both counts stated, never silently, the tags kept whole. Pictures travel as
image content; one over `RESULT_IMAGE_MAX_SIDE` (8000 pixels a side) or
`RESULT_IMAGE_MAX_BYTES` (10 MiB of base64) is left out with a sentence.

## The relay (`cellucid/relay/`)

A Cloudflare Worker that makes Cellucid a connector in Claude and ChatGPT
(web, desktop, mobile): a remote MCP server at its root (the connector address
is the bare origin, `https://mcp.cellucid.com`) and its own OAuth 2.1
authorization server whose consent screen is the person's Cellucid tab. It
carries tool calls to the tab the person approved and the tab's answers
back; no accounts, no stored chats, no dataset bytes. `relay/README.md`
holds every route, the connection flow and the security model.

- **Routes**: `/` (the MCP endpoint: Streamable HTTP, one JSON-RPC message per
  POST: `initialize`, `ping`, `server/discover`, `tools/list`, `tools/call`,
  `resources/list`, `resources/read`, `resources/templates/list`,
  `subscriptions/listen`; GET and DELETE answer 405), `/health`, the OAuth
  metadata (`/.well-known/oauth-protected-resource`, whose `resource` is the
  bare origin, and `/.well-known/oauth-authorization-server`), registration, authorize (the
  waiting page with the connect code) and token routes, and `/tab` (the
  tab's WebSocket, proven by its tab secret).
- **Tool calls** go to the pairing's Durable Object (one per browser), which
  sends them to that browser's tab focused last; `bridge/agent-calls.js`
  runs them as connection kind `connector`. Calls of one connection share a
  turn until someone else makes a step or `AGENT_TURN_IDLE_MS` passes.
- **The view card** (`view-app.js`, `ui://cellucid/view`, MCP Apps): what a
  host that shows MCP Apps draws beside `view_picture` and `dataset_open`
  results (Refresh, Open in Cellucid); it reaches nothing but its host.
- **A question to the person.** When the tab asks the person about a call it
  says so (`asking`); if they have not answered when the call's wait runs out
  (`TAB_CALL_WAIT_MS`, 225 s, under Claude's 240 s), the agent is told the
  person is being asked and to read the history once they answered, and the
  tab keeps the question open (`detach`). MCP elicitation is not used by the
  relay.
- **Reconnects.** A pairing that cannot answer (its Durable Object restarting
  under a new deployment) answers a tool call with a sentence saying so and
  to read the history, a tool list with the bundled one; tabs reconnect by
  themselves with growing waits (`relay-client.js`).
- **The bundle** (`relay/generated/bundle.js`): the catalogue, the
  instructions and the wire protocol, generated from the web sources by
  `node scripts/generate-relay-bundle.mjs` (`--check`;
  `tests/relay-bundle-contract.test.mjs` fails until it is regenerated after
  an action, the instructions, the protocol or the version change). A change
  of the instructions reaches connected agents only once the relay is
  redeployed.
- **Deploy** (the owner's): `wrangler.relay.jsonc`, the secrets
  `RELAY_SIGNING_KEY` and `ALLOWED_ORIGINS`; steps in `relay/README.md`.
  `npm run test:relay-bundle` is the Wrangler dry run.
- **Live check** (after a deploy, not in CI): `node
  scripts/relay-live-check.mjs [--relay <url>]` is an MCP client of the
  deployed relay with a Chromium page as the person's tab (discovery,
  registration, PKCE, approval, tool calls, the view card, an idle pairing,
  token rotation and reuse, Forget); it prints PASS or FAIL per check and no
  token.

## Tests

| Test | Covers |
|---|---|
| `tests/agent-bridge.test.mjs` | The MCP tools are the catalogue's; instructions name only existing tools; content encoding, data tags, refusal marking |
| `tests/agent-data-guard.test.mjs` | Every action that answers with data is classified; cell-level reads, shaped summaries and single-cell statistics refused or left out under the switch |
| `tests/agent-presence.test.mjs` | Permissions (`decide`, modes, grants), the confirmer contract, cards made moot, the activity line, Stop while a card waits |
| `tests/agent-chat.test.mjs`, `tests/agent-providers.test.mjs` | Conversations, the context rule, budgets, keys, export; providers on recorded API shapes (no network) |
| `tests/agent-relay-client.test.mjs`, `tests/agent-connect.test.mjs`, `tests/agent-calls.test.mjs` | The relay's tab end, the Connect panel, the shared call engine |
| `tests/relay-oauth.test.mjs`, `tests/relay-mcp.test.mjs`, `tests/relay-bundle-contract.test.mjs` | The Worker on an in-memory Workers runtime; the bundle is current |
| `tests/browser/agent-chat.spec.mjs`, `agent-chat-features.spec.mjs`, `agent-chat-escaping.spec.mjs`, `agent-chat-interrupted.spec.mjs`, `agent-card-moot.spec.mjs`, `agent-connect.spec.mjs` | The chat bar and the connector on a page, with a scripted model and a played relay |
| `tests/browser/agent-data-guard.spec.mjs` | Every read walked with the switch on (three engines); the notice, the selection chip, a question's choices |
| `tests/browser/agent-matrix.spec.mjs` | Either way in alone, both at once, one failing (the scripted model, a played Anthropic API, a played relay) |
| `tests/browser-specs-hold-tests-contract.test.mjs` | Every browser spec file declares tests (a spec cut to its imports once lost thirteen tests while the run stayed green) |
| `tests/browser/live-agents/` (`npm run test:live-agents`) | Real providers, local only, with the developer's keys (the environment's `ANTHROPIC_API_KEY` / `OPENAI_API_KEY`, else `CELLUCID_AGENT_ENV` or the checkout's `.env`); never in CI. `byok.spec.mjs` in every engine, `every-action.spec.mjs`, `showcase.spec.mjs`, `robustness.spec.mjs`, `undo-own.spec.mjs`, and `injection.spec.mjs` (a category named like an instruction, read as data by a real model) |

Next: {doc}`31_local_agent_connection` (the local way in).
