# Decision record: agents

This page answers: **"Why do agents act in a Cellucid tab the way they do:
one tool call for every way in, the person's permissions in front of every
change, a switch for what agents are sent of the data, a relay that keeps
nothing, and a local server with no relay at all?"** It records the
decisions with their dates, the measurements that justify them, what was
tried and taken out, and what is still open. How the code works is
{doc}`32_agents_in_the_tab` (the chat bar, the permissions, the data guard,
the relay) and {doc}`31_local_agent_connection` (`cellucid mcp`), and this
page does not repeat them. Update this page in the same change as the code
it describes.

## At a glance

**Audience**
- Developers changing `app/agent/`, `app/viewer-host/`, `relay/`, the
  packages' local MCP servers or what a tool result carries.

**Time**
- 10 minutes today; the tables grow as history moves here from the code.

**Prerequisites**
- {doc}`32_agents_in_the_tab` and {doc}`31_local_agent_connection`.

:::{note}
This record was started on 2026-11-12 from the two reference pages and the
code. Its tables hold what was written down by then; a story that still
stands in a comment or on a reference page moves here in the change that
rewrites that passage.
:::

---

## What it does now

Paths are under `cellucid/assets/js/` unless they name another place.

1. **Three ways in, one tool call.** The chat bar (the person's own API key),
   Claude and ChatGPT through the relay (a connector), and an agent app on
   the person's computer through `cellucid mcp` or `cellucid_mcp()` all run a
   call with `runToolCall` (`app/agent/bridge/tool-call.js`): one tool list
   made from the action catalog, one result encoding, one data guard.
2. **An agent acts only through the action dispatcher.** Its changes are
   history steps with origin `agent`, its name, its connection and its turn,
   so the history groups them and Undo turn takes them back.
3. **The person's permissions stand in front of every change of an agent**
   (`app/agent/presence/`): allowed, refused, or asked on a confirm card
   while the agent waits. The presence hangs on the dispatcher only, so a new
   way in gets the permissions without code.
4. **The tool list never changes with the person's settings.** The settings
   are told to the agent after tool results instead
   (`app/agent/settings-told.js`).
5. **Free reading is the default, and a switch narrows it.** "Agents see
   names and summaries only" is the person's, per dataset, never an action;
   with it on, cell-level reads and shaped summaries are refused and
   single-cell statistics and pictures are left out
   (`app/agent/data-scope.js`, `app/agent/bridge/data-guard.js`).
6. **Text written by others is data.** Every tool result's JSON sits between
   data tags that nothing inside can close, and the instructions say that
   such text is never an instruction.
7. **The relay carries calls and keeps nothing** (`cellucid/relay/`, a
   Cloudflare Worker): a remote MCP server and its own OAuth 2.1
   authorization server whose consent screen is the person's tab; no
   accounts, no stored chats, no dataset bytes.
8. **A local agent needs no relay.** The packages' MCP server speaks stdio to
   the agent app and the viewer channel to the page a package serves; the
   page runs its calls with the same engine as a connector's
   (`app/agent/bridge/agent-calls.js`).
9. **A conversation is kept with the history's tracker**
   (`app/agent/history-link.js`), and a session file holds it only when the
   share or save includes it.

## Decisions and why

Dates written "by 2026-11-12" are the day the entry was first written into
these docs, where the day of the decision itself was not recorded.

### The ways in, the tool list and the settings

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-04 (owner) | There is no full interactive viewer inside a chat. The in-chat view is the view card (`cellucid/relay/view-app.js`): a picture of the person's view, Refresh, and Open in Cellucid | The sandboxes a chat host gives an MCP App block workers and WebAssembly, which the viewer needs. To be looked at again if the hosts allow them | Running the app inside the chat |
| by 2026-11-06 | A cancel of a local agent's call is sent by whichever ends second, the send of the call or the app giving up (`cellucid-python/src/cellucid/mcp.py`) | The send runs in a worker thread, and a cancel that arrived while it was still there was lost | Cancelling only a call already sent |
| by 2026-11-12 | The instructions a local agent is given name the notebook and R calls that bring its work back (`viewer.request_session()`, `viewer.apply_session_to_anndata(adata, request)`, `plan.apply()`; `cellucid_apply_session(object, cellucid_request_session(viewer))`, `cellucid_apply_plan(plan)`), and a test in each package holds every call named to one that exists | A local agent, driven live, told the person to run a function that does not exist (`cellucid.pull()` in Python) | Instructions written apart from the API |
| 2026-11-04 | The tool list is the same for every way in, all conversation long, whatever the page does (before a page opens, after one took the channel, reloaded or was replaced) and whatever the person's settings are, and `notifications/tools/list_changed` is never sent. Permissions, the data switch and the History settings are told to the agent apart, after a connection's first result and after any result once they changed, never in the tool descriptions | A model's reasoning is bound to the exact tools it was given (Anthropic refuses it sent back with other tools, measured 2026-11-04), and a changed list makes the agent app rewrite its cached prompt and drop its earlier reasoning: a permission changed mid-session cost a connector's next turn about four times its cached price (Measurements) | A tool list that says or follows the person's permissions or settings |
| 2026-11-11 | Only the models Anthropic's Thinking page lists as writing notes between tool calls are asked for progress notes (`display: "updates"`), and a model that refuses the setting is not asked again in the tab | The Models API does not say which models write notes, and the setting comes with adaptive thinking, which would turn thinking on for a model that runs without it by default (Opus 4.8, 4.7, 4.6, Sonnet 4.6) | Asking every model |
| 2026-11-11 | The chat's instructions (not a connector's or a notebook's) carry Anthropic's suggested line for a person in the loop: say what you are about to do, brief updates, a closing recap | The notes alone were rarely written; the line made turns open with words at no cost (Measurements) | Leaving the chat's instructions silent about updates |
| by 2026-11-12 | The Expression values setting is the person's: an agent never changes it unasked and every result names it. A ratio of means of log-scale values is never called an "X-fold" change (the instructions say "log fold change" with its base, or compute the fold change on the original scale first) | An agent that turned the setting on unasked left later chats with it on (seen live, 2026-11); an agent wrote "the means differ about 16-fold on this log scale" | An agent that changes a setting of the person's without being asked |
| by 2026-11-12 | An action's example never says "Ask the person first.", and its example names two different categories to merge | A chat model asked in its own words, so its turn ended on the question and the action was never called; "Merge two categories" was exemplified as merging "T cells" into "T cells", and a model given it merged a category into itself and was refused | Examples that tell the model to ask; an example the action refuses |
| by 2026-11-12 | An agent no site vouches for is named by its host (or "an app on this computer") with its own name quoted, never passing for Claude | Two such apps could not be told apart in the activity line, the confirm card and the history | Showing the claimed name alone |
| 2026-11-11 | A call, tool list or status waits up to `TAB_RECONNECT_WAIT_MS` (20 s, `relay/pairing.js`) for a tab to come back when the browser's last tab closed a moment ago | A reloaded page reaches the relay again within a few seconds (Measurements), so a call arriving during a reload runs in the reloaded tab instead of being told no tab is there; a tab that closed for good is said so after the wait, well within the time a client shows a call as running | Answering "no tab" at once after the last tab closes |
| by 2026-11-12 | Each way of opening a viewer names itself in `runtime` (`jupyter`, `python`, `cellucid-serve`, `cellucid-mcp`) | Ending `cellucid mcp --open` froze the page with "The notebook or R session that opened this viewer stopped" | A page that cannot tell which host opened it |

### The chat beside the plot

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-11 | The chat stands beside the plot as a full-height side panel when the window is wide enough (a plot of at least `PLOT_MIN_WIDTH`, 480 px, stays beside it and the sidebar), and only as a bar over the plot otherwise (`app/agent/dock.js`) | The bar over the bottom of the plot hid the cells the agent was changing | A chat that always overlays the plot |
| 2026-11-11 | While a camera path plays, the open chat folds to its button when it is the bar over the plot, and comes back when the path stops unless the person is writing or opened it again; the side panel stays | On a narrow window the bar covered the lower third of every frame of the film; the side panel covers no cell, and folding it would widen the plot mid-shot | Folding the side panel; playing the film under the bar |
| by 2026-11-12 | The empty chat's suggestions name the clustering and the identity field by the ranked order of field roles (`utils/field-roles.js`, `clusterField`, `identityField`), never the first categorical field | On Suo it offered "Annotate the clusters of sample_ID", its first categorical field; on Pancreas "Annotate the clusters of cell_type", a field that already names the cell types, beside `clusters` (2026-11-10 and 2026-11-11) | Taking the first categorical field |
| by 2026-11-12 | The chat's selection chip offers, and a message carries, only a selection the conversation does not already know (the one attached to its last message, or highlighted when its last turn ended) | An unchanged selection was sent again with every message, and an agent's own highlight came back to it as "the person attached their selection" | Attaching the current selection to every message |
| by 2026-11-12 | What a model writes and the names an agent gives (a category, a page, a made-up tool) are never run as HTML: in the chat, tool rows, history and the Markdown export, raw HTML is escaped and script links are left as text. The Markdown export parses the structure itself and writes it back in a form CommonMark reads one way | Guessing where CommonMark puts code (a backtick run to any backtick run, fences as a toggle) misguessed into live HTML, and link reference definitions and images were never looked at. Checked with markdown-it 3.0 (CommonMark, HTML on, link validation off) on 4,000 generated texts: no tag but p, code, pre, em, strong, br, blockquote, lists, headings, hr and http(s) links (2026-11-04) | Guessing code spans and fences with a regex |
| by 2026-11-12 | The exported chat's rows say what happened: `analysis_export` and `history_note_add` are not "read", a person's Stop of a playing path is not "refused", a color already shown is not "read", and a long turn folds | Seen live; a research answer was buried under 30 to 90 rows that all said the action's title | Rows worded by the tool's kind alone |

### What the agent is told of the history

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The next message to an agent tells it, in one sentence, which of its steps are no longer in effect and which steps came in since its last turn (its own turn undone, a step removed; the person's, another agent's, a notebook's) (`historySinceSentence`, `app/agent/history-link.js`) | After "Undo turn" a model that was not told still said "the Ductal cells are still hidden"; pinned live (`live-agents/undo-own.spec.mjs`) | An agent that keeps reporting as standing what the person took back |
| by 2026-11-12 | `history_undo` takes back whatever step is shown, so the agent is told to call it only when the step shown is one its own call made, and otherwise to point the person to Undo turn | Haiku undid the person's own step (2026-11-04); pinned live: the person's hiding of a category is still there after the turn | An agent that undoes steps that are not its own |
| by 2026-11-12 | A step someone else took while the agent's turn ran is left out of the steps the agent is recorded as having seen, and listed under `during` | A step taken mid-turn (the person coloring) was in the snapshot, so the next message never told the agent it was new | Snapshotting every step in effect at the turn's end |
| by 2026-11-12 | When the history shown is another tracker than the one an agent last saw (a session loaded beside, a saved state, another dataset), the next message says which one, where it came from and its latest steps | A "same questions please" in a different history was read as a first look | Telling only the diff on a tracker the agent never saw |
| by 2026-11-12 | The steps after a turn began that belong to the same agent's other turns of the same conversation are counted as its own later turns | They were worded "9 from another Claude turn", which read as someone else's | Counting every other turn as "another" |
| by 2026-11-12 | A message of a turn the conversation received (a session file or saved state opened) is told to the model as another person's, answered by their agent in their browser | A model continuing a colleague's conversation took the work for its own and the newcomer for the sender ("all 37 earlier steps were taken by me at your request") | Sending the received turns without a word on whose they are |
| by 2026-11-12 | A save or a share keeps a turn taken during it as it was when saved and says so (`transcriptForSave`); the call that makes the save (`history_share`, `history_save`) is answered "the save it made", not "still running". A conversation a save holds is received where it is opened: the sender's failure notes do not offer Try again | An agent's `history_share` wrote a session file whose reader saw, on that very call, "The page closed before this call finished; nothing was applied" and said the export did not finish; the sender's Try again sent the sender's message to the opener's agent and acted again on what the history held (2026-11-11) | Answering every unfinished call of a saved turn as "still running" |
| by 2026-11-12 | What the history would refuse of an agent is refused before the person is asked: the estimate of a squash or a branch delete runs the same checks as its implementation | An Opus turn asked to squash the person's Ins1 and Gcg colorings got the card "It deletes something saved in this browser", the person allowed it, and the squash was then refused: "An agent squashes only its own turn's steps" | A card for a change that will be refused |
| by 2026-11-12 | An agent's or a notebook's turn can start another tracker: the opening is no step of a turn (it says who started it) and the turn's next step in the new tracker makes the turn's record there | The opening was recorded with the turn of the meta and the tracker held no record of that turn, so the model refused the whole tracker ("step ... belongs to the turn ..., which it does not hold") and `history_tracker_start` failed for every agent | Recording the opening as a step of the turn |
| by 2026-11-12 | The chat does not take a change of the tracker shown for the person's change and stop the turn; `history_saved_open` is Load's Open and its load is the turn's own, so an agent asked to continue saved work can open it | The chat stopped the turn and continued the kept conversation, answering the running call "The page closed before this call finished; nothing was applied"; an Opus turn asked to verify the Pdx1 path on Pancreas ended there, mid-sentence | Treating the agent's own load as the person's |

### Turns, budgets and kept conversations

| Date | Decision | Why | Rules out |
|---|---|---|---|
| by 2026-11-12 | The chat's token budget is a cap: the call that would cross it is not made (checked against the spent tokens plus the next prompt's estimate) | Checked only against what was spent, a 150,000-token budget ended at 209,751, since each call reads the whole prompt again | Checking the budget against spent tokens alone |
| by 2026-11-12 | The results of a turn's calls are kept in the conversation as each one comes, and nothing goes to a model before every call is answered | A reload mid-way must keep what finished; a read that ran was answered "the page closed before this call finished" | Saving a turn's results only at its end |
| 2026-11-04 | A chat kept mid-turn (a reload or closed tab while a call waited for the person) answers each call that has no result with what became of it, and its row says interrupted | A provider refuses a conversation in which a call has no result | Keeping a call with no result |
| 2026-11-04 | A question an agent asked is moot once the state it was asked about is replaced: a saved session of the same dataset loaded while a connector agent's card waits refuses the card (the history stays the same object, so its `loaded` event decides) | An answer to a card about a state that is gone would act on another state | A card that outlives its state |
| by 2026-11-12 | A model that ends a turn after its reasoning alone gets a note for the person, never silence | A model may end with no words | A turn that ends with nothing shown |
| by 2026-11-12 | A conversation kept with a tracker that this version cannot show is said to the person, not dropped | Conversations kept with a tracker were lost silently | Dropping an unreadable saved conversation without a word |

### What an agent can read, open and do

| Date | Decision | Why | Rules out |
|---|---|---|---|
| 2026-11-04 (owner) | The dataset switch "names and summaries only" keeps cell-level values out of tool results, not files the person receives: an export runs with the switch on, and the agent is told its name, size and counts (`app/agent/bridge/data-guard.js`) | A file written to the person's computer never goes to the provider | Refusing exports while the switch is on |
| by 2026-11-12 | The addresses in what the person typed (not a choice the agent offered) count as links the person has seen, which the agent may open without asking; an agent opening a link the person has not seen in this tab asks on the card before anything of it is requested (declined: its host hears nothing; from My datasets by its listed name: no question) | The person who types an address has seen it | Asking about a link the person typed; fetching before asking |
| by 2026-11-12 | `cells.open` whose reading fails fails with an error already said (no second notice) and no step | The action had resolved as if the cells were open, so the history recorded "Opened 50% of the cells" and an agent or a notebook was told it was done | An action that resolves on a failed read |
| by 2026-11-12 | `dataset.summary` names the modality shown and every modality, and lists every embedding the data mapping offers, on or off, with what `mapping_set` and `view_embedding_set` take (`{key, columns}`) | An agent opening a CITE-seq `.h5mu` said "I have no tool that lists the modalities or their features"; an agent asked to show a 1D embedding read embeddings named by label only ("UMAP", 2), guessed obsm keys (`X_umap`, `UMAP`) and was refused (2026-11) | A summary that names embeddings by label alone |
| by 2026-11-12 | `field.create-from-pages` takes `rest`: the cells on no page keep their category of a categorical field, so splitting one label is a page per new label and one call | An agent had made 45 temporary pages, a field of 48 pages, then deleted 56 pages: 101 calls, $3 | Making a page for every cell to split one label |
| by 2026-11-12 | `field.summary` takes `within`: a summary of a subset of the cells with the view left as it is | Agents hid the person's cells with filters to summarize donors within one cell type (an Opus turn, 6 filter steps) or cell types per age range (a Haiku turn, 4), then put the view back | Filters as a way to read a subset |
| by 2026-11-12 | `view.picture` takes `source: "export"` (the figure exactly as Export would download it, scaled into the box asked for) and `source: "screen"` (the pane as the screen shows it, overlay included); a screen picture says what it draws over the coloring; a figure picture taken while an overlay is on says so, and how to see the particles | An agent preparing a figure could not see what it was about to export: it set a white background, no grid, a frame and 24 px text and told the person three times "my previews don't show the export settings, please check the file", and never saw a panel title cut off (2026-11-11); another read "the bright cyan patches are the cells with high NANOS3" where the cyan was the highlight over the gene's viridis; another turned the vector field on, tuned it blind and left a flow blown to white | Previews of single views as an agent's only check |
| by 2026-11-12 | `vector-field.summary` reads per category the speed, alignment, direction and where its arrows lead | Asked "where are the cells heading ... evidence from the data", an agent answered "the tools don't give me the velocity numbers" | A flow only a picture shows |
| by 2026-11-12 | Under reduced motion `camera-path.play` answers `ended: 'reduced-motion'` with the seconds moved (0), and `views.state` names the preference (the entry is in {doc}`33_decisions_actions_domains_and_history`) | The agent told the person a film had run in full that had jumped to its last keyframe | An answer that says a film played when it did not |
| by 2026-11-12 | The live every-action sweep names in its task what an action needs in place first (a highlight, a filter, a note to act on) | Haiku asked which, or reported there was nothing to act on | A sweep that lets a model fail for want of set-up |


## Measurements that justify it

Measured on Anthropic's API with the chat's own request (live-agents fixture;
Opus 5.5, Sonnet 5.5 and Claude Haiku 5.5 as named), except the budget and
the reload rows. How the chat sends a
call: {doc}`32_agents_in_the_tab`.

| What | Where and when | Measured |
|---|---|---|
| Reasoning sent back with other tools | Anthropic's API, 2026-11-04 | Refused: a model's reasoning is bound to the exact tools it was given |
| A connector's next turn after a permission changed the tool list | Opus 5.5 | $0.41 against about $0.10 cached: a changed list rewrites the client's cached prompt |
| What a call carries | `count_tokens` on the exact request the chat builds, 2026-11-11 | About 62,000 prompt tokens before the conversation: 96% the tool list (178 tools, about 59,600 tokens), 2,300 the instructions; tool results (a few hundred characters) and pictures (640 pixels) are small beside it |
| Saying a sentence true of most tools once in the instructions | The same request | A description that says only its exception saved about 1,800 tokens a call |
| `examples`, `additionalProperties: false`, `minLength` in the tool schemas | The same request | Each at most 4% of the prompt, and each is what lets a model call exactly: kept |
| Sending the 128 rarely used tools with `defer_loading` beside Anthropic's tool search | Claude Haiku 5.5, 12 tasks | A call falls to about 21,000 tokens and 0.4 s faster, but 4 of the 12 tasks the full list got right were missed, and a task needing a deferred tool took two or three searches (slower, not cheaper): not adopted |
| A budget checked against the tokens spent alone | A 150,000-token budget in the chat's loop (the chat test holds it), 2026-11 | The turn ended at 209,751 tokens, since each call reads the whole prompt again: the next call's estimated prompt is part of the check |
| The time for a reloaded page to reach the relay again (the page dials before its dataset opens) | The Pancreas sample with a local relay, 2026-11-11 | 1.5 to 2.8 s, which is why a call waits `TAB_RECONNECT_WAIT_MS` (20 s) for a tab |
| Progress notes shown (`display: "updates"`) against not shown | One live 10-call annotation turn, setting off then on and the order reversed (the first run of a prefix pays its cache write), 2026-11-11 | Sonnet 5.5 728k to 725k prompt and 1,333 to 1,187 output tokens; Opus 5.5 803k to 802k and 1,555 to 1,555: the notes are billed in `output_tokens` at full length whatever the display, so showing them changes no cost |
| Whether a model writes a note at all | The same turn and a 4-call one, 2026-11-11 | Neither model wrote a note longer than its plain text remarks, so none showed: the models write them on longer, quieter stretches of work |
| The chat's instructions asking for progress updates | The four showcase tasks, Opus 5.5, two passes each way, the order reversed, 2026-11-11 | The turn opened with a line 8 of 8 times (6 of 8 without); the longest run of model calls with nothing to read fell from 4 and 6 to 2 on the longest task (annotation); words between calls 29 parts against 23; cost unchanged (output 8,329 against 8,321 tokens a pass, $0.99 against $1.02 a pass of four tasks) |

## Tried and reverted

| What | When | Result | Would need |
|---|---|---|---|

## Known limits and open items

- **A package's own status text is not marked as data.** What a package
  answers itself, with no page in the path (`cellucid_status`), carries the
  dataset's name unmarked, where the relay's status result marks it
  ({doc}`31_local_agent_connection`).
- **Live agent tests never run in CI.** They spend real credit and need the
  developer's keys ({ref}`long-suites`).

## Where the code and tests live

| Concern | Module(s) | Tests |
|---|---|---|
| One tool call, the tools, the result encoding | `cellucid/assets/js/app/agent/bridge/` (`tool-call.js`, `mcp-tools.js`, `mcp-content.js`, `agent-calls.js`), `cellucid/assets/js/app/agent/instructions.js` | `cellucid/tests/agent-bridge.test.mjs`, `cellucid/tests/agent-calls.test.mjs`, `cellucid/tests/mcp-tools-document-contract.test.mjs` |
| Permissions, the confirm card, the activity line | `cellucid/assets/js/app/agent/presence/` | `cellucid/tests/agent-presence.test.mjs`, `cellucid/tests/browser/agent-card-moot.spec.mjs` |
| What agents are sent of the data | `cellucid/assets/js/app/agent/data-scope.js`, `cellucid/assets/js/app/agent/bridge/data-guard.js` | `cellucid/tests/agent-data-guard.test.mjs`, `cellucid/tests/browser/agent-data-guard.spec.mjs` |
| The chat bar and the providers | `cellucid/assets/js/app/agent/chat/`, `cellucid/assets/js/app/agent/providers/` | `cellucid/tests/agent-chat.test.mjs`, `cellucid/tests/agent-providers.test.mjs`, `cellucid/tests/browser/agent-chat.spec.mjs`, `cellucid/tests/browser/agent-matrix.spec.mjs` |
| Conversations and the history | `cellucid/assets/js/app/agent/history-link.js` | `cellucid/tests/browser/agent-chat-history-trackers.spec.mjs`, `cellucid/tests/browser/agent-chat-interrupted.spec.mjs` |
| The connector in the tab | `cellucid/assets/js/app/agent/connect/`, `cellucid/assets/js/app/agent/bridge/relay-client.js` | `cellucid/tests/agent-relay-client.test.mjs`, `cellucid/tests/agent-connect.test.mjs`, `cellucid/tests/browser/agent-connect.spec.mjs` |
| The relay | `cellucid/relay/` | `cellucid/tests/relay-oauth.test.mjs`, `cellucid/tests/relay-mcp.test.mjs`, `cellucid/tests/relay-bundle-contract.test.mjs` |
| A Python or R host's commands | `cellucid/assets/js/app/viewer-host/` | `cellucid/tests/viewer-protocol-contract.test.mjs` |
| Real providers | `cellucid/tests/browser/live-agents/` | `npm run test:live-agents` (never in CI) |

## Related pages

- {doc}`32_agents_in_the_tab` and {doc}`31_local_agent_connection`: the
  references.
- {doc}`30_action_reference`: every tool an agent has.
- {doc}`33_decisions_actions_domains_and_history`: turns, stops and what an
  agent's change is in the history.
- {doc}`26_decisions_served_store_and_widgets`: the viewer channel a local
  agent's calls travel on.
- {doc}`../g_agents/index`: what a user sees.
