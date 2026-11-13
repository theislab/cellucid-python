# Action layer, state domains and history

This page documents the one way the open dataset's state changes in the web
app: the **action catalog** (every change and read, defined once), the
**dispatcher** (checks, orders, records and rolls back), the **state domains**
(the units state is captured, compared, undone and saved in), and the
**history** (trackers of step trees over a pool of values kept by content
hash), with the autosave, the exports, and how one catalog reaches Python,
R and agents. It says what each part does and the rule it keeps. Why it is
built this way (the decisions with their dates, what was measured, what was
tried and taken out) is {doc}`33_decisions_actions_domains_and_history`.

**Audience:** contributors adding a control, an action or a piece of state;
maintainers of the Python and R viewer APIs and of the session format

**Time:** 45–75 minutes

**Prerequisites:** {doc}`05_app_architecture_overview` and
{doc}`06_state_datastate_and_events`. The session file that holds a history is
{doc}`10_sessions_persistence_and_serialization`; the user's view of the same
machinery is {doc}`../m_state_history/index`.

The code is the reference: `cellucid/assets/js/app/actions/README.md`,
`app/history/README.md`, `app/history/exports/README.md`,
`app/state/domains/README.md` and `app/state-serializer/README.md` state the
same contracts beside the modules.

---

## The picture

```
 a control           an agent            a Python or R host      the history
 (ui/modules/*)      (agent/bridge)      (viewer-host/)          (delete a step)
   │ origin user       │ origin agent       │ origin host            │ origin replay
   └──────────┬────────┴─────────┬──────────┴────────────┬──────────┘
              ▼                  ▼                       ▼
      app/actions/dispatcher.js   dispatch(type, params, meta)
        1. catalog.check: the parameter schema    → TypeError naming the parameter
        2. a handler provided for the type?        → ActionRefused 'unavailable'
        3. effect read / control / history / dataset / outside: run now, no step, no gate
        4. confirm (agents only), the history's gate, supersede,
           wait for running actions on overlapping domains
        5. recorder.before: capture the touched domains
        6. the handler (the module that owns the state)
              failure → recorder.failed: put the before-values back
        7. recorder.after: capture again; one step if anything changed
              ▼
      app/history/history.js    trackers (tree.js) + value pool (value-pool.js)
              ▼                     undo / redo / jump: put values back
      app/state/domains/*.js    capture() · restore(value) · describe(a, b)
              ▼
      DataState, the viewer, the sidebar controls, the panels
```

A control never calls a mutating state method itself: it dispatches the action
it stands for, and the action's implementation calls the state. Settings
controls and camera gestures (and a few layout and analysis gestures, listed
below) change their owner first and are recorded after, by observation.

---

## The action catalog

### One definition per action

Every action is one call of `defineAction` (`app/actions/catalog.js`) in one
pure module per area under `app/actions/definitions/`, listed in
`definitions/index.js`: `history.js`, `fields.js`, `coloring.js`,
`filters.js`, `mapping.js`, `highlights.js`, `views.js`, `camera.js`,
`settings.js`, `layout.js`, `reads.js`, `analysis.js`, `datasets.js`,
`stored.js`, `exports.js`, `annotation.js`. The
definitions import nothing from the page, so the catalog loads in Node and
every generated surface is made by the code the page runs. The generated
{doc}`30_action_reference` lists every action and counts them by effect
(`state`, `view`, `layout`, `read`, `control`, `history`, `dataset`,
`outside`); no page counts them by hand, so no count goes stale when an
action is added or replaced.

`defineAction` checks a definition and freezes it; an unknown field, a schema
outside the subset, an unknown domain or effect, or a missing `describe` or
`code` throws when the module is imported, so every test that loads the
catalog fails on it.

| Field | Required | What it is |
|---|---|---|
| `type` | yes | Dotted lowercase words, `<group>.<verb>` or deeper (`filter.category.visibility`), matching `ACTION_TYPE_PATTERN`. It is the agent's tool name (`toolNameOf`: dots and hyphens become underscores), the Python method name and the R function name (`cellucid_<tool>`), so it is stable: renaming one is a breaking change for every generated surface and changes the catalog's digest. |
| `title`, `summary`, `group` | yes | One trimmed line each, for a person and an agent. `group` is the heading the action is listed under. |
| `params` | yes | The parameter schema, an `object` schema in the subset below. Targets are named, never positional. |
| `touches` | for every effect but `read`, `control`, `history`, `dataset` and `outside` | The state domains the handler can change, distinct ids of `app/actions/domains.js`. A stepless action touches none; every other action touches at least one. |
| `reads` | yes (empty for none, and for an action that records no step) | Its **result inputs**: what its result is computed from beyond its parameters and the domains it touches (`ACTION_READ_INPUTS`: state domains, and `history` for another branch's tip). A step records each one's value hash (`Step.reads`), so a step whose input changes later is stale ({ref}`declared-inputs`). |
| `captures` | no (default none) | The inputs of `reads` its result holds the value of itself (`camera-path.keyframe.add` captures `camera`: a keyframe is the camera it was saved at). Still recorded, so replay and re-apply know the step read them; a later change of one never makes the step stale (`provenance.js`). In the shared catalog document beside `reads`. |
| `context` | no (default none) | Its **name context**: the state its parameters' names resolve in and its result does not follow (a field's name in `field-overlays`, the view in focus in `views`, a gene's name in `mapping`, a step's id in `history`). Never recorded, so a rename or a focus change stales nothing. |
| `effect` | yes | One of `ACTION_EFFECTS` (below). |
| `confirm` | no (default `null`) | Why the person is asked first when an agent acts: `writes-file`, `deletes-saved`, `large-download`, `publishes` or `pauses-recording` (`ACTION_CONFIRMS`, each a permission kind of `app/agent/presence/permissions.js`; `history.recording.pause` is the one that pauses recording). Anything the history can undo needs none. |
| `opensLinks` | no (default `false`) | The action may fetch an address its parameters name (`dataset.open`: a link, an accession, a server). Its implementation asks, before anything is requested, about an address the person has not seen in this tab (`ctx.confirmLink(address)`, the permission kind `unseen-link`; seen: typed or pasted in the address bar or the chat, open now, in My datasets, `ui/modules/seen-addresses.js`), so an agent misled by text it read cannot make the tab send it to any host. In the catalog document; the MCP tool says what the person's setting makes of it and is `openWorldHint`. An action of effect `dataset` likewise asks, through `ctx.confirmReplacing(detail)` (kind `replaces-unkept`), when the dataset shown holds steps this browser is not keeping (`app/unkept-work-guard.js`, the one guard every replacement passes in `app/main.js`; the person's own opening is asked by a dialog there) |
| `confirmAt` | no (default `start`) | When the person is asked: `start` (before the action runs), or `handler` (its implementation asks through `ctx.confirm(detail)` once it knows the cost, as `dataset.open` does for a whole-file download once the link says its size; never when nothing costs) |
| `coalesce(params)` | no | A key: consecutive steps with the same key merge into one (a slider drag, a camera gesture). |
| `supersede` | no (default `false`) | A newer dispatch of the same type cancels an unfinished one (`color.by`, `cells.open`, `view.embedding.set`, `view.embedding.default`). |
| `keepsWhenStopped` | no (default `false`) | Stopped while it runs (the person's Stop of an agent turn), it keeps what it did and its step instead of being rolled back; its implementation hears the signal and ends where it is (`camera-path.play`: the camera stays where the film stopped). Only for an action that records a step. Not in the shared catalog document: it is the page's own. |
| `rerun` | no (default `false`) | Whether the same parameters dispatched again bring a stale step's own result up to date with its inputs: `analysis.run`, and `highlight.category` and `highlight.range`, which give a group already on the page the cells the view shows now, in place ("Updated the cells of …"). Only a stale step of such an action gets **Run again** in the State panel (`ui/modules/state/step-details.js`); an action that adds another of what it made (`view.keep`, `field.duplicate`) or only sets options does not. Needs `reads` (only a step with inputs is stale). Not in the shared catalog document: it is the page's own. |
| `describe(params, outcome)` | yes | The timeline label and the sentence the methods paragraph uses: one line, past tense. |
| `code(params)` | yes | One line of pseudocode (`color_by("cell_type")`). |
| `result` | for every effect that records no step (`read`, `control`, `history`, `dataset`, `outside`) | The schema of what the action answers. |
| `version` | only for an action whose steps are recorded (effect `state`, `view` or `layout`; default `1`), and refused for any other, whose `version` is `null` | The version of its parameters. A step records it (`Step.v`), so a history recorded with older parameters keeps replaying: `upgradeParams(definition, from, params)` brings them to this version first. A change to `params` (beyond its descriptions) bumps it, adds the upgrader and keeps an example of the new version (`tests/fixtures/action-versions/<type>@<version>.json`, written by `scripts/generate-action-versions.mjs`); `tests/action-versions-contract.test.mjs` fails on a params schema that changed without a bump, and on an older example that no longer upgrades to valid parameters. |
| `upgrade` | with `version` above 1 | `{<from>: params => params}`: one pure function for each earlier version, taking its parameters to the next version's. |

| Effect | Meaning | Records a step |
|---|---|---|
| `state` | The data's meaning: fields, labels, filters, highlights, open cells, results | yes |
| `view` | How it is shown: coloring, camera, render settings, views | yes |
| `layout` | Panels, sections, sidebar, theme | yes |
| `read` | Changes nothing; answers its `result` | no |
| `control` | Stops work under way (`camera-path.stop`, `analysis.cancel`) at once, beside the action still waiting for that work; that action's own step records what it changed (a stopped `camera-path.play` answers ended `stopped`, its step holding the camera where the path stopped). Touching the domain of the work it stops, a stop would wait in order behind it and run only once nothing was left to stop | no |
| `history` | Moves or edits the history itself (`history.undo`, `history.redo`, `history.jump`, branches, groups, notes, the journal, a replay, the guided tour): records no step of its own | no |
| `dataset` | Replaces the open dataset (`dataset.open`), or the history shown with a saved state continued (`history.saved.open`, its dataset first when another is shown): records no step; the next dataset's history starts, or the saved one goes on | no |
| `outside` | Acts outside the open dataset's state and its history: writes a file the person receives (`history.share`, `history.report`, `history.anndata.save`, `figure.export`, `cells.export`, `field.export`), publishes (`annotation.publish`, `dataset.share-link`), or keeps or removes what this browser keeps for later visits (`history.save`, `history.saved.delete`, `my-datasets.forget`, `stored.remove`). Undo cannot take it back, so it names its `confirm` unless it only keeps something new in this browser (`history.save`) | no |

An action whose browser part needs the person's own gesture (a window a
browser opens only from a click, a sign-in, a save dialog) never waits for
one when an agent or a host runs it: it falls back to what needs no gesture
(the report as a downloaded `.html` file or as text, the AnnData copy
downloaded once whole), or is refused with one sentence saying what the
person has to do (`annotation.publish` until the person has connected the
repository and signed in to GitHub).

The actions that set `confirm` are `dataset.open` and `history.saved.open` (`large-download`), the
files (`writes-file`: `history.share`, `history.anndata.save`,
`history.report`, asked by its implementation once it writes a file,
`figure.export`, `cells.export`, `field.export`), the publications
(`publishes`: `annotation.publish`, `dataset.share-link`),
`history.recording.pause` (`pauses-recording`) and the deletions
(`deletes-saved`: `history.branch.delete`, `history.clear`,
`history.journal.empty`, `history.saved.delete`, `history.squash`,
`history.step.delete`, `history.tracker.remove`, `my-datasets.forget`,
`stored.remove`). The agent
presence installs the confirmer (`app/agent/presence/index.js`,
`setConfirmer`); with none installed, an agent's dispatch of a `confirm`
action is refused with code `needs-confirmation`. A person's dispatch
(`origin` other than `agent`) is never asked. Every other change of an agent
passes the confirmer too, with its owner's estimate of the work
({doc}`32_agents_in_the_tab`).

### The parameter schema subset

`app/actions/schema.js` accepts a closed subset of JSON Schema (draft
2020-12), because one schema is read five ways: the dispatch check, the agent
tool's input schema, the Python and R generators, a reference listing, and the
sample parameters of the round-trip tests (`exampleOf`). Supported keywords:
`type` (one or a list), `properties`, `required`, `additionalProperties`
(`false` on every object, so a misspelt parameter is refused by name),
`items` (required on every array), `minItems`, `maxItems`, `uniqueItems`,
`enum`, `const`, `minimum`, `maximum`, `exclusiveMinimum`, `minLength`,
`maxLength`, `pattern`, `oneOf`, `description`, `title`, `default`,
`examples`, and two Cellucid annotations:

- `x-cells: true` on an array of original cell positions: a `Uint32Array` is
  accepted for it, so a binary viewer-protocol frame and a JSON tool call
  land on the same check (the agent tools drop the keyword and say it in
  words);
- `x-values: true` on an array of per-cell numbers (a column a notebook
  pushes): any numeric typed array is accepted.

`schemaProblem(schema, value)` returns one sentence naming the parameter
(`color.by.field.key is required`, `color.by has no parameter "feld" (known:
field)`); `withDefaults` fills `default` values before the handler runs.

The schemas several areas name things with are defined once in
`definitions/params.js`: `FIELD` (`{source: "obs" | "var", key}`, the name as
the app shows it, after a rename the new name), `CATEGORY_FIELD`,
`CATEGORY_LABEL`, `COLOR` (`#rrggbb`), `CELLS` (original positions,
`x-cells`), `PAGE` (a highlight page by name), `VIEW` (`"live"` or
`"View N"`), `EMBEDDING` (`{key, columns}`), with the helpers that word them
in labels and pseudocode (`fieldWords`, `fieldArgument`, `categoriesWords`,
`fieldCoalesceKey`).

### The catalog document

`createActionCatalog(definitions)` refuses a type defined twice and answers
`get`, `has`, `list` (sorted by type) and `check(type, params)`.
`catalogDocument(catalog)` is the language-neutral document every other
surface is generated from: `{schema: 1, actions: [{type, tool, title,
summary, group, version, effect, confirm, opensLinks, touches, reads,
captures, context, params, result, example}]}` (`version` null for an action
that records no step), with no functions, byte-stable for a given catalog.

---

## The dispatcher

`createActionDispatcher({catalog})` (`app/actions/dispatcher.js`); the page
has one, `getActions()` in `app/actions/index.js`, beside `getActionCatalog()`
and `getStateDomains()`.

| Member | What it does |
|---|---|
| `dispatch(type, params, meta)` | Apply (or read) one action; resolves `{outcome, step, unrecorded}` (`step` null for a read, a control, a history move or edit, a dataset opened, a change that changed nothing, or an action applied before the first tracker starts). `unrecorded` says why a change made no step: `unchanged` (it was so already), `retracted` (a merged change put back where its step began), `off` (only domains the History settings leave out changed), `paused`, `elsewhere` (another tracker shown meanwhile), `no-history` (no tracker started yet); null for a step made and for the effects that record none. An agent is told it after the result (`agent/bridge/mcp-content.js` `UNRECORDED_NOTES`); a notebook hears it in its `action_result` (`unrecorded`: the closed set `data/viewer-protocol.js` `VIEWER_ACTION_UNRECORDED`, with the sentence Python and R say), where `ActionRequest.unrecorded` / `request$unrecorded` keep it and collecting the outcome warns it. |
| `provide(type, handler)` | Install the implementation; returns its remover. An unknown type, a second implementation of one type, or a handler that is not a function throws. A handler is `(params, {signal, meta, confirm, confirmLink, heavy}) => outcome` (`confirm(detail)`: for a `confirmAt: 'handler'` action, ask the person now; `confirmLink(address)`: for an `opensLinks` action, ask about an address the person has not seen in this tab, an agent's only; `heavy(estimate)`: heavy work found while it runs, asked of an agent as heavy work unless its start was already judged heavy, declined → `ActionRefused` `declined`). |
| `provideEstimate(type, estimator)` | What an action of `type` would cost, from what this page measured: `{ms, bytes, what}` (each null when not known); the presence judges heavy work by it. Provided for `analysis.run`, `cells.read-genes` and `cells.open`, and by the history (`history/actions.js`) for `history.squash`, `history.branch.delete` and `history.step.delete` (what they would refuse, refused before the person is asked) and for what applies steps again, priced by the time those steps took when recorded: `history.step.delete`, `history.branch.without`, `history.branch.apply`, `history.verify` and `history.replay` (the saved state's steps, read from its first chunks; a replay of a step that summed the cells' totals is heavy before it starts). |
| `recordObserved(type, params, meta, outcome, {after, person})` | Record a change its owner already made (below). |
| `addSettler(settle)`, `settle()` | A gesture recorded only once it settles (a camera drag) registers how to end it now; the history calls `settle()` before it moves, so a move never starts before the gesture is a step. |
| `stop(match)` | Stop the dispatches whose meta matches (the person's Stop of an agent turn), from the moment they are dispatched: one waiting for its estimate, the person's answer or the gate never runs; one running ends as a cancellation and a change is rolled back (one that `keepsWhenStopped` answers and keeps its step). |
| `running()` | The changes running now, oldest first (what a move waits for: "Waiting for …"). |
| `provided(type)`, `unprovided()` | Whether a type has an implementation now; every type that has none. |
| `setRecorder`, `setGate`, `setConfirmer` | Installed by the history (recorder and gate) and by the agent presence (`app/agent/presence/index.js`: the person's permissions). |
| `subscribe(listener)` | Every dispatch, reads included, as it starts and ends: `{kind: 'start' | 'done' | 'failed', type, params, meta, step?, outcome?, error?}`; and `asking`, then `answered` (`allowed`), when the person is asked about an agent's change. |
| `idle()` | Resolves once every dispatch has ended, those made meanwhile included, but for those waiting for the person's go-ahead. A dispatch counts from the moment it is made, so one waiting at the history's gate or for its estimate counts; one waiting for the person's answer does not, from the moment the confirmer asks the person (its `asking()`; a change the settings allow at once is never waiting) until the answer: asked at its start, by its implementation (`ctx.confirm`, `ctx.confirmLink`, `ctx.confirmReplacing`, `ctx.heavy`), or about heavy work an internal dispatch it runs finds. It is judged again whenever a dispatch ends or the person is asked or answers. For a control's caller (tests, a host); never awaited inside a dispatch. |
| `whenApplied()` | Resolves once no change is applying (the changes past the gate). What the history's own moves and the save path wait for from inside their turn, where a change waiting at the gate waits for them. |
| `retire()` | The application's teardown (`main.js`), before the state domains' owners are withdrawn: every dispatch under way ends with one cancellation (kind `closed`) wherever it is, its recording included, so none records a step or puts state back into a withdrawn owner; every later dispatch is refused with it. Synchronous, and final for the page. |

### Who dispatches, and with which origin

`meta` is `{origin, agent?, connection?, host?, turn?, signal?}`, plus, for a
turn, what the history's turn record keeps: `model` (an agent's),
`conversation`, `message` and `turnLabel` (what the turn was asked).
`origin` is one of `ACTION_ORIGINS` and defaults to `user`. An `agent` origin
must name its agent (Claude, ChatGPT, Claude Code); `connection` is `{id,
kind: 'chat' | 'connector' | 'local', client, host}` (`AGENT_CONNECTION_KINDS`),
what a permission "for this chat" or "for this connection" keys on. The
history alone passes `rewrite` (a deleted step's dependant applied again into
its own record), `replayOf` (the step a replay applies) and `internal` (a
pass through its own gate).

| Origin | Dispatched by |
|---|---|
| `user` | Every control (`dispatchFromControl`), and `recordObserved` from a control or gesture |
| `agent` | One tool call (`app/agent/bridge/tool-call.js` `runToolCall`): the chat bar (connection kind `chat`), the relay client for Claude or ChatGPT (`connector`) and a local `cellucid mcp` through the viewer channel (`local`), each `{origin: 'agent', agent, connection, turn, signal}` |
| `host` | A Python or R viewer's `action` command (`app/viewer-host/command-handler.js`), with the host's language, and its `turn` inside a notebook turn |
| `replay` | The history re-applying the steps after a deleted one (`internal: true, rewrite: <step>`), a restore as new, and `history.replay` (`replayOf: <step>`) |
| `startup` | Not dispatched: the origin of a tracker's root step and of a change kept as "Changed outside the history" |
| `session` | `session/recompute-results.js`: once a session file or saved state shared without its results is shown, each analysis whose result was left out is run again through `analysis.run`, an ordinary step. A loaded session otherwise puts values back and does not dispatch. Also the origin of a paused draft that a save, a share or a load settles |

### What one dispatch does

1. **Check.** `catalog.check` refuses parameters outside the schema with a
   `TypeError` whose message names the parameter. Nothing is recorded.
2. **Implementation.** A type nobody provides is refused
   (`ActionRefused`, code `unavailable`: "… is not available on this page
   yet"), so an action whose owner is not on the page (a panel never opened)
   says so instead of doing nothing. An analysis action on a panel that is
   hidden is refused by its implementation with a sentence saying the panel
   is not shown.
3. **Stepless effects.** A `read`, a `control` (a stop of work under way), a
   `history` move, a `dataset` opened and an `outside` act run at once, with no recorder and no gate (the history runs its move exclusively itself; waiting
   at its gate would wait for that move), and is stopped with its caller's
   signal. A read stopped ends at once, wherever its implementation is; the
   others end when their implementation does, with what they did.
4. **Ordering.** For every other action: an agent's change passes the
   confirmer (the person's permissions: allowed, refused, or asked while the
   agent waits); the dispatch waits at the history's gate (unless `internal`);
   a `supersede` action aborts the unfinished dispatch of its type
   (`ActionRefused`, code `superseded`); then it waits for every action in
   flight whose `touches` overlap its own. Two actions on overlapping domains
   run one after the other in dispatch order; actions on disjoint domains run
   side by side.
5. **Before.** The recorder captures the touched domains (the history's
   `recorder.before`), healing first any change made outside an action
   (below). That recorder is the one told how the action ended (steps 6 and
   7), whichever is installed by then: another dataset opened meanwhile has
   a history of its own, which never puts back or records the replaced
   dataset's change; the replaced history, retired, does neither.
6. **Handler.** The implementation runs with an abort signal. If it throws
   after changing something, `recorder.failed` captures the touched domains
   and puts the before-values back, so a failed action leaves no half-state;
   a domain the person changed meanwhile by an act of their own (an observed
   change with `person`) keeps the person's change, its own step.
   An action stopped while it runs (the person's Stop, the page closing)
   ends as a cancellation and is put back the same way even when its
   implementation did not watch the signal and ran to its end: the
   dispatcher checks the signal once the handler settles. One superseded by
   a newer dispatch of its type that still ran to its end keeps its step
   (the newer runs right after it; putting back what it replaces at once
   would read again for nothing), and so does one whose definition says it
   keeps what it did when stopped (`keepsWhenStopped`: `camera-path.play`,
   whose camera stays where the film stopped, as a pause leaves it; its
   implementation hears the signal, ends where it is and answers `ended:
   "stopped"`) (`tests/action-dispatch-stop-contract.test.mjs`).
   The dispatcher therefore waits for a stopped change's handler to end:
   a handler that never ends keeps its dispatch in flight, and
   `dispatcher.whenApplied()` (every history move, a save) waits with it, and `idle()` too. Every
   wait of a handler must end, or hear its signal, and hearing it
   means stopping the work it waits for, not only ceasing to wait: work
   that ran on after the rollback would land on the state put back. The
   shared waits take the signal down to the work
   (`tests/analysis-action-stop-contract.test.mjs`,
   `tests/cells-reading-stop-contract.test.mjs`,
   `tests/action-stop-reaches-work-contract.test.mjs`): an analysis's
   `whenSettled({signal})` and the window manager's
   `whenWindowSettled(id, {signal})` stop the analysis's request (its job
   canceled, an update in its delay dropped) and throw once it has ended;
   `runNow({signal})` never starts a run stopped while it waits for the
   open cells or the gene sample; the open dataset read again
   (`reopenOpenDataset(cellSubset, {signal})`, `reopenKeeping`) is canceled
   as its Cancel cancels it; genes read for more cells
   (`readGenesForMoreCellsSaid(..., {signal})`) stop as Stop on their card
   does; a view's embedding (`setViewRepresentation(id, key, {followsDefault,
   signal})` and its callers) stops waiting for its turn or its coordinates
   and publishes nothing; a field copied (`duplicateField(..., {signal})`)
   stops its read; a report's picture steps (`history.report`) jump with the
   signal and come back to the step shown.
7. **After.** The recorder captures the touched domains again and adds one
   step when any hash changed (`recorder.after`). The handler's outcome is
   kept only as scalars on the step (`summarizeOutcome`: counts and names,
   never arrays of cells).

Before the dataset's first tracker starts (while its starting view is put
back), the recorder answers nothing: an action is applied and not recorded,
and what it did is part of the state the tracker starts from.

### Refusals

`ActionRefused(message, {code, type})` is an action that was not applied, with
one sentence a person, an agent or a replay can act on: it says what there is
instead. Codes in use include `unavailable`, `superseded`, `needs-confirmation`,
`declined`, `no-target`, `no-such-page`, `ambiguous-page`, `no-such-group`,
`ambiguous-group`, `ambiguous-branch`, `no-such-step`, `bad-value` and `busy`.
A refusal whose name resolves to nothing carries `missing` (`{kind: 'field' |
'category' | 'gene', name, field}`), which a replay's report keeps so "Use
instead" can offer another name.
A schema problem is a `TypeError`, not an `ActionRefused`; both reach an agent
as its tool's error and a host as its `command_error`.

`app/actions/control-dispatch.js` is how a control dispatches:
`dispatchFromControl(type, params, {failure, category, signal})` dispatches
with origin `user` and says a failure once in a notification
(`reportActionFailure`: a refusal as its own sentence, any other failure after
what was being done, `"<failure>: <message>"`); a cancellation and a
superseded choice are quiet (`isQuietFailure`) and resolve `null`.
`createLatestDispatch(type)` is for a drag: while one dispatch runs, newer
values replace each other and the last is dispatched next, so a fast drag
never queues one step per frame, and the action's `coalesce` makes the drag one
step.

### Observed changes

`recordObserved(type, params, meta, outcome, {after, person})` records a change its owner has
already made. It must be called synchronously from the owner's change event,
so a change the history itself makes while showing another step is never a
step. Its callers on main:

| Caller | Records |
|---|---|
| `ui/modules/settings-actions.js` | `setting.set` for every `input` or `change` on a settings control |
| `ui/modules/camera-controls.js` | `camera.set` when the camera settles after a gesture (a drag, a wheel, keys, a camera path played from the transport bar); consecutive ones on one view coalesce. A path played by `camera-path.play` (an agent, a notebook) is that action's step: the settles on the way are inside it |
| `ui/modules/view-controls.js` | `view.focus` on a click on a pane |
| `ui/modules/layout-actions.js`, `ui/modules/sidebar-controls.js`, `app/dockable-accordions.js` | sections opened, the sidebar's width, panels docked, floated and moved |
| `app/analysis/actions.js` | a person's changes in an analysis panel or window (view, plot, size and window), recorded once the work they started settles |
| `ui/modules/cinematic-camera/index.js` | `camera-path.playback.set`: the Camera Path's playback settings |
| `ui/modules/figure-export/figure-export-ui.js` | `figure-export.set`: the Figure export panel's settings (the `settings` domain's `figureExport` record) |

**A change that sets off more work** (an analysis's redraw, its run after a
delay) passes that work as `after`: it is under way from the change
(`history.idle` waits for it), registered the moment the person makes the
change and recorded once the work settles, one step holding all of it. A
change made while the history restores is ignored at that moment. A camera
gesture is recorded when the viewer says it settled; an action on the
cameras that starts while it still moves, or a move of the history, ends it
first (`addSettler`), so no camera motion is met as a change made outside an
action. Page chrome is no action and no state: the phone sheet, scroll,
focus, which panel is in front, and input modes (the selection tool, Select
mode); the selection tool's pressed button is carried in a session file's
`chrome` ({doc}`10_sessions_persistence_and_serialization`), never a step.

The history records observed changes one after the other, each against what
it shows when recorded, and never while it is putting values back
(`restoringDepth`). A change observed while a dispatched action on the same
domains is applying belongs to that action's step and is not recorded on its
own — unless its caller says it is the person's own act (`{person: true}`:
the analyses' changes, settings controls, a panel's drag or press, the
sidebar's width): that is its own step, and a rollback of the running action
(`recorder.failed`) leaves the domains it changed as the person made them
(the token's `meanwhile`), so a window closed with × while its run runs
stays closed when the stopped run is rolled back. The events a restore or
`setting.set` itself dispatches on a control are never steps
(`isApplyingControls()` in `state-serializer/ui-controls.js`).

### Outcomes as JSON

`outcomeAsJson` (`app/actions/outcome.js`) is what an agent's tool result and
a host's `action_result` carry: typed arrays become arrays, Maps objects (or
pairs), Sets arrays, BigInts their text, Dates ISO text, a non-finite number
`null`; a picture (an image `Blob`) moves to `images` and the JSON says where
it was. Converting never fails on a value the action made, so an action that
changed the viewer is never answered as if it had not; only an outcome that
refers to itself is refused.

---

## Names, never positions (`app/actions/targets/`)

Every parameter names its target the way the app shows it, so a step replays
on another dataset that has the same names, and a refusal can say what there
is. Every action and read resolves names through these modules, never a copy
of its own:

| Module | Resolves | Refuses with |
|---|---|---|
| `targets/fields.js` | A field or gene by `{source, key}` (`resolveField`), a category by its label, a deleted field by the name it had, a new field name | The fields like it, else the fields there are; for a gene, the gene list's own words (`app/utils/gene-lookup.js`) and a name that differs only in case |
| `targets/pages.js` | A highlight page by its name (`resolvePage`), a group by its label (`resolveGroup`) | `no-such-page`, `ambiguous-page`, `no-such-group`, `ambiguous-group`, listing the names there are |
| `targets/views.js` | A view as the Compare views list names it: `"View N"`, the N-th row (the live view unless hidden, then the kept views in the order they were kept) (`resolveViewName`, `viewNameOf`, `openViewIds`); the name a view was given (`view.rename`) is a label beside that number, never its address (`viewLabelOf`, `viewWords`: `View 2 "Early"`, the name quoted as data; `checkViewLabel`) | `no-target`, saying how many views are open; `bad-name` for a name that is not one trimmed line or is of the form "View N". Kept-view ids (`snap_N`) are the viewer's counters and never leave the page |
| `targets/names.js` | How a refusal lists names: `"a", "b" and "c"`, the first eight and how many more (`listNames`); a category's label (`categoryLabel`) | |

Cells are named by original position (their 0-based `obs` row), the
`x-cells` parameters; `data/cell-identity.js` converts them to the open cells.

---

## The state domains

### The sixteen domains

`app/actions/domains.js` lists every domain, in the order a full restore
applies them (each after what it depends on). Every piece of state a person
can change belongs to exactly one. `domainName(id)` is what a person reads
where a step or a comparison lists the parts it changed.

| Order | Domain | Name | Owner module (`app/state/domains/`) | Value (plain data, typed arrays immutable) |
|---:|---|---|---|---|
| 1 | `cells` | Open cells | `cells.js` | `{open: null or the subsample's specification, genes: null or the count of open cells genes were asked to be read for}` |
| 2 | `mapping` | Data mapping | `mapping.js` | The canonical DataMapping, or null with no dataset |
| 3 | `expression` | Expression values | `expression.js` | `{normalize, log1p}`: the expression transform asked (Data mapping › Expression values); both false by default and in a history without it |
| 4 | `field-overlays` | Fields and labels | `field-overlays.js` | `{renamedFields, renamedCategories, deletedFields, userFields}`; user fields' codes and values as columns by original position |
| 5 | `gene-sets` | Gene sets | `gene-sets.js` | `{sets: [{name, description, genes}]}` |
| 6 | `views` | Views | `views.js` | `{layout, focus, liveHidden, camerasLinked, views: [{embedding, followsDefault, centroids, name?}]}` (`name` only on a view given one, `view.rename`) |
| 7 | `active-field` | Active field | `active-field.js` | `{views: [field ref or null]}` |
| 8 | `filters` | Filters | `filters.js` | `{views: [[{field, off?, hidden?, limits?, missingHidden?, outlier?, outlierOff?, range?}]]}` |
| 9 | `colors` | Colors | `colors.js` | `{views: [[{field, colors?, colormap?, range?, log?, rescale?, zeroGray?}]]}` |
| 10 | `visibility` | Hidden cells | `visibility.js` | `{views: [Uint32Array of original positions hidden by hand]}` |
| 11 | `highlights` | Highlights | `highlights.js` | `{active, pages: [{id, name, color, groups: [{id, type, label, enabled, cells, field?, category?, range?}]}]}` |
| 12 | `camera` | Camera | `camera.js` | `{views: [camera state]}`, navigation mode included |
| 13 | `camera-path` | Camera path | `camera-path.js` | The Camera Path's keyframes and playback settings, or null |
| 14 | `settings` | Settings | `settings.js` | `{controls: {id: record}, figureExport}`: the generic control inventory's settings, and the Figure export panel's record |
| 15 | `analysis` | Analysis | `analysis.js` | `{windows: [{modeId, headerTitle, settings, result}], docked}`: each analysis, its settings and its result |
| 16 | `layout` | Layout | `layout.js` | `{theme, sections, panels: [{id, rect}], analysis: {windows: [{rect, open}], sizes}, sidebar: {open, width}}`: every accordion section open or not (every sidebar section), every panel's and analysis window's place, the sizes remembered per analysis (not which is in front) |

An analysis window is split between the two last domains, as every panel's
content and place are: what it shows is `analysis`'s, where it stands, its
size and whether it is folded are `layout`'s, one place per window in window
order. So History settings › "Panels, windows and theme" covers every
window's geometry and never an analysis: moving,
resizing or folding a window (`analysis.window.place`,
`analysis.window.collapse`) and an expanded view's size (`analysis.size.set`)
touch `layout` only; opening, copying or closing one touches both.
`analysis` is restored first (a window it must open again keeps the place the
window at its position had), then `layout` puts every window at its place.

Each owner module's header gives its value in full. `PER_VIEW_DOMAINS`
(`active-field`, `filters`, `colors`, `visibility`, `camera`) hold one entry
per view, live first, then the kept views by position (`view-list.js`).
`REREADING_DOMAINS` is `cells`: its restore reads the dataset again, so every
other value may come out of it changed, and a restore captures every domain
again after it.

### The owner contract

An owner is `{capture, restore, describe}` (and, optionally, `canonical` and
`absent`), provided once per domain to the
page's registry (`createDomainRegistry`: `provide`, `owner`, `has`,
`missing`) by `installStateDomains` (`app/state/domains/index.js`, from
`main.js` once the state, viewer and interface exist; its `OWNERS` map names
one factory per domain).

- `capture()` returns the value now, and never waits for code still loading:
  before the analysis code has arrived, `analysis` and `layout.analysis` are
  the analyses with nothing opened, the same value the code gives on arrival,
  so the history starts and records while it loads (only a restore that shows
  analyses waits for it). The history awaits each capture; a capture inside
  one synchronous step, as the open cells' in-place growth takes the kept
  views', uses only synchronous ones. A value
  is never changed after capture: arrays the app replaces rather than writes
  (highlight groups, user codes and values) are taken by reference and their
  conversions cached by the array and the open cells; the visibility mask,
  written in place, is copied and cached by its write count. Capturing an
  unchanged domain again returns the same arrays, whose content hashes the
  value pool caches by identity.
- `restore(value, {signal, heavy, target})` puts the value back through the
  state methods the app's own controls use, never by dispatching actions;
  changes only what differs; is idempotent (`restore(v)` then `capture()`
  gives `v`). `target(id)` is the value the same restore puts back for
  another domain (undefined when it puts none back): the data mapping applies
  its matrix in the Expression values it ends on, so no matrix is summed for
  values the restore does not keep. `heavy` asks for heavy work the restore
  finds (the cells' totals of normalized values): the running action's
  `ctx.heavy` (an agent's undo, redo, jump or `history.saved.open`), null for
  the person, a replay, a session, and for a rollback. Per-view
  state of a view that is not focused is changed while that view is focused,
  and the focus is put back (`view-list.js` `inEachView`). An async restore
  loads what it needs first (a field's values, a gene once the gene list is
  read, the analysis code), as choosing it does. A per-view value whose view
  count does not fit the views shown is refused, naming both counts.
  `signal` (a move's: the person's Stop of an agent's undo or jump) reaches
  what the restore waits for: the open cells read again (`cells`), a view's
  coordinates (`views`), a field's or a gene's values (the per-view domains,
  `highlights`); stopped, the restore throws why and the history puts back
  what was shown. A rollback passes none: it must finish.
- `describe(a, b)` returns the lines that say what differs between two values
  ("cell_type: 3 categories hidden"): the State panel's Compare two steps and
  the step details read them.
- `absent()` (optional): the value a saved state that does not hold the
  domain stands for (a file written before the domain, or by hand). Loading
  it onto the dataset open puts that value back with the rest
  (`history.js` `readSaved`), so one file restores one state whatever was
  shown before: the Expression values say "as stored" (`{normalize: false,
  log1p: false}`). Without `absent`, a saved state that does not hold the
  domain keeps the value shown.

### Invariants

- Values name things: a field by `{source, key}` with the key the dataset gave
  it (or `{source, user}` for a user field's id), so a rename changes
  `field-overlays` only; a category by its label; a view by its position;
  cells by original position (`cell-positions.js`: a cell list is an
  ascending `Uint32Array`; a column is `{cells, values}`), so a value applies
  to whatever cells are open: cells not open are left out, and a newly open
  cell gets a column's missing marker.
- A value records only what differs from the dataset's opening state where
  opening fills defaults in lazily (filters, colors): opening a field changes
  no value. A color range is measured against the data, never the filter, so
  each change belongs to one domain.
- The sidebar's generic control inventory is split between domains by one map,
  `controlDomainOf(key, type)` (`app/state-serializer/ui-controls.js`): the
  theme and the accordion sections are `layout`'s, the centroid toggles
  `views'` (they show the focused view's own flags), every other control
  `settings'`. A control added to the markup is a setting with no code.
- The state fingerprint (`fingerprint.js`: `stateFingerprint`,
  `captureEveryDomain`, `fingerprintOf`, `differingDomains`) is the hash of
  every domain's content hash; a mismatch names the domain. It is what the
  round-trip specs compare.

---

## The history

The history's records are `app/history/model.js`'s
(`Tracker`, `Branch`, `Group`, `Step`, `Squash`, `By`, `Turn`, `Note`,
`JournalEntry`, `Start`, `Scope`), each with exact fields, rules,
constructors and a validator; `historyModelDocument()` is what
`history_model.schema.json`, the Python and R reading model, is generated
from (`scripts/generate-history-model.mjs`). `app/history/README.md` is the
full map of the core.

### Trackers, branches and steps (`model.js`, `tree.js`)

A **tracker** belongs to one dataset and holds a tree of steps: its root
(`s0`, kind `root`) holds every domain's value hash when the tracker
started; each step below records one change. A dataset can have several
trackers; one is in use. A tracker is `t-<uuid>`; its first branch is `b0`
"Main". A step's id is always greater than its parent's, so id order is
parents-first order, and every list (steps, branches, journal) is stored in
id order. `tree.js` holds a tracker live (its records by id in maps) with
each step's depth, so a common ancestor, a jump's path and redo cost the
steps between, never the depth; walks are loops, never recursion.

| Step field | Holds |
|---|---|
| `id`, `parent`, `children`, `branch` | `s<n>`; the branch it was made on |
| `kind` | `root`, `action`, `outside` (a change made outside an action, "Changed outside the history") or `squashed` |
| `type`, `params`, `v`, `code` | The action (`history.outside-change` for an outside change), its parameters by name (cells by original position), the version of its parameters (`Step.v`) and its pseudocode line; null for the root and a squashed step |
| `by`, `at`, `until`, `merged`, `ms` | Who acted (`By`: origin, agent and model, host, turn, the person's name when set), when, and how long (a coalesced gesture's span and count) |
| `label` | `describe(params, outcome)` |
| `touches`, `before`, `after` | The domains whose hash changed, and their value hashes before and after |
| `reads` | The value hash of each result input the action declares (`reads`): what decides staleness |
| `status`, `reason` | `applied`, or `unapplicable` with the reason; a step is never dropped |
| `outcome`, `coalesceKey`, `replayOf` | The outcome's scalars; the coalescing key; the step a replay applied |
| `squashed` | For a squashed step, what it stands for: `{cause: 'paused' \| 'steps' \| 'clear' \| 'share', count, by, carried}` |
| `notes`, `bookmark`, `group` | Notes on the step, its bookmark name, the group it belongs to |

The value of domain D at step S is the after-value of the nearest step on the
path root → S that touched D, or the root's (`valuesAt`; `valuesFrom` moves
values from one step to another through their common ancestor), so showing
any step needs no recomputation.

- **Branches** (`branches.js`) are records (`Branch`: `id`, a `name` unique
  by key, `base` (the fork step), `tip`, `lane` (its color), `locked`,
  `archived`, `notes`); the tracker names its `activeBranch` (what redo
  follows and the panel shows) and `mainBranch` (what sharing, replay and
  load default to). A new step goes on a new
  branch when the head is not its branch's tip, the branch is locked, or
  "Branch here" is pending; an agent's new branch is named after it.
  Deleting a branch removes its own steps after the last fork; shared ones
  join the oldest branch forking there.
- **Groups** (`Group`, `g<n>`): a named run of consecutive steps of one
  branch, metadata only. The step names its group and the group lists no
  step, so a group costs the same whatever it spans; it holds whole turns
  and no group.
- **Squash** (`history.squash`): consecutive steps become one squashed step
  (cause `steps`) that takes the id of the range's last step, the one
  exception to ids never being reused, so a step or a branch after the range
  keeps its parent; it carries the notes and bookmarks of the steps it
  replaces. `history.clear` keeps the state as one squashed step on a fresh
  branch (cause `clear`), resuming after a pause adds one (cause `paused`),
  and a state shared without its history is one (cause `share`). A squashed
  step restores exactly on its own dataset and replays elsewhere as values
  only, flagged.
- **Names** (branches, trackers, bookmarks): stored NFC with white space
  trimmed and inner runs one space; compared lower-cased; a name of the form
  of an id (`b7`, `t-<uuid>`) is refused, so every parameter takes a name or
  an id without ambiguity.

### The value pool (`value-pool.js`)

Values are kept once, by the SHA-256 of their canonical form, so a value that
comes back (a recolor undone, a filter toggled twice) is stored once, and two
steps show the same state exactly when their hashes agree.

- The canonical text is JSON with object keys sorted, `undefined` properties
  left out, and what JSON cannot hold written as the stored skeleton writes it
  (`value-parts.js`): a typed array as `{"$array": "<digest>"}`, its digest
  `<TypeName>:<length>:<sha256 of its bytes>` (`typedArrayDigest`, computed
  once per array); a number JSON cannot hold as `{"$number": "NaN" |
  "Infinity" | "-Infinity" | "-0"}`. So NaN never hashes like the text
  "NaN", and -0 never like 0.
- A value is plain data: an object of another prototype, a function or a
  symbol is refused where it is captured, and an object of the value's own
  may not use the key `$array`, `$number` or `$part` (it would read as one of
  them; `$part` names a sub-value a session file stores once).
- The pool keeps only what something refers to: a capture nothing took is
  released once nothing shows it, a gesture merged into one step keeps one
  value (`release`, `hold`), and its byte count is a counter kept as values
  come and go.
- SHA-256 is `utils/sha256.js`, so the history works outside a secure
  context (the Python package's plain-HTTP LAN address).
- The pool's typed arrays are part of the memory working set
  (`data/memory-budget.js`, part `'history'`): the person's work, never
  evicted, held while the dataset is shown.

cellucid-python (`src/cellucid/session.py`) and cellucid-r (`R/json.R`,
`src/number_text.c`) compute the same hashes when they read a session file.

### Recording

`createHistory({dispatcher, domains, dataset, newTrackerId})` installs itself
as the dispatcher's recorder and gate. A step is added when an action changes
at least one hash; a no-op is no step. It goes on the branch the head is on,
with who acted, how long, the version of its parameters, its pseudocode and
the hashes of its result inputs. A step whose `coalesceKey` equals the step
just recorded merges into it (after-values and parameters become the newer
ones, before-values stay) when that step has no children and has the same
origin and turn. A domain that follows another (`DERIVED_DOMAINS`) is
captured with it.

A capture that disagrees with what the history shows is a change made outside
an action: it is kept as its own step (kind `outside`, "Changed outside the
history"), so undo still lands where it was, and warned about in the
console; tests treat it as a failure. It is healed before a dispatched
action's before-capture (for the domains it touches) and by `reconcile()`
(every domain), which a save runs first.

**The recording switch** (`Tracker.paused`: null while recording, else
`{at, by, draft}`; `history.recording.pause`, `history.recording.resume`).
While paused, changes apply and make no step; the draft of what differs from
the head is kept on the tracker (and by the autosave, so a crash keeps it).
Resuming, a move (the person's choice: keep, the default, or discard), a
save, a share or a load (`settlePaused`) make it one squashed step. A file
keeps `draft: null`.

**History settings** (`history-settings.js`, this browser's `localStorage`
`cellucid:history-settings`, never part of a session or a step):
`recordCamera` and `recordLayout` (the two detail switches: off, that domain
makes no step and the history never moves it, while steps still hold its
value, `offDomainsOf`), `recordNewDatasets` (off: a new dataset's first
tracker starts paused), `askWhilePaused`, `personName` (`By.person` in what
the person shares; unset by default), `autosave` and `reportPictures`.
Nothing else can be switched off: everything else is the analysis.

### Moves put values back

`undo()`, `redo()`, `jump(stepId)` and a branch or tracker switch show
another step by restoring only the domains whose value differs from what is
shown, in domain order, then capturing what is shown again, so a value its
owner normalizes (a camera's last digits) does not read as a change later;
after `cells`, every domain is captured again. An analysis result, a merged
field or a camera comes back as it was, without being computed again. If a
restore fails, what was shown before is put back, and a rollback failure is
attached to the error.

**Every move and every journaled edit runs alone** (`exclusive`): gestures
still settling are ended (`dispatcher.settle()`), actions wait at the gate,
the actions and observed changes under way finish first, and a change
observed meanwhile is recorded after it against what is shown then. The
inner operations refuse to run outside one. A history retired (its dataset
replaced) puts nothing more back. Each is given the signal of the action
that asked for it (the person's Stop of an agent turn reaches it): stopped
while it waits for its turn it never runs and its caller hears so at once;
a move stopped while it puts a state back puts back what was shown and keeps
the head (`tests/history-move-stop-contract.test.mjs`). The waits before the
work changes anything hear the Stop the same way: a replay's read of the saved
state (`readReplayPlan`, the one read Replay… and `history.replay` share), and
a saved file's wait for its turn among the session's operations and for the
actions running (`session-serializer.js` `_exclusive`, `writeShare`;
`tests/history-replay-stop.test.mjs`).

### The journal and Recently deleted (`journal.js`, `delete-restore.js`)

Every structural edit of the history (a branch or step deleted, a squash, a
clear, a group or note edit, a rename, a tracker removed, a file loaded
beside) is a `JournalEntry` in its tracker's own journal (`j<n>`), holding
the exact records it replaced and made (a group edit only which steps each
group held, `memberships`). The State panel's Recently deleted lists them.

- **Restore** (`history.journal.undo`) puts them back exactly, as an entry
  of its own (undoable; tried on the records and checked whole, rolled back
  when the result would not be a valid tracker; a step added below one since,
  or a move since, keeps it exact). When later edits changed what the entry
  made, the history **restores as new**: the steps (a squash's run in full)
  are applied again by name as branches, each line on its restored parent,
  and the restore is journaled itself. Its answer is `{entry, restored:
  'exact' | 'as-new', journal, branch, report}`: `branch` the first branch a
  restore as new made (or null), `report` each step applied again and what
  became of it (the replay entries below).
- Undoing a rename whose old name was taken since is refused, naming the
  name to use instead.
- **Delete a step** applies again only the steps below it whose action writes
  or reads a domain whose value changed there (by touched domains and result
  inputs, ignoring generated ids), each on its new parent's values (origin
  `replay`, `rewrite: <step>`); a step that no longer applies is kept,
  `unapplicable`, with the refusal as its reason. Cancel puts everything
  back exactly. Values are compared in their owners' canonical form
  (`canonical`: generated ids renamed, times dropped), and every user
  field's generated id a value names (its colors, a filter of it) is
  renamed as the `field-overlays` value of the same state names it
  (`user_1`, …): a field made again under a new id is the same field.
- **Removing a tracker** is soft (`Tracker.removed = {at, by}` and an entry
  in its own journal); restore is `history.journal.undo`; emptying the
  journal (`history.journal.empty`) drops the patches and so purges.

### Turns (`turns.js`)

An agent's dispatches carry its `turn`, and so do a Python or R host's
inside a notebook turn (`viewer.turn()`, `cellucid_turn()`). A `Turn` record
keeps the agent, its model, the conversation, the message, the label and the
branch. Edits inside a turn are exclusive to it: the turn's binding to its
branch refuses only when the person moved meanwhile (`moved`).
`turnUndoPlan` says what undoing it would do: back to before it when it is
alone at the end of its branch, or only its own steps when others
interleaved (the person chooses: "Undo turn" with two choices).
`history.turn.undo` and the presence's Undo turn (`app/agent/history-link.js`
`undoAgentTurn`) take this path.

(declared-inputs)=
### Declared inputs: reads and context

A step records the value hash of each result input its action declares
(`Step.reads`). `provenance.js` finds which steps are stale at the head: a
step of any action whose result input differs later, with each changed
input's own cause (`causes`; `staleReasons` words them, one sentence per
input, for `history.stale`, `history.step` and the State panel). An input
held per view is compared at the view the step acted in (the domain owner's
`readPart(value, views)`, `view-list.js` `viewReadPart`: the entry at the
`views` value's `focus` when the step ran), so keeping a view or changing
another view's filters marks nothing stale. With several views, a
step that changed some per-view domain in some views, not all, names them in
its label (`history.js` `labelWithViews`, `viewNameInValue`), which the
methods text and the pseudocode report repeat. The state an action's parameters' names resolve in is its
`context`, never recorded, so a rename or a focus change never marks a step
stale. Result inputs and name context together decide which dependants a
deleted step applies again. `tests/action-reads-contract.test.mjs` checks the
declarations (distinct, known, never a touched domain, never both, empty for
a stepless action, every name a parameter gives resolved in one of them);
`tests/browser/action-reads.spec.mjs` checks the property on the real page:
an action applied in two states that differ only in domains it neither
touches, reads nor resolves names in gives the same after-values.
A stale step offers **Run again** only when its action is marked `rerun` in
the catalog: the same parameters dispatched again as a new step, by
`ops.runAgain`.

### Replay on another dataset (`replay.js`)

`history.replay` ("run my workflow on my next sample") reads a saved history
(a session file, a saved state, of this dataset or another) whole, checks it
as a file, and takes a branch (main by default). It starts a new tracker on
the state shown whose start names where its steps come from (`start.kind
'replay'`, a `Link`), and dispatches each step again by name (origin
`replay`, `replayOf: <step>`, its parameters brought to this catalog's
version by `upgradeParams`) while the person's actions wait at the gate.
Names map by the parameter schema (`mapParams`: fields, genes, categories),
and the person can map a missing name to another. The source's domain values
are never put back.

Each step's report entry (`ReplayEntry`) is `{source, label, outcome,
step, reason, missing}`: `outcome` `applied`, `unchanged`, `unmatched`,
`skipped` or `stopped`; `missing` the name that resolved to nothing here
(`{kind: 'field' | 'category' | 'gene', name, field}`, from the refusal's
`missing`) or null. A step that does not apply here changes nothing and is
kept, flagged with the reason; a step flagged where it was recorded is kept
flagged; a step of an unknown action, or recorded with a newer version, is
kept flagged ("Recorded with a newer Cellucid"); a squashed step is values
only, flagged; a change made outside the steps is skipped; a stop keeps the
steps applied before it. A file of steps only (a share of `content:
'steps'`) opens as a replay plan.

### Loading, conversations and the rest

- **Load beside** (`importTrackers`): a session file chosen with Load is added
  beside the person's trackers; a tracker whose id is taken here gets another
  (its conversations, Recently deleted and late parameters follow it), a
  name taken is suffixed "(shared)", and each tracker added is journaled
  (`load`): the notice's **Undo** removes them to Recently deleted.
- **Load** (`load(saved, {values, rest, signal})`): a saved state of this
  browser replaces the history, all or nothing
  ({doc}`10_sessions_persistence_and_serialization`).
- **Conversations** kept with a tracker: `setChat(tracker, conversation,
  transcript)` and `chatsOf(tracker)`, plain JSON. The chat bar keeps its
  conversation after every message (`app/agent/chat/conversation.js`
  `transcriptOf`: no providers' reasoning, each picture replaced by a line
  saying it was there). A session file holds them only when the share or
  save includes them (`chat/<n>/<m>` chunks).
- **Before the first tracker**: `showValues(values)` puts a starting view back
  without a step; the first tracker then starts from it.
- **Diff** (`diff(a, b)`): for each domain whose value differs between two
  steps, the owner's `describe` lines, and their common ancestor.
- **Search** (`search.js`, `history.search`): steps by words and filters.
- **Verify** (`history.verify`): whether the steps still give the values
  they recorded, comparing values in their owners' canonical form.

### The history's own actions

`app/history/actions.js` provides every `history.*` action but the methods,
the pseudocode, share, save and the saved states' list and delete, which the
persistence provides (`persistence-actions.js`), and `history.anndata.save`,
which `anndata-copy-save.js` provides ({ref}`history-anndata-save-action`). Reads (effect `read`): `history.get` (the
tracker in use, the trackers, the head, the branches and the steps on the
path to the step shown, or every branch with `all_branches`, each with who
took it and its pseudocode line, a page at a time: `limit` steps from
`cursor`, answered with `total` and the `next` cursor; `describe.js`),
`history.summary`, `history.step`, `history.diff`, `history.search`,
`history.stale`, `history.methods`, `history.pseudocode`,
`history.saved.list`. Deleting a saved state of this browser
(`history.saved.delete`, effect `outside`, `deletes-saved`): by its id from
`history.saved.list`, through `autosave.js` `removeSavedStates` (the one way
saved states are removed, which resumes an autosave paused for room); the
history shown stays open. Opening one to continue it (`history.saved.open`,
effect `dataset`, asked as `large-download` from its handler): Load…'s Open
(`SavedHistories.open`, the session controls' `openSavedState`: its dataset
first, with its cells, then its history), run as an agent's opening by
`ui/modules/dataset-open-action.js`, so a whole-file download of its
dataset is asked about as `dataset.open`'s is; refused before anything opens
for an unknown id or a state this version cannot open, and, when it did not
open, with why (`needs-person` when the page offers the person to reopen a
file of this computer, `canceled` when the person stopped the loading, the
load's failure, or the opening's); the history
shown already answers `already_shown`. A chat turn that opens one goes on
with the history it loaded (`agent/chat/index.js` `followTracker`: a load
during the turn's own dataset action). Moves (effect `history`): `history.undo` and
`history.redo` answer `{moved, head}` (`moved` false when there was nothing
to undo or redo), `history.jump` answers `{head}`; none records a step. Edits
(effect `history`): branches (create, rename, switch, delete, archive, lock,
set main, cancel, apply, without), groups, notes, bookmarks, trackers
(start, use, rename, remove), the journal (`history.journal.undo`,
`history.journal.empty`), a step deleted, a turn undone, `history.replay`,
`history.squash`, `history.clear`, `history.verify` and the recording
switch. Writing it out is effect `outside`: `history.share`, `history.save`,
`history.anndata.save` and `history.report`; the State panel provides
`history.report` and the guided tour (`history.tour.start`,
`history.tour.stop`, effect `history`) from `ui/modules/state/index.js`.
`history.step.revert` is a `state` step. A jump to a step the tracker does not hold is refused (`no-such-step`)
naming the tracker. Before the dataset's first tracker starts each is
refused with "The history starts once a dataset is open." The State panel,
agents and Python or R hosts all act through these actions.

### A dataset's history from opening to replacement (`dataset-history.js`)

`createDatasetHistories` (made in `main.js`; `window._cellucidHistories` for tests):

- **open**, when a dataset is published: its history is made with no tracker;
  a starting view the catalog names is put back into it without a step
  (`restorePublishedDefaultState` → `showValues`);
- **begin**, once the starting view is applied (or there is none): the first
  tracker starts from what is shown (paused when History settings say not to
  record new datasets), the autosave starts, and the dataset's latest
  autosave of an earlier visit is offered;
- **retire**, before another dataset (or none) is published: the history
  stops recording at once and the autosave writes what changed.

The same dataset read again on other cells is not another dataset: its history
is kept (`app/dataset/publication.js` `keepHistory`). **Cells to open** is the action `cells.open`
(provided by `ui/modules/subsample-control.js`), a step that touches `cells`
and every domain a cell can be named in; undoing it puts the `cells` value
back, which reads the cells shown before again (`state/domains/cells.js`), and
redo reads the share again, on the same tracker
(`tests/browser/history-open-cells.spec.mjs`). The reading keeps what was made
by original cell position (`session/view-carry.js`, `app/cell-subset-carry.js`).

(history-autosave)=
### Autosave and saved states (`autosave.js`, `autosave-store.js`)

- **Every change as it is made.** The write of each change of the history
  starts within a task of it, to IndexedDB `cellucid-history`: only what
  changed (the steps, turns, journal entries and messages an event names, the
  tracker's own record), each value once by its hash, each sub-value once by
  name (parts) and each array once by its digest (cell sets in the cell-set
  encoding), in one transaction, so a crash leaves all of a write or none. A
  history is kept once it holds anything of the person's (a step, a
  conversation, a note, a paused draft).
- **All of a write or none, whatever becomes of the page.** A transaction of
  `utils/indexed-db-database.js` reads first and then makes every write at
  once, none awaited, and a read after a write throws: the writes are all in
  the transaction before any has run, so the browser runs them all or keeps
  none. It is then told to commit (`IDBTransaction.commit`), so a write handed
  to the browser is finished whether or not the page is still there.
  `autosave-store.js` plans each write from its reads (`planWrite`,
  `planCollect`) and makes it in one go; `tests/indexed-db-database.test.mjs`
  holds the rule, `tests/browser/saved-state-close.spec.mjs` closes a tab
  during a write in each engine (what the engines did before, and the counts,
  are in {doc}`33_decisions_actions_domains_and_history`).
- **What a closed tab can lose.** The change of its last moment: a write is
  handed to the browser some tens of milliseconds after the change (a task,
  the session header, the digests of new values), and a page that goes
  inside that time has not handed it over. A page that is left retires its
  autosave, which writes what is left if the browser gives it the time;
  nothing can make it.
- **One writer.** One page writes a saved state (a Web Lock, and the store's
  writer check): a second tab continuing it opens a copy.
- **Continue.** When a dataset opens whose autosave holds steps, a
  notification offers Continue or Start fresh; it is never loaded by itself.
  My datasets' Continue opens the dataset, then its latest saved state.
- **Named saves** (`history.save`: the history, the state only or the steps
  only) are kept beside the autosave and do not change afterwards. The Load
  dialog lists, opens, replays here, renames and deletes them.
- **Room.** A write the browser refuses for room pauses the autosave and says
  so once; any removal of saved states, in any tab, resumes it. A saved state
  removed while it is written is written whole again.
- A record of another layout (`SAVED_STATE_SCHEMA`) is listed, never
  opened. A browser with no IndexedDB keeps no saved state.

### The publication (`current.js`)

The page's history and this browser's saved histories are published once for
everything that shows or reads them: `setCurrentHistory`,
`getCurrentHistory`, `setSavedHistories`, `getSavedHistories`,
`subscribeCurrentHistory`. The dataset lifecycle publishes the history it
creates (and `null` when the dataset goes); the autosave publishes the store.
The State panel (`ui/modules/state/`), the keyboard's undo keys and the agent's
history reads subscribe; none of them keeps a history of its own.

### The history in words (`app/history/exports/`)

Pure modules (no DOM, importable in Node) that take the history's records as
they are (`app/history/model.js`): a tracker, one of its branches (by id or name;
the main branch when left out) and optionally the last step to cover (`to`, the
step shown after an undo), the dataset's facts and the version. Any branch can
be exported without showing it. They read each step's `kind`, `by`, `squashed`,
`group`, `bookmark` and `notes`, and the tracker's `groups` and `notes`; never
its recording switch. The State panel (`ui/modules/state/share.js`) calls them;
Python (`cellucid.history`) and R (`R/history.R`) write the same methods text
and pseudocode, pinned by the golden history's `expected.methods` and
`expected.pseudocode`.

| Module | Entry points | Writes |
|---|---|---|
| `step-text.js` | `originOf(by)`, `originPhrase(by)`, `actorsOf`, `squashPhrase`, `bookmarksOf`, `parameterFacts`, `branchPath`, `groupRuns`, `runHeadings`, `pathOf` | How a step is said: who acted (a `By`), what a squashed step stands for, its parameters as facts; the path an export covers and its runs of grouped steps with their headings |
| `pseudocode.js` | `historyPseudocode({tracker, branch, to, dataset, version})` | Each step's stored `code` line; agents and notebooks as trailing comments; `# --- QC ---` … `# --- end of QC ---` around a group; `# Bookmark: <name>` (a squashed step's carried bookmarks marked `(on a step merged into this one)`); a flagged step, an outside change and a squashed step (values only) as comments |
| `methods.js` | `methodsParagraph(…)`, `methodsParts(…)` (same options plus `cellucidCitation`) | The tracker's first note (its abstract) as the opening paragraph; software and version, the dataset with its size, citation and cells opened; the steps that changed the data in order, under group headings when a listed step is in a group; a squashed step said as what it stands for, with "only the state it reached was kept"; presentation steps counted, not listed |
| `report.js` | `reportDocument({tracker, branch, to, catalog, generatedAt, …, pictures})`, `reportEntries`, `reportPictureDefaults`, `formatWhen` | One self-contained, escaped HTML document: the abstract and the other tracker and branch notes, the dataset, the methods text, every step by group (its bookmarks and notes, those a squashed step carried from the steps it replaced marked, who acted, its parameters, a flag for a flagged, outside or squashed step, its notes, and a picture when the caller passes one), the pseudocode. No picture is stored in a history; the caller draws the picked steps from the live view (`reportPictureDefaults`: the bookmarks and the last step of each group) |

(history-anndata-save-action)=
### Save AnnData with this state (`anndata-copy-save.js`, `exports/anndata-copy.js`)

`history.anndata.save` (effect `outside`, `confirm: 'writes-file'`, provided
by `provideAnnDataCopyAction` in `main.js`) is one path for the History ⋯
menu's **Save AnnData with this state…** and for agents, Python and R
(`save_anndata(...)` in its pseudocode line). It writes a **new**
Cellucid-ready AnnData `.zarr.zip` of the open dataset with the state shown
applied, never over the dataset's own file:

1. **The destination first**, while the person's click still lets the page
   open a dialog: `showSaveFilePicker` (origin `user`, user activation, a
   browser that has it), whose chosen file is refused when `isSameEntry` with
   the dataset's own file handle; otherwise a spool in the origin-private
   file system (`exports/cell-exports.js` `createSpool`), downloaded once
   whole. Agents and hosts always get the spool. The name is
   `copyFileName` ("<dataset> with this state.zarr.zip") or `file`
   (`copyFileNamed`: `.zarr.zip` added unless the input-format registry reads
   the name as a zipped AnnData or MuData Zarr).
2. **The plan.** The paused changes become a step as Share writes them
   (`writeShare` with `STATE_DEFAULTS`, whose blob's SHA-256 is the
   provenance's `session_sha256`). `planAnnDataCopy` (pure) plans the step
   shown against the dataset as stored (`BaseAnnDataAdapter.readCopyTarget`:
   obs columns, the columns a row reads whole, `uns['cellucid']`, the var
   columns of the mapping's modality), with the open cells as the cells the
   state is known on. Its rows, statuses and var columns are the packages'
   (`expected.apply.sequence[].rows` and `.var_columns` of the shared golden
   history, `tests/history-anndata-copy.test.mjs`); its provenance has
   `package: 'web'`. Gene sets get `cellucid_gene_set__<name>` booleans on
   the var of the mapping's genes (`geneSetColumns`); a mapping reading
   `raw` has none (the sets stay in `uns`, said in the left-out lines).
   Each written, renamed or relabeled column's colors go with it as
   `CopyChanges.uns` (`colorChanges`: `uns['<column>_colors']`, a
   categorical field's colors as the view in focus shows them,
   `shownCategoryColors` -- its recolored categories (`focusedRecolors`, the
   `colors` value), else its `palette`, else `getCategoryColor` -- or none,
   carried by a rename, a merged category taking its
   first category's color; the target reads the root `uns` keys and the
   colors it carries, `readCopyTarget(..., {colors})`), the packages' rule
   ({ref}`apply-column-colors`). A dataset column the view in focus
   recolored (`focusedColumnRecolors`, by the key the dataset gives it) is
   one `column-colors` row (domain `colors`, before the colors presentation
   row) from `planColumnColors`: its categories, by label once the apply's
   renames are applied, in the colors given, the rest in the colors it opens
   with (`openingPalette`: the file's colors as drawn when there is one a
   category and each is hex, a merged category taking its first source's,
   else `getCategoryColor` by position); `write`, `unchanged` (the colors it
   would hold after the apply's other changes), or `not-applied` (absent, a
   modality's column, not categorical, none of its labels there). The
   presentation row's note is `COLORS_NOTE` when a field's recolors were
   written or a `column-colors` row is not `not-applied`.
3. **The checks**, before a byte: `unwritableName` (a gene set or field
   name `requireElementKey` refuses), then `requireCopyChanges`
   (`data/anndata/write/changes.js`: every name, column and `uns` value;
   `PRIORITY_COLUMN_REFUSAL` for a change that writes, renames or replaces
   `obs['cellucid_original_position']`). A dataset whose matrix still
   streams in (`AnnDataReader.storedRoot`) and a cell count other than the
   open dataset's are refused too.
4. **The copy**: `BaseAnnDataAdapter.writeAnnDataCopy` into the
   destination, as one `CancellableJob` (kind `cells` reads, Cancel on its
   card); the outcome is told once with the left-out lines (the plan's rows
   not applied, the copy's notes) in its details. On a failure or a cancel
   the chosen file's writable is aborted (never closed, so no partial file
   is left under its name) and the spool's files are removed; nothing is
   downloaded.

The result (`file`, `to: 'file' | 'download'`, `bytes`, `layout`, `cells`,
the tracker, branch and step, `applied`, `rows`, `left_out`) is the action's.
Tests: `tests/history-anndata-copy.test.mjs` (the plan against the golden
history), `tests/history-anndata-copy-save.test.mjs` (destinations, refusals,
cancel), `tests/anndata-copy-contract.test.mjs` (the writer's checks) and
`tests/browser/history-anndata-copy.spec.mjs` (saved from the menu and read
back by the app and by anndata and cellucid-python, a cancel mid-write, a
MuData). The file's profile:
{ref}`browser copy <browser-copy-profile>`.

---

## One catalog, every surface

| Surface | How it gets the catalog | Regenerated by | Fails when stale |
|---|---|---|---|
| The page | `ACTION_DEFINITIONS` → `createActionCatalog` (`getActionCatalog()`) | — | `defineAction` at import; `tests/action-catalog-contract.test.mjs` |
| The shared document | `tests/fixtures/action-catalog.json` (cellucid), `src/cellucid/action_catalog.json` (cellucid-python), `inst/action_catalog.json` (cellucid-r): `catalogDocument` as canonical JSON (keys sorted, two-space indent, one trailing newline), byte for byte the same: the web repository's file is its home, and each package's is a copy of it ({ref}`shared-files`) | `node scripts/sync-shared-files.mjs` (it runs `scripts/generate-action-catalog.mjs` and writes every copy; finds the siblings beside the checkout or at `CELLUCID_PYTHON_ROOT` / `CELLUCID_R_ROOT`; `--check` writes nothing) | `tests/action-catalog-document-contract.test.mjs` and `tests/shared-files-contract.test.mjs`, cellucid-python `tests/test_action_catalog.py`, cellucid-r `tests/testthat/test-actions.R` |
| Python | One method of `cellucid.jupyter.AnnDataViewer` (and the notebook widget) per action, `src/cellucid/jupyter/_actions_generated.py`, each calling `act(type, …)`, which converts parameters by the catalog (`_actions.py`) | `python scripts/generate_actions.py` (`--check`) | `tests/test_action_catalog.py` |
| R | One exported `cellucid_<tool>(viewer, …)` per action (`R/actions-generated.R`), its `man/cellucid_<tool>.Rd` page, the NAMESPACE export block and `_pkgdown.yml`'s actions group; each calls `.act()` (`R/actions.R`) | `Rscript tools/generate-actions.R` | `tests/testthat/test-actions.R` |
| The viewer channel | A host sends one `action` command (`action`: the dotted type, `params`; `turn` and `turnLabel`: the host's turn, both text or both null), answered by `action_result` (its outcome, the step it made, or why a change made none: `unrecorded`) or `command_error`; it pushes a session file with `load_session` (answered by `session_loaded`) and hears every change of the history as `history_changed`; viewer protocol schema 5 (`data/viewer-protocol.js`, mirrored as `viewer_protocol.json`; schema 5 added the local agent messages, {doc}`31_local_agent_connection`) | — | `tests/viewer-protocol-contract.test.mjs` (the protocol's copy of `ACTION_TYPE_PATTERN` is the catalog's) |
| Agents | Every action is one MCP tool (`app/agent/bridge/mcp-tools.js` `mcpToolsOf`: its tool name, input schema, title, and hints from `effect` and `confirm`; the same list whatever the person's settings, which are told after tool results instead, `settings-told.js`); the relay's copy for when no tab is connected is generated, and so is `mcp-tools.json` (the tools, the local instructions and status tool) for the packages' `cellucid mcp` | `node scripts/generate-relay-bundle.mjs` (`--check`); `node scripts/generate-action-catalog.mjs` writes `mcp-tools.json` into the three repositories | `tests/relay-bundle-contract.test.mjs`; `tests/mcp-tools-document-contract.test.mjs`; `tests/agent-bridge.test.mjs` (the agent instructions name only tools that exist) |
| The history model | `history_model.schema.json` (`historyModelDocument()`), in the three repositories: the web repository's `tests/fixtures/history-model.schema.json` is its home, each package's a copy ({ref}`shared-files`) | `node scripts/sync-shared-files.mjs` (it runs `scripts/generate-history-model.mjs`), then `python scripts/generate_history_model.py` and `Rscript tools/generate-history-model.R` | `tests/history-model-schema.test.mjs` and `tests/shared-files-contract.test.mjs`, cellucid-python `tests/test_history_model.py`, cellucid-r `tests/testthat/test-history-model.R` |
| The docs | {doc}`30_action_reference`: every action with its parameters, effect, label and pseudocode, every read with what it returns | `node scripts/generate-action-reference.mjs` (`--check`; writes into `CELLUCID_PYTHON_DOCS`, else the sibling checkout) | `tests/action-reference-contract.test.mjs` |

---

## How to add an action

1. **Define it** in `app/actions/definitions/<area>.js` with `defineAction`
   (a new area is one module plus one import and spread in
   `definitions/index.js`). Name targets with the shared schemas of
   `definitions/params.js`; every object schema closes its properties; give
   `examples` where the derived example would not be meaningful.
2. **`touches`: exactly the domains the handler can change.** Not fewer (the
   history would not capture what changed, so undo would not put it back and
   the next capture would record it as a change made outside an action), and
   not more than it can change. A rename that reaches every reference of a
   field touches every domain that names fields (`field.rename`).
3. **`effect`**, and `coalesce` for a value dragged or typed, `supersede` for a
   choice a newer one replaces while it is still read, `confirm` only for
   what the history cannot take back. Declare `reads` (its result inputs)
   and, where its names resolve in other state, `context`
   ({ref}`declared-inputs`). An action that answers with data is classified
   for the agent data guard (`app/agent/bridge/data-guard.js`,
   {doc}`32_agents_in_the_tab`); one nobody classified is refused to agents
   while a dataset's switch "names and summaries only" is on. A later change
   to its `params` bumps `version` and adds an upgrader.
4. **`describe(params, outcome)` and `code(params)`**: one line each, from the
   parameters (names, never positions), so a shared history reads the same on
   any machine.
5. **Provide the implementation from the module that owns the state**:
   `getActions().provide(type, handler)` when the owner starts, and call the
   returned remover when it ends (`provided` and `unprovided` say what is
   installed). The state's own actions are in `app/state/actions/field-actions.js`;
   a control's orchestration stays with the control (`ui/modules/*.js`). The
   handler resolves every name through `app/actions/targets/`, refuses what
   does not resolve with `ActionRefused` and one sentence saying what there
   is, calls the `DataState` (or viewer) method the control used to call, and
   returns a small outcome (names, counts, ids) or `null` when nothing
   changed, never a per-cell array.
6. **The control dispatches**: `dispatchFromControl(type, params, {failure,
   category})`, or `createLatestDispatch(type)` for a drag; it calls no
   mutating state method. A change the owner must make first (a gesture, a
   settings control) is recorded with `recordObserved`, synchronously from its
   event.
7. **Regenerate every surface**: `node scripts/sync-shared-files.mjs` (it
   writes the catalog and `mcp-tools.json` and copies them into both
   packages; {ref}`shared-files`), `python scripts/generate_actions.py` in
   cellucid-python, `Rscript tools/generate-actions.R` in cellucid-r, and run
   `node scripts/generate-relay-bundle.mjs`, `node
   scripts/generate-action-versions.mjs` (a new action's or version's
   examples) and `node scripts/generate-action-reference.mjs` (the docs'
   reference page). Nothing else is updated by hand: each suite compares its
   copy with the web repository's file.
8. **Tests.** `npm test` runs `tests/action-catalog-contract.test.mjs` (the
   example, label and pseudocode, the domains, a result schema for a read)
   and the document and relay checks; add the action to its area's tests
   (`action-fields-contract`, `action-highlights`, `action-camera-path`,
   `action-views-camera`, `actions-settings-layout-analysis`) and to its
   area's browser spec, which dispatches it on the page's own history and
   proves `touches` with `applyUndoRedo`
   (`tests/browser/helpers/state-history.mjs`): the domains it changed must be
   among those it names, undo must put every domain back (the state
   fingerprint) and redo must give them again. Run `python -m pytest
   tests/test_action_catalog.py` in cellucid-python and
   `Rscript -e 'testthat::test_local(filter = "actions")'` in cellucid-r.

## How to add a piece of state

**A part of an existing domain** (a new option of a filter, a new setting of a
panel):

1. Extend the owner module in `app/state/domains/`: `capture` names it by name
   (cells by original position, typed arrays never written after capture),
   records only what differs from the dataset's opening state when opening
   fills it lazily, and `restore` puts it back through the state methods the
   control uses, changing only what differs, idempotently; `describe` says
   it. Update the module's header, which gives the value in full.
2. Every action that can change it lists the domain in `touches`, and the
   control dispatches one of those actions. A new sidebar `input`, `select` or
   pressed-button group with a stable id outside a skipped subtree is already
   a setting (`settings` domain, `setting.set`, recorded by observation) with
   no code.
3. The session file stores it as it is: there is no encoder per feature. If
   the domain is one the packages read (`cells`, `field-overlays`,
   `highlights`, `gene-sets`: cellucid-python `anndata_session/`, cellucid-r
   `R/session-apply.R`), the shape is a contract between three repositories:
   change their readers in the same wave and regenerate the shared format-4
   fixture ({doc}`10_sessions_persistence_and_serialization`).
4. Cover it in `tests/state-domains-contract.test.mjs` (value shape, exact and
   idempotent restore, describe lines) and in
   `tests/browser/state-domains-roundtrip.spec.mjs` (changed through the
   interface, captured, changed again, restored: value, visible state and
   fingerprint back).

**A new domain**: add `[id, name, description]` to `DOMAINS` in
`app/actions/domains.js` at its place in restore order (after what it depends
on; a per-view domain after `views`, and in `PER_VIEW_DOMAINS`), one owner
module in `app/state/domains/` and one line in `OWNERS`
(`app/state/domains/index.js`); then everything above. A session or saved
state written before it does not hold the new domain: loading it keeps the
value shown now for that domain. A history holding a domain this version does
not know is refused, naming it. A domain removed or renamed is a format
change: update every producer, consumer and fixture together.

---

## Drift guards

The guards that fail when the layers drift apart (the catalog and its
generated surfaces, the domains and their owners, the history's semantics,
the control walk) are listed with the other test rules, each with what makes
it fail: {ref}`drift-guards`; no control calls a mutating state method
itself, which the static mutation scan holds at the source.

Next: {doc}`10_sessions_persistence_and_serialization` (the file a history is
saved in).
