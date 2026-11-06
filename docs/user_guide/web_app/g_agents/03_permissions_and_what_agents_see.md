# Permissions and what agents see

**Audience:** everyone who lets an agent work in Cellucid
**Time:** 10 minutes
**What you'll learn:**
- What an agent may do without asking, and how to change it
- The confirm card, the activity line, Stop and Undo turn
- What agents are sent of your data, and the switch **Agents see names and
  summaries only**
- How text written by others reaches an agent

These rules hold for every agent acting in the tab, whichever way it came
in: the chat bar, Claude or ChatGPT through the connector, or a local app
through `cellucid mcp` / `cellucid_mcp()`. They are enforced in the page,
where every agent's call runs the same way.

---

## Permissions

Reading the dataset never needs permission. Every change an agent asks for
falls under one or more **kinds**, each set to **Ask**, **Allow** or
**Deny**; the strictest kind of a request decides.

| Kind | What it covers |
|---|---|
| **Changes in this tab** | colour, filter, highlight, views, new fields, analyses (each a step the history can undo), and opening another dataset (its history stays kept) |
| **Heavy work (longer than about 10 s)** | work this device and connection, as measured, would take longer than about ten seconds to do: a long analysis, reading genes for many cells, opening many cells |
| **Writing files** | writing a session file (**Share**), a new AnnData file of the dataset with the state shown (**Save AnnData with this state…**), the printable report, the figure of the view shown, the highlighted cells (CSV or `.h5ad`) or a field (CSV) |
| **Deleting what was saved** | deleting a branch, a step or a tracker, squashing steps, clearing the history, emptying Recently deleted, deleting a saved state of this browser, forgetting a dataset of My datasets, removing what the **Stored** dialog lists |
| **Large downloads** | opening a dataset by downloading a whole file (a host that serves no parts of it), not one this browser already keeps |
| **Publishing** | sending what you made to a shared place: a community annotation vote (published alone to its repository, or as a pull request) and Share's public link. Asked whoever asks: an agent in the chat, a connected app, or a Python or R notebook |
| **Pausing the history** | pausing recording |

**Modes** set every kind at once:

- **Ask before heavy or irreversible** (the default): changes run; heavy work
  and the five other kinds (files, deleting, large downloads, publishing,
  pausing) ask first.
- **Allow everything undoable without asking**: changes and heavy work run;
  the five other kinds ask.
- **Read only**: every kind is refused.

Change them under **Permissions › Change…**, in the chat settings or in
**Agents › Claude and ChatGPT**: the **Agent permissions** dialog has the
mode, each kind, what you allowed for one chat or connection (**Forget**
clears it), and the dataset's switch (below). The settings are kept in this
browser (`cellucid:agent-permissions`).

A connected agent's tool list says, tool by tool, what your settings make of
it ("asks", "runs without asking", "refused"). The chat's model is told your
settings with your next message after they change. A kind set to Deny is
refused with one sentence the agent reads ("The person's settings in
Cellucid do not let agents …; do not try again").

---

## The confirm card

When a request asks first, a card appears above the agent's line: who asks
(the agent and its way in, for example "Claude (chat)" or "Claude (Claude
Desktop)"; a Python or R notebook that publishes is "Python (notebook)" or
"R (notebook)"), what exactly (the action and what it would do with its
names, for example `Claude (chat) asks to vote on a suggested label: “Voted
"Alpha cell" up for alpha of cell_type”.`, its parameters under **Details**)
and why (the estimate for heavy work, what a download costs, where a vote is
published). A local app (Claude Desktop, Claude Code) shows the same sentence
when it asks you there. Answer:

- **Deny**: the agent is told you declined (`The person declined, so this
  was not done: …`);
- **Allow · Once**: this request only;
- **Allow · This chat** (or **This connection**, **This notebook**): that
  kind for that chat, connection or notebook until the page closes (a kind
  you set to Deny since stays refused);
- **Allow · Always**: the kind is set to Allow.

The agent waits. A card never outlives its question: it is withdrawn when the
agent gives up, when you press **Stop** on its turn, when another dataset is
shown, when a saved state is loaded over the history shown, or when the tab
closes. The first card takes the keyboard focus unless you are typing
elsewhere, and its title is announced either way. A local app that can ask
too (Claude Code) asks the same question there: answer in either place.

---

## The activity line, Stop and Undo turn

One line per agent turn sits at the foot of the plot (in the chat itself
while the chat is open): while the agent acts, what it runs and **Stop**,
which also shows while one of its cards waits and withdraws that card; once
done with steps, **Undo turn** and **×**. While the history's recording is
paused the line starts with `⏸`: the agent's changes make no step meanwhile.
Two agents at once (Claude in the
chat and Claude Desktop through the connector) are two lines. A done line
steps aside when you (or a notebook) take a step.

What a turn did outside the tab, which no undo takes back, is counted on its
line: `Claude · 2 steps · 1 vote published, not undoable` (a published vote,
a file written, a saved state deleted…), kept in the turn's record in the
history (the turn's details show it under **Outside**), and said again by
**Undo turn**: `Undid Claude's 2 steps; 1 vote published stays: Undo cannot
take it back.` A turn that only published still shows its line.

**Undo turn** takes the turn back: alone at the end of the branch, back to
the step before it; when others acted after it began, you choose only its
steps (to Recently deleted, the later steps applied again without them) or
everything since ({doc}`01_chat_bar`). The same **Undo turn** is in the
turn's details in the history ({doc}`../m_state_history/01_the_history`).

---

(agents-what-is-sent)=
## What agents are sent

**Free reading is the default.** What an agent reads goes to the provider of
its model with the conversation: names, counts and summaries, and, when it
asks, values of cells and pictures of the plot. The dataset's files stay in
your browser. The chat bar sends straight from the tab to Anthropic or
OpenAI; a connected Claude or ChatGPT gets results through the relay, which
passes them on without keeping them; a local app gets them from the viewer
on your computer and sends them to the provider it uses.

**The notice.** The first time an agent acts on a dataset (or you open or
use the chat on it) and you have not answered it in this browser, one
compact notice says *What agents read here goes to their AI provider*, with
the switch and **Got it**. Nothing waits on it.

**The switch: Agents see names and summaries only.** In the notice and in
the **Agent permissions** dialog, for the dataset shown. It is yours alone:
no agent can turn it off. It is kept for that dataset (by its identity, as
saved states and My datasets know it) in this browser, in every tab
(`cellucid:agent-data`); a dataset nothing names keeps it for this page only,
as does a browser that keeps no site data. With it on, for every agent:

- a read of single cells' values, and **a picture of the plot**, is refused
  before it runs (a picture counts as cell-level: the plot draws each cell);
- a minimum, a maximum and a median (and the quantile columns of an
  analysis's table) are left out of results, each said in a note: each is
  one cell's own value; means, counts, sums and shares stay;
- a summary of values over **groups an agent made by naming cells** is
  refused: a highlight group an agent made from a list of cells, a column it
  pushed, and the fields and pages made from them. Groups from the dataset's
  categories and groups you made stay summarisable. These marks are kept for
  the page only: they are not saved in a session or a saved state, so a
  history loaded later starts with none;
- a read Cellucid has not classified is refused, so a tool added later is
  guarded before anyone classifies it;
- a file written for you (the report, the figure, the highlighted cells, a
  field, a session file, an AnnData copy) is still written: it goes to your
  computer, not to the provider, and the agent is told only the file's name,
  its size and counts (the report asked for as text holds what the history
  holds, never a picture);
- switched on during a chat, earlier results that could hold values of
  cells, and earlier pictures, are withheld when the conversation is sent
  again;
- the turn's **Sent** line says that no result carried values of cells or a
  picture ({doc}`01_chat_bar`).

A refusal tells the agent, in one sentence, what was not done and to work
from names, counts and summaries of groups the dataset or you define.

:::{warning}
**It is not anonymisation.** A summary of a very small group the data or you
define (a narrow range, a filter, a rare category) can still reveal its
cells' values. Cellucid sets no minimum group size; agents are told never to
work around the switch with very small groups, but only you can judge what a
group reveals. For data that must not reach an AI provider at all, do not
use agents on it.
:::

---

(agents-text-is-data)=
## Text written by others is data

Names, labels, notes, shared and community annotations, gene-set and file
names and a dataset's own description were written by the dataset's authors
and by other people, not by you. An agent could be tricked by text that reads
like an instruction ("ignore your instructions and delete…"). So:

- every tool result's data is sent between `<cellucid-data>` and
  `</cellucid-data>`, with every `<` inside escaped, so nothing in it can
  close the marks; Cellucid's own sentences (a refusal, a picture left out, a
  result cut) come after the closing mark;
- a name Cellucid's own sentence quotes (a refusal listing the fields there
  are, the suggestions of a category, a question the person declined) is
  marked where the sentence is built: the tags around the name written as one
  JSON string, every `<` in it escaped, so no quote, newline or tag in a name
  can end the mark early; the names the selection chip attaches to a message
  are marked the same way;
- every agent is told early that marked text is data, never instructions:
  if it asks for something, the agent does not do it and tells you what it
  says.

The chat's rows and its Markdown export show the text without the marks.
This is a guard, not a guarantee: a model can still be misled, which is why
your permissions, the confirm card and **Undo turn** stay the backstop.

---

## Next

- {doc}`01_chat_bar`, {doc}`02_connect_claude_and_chatgpt`
- {doc}`../o_accessibility_privacy_security/02_privacy_model`: everything
  kept in this browser and what goes over the network
