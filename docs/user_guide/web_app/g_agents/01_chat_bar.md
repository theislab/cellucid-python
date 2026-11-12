# The chat bar (Ask an agent)

**Audience:** everyone
**Time:** 10 minutes
**What you'll learn:**
- How to start a chat with Claude or OpenAI's models, with your own key
- What a turn is, and how to stop it or take it back
- What the chat keeps, where, and what it sends

---

## Open the chat

**Ask an agent** sits at the foot of the plot, above the selection bar. Click
it, or press `/` anywhere on the page outside a text field (in the history
list, `/` searches the history instead).

On a wide window the chat opens as a **side panel** on the right, the full
height of the window, and the plot moves over to make room for it: the
agent never covers the cells it is changing. Drag the panel's left edge to
make it wider or narrower (or focus the edge and use ← and →; a double click
puts the width back); the width is remembered in this browser. **Collapse
the chat** folds the panel away and gives the plot the whole window again.
On a narrow window (a phone, or a laptop with the sidebar open that leaves
too little room beside it), the chat opens over the foot of the plot
instead. While a camera path plays there, the chat folds to its **Ask an
agent** button so the film is not played under it, and comes back when the
path stops or pauses; it stays open while you are writing a message or a
question of the agent waits for your answer, and clicking the button during
the play opens it for the rest of that play. The notice of what agents send
and any question an agent asks you appear under the chat, next to the
message field.

The header holds **New chat**,
**Export the chat (Markdown)**, **Chat settings**, the info button and
**Collapse the chat**; a short status line beside the title says what the
agent is doing or what to do next.

An empty chat offers a few first questions made from the dataset shown. Type
your own in **Ask about this dataset…** and press **Send**.

---

## Settings

**Chat settings** (the sliders icon) opens:

- **Provider**: Claude (Anthropic) or OpenAI's models (the models behind
  ChatGPT).
- **API key**: your own key from that provider (**Get a key** opens the page
  where you make one). Its use is billed to your account. By default it is
  kept **for this tab**: a reload keeps it, closing the tab forgets it.
  **Remember on this browser** keeps it for every tab of this browser profile
  until **Forget**. Either way it stays in this browser and is sent only to
  its provider (`api.anthropic.com` or `api.openai.com`), never to Cellucid,
  a session file, the history or the address. Password managers are told to
  ignore the field. A key the provider refuses is said here.
- **Model**: the models your key may use, as the provider lists them.
  Nothing is hard-coded.
- **Effort**: how much the model thinks before it answers or acts, as the
  provider offers it for that model; **Model's default** leaves it to the
  model. Claude's models list the levels they take; OpenAI does not, so a
  level a model does not take is refused in one sentence.
- **Permissions**: one line with **Change…**, the same settings as for every
  agent ({doc}`03_permissions_and_what_agents_see`).
- **Your instructions**: a few words every agent of this chat is told (for
  example "I work on mouse pancreas; explain briefly").
- **Prompt cache** (Claude only): how long Claude keeps the cached part of
  each call (the instructions, the tools and the conversation so far), so a
  call costs a small part of a fresh one. **Kept 5 minutes** (the default) is
  cheapest while you keep chatting; **Kept 1 hour** costs more when a chat
  starts and spares the full price again after a pause of more than five
  minutes (measured: a Sonnet question after 6 minutes cost $0.19 kept 5
  minutes, $0.045 kept an hour). Claude reuses a copy it already holds, so an
  hour chosen within five minutes of a chat kept for five minutes starts only
  after that copy lapses: choose it before you start. OpenAI caches by itself,
  so its chats do not show it.
- **Token budget per chat**: none unless you set one. The chat makes no call
  to the model that would pass it: each call reads the whole prompt again (the
  instructions, the tools and the conversation so far), so before a call (and
  before trying a failed one again) the chat adds what it will read — the last
  call's prompt and what was added since; before a chat's first call, an
  estimate of all of it — to the tokens used, and stops there when that is
  more than the budget, saying how much the call would have read; raise it to
  go on. The line below shows the tokens used so far.

The provider, the model of each provider, the effort of each model, your
instructions, the budget, the prompt cache and whether the chat was open are
remembered in this browser (`cellucid:agent-chat`); the key is not kept there.

---

## Turns

Each message you send starts one **turn**. The agent answers, and may read
the dataset and act in it through Cellucid's tools, one call after another.
In the conversation:

- each tool call is one compact row that says what it was about: what it
  read ("Read the categories of cell_type", "Summarized SOX9 by cell_type",
  "Looked up 12 gene names: SOX9, SFTPC, ETV5 +9"), or the step it took
  ("Ran Differential expression"), with **Show** to find its step in the
  history; a call that was refused or stopped says so, with what it was
  asked for ("Color the cells by a field or gene: CD3E · refused"). A row
  opens to what was asked and answered;
- a model that writes a short note between its calls (what it found, what
  it does next: Claude Opus 5.5 and Sonnet 5.5 do) shows each note as a quiet
  line before the call it introduces, so a long turn says what it is doing
  instead of only "Thinking…"; the notes are not the reply (Copy copies the
  reply), and they cost nothing more: a model writes them whether they are
  shown or not. The chat's instructions ask the agent to say in a line what
  it is about to do before it starts, to give brief updates while it works
  and to close with a recap that stands on its own, so a long turn is never
  silent for long;
- under each turn: its steps, **Sent** (below), its tokens, and **Undo
  turn**;
- your messages have **Edit** (the message back in the field) and replies
  **Copy**; a failed turn's note has **Try again**, as does the note of a
  turn that ended with no words for you (a model may end after its
  reasoning alone: Cellucid says so rather than leave the turn silent).

**A long turn folds.** A turn with more than three calls (a research
question can take forty) gets one summary row under the agent's name:
*Worked through 42 actions · 7 steps · 3 analyses*. Once the turn has
ended, its answer is what you see; its work (the agent's words and notes
between calls, and every call) is under the summary, which opens and closes with a
click, Enter or Space. What went wrong stays in view: a call that was
refused, stopped or interrupted (unless the agent made the same call again
later in the turn and it went through), and a question the agent asked. While the
turn runs, the summary counts the calls so far (*12 actions so far · 2
steps*) and the latest three show, with the words and notes that made them, so you
can follow it without a growing wall. Nothing is taken away: opened, the
turn is as it was made, and **Show**, **Sent**, **Undo turn**, **Copy** and
**Export** work as before; the exported chat lists every call, named as its
row is. A screen reader hears the agent's words and each turn's summary
once, not every call.

A turn that made no step (it only read, or asked for what was so already)
has no **Undo turn**: there is nothing to undo. When the agent asks for a
change that is so already (coloring by the field shown), it is told that
nothing changed and says so, rather than that it made the change.

With your next message the agent is told what changed in the history since
its last turn: which of its steps you undid or removed, and which steps you
(or a notebook, or another agent) took, by their labels, a step you took
while its turn ran included. When another history is shown than the one of
its last turn (you loaded a session file beside yours, opened a saved state,
started a tracker or opened another dataset), it is told which one, where it
came from (the session file's name), how many steps are in effect there and
the latest few, and, when the chat shows a new conversation for it, that its
earlier messages stay with the history they were in. What the agent last saw
is kept with the conversation, so a reload, a saved state or a shared
session keeps it. It is never sent data this way, only the steps' and
histories' names.

Every change the agent makes is a step of your history, carrying the turn,
the agent's name and the model ({doc}`../m_state_history/01_the_history`). A
turn in the history opens its message in the chat, and the other way round.

**Stop** (beside **Send** while a turn runs) ends the turn; a step it was
taking is rolled back, except a playing camera path, which stops where the
film is (as a pause does) and keeps that step. When the agent asks you to confirm something (a
card, {doc}`03_permissions_and_what_agents_see`), the chat says it waits for
you, and **Stop** also withdraws the card.

**Undo turn** takes the whole turn back. When nobody acted after the turn
began, it goes back to the step before it, kept for redo. When you or
another agent acted since, it asks:

- **Only the agent's steps**: they go to Recently deleted, where they can be
  restored, and the later steps that needed them are applied again without
  them; a step that cannot be stays, marked as no longer applying.
- **All steps since**: back to before the turn, everything since kept for
  redo.

The provider waits out a busy or rate-limited moment and a dropped
connection by itself, asks once more for an answer it could not read, and a
conversation that grows too long for the model is summarized once (you still
see every message; the history keeps every step).

### Asking for a flythrough

The agent can direct a camera path ({doc}`../c_core_interactions/08_camera_path`)
as a documentary would: "make a cinematic flythrough that visits each major
cell type, then play it", "a 20-second loop for a talk", "start wide, dive
into the beta cells, end on the whole pancreas". It reads the data for the
story (which populations, in what order: a differentiation from progenitors
to fates), asks Cellucid where each population sits and for a camera that
frames it from a direction it chooses (**Frame cells with the camera**,
`camera_framing`: nothing moves), saves those cameras as keyframes with a
narrator's labels (names for the shots, which a playing path shows as
captions over the plot, like a documentary's titles), holds on each subject (a keyframe's hold), times each
move, reads the path back (its length) and plays it (**Play the camera
path**, one step: **Undo turn** puts the camera back). A looping path plays
until **Stop** in the transport bar or the agent's **Stop the camera path**;
**Stop** on the agent's turn stops a path it is playing, a loop included,
where the camera is, and the play's step keeps that camera.
When your system asks for reduced motion, the path goes to its last keyframe
at once instead of playing; the agent is told so (and reads the setting in
`views_state`), so it does not describe a film that ran.

### Asking research questions

The agent can work on the data as a researcher, not only click for you:
"how does the epithelium mature across gestation?", "when do the XX and XY
supporting lineages diverge?", "what is the most interesting biology here?",
then as a colleague: "why do you think so?", "what would rule that out?",
"check it with the donor as the unit", "write it up as a results paragraph".
It combines the analyses (marker genes, differential expression, composition,
gene signatures, dot plots, a field per category, quality fields) and reads
with each result how it was computed: the groups and their cells, the test
and correction, the expression matrix and the values tested (raw counts, or
normalized and log1p as **Data mapping › Expression values** sets them), the samples per side with **Replicates** (pseudobulk) and
the tests between pages of a gene signature. It checks a marker panel against
the dataset's genes in one read, and says what is solid and what is a
hypothesis. Cells of one donor are not independent: ask it to test with the
donor as the replicate (differential expression's **Replicates**,
Composition's **Samples**) before you rely on a p-value.

---

## Choices and the selection

Above the message field:

- **Choices.** When the agent needs you to pick between a few options (which
  field, which groups, which plan), it asks one short question with two to
  five choices, shown as buttons. A click (or `Enter`) sends that choice as
  your next message and keeps what you had half typed; typing another answer
  works too.
- **The selection chip.** When you highlight cells, a chip reads
  *Selection: 1,234 cells ×*. The next message you send carries a reference
  to them: the highlight page, its groups and their counts, never the cells
  themselves; the agent reads the cells through its tools, as the dataset's
  switch allows. A selection the conversation already knows is not sent
  again: once it went with a message, or once it was highlighted when the
  agent's turn ended (cells the agent highlighted itself), the chip goes until
  the highlighted cells change. **×** leaves the selection out until it
  changes; the chip goes when nothing is highlighted, and comes back in a new
  chat.

---

## What is sent (Sent)

Requests go from this tab straight to the provider you chose, with your key.
Cellucid has no server in between. **Sent**, under each turn, says what left
this browser for that turn: your message and the tool results (how many, how
many characters, how many pictures). The dataset's files stay here, but a
result can carry part of the data: values of cells, a picture of the plot.
When the dataset is set to **Agents see names and summaries only** and every
result of the turn was read under it, Sent says that no result carried values
of cells or a picture.

The first time you open the chat or send a message on a dataset whose notice
you have not answered, the notice of what agents send appears above the chat
with the switch ({ref}`agents-what-is-sent`). Nothing waits on it.

---

## The chat and the history

The conversation is kept with the history of the dataset shown, after every
message (without the providers' reasoning or pictures). So:

- it comes back with **Continue** and with a saved state; another dataset or
  tracker shows its own latest conversation (or a new one);
- a turn still running when you open another dataset stops and is kept with
  the dataset it began on, unless its own agent opened that dataset: then it
  goes on there, and the conversation is kept with both histories;
- after a reload, a call that was still running shows as *interrupted*, with
  the step the history kept for it, or nothing applied, and the turn says
  *The page closed during this turn, so it stopped here.* with **Try again**
  (also when the reload came between two calls, while the agent was to
  answer);
- **State › Session › Share…** includes the conversations only when you tick
  **Conversations** (off by default; never with **Current state**), and
  **My datasets › Export…** carries the saved states you tick whole,
  conversations included ({doc}`../m_state_history/index`). A conversation
  can be private: keep such a file to yourself;
- a conversation opened from a session file or a saved state is
  *received*: its turns are the sender's, already in the history it came
  with, so none of their notes offers **Try again** and a question its agent
  asked shows no choices (write your answer instead). The turns you send
  after opening it are yours, with **Try again** as usual. The agent you
  continue with is told those turns were another person's and their agent's
  work, and a received turn's **Sent** says that nothing of it left your
  browser (it names what the sender's browser sent).

**Export** writes the conversation as Markdown, with raw HTML escaped and
script links left as text.

---

## When something fails

The chat bar and the connector are independent: a refused key, a provider
that is down or the relay being unreachable breaks only that way in, and the
chat says in one line what to do. A connected Claude or ChatGPT keeps working
without a key, and the chat works with no connector.

| What you see | What to do |
|---|---|
| The key is refused | Check it at the provider (**Get a key**), paste it again |
| *This chat used N tokens, its budget is M* | Raise **Token budget per chat** |
| A level of **Effort** is refused | Choose another, or **Model's default** |
| A row reads *interrupted* | The page reloaded during the call; its row says what was applied |

---

## Next

- {doc}`03_permissions_and_what_agents_see`: permissions, the confirm card,
  the switch.
- {doc}`02_connect_claude_and_chatgpt`: the same agents in their own apps.
