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
list, `/` searches the history instead). The header holds **New chat**,
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
- **Token budget per chat**: none unless you set one. Once a chat has used
  it, the chat stops before the next call to the model (and before trying a
  failed call again) and says so; raise it to go on. The line below shows the
  tokens used so far.

The provider, the model of each provider, the effort of each model, your
instructions, the budget and whether the chat was open are remembered in this
browser (`cellucid:agent-chat`); the key is not kept there.

---

## Turns

Each message you send starts one **turn**. The agent answers, and may read
the dataset and act in it through Cellucid's tools, one call after another.
In the conversation:

- each tool call is one compact row (what it read, or the step it took), with
  **Show** to find its step in the history; a call that was refused or
  stopped says so;
- under each turn: its steps, **Sent** (below), its tokens, and **Undo
  turn**;
- your messages have **Edit** (the message back in the field) and replies
  **Copy**; a failed turn's note has **Try again**.

Every change the agent makes is a step of your history, carrying the turn,
the agent's name and the model ({doc}`../m_state_history/01_the_history`). A
turn in the history opens its message in the chat, and the other way round.

**Stop** (beside **Send** while a turn runs) ends the turn; a step it was
taking is rolled back. When the agent asks you to confirm something (a
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
conversation that grows too long for the model is summarised once (you still
see every message; the history keeps every step).

---

## Choices and the selection

Above the message field:

- **Choices.** When the agent needs you to pick between a few options (which
  field, which groups, which plan), it asks one short question with two to
  five choices, shown as buttons. A click (or `Enter`) sends that choice as
  your next message and keeps what you had half typed; typing another answer
  works too.
- **The selection chip.** While cells are highlighted, a chip reads
  *Selection: 1,234 cells ×*. Each message you send then carries a reference
  to them: the highlight page, its groups and their counts, never the cells
  themselves; the agent reads the cells through its tools, as the dataset's
  switch allows. **×** leaves the selection out until it changes; the chip
  goes when nothing is highlighted.

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
  the step the history kept for it, or nothing applied;
- **State › Session › Share…** includes the conversations only when you tick
  **Conversations** (off by default; never with **Current state**), and
  **My datasets › Export…** carries the saved states you tick whole,
  conversations included ({doc}`../m_state_history/index`). A conversation
  can be private: keep such a file to yourself.

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
