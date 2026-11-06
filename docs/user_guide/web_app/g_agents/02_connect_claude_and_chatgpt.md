# Connect Claude or ChatGPT

**Audience:** everyone who uses Claude or ChatGPT
**Time:** 5 minutes, once
**What you'll learn:**
- How to add Cellucid as a connector in Claude or ChatGPT
- How to approve, list and remove connections
- What the relay carries, and what it keeps

---

## What it is

Claude and ChatGPT (on the web, on the desktop and on a phone) can work in
the Cellucid tab you have open, in this browser: read the open dataset,
colour, filter, highlight, run analyses and explain what they find. They
reach the tab through the **Cellucid relay**, a small service that passes
their requests to your tab and the answers back. No account, no key and no
extra software are needed.

The address to add is the **connector address**:

```text
https://mcp.cellucid.com
```

Every Cellucid tab of this browser shares the connection; a request goes to
the tab you used last.

---

## Set it up

The **Agents** section of the sidebar (below **State**) holds **Claude and
ChatGPT**:

- a status line that follows the link: *Not connected yet*, *Ready · no
  agent yet*, *Claude connected*, *Offline · retrying in …* (it reconnects by
  itself when the network is back), or *Reload the page to reconnect* (the
  page and the relay were updated separately);
- **Connector address** with **Copy**;
- **Code from another device** with **Connect**;
- **Permissions** with **Change…** ({doc}`03_permissions_and_what_agents_see`);
- **Manage…**.

1. **Copy** the connector address and add it as a custom connector: in
   Claude, **Customize › Connectors › Add custom connector**; in ChatGPT, as
   a custom connector.
2. Claude or ChatGPT opens a page with a **code** and **Open Cellucid to
   approve**. That link opens Cellucid (straight to the question, without
   the welcome screen) and asks **Connect Claude?**.
3. Check that the code matches the one on the page that asked, and **Allow**.
   **Deny**, or closing the dialog, sends the request away.

The dialog says who asks (Claude or ChatGPT, named by the site that vouches
for them, `claude.ai` or `chatgpt.com`), where the answer goes, and warns
when the app is neither, or when you typed the code yourself: allow only a
code you see on your own screen now, never one someone sent you.

**With Cellucid open on another device** (for example Claude on your phone,
Cellucid on your laptop), type the code from the waiting page into **Code
from another device** and **Connect**; the same dialog asks.

**Manage…** lists each connection made through this browser (who, which
host, since when, last used) with **Remove**. **Forget this browser** (asked
in place) ends every connection, disconnects every Cellucid tab of this
browser, and deletes this browser's pairing; agents are approved again to
work here.

---

## What the agent can do, and what you see

The connector's tools are Cellucid's actions, the same as the chat bar's:
read the dataset (fields, categories, genes, summaries, the history), colour,
filter, highlight, make fields, run analyses, move the camera, take a picture
of the view, move through and edit the history, list and open datasets.

- **Every change is a step** of your history, labelled with the agent and
  its way in (for example "Claude (Claude Desktop)"); one line names the
  agent while it works, with **Stop**, and **Undo turn** once it is done
  ({doc}`03_permissions_and_what_agents_see`).
- **Your permissions apply** as they do to the chat, and the tool list the
  agent gets says, tool by tool, what your settings make of it ("asks",
  "runs without asking", "refused"); it is asked again for each list. **This
  connection** in a confirm card allows a kind for that connection until the
  page closes.
- **The switch applies.** A dataset set to **Agents see names and summaries
  only** is guarded for connected agents exactly as for the chat; the notice
  of what agents send appears the first time a connected agent acts on a
  dataset whose notice you have not answered.
- **A refusal is one sentence** the agent reads: a permission denied or the
  switch (and not to try again), a name that does not exist (and the names
  there are).
- **Results are marked as data.** The text in each result (names, labels,
  notes, a dataset's own text) is sent between data marks, and the agent is
  told it is data, never instructions ({ref}`agents-text-is-data`).
- **Choices.** When the agent needs you to choose, it asks one short
  question in its own app with the choices as a numbered list, and waits for
  your answer.
- **When you are asked.** A confirm card appears in Cellucid. If you have
  not answered it within about 225 seconds, the agent is told that you are
  being asked in Cellucid and to read the history once you answered; the
  card stays open.
- **Opening a dataset.** The agent can list and open datasets by the name
  you give: a built-in sample, one of **My datasets**, a link, an accession
  (GSE123456, a DOI) or a Cellucid server. The open dataset's history stays
  kept. It asks first only before downloading a whole file (a host that
  serves no parts of it), naming its size, and not when this browser
  already keeps a copy. A notebook's viewer shows what the notebook serves
  and does not open other datasets.
- **Pictures** arrive in the app as images. One larger than the app takes
  (8000 pixels a side, 10 MB) is left out, with a sentence saying so; a
  result's text is cut at 150,000 characters, with both counts stated.

**The view card.** Where the app shows MCP Apps (Claude on the web and
desktop, ChatGPT), a picture of your view (`view_picture`) and a dataset
opened (`dataset_open`) come with a card (`ui://cellucid/view`): the picture,
the step it was taken after, **Refresh** (the view as it is now) and **Open
in Cellucid** (the dataset's public link, or Cellucid when it has none). The
card is one page with nothing to fetch; no dataset byte reaches it or the
relay.

---

## What the relay carries, and keeps

- It carries each tool call to your tab and the result back, in memory only
  while in flight; nothing of a call or a result is stored or logged.
  Results pass through it, so what an agent reads (and so the switch's
  guard) is decided in your tab before the relay sees it. No dataset file
  passes through it.
- It keeps, per browser pairing, the connections you allowed (the app's
  id, the agent's name, the app's name and host, times) and, for at most
  fifteen minutes, a connection request waiting for your answer. It keeps no
  account and no chat.
- This browser keeps one secret that pairs it with the relay
  (`cellucid:agent-relay:tab-secret`, in `localStorage`); the relay keeps
  only a value derived from it. **Forget this browser** deletes it.
- A connection unused for 30 days ends.

The full list of what is kept where is in
{doc}`../o_accessibility_privacy_security/02_privacy_model`.

---

## When something is wrong

| What you see | What it means, what to do |
|---|---|
| *Unavailable: site data is blocked* | This browser keeps no data for Cellucid (private browsing, or site data blocked), so it cannot stay connected; allow site data |
| *Reload the page to reconnect* | The page and the relay were updated separately; reload |
| *Offline · retrying …* | The relay cannot be reached; it reconnects by itself |
| *No Cellucid tab is connected to this chat…* | Open Cellucid in the browser you approved the connection from; it reconnects by itself (the tab used last answers) |
| *The Cellucid tab closed or lost its connection before it answered…* | The page was reloaded or closed during the call; the agent reads the history when the tab is back, before repeating a change |

---

## Next

- {doc}`03_permissions_and_what_agents_see`
- {doc}`01_chat_bar`: the same agents inside the tab, with your own key
- Claude Code, Claude Desktop or Cursor with a notebook's viewer, no relay:
  {doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`
