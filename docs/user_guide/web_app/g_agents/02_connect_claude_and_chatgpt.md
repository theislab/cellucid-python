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
color, filter, highlight, run analyses and explain what they find. They
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
  itself when the network is back), or *Reload the page to reconnect.* (the
  page and the relay were updated separately);
- **Connector address** with **Copy**;
- **Connect code** with **Connect**;
- **Permissions** with **Change…** ({doc}`03_permissions_and_what_agents_see`);
- **Manage…**.

1. **Copy** the connector address and add it as a custom connector: in
   Claude, **Customize › Connectors › Add custom connector**; in ChatGPT, as
   a custom connector.
2. Claude or ChatGPT opens a page with a **code** and **Open Cellucid to
   approve**, in your browser. That link opens Cellucid in the same browser
   (straight to the question, without the welcome screen) and asks
   **Connect Claude?**.
3. Check that the code matches the one on the page that asked, and **Allow**.
   **Deny**, or closing the dialog, sends the request away.

**Only the browser that shows the waiting page can approve it.** The page
leaves a short-lived cookie for the relay in that browser, and Cellucid's
question to the relay carries it; a code or link that someone else started
in their own browser and sent you is refused in your tab (*This connect code
was started in another browser…*), so nobody can connect *their* Claude to
your Cellucid by sending you a link. A Cellucid tab already open in the same
browser can take the code too: type it into **Connect code** and
**Connect**. Another browser or another device cannot: open Cellucid where
the waiting page is (on a phone, in the browser that showed the page), or
start connecting again from Claude or ChatGPT in the browser you use for
Cellucid.

The dialog says who asks and where the answer goes. Claude and ChatGPT are
named by the site that vouches for them (`claude.ai`, `chatgpt.com`); any
other app is named by its own site (or as *an app on this computer*), with
the name it gives itself only in quotes: an app at `evil.example` that calls
itself "Claude" asks as *evil.example, which calls itself "Claude"*, and is
*evil.example (calls itself "Claude")* in **Manage…**, its activity line, its
confirm cards and the history's steps, so two such apps are told apart.
The dialog warns when the app is neither Claude nor ChatGPT; when you came by
a link (allow only right after you asked to connect, never for a link
someone sent you); and when you typed the code yourself (allow only a code
you see on your own screen now). A Cellucid page shown inside another site's
frame never asks: the link's code is dropped.

**Manage…** lists each connection made through this browser (who, which
host, since when, last used) with **Remove**. **Forget this browser** (asked
in place) ends every connection, disconnects every Cellucid tab of this
browser, and deletes this browser's pairing; agents are approved again to
work here. Removing a connection (or forgetting the browser) while its agent
is working stops the call it is waiting on, in whichever tab runs it: the
agent reads that the connection was removed during the call, and nothing it
asked for goes on running or reaches it.

---

## What the agent can do, and what you see

The connector's tools are Cellucid's actions, the same as the chat bar's:
read the dataset (fields, categories, genes, summaries, the history), color,
filter, highlight, make fields, run analyses, move the camera, take a picture
of the view, move through and edit the history, list and open datasets.

- **Every change is a step** of your history, labeled with the agent and
  its way in (for example "Claude (Claude Desktop)"); one line names the
  agent while it works, with **Stop**, and **Undo turn** once it is done
  ({doc}`03_permissions_and_what_agents_see`).
- **Your permissions apply** as they do to the chat. The agent is told your
  settings (what runs without asking, what is asked first, what is refused,
  the dataset's switch) after its first tool result, and again after a
  result once you changed them; its tool list stays the same, so changing a
  setting never makes the app rebuild its cached prompt. **This
  connection** in a confirm card allows a kind for that connection until the
  page closes.
- **What you did meanwhile is told.** With the first result after you acted
  (you undid its step, colored by another field, loaded a session file), the
  agent is told what changed in the history since its last call, as the chat
  bar's agent is with your next message ({doc}`01_chat_bar`): steps by their
  labels, another history by its name. The tab keeps what each connection
  last saw, so a reload keeps it too.
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
  kept. A link, accession or server you have not seen in this tab (typed or
  pasted into the address bar or the chat, open now, or in **My datasets**)
  is asked about first, before anything of it is requested: an agent misled
  by text it read could otherwise make your tab send what it read to any
  host ({ref}`agents-unseen-links`). It also asks before downloading a whole
  file (a host that serves no parts of it), naming its size, and not when
  this browser already keeps a copy. A notebook's viewer shows what the notebook serves
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
- It keeps, per browser pairing, the connections you allowed (a digest of
  the app's id, the agent's name, the app's name and host, times) and, for at most
  fifteen minutes, a connection request waiting for your answer (with
  digests of the waiting page's secret and of its browser's cookie). It
  keeps no account and no chat.
- The waiting page leaves one cookie in the browser that shows it
  (`__Host-cellucid-approver`, for `mcp.cellucid.com` only, unreadable by
  scripts, sent only from Cellucid's own site), for ten minutes: it is what
  lets only that browser approve.
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
| *Unavailable: site data is blocked.* | This browser keeps no data for Cellucid (private browsing, or site data blocked), so it cannot stay connected; allow site data |
| *Reload the page to reconnect.* | The page and the relay were updated separately; reload |
| *This connect code was started in another browser…* | The waiting page is in another browser (or someone sent you the link): approve in the browser that shows the page, or start connecting again from Claude or ChatGPT in this browser |
| *Offline · retrying …* | The relay cannot be reached; it reconnects by itself |
| *No Cellucid tab is connected to this chat…* | Open Cellucid in the browser you approved the connection from; it reconnects by itself (the tab used last answers). While a tab is being reloaded, the agent's next call waits for it to be back instead (up to 20 seconds) |
| *The Cellucid tab closed, was reloaded or lost its connection before it answered…* | The page was reloaded or closed during the call, so whether it was done is not known; the agent reads the history (its call waits for a reloaded tab) before repeating a change |

---

## Next

- {doc}`03_permissions_and_what_agents_see`
- {doc}`01_chat_bar`: the same agents inside the tab, with your own key
- Claude Code, Claude Desktop or Cursor with a notebook's viewer, no relay:
  {doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`
