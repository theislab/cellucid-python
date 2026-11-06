# Agents

An AI agent can work in Cellucid with you: read the dataset shown, colour,
filter, highlight, make annotation fields, run analyses, move through the
history and explain what it finds. It acts only through the same actions the
interface uses, so **every change it makes is a step of your history**,
labelled with its name, and one click takes its whole turn back.

It can also do what lies outside the history: play and stop the guided tour,
make the printable report (as a file or as text for itself), export the
figure, the highlighted cells (CSV or `.h5ad`) or a field (CSV), list or
remove what this browser keeps (the **Stored** dialog) and forget a dataset
of **My datasets**, give Share's public link of the dataset, and read and
vote on community annotation suggestions. None of these is a step, and
Undo does not take it back, so writing a file, removing what is kept and
publishing ask you first by default
({doc}`03_permissions_and_what_agents_see`). Voting publishes to the
repository, so it needs you to have connected the repository and signed in
to GitHub first: both need your click.

There are three ways in. All three use the same tools, the same permissions
and the same rules for what is sent:

| Way in | What it is | Read |
|---|---|---|
| **The chat bar** (**Ask an agent**, at the foot of the plot) | Cellucid's own chat with Claude or OpenAI's models, with your own API key, inside the tab | {doc}`01_chat_bar` |
| **The connector** (**Agents › Claude and ChatGPT**) | Claude or ChatGPT, in their own apps, reaching your open Cellucid tab through the Cellucid relay | {doc}`02_connect_claude_and_chatgpt` |
| **A local agent app** (`cellucid mcp`, `cellucid_mcp()`) | Claude Code, Claude Desktop or Cursor on your computer, driving a viewer the Python or R package serves | {doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`, {doc}`../../r_package/d_viewing_loading/09_local_agents_mcp` |

Whichever way an agent came in, what you see and decide is the same
({doc}`03_permissions_and_what_agents_see`):

- **Your permissions decide** what it may do without asking: changes, heavy
  work and what the history cannot take back are each Ask, Allow or Deny.
- **One line says which agent is acting**, with **Stop**, and once it is done
  **Undo turn**.
- **What agents read goes to their AI provider.** Once per dataset Cellucid
  says so, with the switch **Agents see names and summaries only**.
- **Text written by others is data.** Names, labels, notes and a dataset's
  own text reach an agent marked as data, never as instructions.

:::{important}
The dataset's files never leave your browser through an agent. What an agent
reads (names, counts, summaries and, when it asks, values of cells or a
picture of the plot) goes to the provider of its model with the
conversation. See {doc}`03_permissions_and_what_agents_see` and
{doc}`../o_accessibility_privacy_security/02_privacy_model`.
:::

---

## Fast path

| You want to… | Do this |
|---|---|
| Ask a question about the dataset shown | **Ask an agent** (or press `/`), add your API key once in the settings, type |
| Have Claude or ChatGPT work in Cellucid from their own app | **Agents › Claude and ChatGPT › Copy** the connector address and add it in the app as a custom connector |
| Use Claude Code or Cursor on a notebook's viewer | `claude mcp add cellucid -- cellucid mcp` ({doc}`../../python_package/d_viewing_apis/18_local_agents_mcp`) |
| Take back what an agent did | **Undo turn**, under its turn in the chat, on its line, or in the history |
| Stop an agent now | **Stop** on its line (or in the chat) |
| Keep cell values and pictures from AI providers | Switch on **Agents see names and summaries only** (the notice, or **Permissions › Change…**) |
| Let it work without asking each time | **Permissions › Change…**: a mode, or Allow per kind |

```{toctree}
:maxdepth: 1
:hidden:

01_chat_bar
02_connect_claude_and_chatgpt
03_permissions_and_what_agents_see
```
