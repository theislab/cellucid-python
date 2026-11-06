"""
``cellucid mcp``: a local MCP server that lets an agent app on this computer
-- Claude Code, Claude Desktop, Cursor -- work in a Cellucid viewer.

It speaks MCP over stdio (the official ``mcp`` SDK: ``pip install
'cellucid[mcp]'``) and drives a viewer the cellucid Python or R package
serves on this computer, with no relay and no deployed site in between:

- ``cellucid mcp --open data.h5ad`` serves the file and opens its viewer in
  the browser itself, and stops it when the app ends the server;
- ``cellucid mcp`` drives a viewer opened elsewhere -- in a notebook
  (``cellucid.show_anndata(adata)``), a script (``cellucid.AnnDataViewer``)
  or R (``cellucid::cellucid_view(obj)``) -- found among those published on
  this computer (:mod:`cellucid._agent_link`); ``--viewer`` names one.

The tools, their descriptions and the instructions are the web app's (the
web connector's tools, generated into ``mcp_tools.json``); the viewer page
runs each call as the app's agent ("Claude Code (local)") under the
person's permission settings: a change runs or waits for the person's
answer in the viewer's confirm card -- asked in the app too, as an MCP
elicitation, when the app can ask -- or is refused with one sentence. Every
step is in the person's history, and the person can Stop the agent's turn or
undo it whole.
"""

from __future__ import annotations

import asyncio
import contextlib
import logging
import sys
import threading
import webbrowser
from pathlib import Path
from typing import TYPE_CHECKING, Any

from ._agent_link import STOPPED, CallListener, ViewerLink, ViewerUnavailableError
from ._mcp_core import SERVER_NAME, SERVER_TITLE, LocalAgentSession, mcp_tools_document

if TYPE_CHECKING:
    from .jupyter import AnnDataViewer

logger = logging.getLogger("cellucid.mcp")

#: What the person is asked in the agent's app: allow (accept) or deny (decline).
_ELICIT_SCHEMA: dict[str, Any] = {"type": "object", "properties": {}}

SDK_MISSING = (
    "cellucid mcp needs the MCP SDK, which this environment lacks. Install it with: "
    "pip install 'cellucid[mcp]'"
)


def _import_sdk() -> tuple[Any, Any, Any]:
    try:
        from mcp.server.lowlevel.server import NotificationOptions, Server
        from mcp.server.stdio import stdio_server
    except ImportError as error:
        raise RuntimeError(SDK_MISSING) from error
    return Server, NotificationOptions, stdio_server


def _client_supports_elicitation(ctx: Any) -> bool:
    """Whether the app takes a form elicitation: it declared ``elicitation``,
    and not for URLs alone (an empty declaration means forms, as before
    modes existed)."""
    with contextlib.suppress(Exception):
        capabilities = ctx.session.client_capabilities
        elicitation = None if capabilities is None else getattr(capabilities, "elicitation", None)
        if elicitation is None:
            return False
        form = getattr(elicitation, "form", None)
        url = getattr(elicitation, "url", None)
        return form is not None or url is None
    return False


def _adopt_client(session: LocalAgentSession, ctx: Any) -> None:
    with contextlib.suppress(Exception):
        params = ctx.session.client_params
        info = None if params is None else params.client_info
        if info is not None:
            session.adopt_client(info.name, getattr(info, "title", None))


async def _ask_in_app(ctx: Any, session: LocalAgentSession, link: Any, request_id: str, message: str) -> None:
    """The person's question, asked in the agent's app too (an MCP form
    elicitation): Accept allows, Decline denies; dismissing it leaves the
    viewer's card to answer. The model never answers it: the app shows it."""
    logger.debug("call %s: asking in the app", request_id)
    try:
        answer = await ctx.session.elicit_form(
            f"{message} (You can also answer in the Cellucid viewer.)",
            _ELICIT_SCHEMA,
            related_request_id=ctx.request_id,
        )
    except asyncio.CancelledError:
        raise
    except Exception as error:  # an app that cannot ask after all: the card asks alone
        logger.debug("The app could not ask: %s", error)
        return
    action = getattr(answer, "action", None)
    logger.debug("call %s: the app answered %s", request_id, action)
    if action in {"accept", "decline"}:
        await asyncio.to_thread(link.answer, request_id, action == "accept")


async def _progress(ctx: Any, step: int, message: str) -> None:
    """Tell an app that asked for progress what the call waits on (``step``
    grows with every notification of the call, as progress must)."""
    meta = ctx.meta if isinstance(ctx.meta, dict) else {}
    token = meta.get("progress_token")
    if token is None:
        return
    with contextlib.suppress(Exception):
        await ctx.session.send_progress_notification(token, step, None, message, related_request_id=ctx.request_id)


class _ViewerCall:
    """One call sent to the viewer, stopped there whenever the app gives up on
    it -- also while it is still being sent.

    The send runs in a worker thread, and a cancelled task stops waiting for
    it at once while the thread goes on: the page may already run the call
    when the app's ``notifications/cancelled`` arrives, before the link that
    sent it is back. Whichever of the two ends second -- the send, or the
    app giving up -- sends ``agent_cancel``, so a cancel is never lost to how
    the two interleave (``tests/test_mcp.py`` holds a send until the server
    took the cancel, which lost it before).
    """

    def __init__(self, session: LocalAgentSession, request_id: str) -> None:
        self._session = session
        self.request_id = request_id
        self._lock = threading.Lock()
        self._link: ViewerLink | None = None
        self._abandoned = False

    def send(self, tool: str, args: dict[str, Any], listener: CallListener) -> ViewerLink:
        """Send the call (in a worker thread); stop it at once when the app gave up meanwhile."""
        link = self._session.start_call(self.request_id, tool, args, listener)
        with self._lock:
            self._link = link
            abandoned = self._abandoned
        if abandoned:
            self._cancel(link)
        return link

    def abandon(self) -> None:
        """The app gave up on the call: stop it in the viewer once it is sent
        (a change rolled back, a question withdrawn). From a thread of its
        own -- a cancelled task may not await anything more -- that the
        process waits for when it ends, so a cancel sent as the app quits
        still reaches the viewer."""
        with self._lock:
            self._abandoned = True
            link = self._link
        if link is not None:
            threading.Thread(target=self._cancel, args=(link,), name="cellucid-mcp-cancel").start()

    def _cancel(self, link: ViewerLink) -> None:
        try:
            link.cancel(self.request_id)
        except Exception as error:  # the viewer went meanwhile: nothing runs there to stop
            logger.debug("call %s: not stopped in the viewer: %s", self.request_id, error)


def _build_server(session: LocalAgentSession, holder: dict[str, Any]) -> tuple[Any, Any]:
    Server, NotificationOptions, _stdio = _import_sdk()
    from . import __version__

    async def list_tools(ctx: Any, params: Any) -> dict[str, Any]:
        _adopt_client(session, ctx)
        holder.setdefault("session", ctx.session)
        return {"tools": await asyncio.to_thread(session.tools)}

    async def call_tool(ctx: Any, params: Any) -> dict[str, Any]:
        _adopt_client(session, ctx)
        holder.setdefault("session", ctx.session)
        name = params.name
        args = dict(params.arguments or {})
        if name == mcp_tools_document()["statusTool"]["name"]:
            return await asyncio.to_thread(session.status_result)
        loop = asyncio.get_running_loop()
        events: asyncio.Queue[tuple[str, dict[str, Any]]] = asyncio.Queue()

        def listener(kind: str, data: dict[str, Any]) -> None:
            loop.call_soon_threadsafe(events.put_nowait, (kind, data))

        call = _ViewerCall(session, session.next_request_id())
        asking: asyncio.Task[None] | None = None
        progress = 0
        try:
            try:
                link = await asyncio.to_thread(call.send, name, args, listener)
            except (ViewerUnavailableError, TypeError, ValueError) as error:
                # No viewer can take it, or the call is not one a viewer reads
                # (a tool name no list gives): one sentence, as a tool error.
                return session.error_result(str(error))
            while True:
                kind, data = await events.get()
                logger.debug("call %s (%s): %s", call.request_id, name, kind)
                if kind == "asking":
                    progress += 1
                    await _progress(ctx, progress, data["message"])
                    if asking is None and _client_supports_elicitation(ctx):
                        asking = asyncio.create_task(_ask_in_app(ctx, session, link, call.request_id, data["message"]))
                elif kind == "answered":
                    if asking is not None:
                        asking.cancel()
                        asking = None
                elif kind == "result":
                    return dict(data["result"])
                elif kind == "failed":
                    # The page refused it, or went before it answered: the sentence says which.
                    return session.error_result(str(data["reason"]))
                else:
                    return session.error_result(STOPPED.format(reason=data["reason"]))
        except asyncio.CancelledError:
            # The app gave up on the call, while it was being sent or after.
            call.abandon()
            raise
        finally:
            if asking is not None:
                asking.cancel()

    server = Server(
        SERVER_NAME,
        version=__version__,
        title=SERVER_TITLE,
        instructions=session.instructions,
        on_list_tools=list_tools,
        on_call_tool=call_tool,
    )
    return server, NotificationOptions(tools_changed=True)


def _open_viewer(path: str, *, port: int | None, open_browser: bool, web_source_url: str | None, web_cache_dir: str | None) -> AnnDataViewer:
    """Serve ``path`` in a viewer of this process (``--open``) and open its page."""
    from . import jupyter

    options: dict[str, Any] = {}
    if web_source_url is not None:
        options["web_source_url"] = web_source_url
    if web_cache_dir is not None:
        options["web_cache_dir"] = web_cache_dir
    # Started once the MCP transport holds stdout (it sends what anything
    # else writes there to stderr), so nothing it prints reaches the channel.
    viewer = jupyter.AnnDataViewer(str(Path(path).expanduser()), port=port, auto_open=False, **options)
    print(f"Cellucid viewer: {viewer.viewer_url}", file=sys.stderr, flush=True)
    if open_browser and webbrowser.open(viewer.viewer_url) is not True:
        print("Open the address above in a browser: no browser could be opened from here.", file=sys.stderr, flush=True)
    return viewer


def run_mcp(
    open: str | Path | None = None,  # noqa: A002 - the command's own word
    *,
    viewer: str | None = None,
    port: int | None = None,
    open_browser: bool = True,
    web_source_url: str | None = None,
    web_cache_dir: str | Path | None = None,
) -> None:
    """Run ``cellucid mcp`` on this process's stdin and stdout until the app
    closes them.

    Args:
        open: A file or folder to serve and open in a viewer of this process
            (any format Cellucid reads); its page opens in the browser, and
            the viewer stops when the server ends. None drives a viewer
            opened elsewhere on this computer.
        viewer: Which published viewer to drive (its viewer id, dataset id or
            dataset name); None takes the most recently opened one whose page
            is open. Ignored with ``open``.
        port: The port of the viewer ``open`` starts (None: 8765, or the next
            free one).
        open_browser: Whether ``open`` opens the viewer's page in the browser
            (False: its address is printed on stderr).
        web_source_url: Where the viewer's page is loaded from (default:
            https://www.cellucid.com).
        web_cache_dir: Where the viewer keeps its copy of the page.
    """
    _server_class, _options_class, stdio_server = _import_sdk()
    if open is not None and viewer is not None:
        raise ValueError("Give --open or --viewer, not both: --open drives the viewer it opens.")
    session = LocalAgentSession(choice=viewer)
    holder: dict[str, Any] = {}
    opened: list[AnnDataViewer] = []
    server, notification_options = _build_server(session, holder)

    async def main() -> None:
        loop = asyncio.get_running_loop()

        def tools_changed() -> None:
            mcp_session = holder.get("session")
            if mcp_session is None:
                return

            async def send() -> None:
                with contextlib.suppress(Exception):
                    await mcp_session.send_tool_list_changed()

            loop.call_soon_threadsafe(lambda: asyncio.ensure_future(send()))

        session.add_tools_listener(tools_changed)

        async with stdio_server() as (read_stream, write_stream):
            if open is not None:
                target = str(open)
                session.set_opened(None, opening=f"The Cellucid viewer of {target} is still opening; try again in a moment.")

                def start() -> None:
                    try:
                        started = _open_viewer(
                            target,
                            port=port,
                            open_browser=open_browser,
                            web_source_url=web_source_url,
                            web_cache_dir=None if web_cache_dir is None else str(web_cache_dir),
                        )
                    except Exception as error:  # said to the agent and on stderr; the server stays up
                        message = f"The Cellucid viewer of {target} could not open: {error}"
                        print(message, file=sys.stderr, flush=True)
                        session.set_opened(None, opening=message)
                        return
                    opened.append(started)
                    session.set_opened(started._agent_port)

                threading.Thread(target=start, name="cellucid-mcp-open", daemon=True).start()
            elif await asyncio.to_thread(session.link) is None:
                print(f"cellucid mcp: {session.status()['advice']} The server waits; its tools work once one is open.", file=sys.stderr, flush=True)
            await server.run(read_stream, write_stream, server.create_initialization_options(notification_options))

    try:
        asyncio.run(main())
    finally:
        session.close()
        for started in opened:
            with contextlib.suppress(Exception), contextlib.redirect_stdout(sys.stderr):
                started.stop()
