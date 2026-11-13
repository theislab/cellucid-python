"""``cellucid mcp``: the local MCP server for agent apps on this computer.

- ``src/cellucid/mcp_tools.json`` is the web connector's tool list for this
  package's own ``action_catalog.json``, byte for byte the web repository's
  ``tests/fixtures/mcp-tools.json`` (written by its
  ``node scripts/generate-action-catalog.mjs``; cellucid-r holds the same
  copy).
- ``tests/fixtures/local_agents.json`` is shared with cellucid-r byte for
  byte (its home is here; cellucid-r copies it): the agents folder, the
  record, the endpoint, the app names and the sentences an agent reads, so
  either package's ``cellucid mcp`` drives either package's viewer and says
  the same things.
- The agent commands and events of the viewer channel, a viewer's agent
  port, its publication on this computer, and the MCP server itself over
  stdio (initialize, tools/list, tools/call, errors, cancellation, the
  person asked in the app, a tool list that never changes) against a fake page.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import queue
import re
import subprocess
import sys
import threading
import time
import urllib.parse as urllib_parse
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from _serverless_viewer import serverless_viewer
from _viewer_channel_client import TEXT, FakeViewer, binary_event
from _web_root import SYNC_SHARED_FILES, assert_shared_file, rendered_by_web

import cellucid._agent_link as agent_link
import cellucid._mcp_core as core
import cellucid.jupyter._wire as wire
from cellucid.anndata_server import AnnDataServer

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).parent / "fixtures" / "local_agents.json"


def _fixture() -> dict[str, Any]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _has_sdk() -> bool:
    try:
        import mcp.server.lowlevel.server  # noqa: F401
    except ImportError:
        return False
    return True


needs_sdk = pytest.mark.skipif(not _has_sdk(), reason="the MCP SDK (cellucid[mcp]) is not installed")


# =============================================================================
# The generated tools and the shared fixture
# =============================================================================


@pytest.mark.core
def test_mcp_tools_json_is_the_web_apps_byte_for_byte() -> None:
    assert_shared_file("src/cellucid/mcp_tools.json", "tests/fixtures/mcp-tools.json")


@pytest.mark.core
def test_mcp_tools_json_is_the_list_the_web_apps_source_makes_now() -> None:
    """Against the web app's source (``renderMcpTools``), not only its committed file."""
    shipped = (ROOT / "src" / "cellucid" / "mcp_tools.json").read_text(encoding="utf-8")
    assert rendered_by_web("scripts/generate-action-catalog.mjs", "renderMcpTools") == shipped, (
        f"src/cellucid/mcp_tools.json is not what the web app's source makes: run `{SYNC_SHARED_FILES}` "
        "in the web repository (cellucid)."
    )


@pytest.mark.core
def test_mcp_tools_json_is_the_web_connectors_list_for_this_packages_catalogue() -> None:
    document = core.mcp_tools_document()
    catalog = json.loads((ROOT / "src" / "cellucid" / "action_catalog.json").read_text(encoding="utf-8"))
    assert document["catalogDigest"] == hashlib.sha256((ROOT / "src" / "cellucid" / "action_catalog.json").read_bytes()).hexdigest()
    # One tool per action, named, titled and annotated from it, in the catalogue's order.
    assert [tool["name"] for tool in document["tools"]] == [action["tool"] for action in catalog["actions"]]
    for tool, action in zip(document["tools"], catalog["actions"], strict=True):
        assert tool["title"] == action["title"] == tool["annotations"]["title"]
        assert tool["description"].startswith(action["summary"])
        assert tool["annotations"]["readOnlyHint"] is (action["effect"] == "read")
        assert tool["annotations"]["destructiveHint"] is (action["effect"] != "read" and action["confirm"] is not None)
        assert '"x-' not in json.dumps(tool["inputSchema"]), "Cellucid's own keywords are said in words"
    assert document["statusTool"]["name"] == _fixture()["server"]["statusTool"]
    assert document["statusTool"]["name"] not in {tool["name"] for tool in document["tools"]}
    assert "a Cellucid viewer on this computer" in document["instructions"]
    assert "https://www.cellucid.com" not in document["instructions"]
    # Every tool name the instructions use exists.
    names = {tool["name"] for tool in document["tools"]} | {document["statusTool"]["name"]}
    for name in re.findall(r"`([a-z][a-z0-9_]*)`", document["instructions"]):
        assert name in names, name


@pytest.mark.core
def test_the_local_instructions_name_the_notebook_calls_that_bring_the_agents_work_back() -> None:
    """The instructions name the calls that apply the viewer's session, and they exist: an
    agent told of a call that does not exist tells the person to run it."""
    from cellucid.jupyter import AnnDataViewer

    text = core.mcp_tools_document()["instructions"]
    for call in ("viewer.request_session()", "viewer.apply_session_to_anndata(adata, request)", "plan.apply()"):
        assert call in text, call
    assert callable(AnnDataViewer.request_session) and callable(AnnDataViewer.apply_session_to_anndata)
    from cellucid.anndata_session._apply import ApplyPlan

    assert callable(ApplyPlan.apply)
    # The R calls it names are cellucid-r's (checked there against its own functions).
    for call in ("cellucid_apply_session(object, cellucid_request_session(viewer))", "cellucid_apply_plan(plan)"):
        assert call in text, call


@pytest.mark.core
def test_a_stale_tools_document_is_refused_by_name(monkeypatch: pytest.MonkeyPatch) -> None:
    core.mcp_tools_document.cache_clear()
    real = core.resources.files

    class Package:
        def __init__(self, base: Any) -> None:
            self.base = base

        def joinpath(self, name: str) -> Any:
            if name == "action_catalog.json":
                return type("Stale", (), {"read_bytes": staticmethod(lambda: b"{}\n")})()
            return self.base.joinpath(name)

    monkeypatch.setattr(core.resources, "files", lambda package: Package(real(package)))
    try:
        with pytest.raises(core.McpToolsDocumentError, match="another action catalog"):
            core.mcp_tools_document()
    finally:
        core.mcp_tools_document.cache_clear()


@pytest.mark.core
def test_cellucid_r_holds_the_same_local_agents_fixture() -> None:
    assert_shared_file(
        "tests/fixtures/local_agents.json", "tests/testthat/fixtures/local_agents.json", repository="cellucid-r"
    )


@pytest.mark.core
def test_the_shared_local_agents_fixture_is_this_packages_behaviour() -> None:
    fixture = _fixture()
    assert fixture["server"]["name"] == core.SERVER_NAME
    assert fixture["server"]["title"] == core.SERVER_TITLE
    assert fixture["record"]["schema"] == agent_link.AGENT_RECORD_SCHEMA
    assert set(fixture["record"]["fields"]) == agent_link.RECORD_FIELDS
    endpoint = fixture["endpoint"]
    assert endpoint["route"] == agent_link.AGENT_ROUTE
    assert endpoint["eventWaitLimitSeconds"] == agent_link.EVENT_WAIT_LIMIT_SECONDS
    assert endpoint["consumerIdleSeconds"] == agent_link.CONSUMER_IDLE_SECONDS
    assert endpoint["commandBodyLimitBytes"] == agent_link.COMMAND_BODY_LIMIT_BYTES
    for case in fixture["agentNames"]:
        assert core.agent_name_of(case["name"], case["title"]) == case["agent"], case
    sentences = fixture["sentences"]
    assert sentences["noViewer"] == core.NO_VIEWER
    assert sentences["notFound"] == core.NOT_FOUND
    assert sentences["notOpen"] == agent_link.not_open_sentence({"dataset": {"name": "PBMC"}, "url": "http://127.0.0.1:8765"})
    assert sentences["pageGone"] == agent_link.PAGE_GONE
    assert sentences["notTaken"] == agent_link.NOT_TAKEN
    assert sentences["stopped"] == agent_link.STOPPED
    assert sentences["folderUnsafe"] == agent_link.FOLDER_UNSAFE
    assert fixture["optOut"]["variable"] == agent_link.OPT_OUT_VARIABLE
    assert sentences["notOpenWidget"] == agent_link.not_open_sentence({"dataset": {"name": "PBMC"}, "url": None})
    import cellucid.mcp as mcp_module

    assert sentences["sdkMissing"] == mcp_module.SDK_MISSING
    session = core.LocalAgentSession(discover=list)
    assert re.match(fixture["server"]["connectionIdPattern"], session.connection_id)


@pytest.mark.core
def test_the_agents_folder_follows_the_user_cache_and_its_override(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("CELLUCID_AGENTS_DIR", str(tmp_path))
    assert agent_link.agents_dir() == tmp_path
    monkeypatch.setenv("CELLUCID_AGENTS_DIR", "relative/folder")
    with pytest.raises(agent_link.AgentsFolderError, match="absolute"):
        agent_link.agents_dir()
    # Refused as a shared folder is: nothing found, and the agent reads why
    # (cellucid-r's sentence), where a ValueError ended cellucid mcp at start.
    with pytest.raises(agent_link.AgentsFolderError, match="absolute"):
        agent_link.discover_viewers()
    session = core.LocalAgentSession()
    assert session.link() is None
    assert session.status()["advice"] == "CELLUCID_AGENTS_DIR must be an absolute path"
    monkeypatch.delenv("CELLUCID_AGENTS_DIR")
    folder = agent_link.agents_dir()
    assert folder.parts[-2:] == ("cellucid", "agents")
    if sys.platform == "darwin":
        assert folder == Path.home() / "Library" / "Caches" / "cellucid" / "agents"


# =============================================================================
# The agent commands and events on the viewer channel
# =============================================================================


@pytest.mark.core
def test_agent_commands_serialize_as_the_viewer_protocol_has_them() -> None:
    protocol = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))
    for name in ("agent_list", "agent_call", "agent_cancel", "agent_answer"):
        assert wire._COMMAND_MESSAGE_FIELDS[name] == frozenset(protocol["commands"][name])
    for name in ("agent_asking", "agent_answered"):
        assert wire._INBOUND_EVENT_FIELDS[name] == frozenset(protocol["events"][name])
    # Exactly the agent events the protocol has: nothing says the tool list changed.
    agent_events = {name for name in [*protocol["events"], *protocol["binaryEvents"]] if name.startswith("agent_")}
    assert agent_events == wire._AGENT_EVENT_TYPES == {"agent_tools", "agent_result", "agent_asking", "agent_answered"}
    for name in ("agent_tools", "agent_result"):
        assert wire._INBOUND_EVENT_FIELDS[name] == frozenset(protocol["binaryEvents"][name]) | {"body"}
    assert frozenset(protocol["agent"]["connectionFields"]) == wire._AGENT_CONNECTION_FIELDS
    assert wire._AGENT_CONNECTION_ID_PATTERN.pattern == protocol["agent"]["connectionIdPattern"]
    for example in protocol["examples"]["commands"]:
        if example["type"].startswith("agent_"):
            assert json.loads(wire._require_frontend_message(example)) == example
    connection = {"id": "local-1", "client": None}
    with pytest.raises(ValueError, match="connection id"):
        wire._require_frontend_message({"type": "agent_list", "requestId": "a", "connection": {"id": "has space", "client": None}})
    with pytest.raises(TypeError, match="args must be one object"):
        wire._require_frontend_message({"type": "agent_call", "requestId": "a", "tool": "color_by", "args": [], "agent": "A", "connection": connection})
    with pytest.raises(ValueError, match="names its tool"):
        wire._require_frontend_message({"type": "agent_call", "requestId": "a", "tool": "Color By", "args": {}, "agent": "A", "connection": connection})
    with pytest.raises(TypeError, match="allowed"):
        wire._require_frontend_message({"type": "agent_answer", "requestId": "a", "allowed": "yes"})


@pytest.mark.core
def test_agent_events_are_read_as_the_viewer_protocol_has_them() -> None:
    protocol = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))
    for example in protocol["examples"]["events"]:
        if example["type"].startswith("agent_"):
            event, data = wire._require_inbound_jupyter_event(example)
            assert {"type": event, **data} == example
    for example in protocol["examples"]["binaryEvents"]:
        header = example.get("header", {})
        if header.get("type") not in {"agent_tools", "agent_result"}:
            continue
        body = json.dumps(example["body"]).encode("utf-8")
        event, data = wire._require_inbound_jupyter_event({**header, "body": body})
        key = "tools" if event == "agent_tools" else "result"
        assert data == {"requestId": header["requestId"], key: example["body"]}
    with pytest.raises(TypeError, match="MCP tool result"):
        wire._require_inbound_jupyter_event({"type": "agent_result", "requestId": "a", "body": b'{"content": "no"}'})
    with pytest.raises(ValueError, match="not UTF-8 JSON"):
        wire._require_inbound_jupyter_event({"type": "agent_tools", "requestId": "a", "body": b"\xff"})


# =============================================================================
# A viewer's agent port, its publication, and a fake page
# =============================================================================


@contextlib.contextmanager
def _published_viewer() -> Iterator[tuple[Any, int]]:
    """A registered viewer published to agents, a running data server, no browser."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.AnnData(
        X=np.zeros((2, 1), np.float32),
        obs=pd.DataFrame(index=["a", "b"]),
        var=pd.DataFrame(index=["g"]),
        obsm={"X_umap_2d": np.array([[0, 0], [1, 1]], np.float32)},
    )
    server = AnnDataServer(adata, host="127.0.0.1", port=0, quiet=True, serve_web_ui=False, dataset_name="Agents", dataset_id="agents")
    server.start_background()
    assert server._server is not None
    port = server._server.server_address[1]
    viewer = serverless_viewer(client_server_url=f"http://127.0.0.1:{port}")
    viewer._server = server
    viewer._activate()
    try:
        yield viewer, port
    finally:
        viewer._server = None
        viewer.stop()
        server.stop()


class FakePage:
    """A page on the viewer channel that answers agent commands as the web app does."""

    def __init__(self, port: int, viewer: Any, *, tools: list[dict[str, Any]] | None = None) -> None:
        self.socket = FakeViewer(port, timeout=30)
        assert self.socket.hello(viewer._viewer_id, viewer._viewer_token) == {"type": "welcome"}
        self.tools = tools if tools is not None else [{"name": "color_by", "title": "Color", "description": "From the page.", "inputSchema": {"type": "object"}, "annotations": {}}]
        self.commands: queue.Queue[dict[str, Any]] = queue.Queue()
        self.ask: set[str] = set()
        self.hold: set[str] = set()
        self._held: dict[str, dict[str, Any]] = {}
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                opcode, payload = self.socket.receive()
            except (OSError, AssertionError, ConnectionError):
                return
            if opcode != TEXT:
                continue
            command = json.loads(payload.decode("utf-8"))
            self.commands.put(command)
            kind = command["type"]
            if kind == "agent_list":
                self.socket.send_binary(binary_event({"type": "agent_tools", "requestId": command["requestId"]}, json.dumps(self.tools).encode()))
            elif kind == "agent_call":
                if command["tool"] in self.ask:
                    self._held[command["requestId"]] = command
                    self.socket.send_json({"type": "agent_asking", "requestId": command["requestId"], "message": "Claude Code (local) asks to save a copy. It writes a file.", "detail": None})
                elif command["tool"] in self.hold:
                    self._held[command["requestId"]] = command
                else:
                    self.result(command["requestId"], {"content": [{"type": "text", "text": f"ran {command['tool']} as {command['agent']}"}], "isError": False})
            elif kind == "agent_answer":
                self.socket.send_json({"type": "agent_answered", "requestId": command["requestId"], "allowed": command["allowed"]})
                text = "saved" if command["allowed"] else "The person declined: save a copy."
                self.result(command["requestId"], {"content": [{"type": "text", "text": text}], "isError": not command["allowed"]})
            elif kind == "agent_cancel":
                self._held.pop(command["requestId"], None)
                self.result(command["requestId"], {"content": [{"type": "text", "text": "Stopped before it finished: the agent stopped waiting"}], "isError": True})

    def result(self, request_id: str, result: dict[str, Any]) -> None:
        self.socket.send_binary(binary_event({"type": "agent_result", "requestId": request_id}, json.dumps(result).encode()))

    def next_command(self, kind: str, timeout: float = 10.0) -> dict[str, Any]:
        deadline = time.monotonic() + timeout
        while True:
            command = self.commands.get(timeout=max(0.01, deadline - time.monotonic()))
            if command["type"] == kind:
                return command

    def close(self) -> None:
        self._stop.set()
        self.socket.close()


def _listen() -> tuple[queue.Queue[tuple[str, dict[str, Any]]], Any]:
    heard: queue.Queue[tuple[str, dict[str, Any]]] = queue.Queue()
    return heard, lambda kind, data: heard.put((kind, data))


def _page_gone(viewer: Any, monkeypatch: pytest.MonkeyPatch) -> threading.Event:
    """Set once ``viewer`` heard its page go (its channel let the socket go
    first, so from then on it says no page is open): what a test waits on
    after closing the fake page, instead of asking again and again."""
    gone = threading.Event()
    port = viewer._agent_port
    page_gone = port.page_gone

    def heard() -> None:
        page_gone()
        gone.set()

    monkeypatch.setattr(port, "page_gone", heard)
    return gone


@pytest.mark.core
def test_a_viewer_is_published_and_its_port_carries_calls_both_ways() -> None:
    with _published_viewer() as (viewer, port):
        port_ = viewer._agent_port
        assert port_.connected() is False
        with pytest.raises(agent_link.ViewerUnavailableError, match=r"No page shows the Cellucid viewer of Agents \(served at http://127\.0\.0\.1:\d+\) now\. Ask the person to open its page again"):
            port_.start_call("r1", "color_by", {}, "Claude Code", {"id": "local-1", "client": None}, lambda *_: None)
        records = agent_link.discover_viewers()
        assert [record["viewer"] for record in records] == [viewer._viewer_id]
        record = records[0]
        assert set(record) == agent_link.RECORD_FIELDS
        assert record["dataset"] == {"id": "agents", "name": "Agents"}
        assert record["pid"] == os.getpid()
        assert record["host"] == agent_link._this_host()
        # Where it is served, never the page address with the viewer's token.
        assert record["url"] == f"http://127.0.0.1:{port}"
        assert viewer._viewer_token not in json.dumps(record["url"])
        if sys.platform != "win32":
            assert (Path(os.environ["CELLUCID_AGENTS_DIR"]) / f"{viewer._viewer_id}.json").stat().st_mode & 0o077 == 0
        page = FakePage(port, viewer)
        try:
            assert port_.connected() is True
            connection = {"id": "local-abc", "client": "claude-code"}
            assert port_.list_tools(connection, timeout=10) == page.tools
            heard, listener = _listen()
            port_.start_call("local-abc-1", "color_by", {"field": "x"}, "Claude Code", connection, listener)
            command = page.next_command("agent_call")
            assert command == {"type": "agent_call", "requestId": "local-abc-1", "tool": "color_by", "args": {"field": "x"}, "agent": "Claude Code", "connection": connection}
            assert heard.get(timeout=10) == ("result", {"result": {"content": [{"type": "text", "text": "ran color_by as Claude Code"}], "isError": False}})
            # The person asked: heard, answered from the app, settled.
            page.ask.add("save_copy")
            heard, listener = _listen()
            port_.start_call("local-abc-2", "save_copy", {}, "Claude Code", connection, listener)
            kind, data = heard.get(timeout=10)
            assert kind == "asking" and data["message"].startswith("Claude Code (local) asks")
            port_.answer("local-abc-2", False)
            assert heard.get(timeout=10) == ("answered", {"allowed": False})
            kind, data = heard.get(timeout=10)
            assert kind == "result" and data["result"]["isError"] is True
            # Stopping the viewer ends what waits.
            page.hold.add("slow")
            heard, listener = _listen()
            port_.start_call("local-abc-3", "slow", {}, "Claude Code", connection, listener)
            page.next_command("agent_call")
            # The viewer stops while its page is still open: the call hears it stopped.
            viewer._agent_port.stopped("the viewer was stopped")
            assert heard.get(timeout=10) == ("stopped", {"reason": "the viewer was stopped"})
        finally:
            page.close()
    assert agent_link.discover_viewers() == []


def test_a_command_error_of_an_agent_command_reaches_its_caller_not_the_notebook(capsys: pytest.CaptureFixture[str]) -> None:
    viewer = serverless_viewer()
    heard, listener = _listen()
    viewer._agent_port._calls["r1"] = listener
    viewer._handle_frontend_message({"type": "command_error", "command": "agent_call", "requestId": "r1", "reason": "This viewer does not know that host command; reload the page."})
    reason = "This viewer does not know that host command; reload the page."
    assert heard.get(timeout=1) == ("failed", {"reason": agent_link.NOT_TAKEN.format(reason=reason)})
    assert capsys.readouterr().err == ""
    assert viewer.state.command_error is None


def test_the_endpoint_refuses_browsers_other_hosts_and_wrong_tokens() -> None:
    with _published_viewer() as (viewer, _port):
        record = agent_link.discover_viewers()[0]
        endpoint = record["endpoint"]

        def ask(headers: dict[str, str], path: str = "") -> int:
            import urllib.error
            import urllib.request

            request = urllib.request.Request(endpoint + path, headers=headers)
            try:
                with agent_link._NO_PROXY_OPENER.open(request, timeout=5) as response:
                    return response.status
            except urllib.error.HTTPError as error:
                return error.code

        good = {"Authorization": f"Bearer {record['token']}"}
        assert ask(good) == 200
        assert ask({"Authorization": "Bearer nope"}) == 401
        assert ask({**good, "Origin": "https://evil.example"}) == 403
        assert ask({**good, "Host": "evil.example"}) == 403
        assert ask(good, "/nothing") == 404
        assert ask(good, "/events") == 400
        status, payload = agent_link._endpoint_request(record, "POST", "/command", {"type": "agent_call", "requestId": "r", "tool": "x", "args": {}, "agent": "A", "connection": {"id": "c", "client": None}})
        assert status == 409 and "No page shows the Cellucid viewer of Agents (served at" in payload["error"]
        status, payload = agent_link._endpoint_request(record, "POST", "/command", {"type": "nothing", "requestId": "r"})
        assert status == 400
        assert viewer._viewer_id in endpoint


def test_a_remote_link_drives_a_published_viewer_and_hears_it_stop() -> None:
    with _published_viewer() as (viewer, port):
        record = agent_link.discover_viewers()[0]
        page = FakePage(port, viewer)
        link = agent_link.RemoteViewerLink(record, "local-remote")
        try:
            assert link.describe()["connected"] is True
            connection = {"id": "local-remote", "client": None}
            assert link.list_tools(connection, timeout=10) == page.tools
            heard, listener = _listen()
            link.start_call("local-remote-1", "color_by", {}, "Cursor", connection, listener)
            assert heard.get(timeout=10)[0] == "result"
            page.hold.add("slow")
            heard, listener = _listen()
            link.start_call("local-remote-2", "slow", {}, "Cursor", connection, listener)
            page.next_command("agent_call")
            link.cancel("local-remote-2")
            assert page.next_command("agent_cancel")["requestId"] == "local-remote-2"
            heard, listener = _listen()
            link.start_call("local-remote-3", "slow", {}, "Cursor", connection, listener)
            page.next_command("agent_call")
            # The viewer stops while its page is still open: the call hears it stopped.
            viewer._agent_port.stopped("the viewer was stopped")
            agent_link.withdraw_viewer(viewer._viewer_id)
            assert heard.get(timeout=30) == ("stopped", {"reason": "the viewer was stopped"})
        finally:
            page.close()
        assert link.stopped == "the viewer was stopped"
        link.close()


def test_the_session_chooses_the_newest_viewer_whose_page_is_open_and_names_its_agent() -> None:
    session = core.LocalAgentSession(discover=list)
    assert session.link() is None
    status = session.status()
    assert status == {"agent": "Local agent", "connected": False, "viewer": None, "viewers": [], "advice": core.NO_VIEWER}
    session.adopt_client("claude-code")
    assert session.agent == "Claude Code"
    assert session.connection == {"id": session.connection_id, "client": "claude-code"}
    names = [tool["name"] for tool in session.tools()]
    assert names[-1] == "cellucid_status"
    assert names[:-1] == [tool["name"] for tool in core.mcp_tools_document()["tools"]]
    with pytest.raises(agent_link.ViewerUnavailableError, match="No Cellucid viewer is open on this computer"):
        session.start_call("r", "color_by", {}, lambda *_: None)
    chosen = core.LocalAgentSession(discover=list, choice="pbmc")
    assert chosen.status()["advice"] == core.NOT_FOUND.format(choice='"pbmc"')
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        try:
            session = core.LocalAgentSession()
            session.adopt_client("cursor-vscode")
            assert session.tools()[:-1] == page.tools
            assert page.next_command("agent_list")["connection"] == {"id": session.connection_id, "client": "cursor-vscode"}
            status = session.status()
            assert status["connected"] is True and status["viewer"]["viewer"] == viewer._viewer_id
            assert status["viewers"][0]["driven"] is True
            heard, listener = _listen()
            session.start_call(session.next_request_id(), "color_by", {}, listener)
            assert page.next_command("agent_call")["agent"] == "Cursor"
            assert heard.get(timeout=10)[0] == "result"
            session.close()
        finally:
            page.close()


# =============================================================================
# The MCP server over stdio
# =============================================================================


class StdioClient:
    """An MCP client app: JSON-RPC 2.0 over the server's stdin and stdout.

    Every wait ends on what it waits for or on the server ending (its stdout
    closing), which fails it at once with the server's stderr. Starting and
    ending the server are waited for on those events alone, never on a
    clock: they are its imports and its teardown, as long as the computer
    makes them (an antivirus starved by CPU load stretched a start past a
    minute on a developer's Mac), and the server answers or ends either way.
    The other waits are bounded by what the test sends and the fake page
    answers; their deadline only turns a hang into a failure.
    """

    def __init__(self, env: dict[str, str], *args: str) -> None:
        self.process = subprocess.Popen(
            [sys.executable, "-m", "cellucid.cli", "mcp", *args],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # stdout is the MCP channel, UTF-8 on every platform; the pipes'
            # default encoding is the Windows code page, which garbles it.
            encoding="utf-8",
            errors="replace",
            env=env,
        )
        #: The server's messages; None once its stdout closed (it ended).
        self.messages: queue.Queue[dict[str, Any] | None] = queue.Queue()
        self.stderr: list[str] = []
        self._stderr_changed = threading.Condition()
        self._stderr_ended = False
        self.garbage: list[str] = []
        threading.Thread(target=self._read, daemon=True).start()
        threading.Thread(target=self._read_stderr, daemon=True).start()
        self._held: list[dict[str, Any]] = []

    def _read(self) -> None:
        assert self.process.stdout is not None
        try:
            for line in self.process.stdout:
                try:
                    self.messages.put(json.loads(line))
                except ValueError:
                    # stdout is the MCP channel: anything else there is a defect.
                    self.garbage.append(line)
        finally:
            self.messages.put(None)

    def _read_stderr(self) -> None:
        assert self.process.stderr is not None
        try:
            for line in self.process.stderr:
                with self._stderr_changed:
                    self.stderr.append(line)
                    self._stderr_changed.notify_all()
        finally:
            with self._stderr_changed:
                self._stderr_ended = True
                self._stderr_changed.notify_all()

    def send(self, message: dict[str, Any]) -> None:
        assert self.process.stdin is not None
        self.process.stdin.write(json.dumps({"jsonrpc": "2.0", **message}) + "\n")
        self.process.stdin.flush()

    def next(self, match: Any, timeout: float | None = 30.0) -> dict[str, Any]:
        """The next message ``match`` takes (``timeout`` None: until the server answers or ends)."""
        for index, message in enumerate(self._held):
            if match(message):
                return self._held.pop(index)
        deadline = None if timeout is None else time.monotonic() + timeout
        while True:
            message = self.messages.get(timeout=None if deadline is None else max(0.01, deadline - time.monotonic()))
            if message is None:
                self.messages.put(None)  # every later wait hears it ended too
                with self._stderr_changed:
                    self._stderr_changed.wait_for(lambda: self._stderr_ended)
                    stderr = "".join(self.stderr)
                raise AssertionError(f"The server ended (exit status {self.process.wait()}) before it sent the message waited for:\n{stderr}")
            if match(message):
                return message
            self._held.append(message)

    def reply(self, request_id: int, timeout: float | None = 30.0) -> dict[str, Any]:
        return self.next(lambda message: message.get("id") == request_id and "method" not in message, timeout)

    def wait_stderr(self, found: Any) -> list[str]:
        """The server's stderr lines once ``found(lines)`` holds; fails when it ends first."""
        with self._stderr_changed:
            self._stderr_changed.wait_for(lambda: found(self.stderr) or self._stderr_ended)
            assert found(self.stderr), "The server ended first:\n" + "".join(self.stderr)
            return list(self.stderr)

    def initialize(self, capabilities: dict[str, Any] | None = None, name: str = "claude-code") -> dict[str, Any]:
        self.send({"id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18", "capabilities": capabilities or {}, "clientInfo": {"name": name, "version": "1"}}})
        # The server's start: until it answers or ends.
        answer = self.reply(1, timeout=None)
        self.send({"method": "notifications/initialized"})
        return answer

    def close(self) -> int:
        """Close the server's stdin: its exit status once it ended (it ends by itself)."""
        with contextlib.suppress(OSError):
            assert self.process.stdin is not None
            self.process.stdin.close()
        return self.process.wait()


def _server_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT / "src"), env.get("PYTHONPATH", "")])
    return env


@needs_sdk
def test_the_server_without_a_viewer_lists_the_generated_tools_and_says_what_is_missing() -> None:
    client = StdioClient(_server_env())
    try:
        answer = client.initialize()["result"]
        assert answer["serverInfo"]["name"] == "cellucid" and answer["serverInfo"]["title"] == "Cellucid"
        # The list never changes while the server runs, and it says so.
        assert answer["capabilities"]["tools"] == {"listChanged": False}
        assert answer["instructions"] == core.mcp_tools_document()["instructions"]
        client.send({"id": 2, "method": "tools/list"})
        tools = client.reply(2)["result"]["tools"]
        assert tools == [*core.mcp_tools_document()["tools"], core.mcp_tools_document()["statusTool"]]
        client.send({"id": 3, "method": "tools/call", "params": {"name": "color_by", "arguments": {"field": {"source": "obs", "key": "x"}}}})
        assert client.reply(3)["result"] == {"content": [{"type": "text", "text": core.refusal_text(core.NO_VIEWER)}], "isError": True}
        client.send({"id": 4, "method": "tools/call", "params": {"name": "cellucid_status", "arguments": {}}})
        status_text = client.reply(4)["result"]["content"][0]["text"]
        tags = core.mcp_tools_document()["resultData"]["tags"]
        assert status_text.startswith(tags["open"]) and status_text.endswith(tags["close"])
        status = json.loads(status_text[len(tags["open"]) : -len(tags["close"])])
        assert status == {"agent": "Claude Code", "connected": False, "viewer": None, "viewers": [], "advice": core.NO_VIEWER}
        client.send({"id": 5, "method": "no/such/method"})
        assert client.reply(5)["error"]["code"] == -32601
        client.send({"id": 6, "method": "ping"})
        assert client.reply(6)["result"] == {}
        client.send({"id": 7, "method": "tools/call", "params": {"name": "color_by", "arguments": "not an object"}})
        assert client.reply(7)["error"]["code"] == -32602
    finally:
        assert client.close() == 0
    assert any("No Cellucid viewer is open on this computer" in line for line in client.stderr)
    assert not any(line.startswith("{") for line in client.stderr)


@needs_sdk
def test_the_tool_list_is_the_same_whatever_the_page_does_and_is_never_said_to_change() -> None:
    """The person's settings are told to the agent after its tool results,
    never in its tool list: a list change makes the app rewrite its
    cached prompt and drop the earlier reasoning. So the list an agent gets is
    the generated one before a page opens, the same once a page (whatever its
    settings) took the channel, reloaded or was replaced, and no
    ``notifications/tools/list_changed`` is ever sent."""
    document = core.mcp_tools_document()
    generated = [*document["tools"], document["statusTool"]]
    with _published_viewer() as (viewer, port):
        client = StdioClient(_server_env())
        pages: list[FakePage] = []
        try:
            assert client.initialize()["result"]["capabilities"]["tools"] == {"listChanged": False}
            client.send({"id": 2, "method": "tools/list"})
            assert client.reply(2)["result"]["tools"] == generated
            # A page takes the channel, reloads, and is replaced by another: a
            # page lists the same tools whatever the person's settings in it.
            pages.append(FakePage(port, viewer, tools=document["tools"]))
            pages[-1].close()
            pages.append(FakePage(port, viewer, tools=document["tools"]))
            pages.append(FakePage(port, viewer, tools=document["tools"]))
            client.send({"id": 3, "method": "tools/list"})
            assert client.reply(3)["result"]["tools"] == generated
            pages[-1].next_command("agent_list")
            with pytest.raises(queue.Empty):
                client.next(lambda message: message.get("method") == "notifications/tools/list_changed", timeout=2)
        finally:
            for page in pages:
                page.close()
            assert client.close() == 0


@needs_sdk
def test_the_server_drives_a_published_viewer_cancels_asks_and_asks_in_the_app(monkeypatch: pytest.MonkeyPatch) -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page_gone = _page_gone(viewer, monkeypatch)
        client = StdioClient(_server_env())
        try:
            client.initialize(capabilities={"elicitation": {}}, name="claude-code")
            client.send({"id": 2, "method": "tools/list"})
            assert [tool["name"] for tool in client.reply(2)["result"]["tools"]] == ["color_by", "cellucid_status"]
            connection = page.next_command("agent_list")["connection"]
            assert connection["client"] == "claude-code"
            # A call runs in the page as the app's agent.
            # A tool name no list gives is a tool error of one sentence; nothing reaches the page.
            client.send({"id": 9, "method": "tools/call", "params": {"name": "Color By", "arguments": {}}})
            assert client.reply(9)["result"] == {"content": [{"type": "text", "text": core.refusal_text("agent_call names its tool as the tool list does, such as 'color_by'")}], "isError": True}
            client.send({"id": 3, "method": "tools/call", "params": {"name": "color_by", "arguments": {"field": {"source": "obs", "key": "leiden"}}}})
            assert client.reply(3)["result"] == {"content": [{"type": "text", "text": "ran color_by as Claude Code"}], "isError": False}
            call = page.next_command("agent_call")
            assert call["args"] == {"field": {"source": "obs", "key": "leiden"}} and call["connection"] == connection
            # The app gives up on a call: it is stopped in the page.
            page.hold.add("slow")
            client.send({"id": 4, "method": "tools/call", "params": {"name": "slow", "arguments": {}}})
            request_id = page.next_command("agent_call")["requestId"]
            client.send({"method": "notifications/cancelled", "params": {"requestId": 4, "reason": "The person pressed Esc"}})
            assert page.next_command("agent_cancel")["requestId"] == request_id
            # The person is asked: in the viewer's card, and in the app (an elicitation); the app's answer settles it.
            page.ask.add("save_copy")
            client.send({"id": 5, "method": "tools/call", "params": {"name": "save_copy", "arguments": {}, "_meta": {"progressToken": "p5"}}})
            # An app that asked for progress hears what the call waits on.
            progress = client.next(lambda message: message.get("method") == "notifications/progress")
            assert progress["params"]["progressToken"] == "p5"
            assert progress["params"]["message"].startswith("Claude Code (local) asks to save a copy.")
            asked = client.next(lambda message: message.get("method") == "elicitation/create")
            assert asked["params"]["message"].startswith("Claude Code (local) asks to save a copy.")
            client.send({"id": asked["id"], "result": {"action": "accept", "content": {}}})
            assert page.next_command("agent_answer")["allowed"] is True
            assert client.reply(5)["result"] == {"content": [{"type": "text", "text": "saved"}], "isError": False}
            # The next question is asked in the app too, and the app's Decline is the person's.
            client.send({"id": 7, "method": "tools/call", "params": {"name": "save_copy", "arguments": {}}})
            asked_again = client.next(lambda message: message.get("method") == "elicitation/create")
            client.send({"id": asked_again["id"], "result": {"action": "decline"}})
            assert page.next_command("agent_answer")["allowed"] is False
            assert client.reply(7)["result"]["isError"] is True
            # The page goes: a call says so in one sentence.
            page.close()
            assert page_gone.wait(30)
            client.send({"id": 6, "method": "tools/call", "params": {"name": "color_by", "arguments": {}}})
            text = client.reply(6)["result"]["content"][0]["text"]
            assert text.startswith("No page shows the Cellucid viewer of Agents (served at")
        finally:
            page.close()
            assert client.close() == 0


@needs_sdk
def test_a_call_the_app_cancels_while_it_is_still_being_sent_is_stopped_in_the_page(monkeypatch: pytest.MonkeyPatch) -> None:
    """The page can run a call before the server's send of it is back, and
    the app's cancel can come in between: the call is stopped in the page all
    the same, once the send is back. The viewer's answer to the send is held
    until the server has taken the cancel -- its answer to a ping sent after
    the cancel says so, as it reads its input in order -- so the cancel lands
    while the send is under way, every time."""
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.hold.add("slow")
        port_ = viewer._agent_port
        start_call = port_.start_call
        sent_back = threading.Event()

        def start_call_answered_later(*args: Any, **kwargs: Any) -> None:
            start_call(*args, **kwargs)
            sent_back.wait()

        monkeypatch.setattr(port_, "start_call", start_call_answered_later)
        client = StdioClient(_server_env())
        try:
            client.initialize()
            client.send({"id": 2, "method": "tools/call", "params": {"name": "slow", "arguments": {}}})
            request_id = page.next_command("agent_call")["requestId"]
            client.send({"method": "notifications/cancelled", "params": {"requestId": 2, "reason": "The person pressed Esc"}})
            client.send({"id": 3, "method": "ping"})
            assert client.reply(3)["result"] == {}
            sent_back.set()
            assert page.next_command("agent_cancel")["requestId"] == request_id
            # Canceled: the app is not answered (a ping after it is, and nothing before it).
            client.send({"id": 4, "method": "ping"})
            client.reply(4)
            assert not any(message.get("id") == 2 for message in client._held)
        finally:
            sent_back.set()
            page.close()
            assert client.close() == 0


@needs_sdk
def test_the_server_reports_a_viewer_that_stops_while_a_call_waits() -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.hold.add("slow")
        client = StdioClient(_server_env())
        try:
            client.initialize()
            client.send({"id": 2, "method": "tools/call", "params": {"name": "slow", "arguments": {}}})
            page.next_command("agent_call")
            viewer._agent_port.stopped("the viewer was stopped")
            agent_link.withdraw_viewer(viewer._viewer_id)
            text = client.reply(2)["result"]["content"][0]["text"]
            assert text == "The Cellucid viewer stopped (the viewer was stopped); nothing more was done."
        finally:
            page.close()
            assert client.close() == 0


def test_cellucid_mcp_help_and_its_refusals() -> None:
    from cellucid import cli

    parser = cli.create_parser()
    help_text = parser._subparsers._group_actions[0].choices["mcp"].format_help()  # type: ignore[union-attr]
    for phrase in ("--open PATH", "--viewer NAME", "pip install 'cellucid[mcp]'", "claude mcp add cellucid -- cellucid mcp", '"mcpServers"'):
        assert phrase in help_text
    assert "cellucid mcp --open" in parser.format_help()
    from cellucid.mcp import run_mcp

    if _has_sdk():
        with pytest.raises(ValueError, match="not both"):
            run_mcp("data.h5ad", viewer="pbmc")


def test_without_the_sdk_cellucid_mcp_says_how_to_install_it(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    import builtins

    real_import = builtins.__import__

    def without_mcp(name: str, globals: Any = None, locals: Any = None, fromlist: Any = (), level: int = 0) -> Any:  # noqa: A002 - __import__'s own names
        if level == 0 and (name == "mcp" or name.startswith("mcp.")):
            raise ImportError(f"No module named {name!r}", name=name)
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", without_mcp)
    from cellucid import cli

    assert cli.main(["mcp"]) == 1
    assert capsys.readouterr().err.strip().splitlines()[-1] == f"Error: {_fixture()['sentences']['sdkMissing']}"


# =============================================================================
# The page going, a viewer of cellucid serve, the agents folder's safety
# =============================================================================


def test_a_page_that_goes_mid_call_fails_the_call() -> None:
    """A reload or a closed tab: the call it was running gets no answer, so it
    fails at once with the page-gone sentence (never a hang)."""
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.hold.add("slow")
        heard, listener = _listen()
        connection = {"id": "local-x", "client": None}
        viewer._agent_port.start_call("local-x-1", "slow", {}, "A", connection, listener)
        page.next_command("agent_call")
        page.close()
        assert heard.get(timeout=10) == ("failed", {"reason": agent_link.PAGE_GONE})
        page = FakePage(port, viewer)
        page.hold.add("slow")
        try:
            # A page replacing the one holding the channel (opened again elsewhere) fails its calls too.
            heard, listener = _listen()
            viewer._agent_port.start_call("local-x-2", "slow", {}, "A", connection, listener)
            page.next_command("agent_call")
            newer = FakePage(port, viewer)
            assert heard.get(timeout=10) == ("failed", {"reason": agent_link.PAGE_GONE})
            newer.close()
        finally:
            page.close()


@needs_sdk
def test_a_page_reload_mid_call_is_answered_through_cellucid_mcp() -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.hold.add("slow")
        client = StdioClient(_server_env())
        try:
            client.initialize()
            client.send({"id": 2, "method": "tools/call", "params": {"name": "slow", "arguments": {}}})
            page.next_command("agent_call")
            page.close()
            page = FakePage(port, viewer)
            assert client.reply(2, timeout=30)["result"] == {"content": [{"type": "text", "text": core.refusal_text(agent_link.PAGE_GONE)}], "isError": True}
        finally:
            page.close()
            assert client.close() == 0


@needs_sdk
def test_progress_grows_and_a_url_only_app_is_not_asked_a_form() -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.ask.add("save_copy")
        client = StdioClient(_server_env())
        try:
            client.initialize(capabilities={"elicitation": {"url": {}}})
            client.send({"id": 2, "method": "tools/call", "params": {"name": "save_copy", "arguments": {}, "_meta": {"progressToken": "p"}}})
            progress = client.next(lambda message: message.get("method") == "notifications/progress")
            assert progress["params"]["progress"] == 1
            request_id = page.next_command("agent_call")["requestId"]
            with pytest.raises(queue.Empty):
                client.next(lambda message: message.get("method") == "elicitation/create", timeout=2)
            page.socket.send_json({"type": "agent_answered", "requestId": request_id, "allowed": False})
            page.result(request_id, {"content": [{"type": "text", "text": "The person declined"}], "isError": True})
            assert client.reply(2)["result"]["isError"] is True
        finally:
            page.close()
            assert client.close() == 0


@needs_sdk
def test_a_question_answered_in_the_page_is_withdrawn_in_the_app() -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.ask.add("save_copy")
        client = StdioClient(_server_env())
        try:
            client.initialize(capabilities={"elicitation": {}})
            client.send({"id": 5, "method": "tools/call", "params": {"name": "save_copy", "arguments": {}}})
            asked = client.next(lambda message: message.get("method") == "elicitation/create")
            request_id = page.next_command("agent_call")["requestId"]
            page.socket.send_json({"type": "agent_answered", "requestId": request_id, "allowed": True})
            page.result(request_id, {"content": [{"type": "text", "text": "saved"}], "isError": False})
            assert client.reply(5)["result"]["content"][0]["text"] == "saved"
            withdrawn = client.next(lambda message: message.get("method") == "notifications/cancelled")
            assert withdrawn["params"]["requestId"] == asked["id"]
            # The app's late answer is not passed on.
            client.send({"id": asked["id"], "result": {"action": "accept", "content": {}}})
            with pytest.raises(queue.Empty):
                page.next_command("agent_answer", timeout=2)
        finally:
            page.close()
            assert client.close() == 0


def test_a_silent_consumer_is_swept_its_calls_cancelled_and_its_events_let_go(monkeypatch: pytest.MonkeyPatch) -> None:
    with _published_viewer() as (viewer, port):
        page = FakePage(port, viewer)
        page.hold.add("slow")
        try:
            record = agent_link.discover_viewers()[0]
            command = {"type": "agent_call", "requestId": "gone-1", "tool": "slow", "args": {}, "agent": "A", "connection": {"id": "gone", "client": None}}
            assert agent_link._endpoint_request(record, "POST", "/command", command)[0] == 202
            page.next_command("agent_call")
            published = agent_link._endpoint.viewers[viewer._viewer_id]  # type: ignore[union-attr]
            with published.lock:
                published.consumers["gone"].events.append({"type": "agent_failed", "requestId": "gone-0", "reason": "not read"})
                published.consumers["gone"].seen -= agent_link.CONSUMER_IDLE_SECONDS + 1
                published.sweep()
                assert "gone" not in published.consumers
            assert page.next_command("agent_cancel")["requestId"] == "gone-1"
        finally:
            page.close()


def test_the_session_leaves_a_viewer_whose_page_closed_for_one_whose_page_is_open(monkeypatch: pytest.MonkeyPatch) -> None:
    with _published_viewer() as (first, first_port):
        first_page = FakePage(first_port, first)
        first_gone = _page_gone(first, monkeypatch)
        session = core.LocalAgentSession()
        assert session.link().record["viewer"] == first._viewer_id  # type: ignore[union-attr]
        with _published_viewer() as (second, second_port):
            second_page = FakePage(second_port, second)
            try:
                # Still open: the driven viewer stays.
                assert session.link().record["viewer"] == first._viewer_id  # type: ignore[union-attr]
                first_page.close()
                assert first_gone.wait(30)
                assert session.link().record["viewer"] == second._viewer_id  # type: ignore[union-attr]
            finally:
                second_page.close()
                session.close()


def test_a_viewer_of_cellucid_serve_holds_a_channel_and_is_published() -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.AnnData(X=np.zeros((2, 1), np.float32), obs=pd.DataFrame(index=["a", "b"]), var=pd.DataFrame(index=["g"]),
                       obsm={"X_umap": np.array([[0, 0], [1, 1]], np.float32)})
    server = AnnDataServer(adata, host="127.0.0.1", port=0, quiet=True, serve_web_ui=False, dataset_name="Served", dataset_id="served")
    server.viewer_channel = True
    server.start_background()
    try:
        host_viewer = server._host_viewer
        assert host_viewer is not None
        assert server.viewer_url == (
            f"{server.url}/?viewer=true&viewerId={host_viewer.viewer_id}&viewerToken={host_viewer.viewer_token}&anndata=true"
        )
        [record] = agent_link.discover_viewers()
        assert record["viewer"] == host_viewer.viewer_id and record["dataset"] == {"id": "served", "name": "Served"}
        assert record["url"] == server.url
        page = FakeViewer(server.port, timeout=30)
        assert page.hello(host_viewer.viewer_id, host_viewer.viewer_token) == {"type": "welcome"}
        heard, listener = _listen()
        host_viewer.port.start_call("local-s-1", "color_by", {}, "Cursor", {"id": "local-s", "client": None}, listener)
        command = json.loads(page.receive()[1].decode("utf-8"))
        assert command["type"] == "agent_call" and command["agent"] == "Cursor"
        # Its other events (a hover) reach no notebook, and are no error.
        page.send_json({"type": "hover", "cell": None, "position": None})
        page.send_binary(binary_event({"type": "agent_result", "requestId": "local-s-1"}, b'{"content": [], "isError": false}'))
        assert heard.get(timeout=10) == ("result", {"result": {"content": [], "isError": False}})
        page.close()
    finally:
        server.stop()
    assert agent_link.discover_viewers() == []


def test_serve_anndata_asks_for_the_viewer_channel(monkeypatch: pytest.MonkeyPatch) -> None:
    import cellucid.anndata_server as served

    started: list[bool] = []
    monkeypatch.setattr(served.AnnDataServer, "start", lambda self, blocking=True: started.append(self.viewer_channel))
    served.serve_anndata(Path(__file__).parent / "fixtures" / "older_anndata" / "legacy-0.6-twin-current.h5ad", open_browser=False, quiet=True, serve_web_ui=False)
    assert started == [True]


def _record(folder: Path, viewer: str, **changes: Any) -> Path:
    record = {
        "schema": agent_link.AGENT_RECORD_SCHEMA, "viewer": viewer, "runtime": "python", "host": agent_link._this_host(),
        "pid": os.getpid(), "started": "2026-11-04T00:00:00.000000Z", "dataset": {"id": "d", "name": "D"},
        "url": "http://127.0.0.1:9", "endpoint": f"http://127.0.0.1:9/_cellucid/agent/{viewer}", "token": "t",
    }
    record.update(changes)
    path = folder / f"{viewer}.json"
    path.write_text(json.dumps(record), encoding="utf-8")
    return path


def test_records_of_another_host_malformed_or_off_loopback_are_not_read(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    folder = tmp_path / "agents"
    folder.mkdir(mode=0o700)
    monkeypatch.setenv("CELLUCID_AGENTS_DIR", str(folder))
    good = _record(folder, "good")
    other_host = _record(folder, "elsewhere", host="another-machine", pid=999_999_999)
    _record(folder, "remote", endpoint="http://10.0.0.5:8000/_cellucid/agent/remote")
    _record(folder, "typed", pid="1")
    _record(folder, "dataset", dataset={"id": 3})
    (folder / "broken.json").write_text("{", encoding="utf-8")
    assert [record["viewer"] for record in agent_link.discover_viewers()] == ["good"]
    assert other_host.exists(), "another computer's record is left alone"
    dead = _record(folder, "dead", pid=999_999_999)
    agent_link.discover_viewers()
    assert not dead.exists() and good.exists()


@pytest.mark.skipif(os.name != "posix", reason="POSIX ownership and modes")
def test_an_agents_folder_others_may_write_is_refused_with_one_sentence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    folder = tmp_path / "agents"
    folder.mkdir(mode=0o700)
    os.chmod(folder, 0o777)
    monkeypatch.setenv("CELLUCID_AGENTS_DIR", str(folder))
    with pytest.raises(agent_link.AgentsFolderError) as refusal:
        agent_link.discover_viewers()
    assert str(refusal.value) == agent_link.FOLDER_UNSAFE.format(folder=folder, why="others may write in it")
    session = core.LocalAgentSession()
    assert session.status()["advice"] == str(refusal.value)
    # Publication skips it, said in the log; the viewer works on.
    with _published_viewer() as (viewer, _port):
        assert list(folder.glob("*.json")) == []
        assert viewer._viewer_id not in agent_link._endpoint.viewers  # type: ignore[union-attr]
    monkeypatch.setattr(os, "getuid", lambda: os.stat(folder).st_uid + 1)
    os.chmod(folder, 0o700)
    with pytest.raises(agent_link.AgentsFolderError, match="another user owns it"):
        agent_link.discover_viewers()


def test_cellucid_agents_off_publishes_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CELLUCID_AGENTS", "Off ")
    with _published_viewer() as (viewer, _port):
        assert agent_link.discover_viewers() == []
        assert agent_link._endpoint is None or viewer._viewer_id not in agent_link._endpoint.viewers


def test_a_command_with_a_length_that_is_no_number_is_refused() -> None:
    import http.client

    with _published_viewer() as (_viewer, _port):
        record = agent_link.discover_viewers()[0]
        parts = urllib_parse.urlsplit(record["endpoint"])
        connection = http.client.HTTPConnection(parts.hostname, parts.port, timeout=5)
        connection.putrequest("POST", parts.path + "/command")
        connection.putheader("Authorization", f"Bearer {record['token']}")
        connection.putheader("Content-Length", "many")
        connection.endheaders()
        assert connection.getresponse().status == 400


def test_a_command_nested_past_the_interpreter_is_refused_as_no_command() -> None:
    import http.client

    with _published_viewer() as (_viewer, _port):
        record = agent_link.discover_viewers()[0]
        parts = urllib_parse.urlsplit(record["endpoint"])
        body = b"[" * 200_000
        connection = http.client.HTTPConnection(parts.hostname, parts.port, timeout=5)
        connection.request("POST", parts.path + "/command", body=body, headers={"Authorization": f"Bearer {record['token']}"})
        response = connection.getresponse()
        assert response.status == 400
        assert json.loads(response.read()) == {"error": "A command is one JSON object."}


def test_a_connection_to_the_endpoint_that_sends_nothing_is_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    # Anyone on this computer may connect without the token: an idle
    # connection held its thread for as long as it stayed open.
    import socket

    monkeypatch.setattr(agent_link._EndpointHandler, "timeout", 0.3)
    with _published_viewer() as (_viewer, _port):
        parts = urllib_parse.urlsplit(agent_link.discover_viewers()[0]["endpoint"])
        with socket.create_connection((parts.hostname, parts.port), timeout=10) as idle:
            started = time.monotonic()
            assert idle.recv(1) == b""
            assert time.monotonic() - started < 5
    assert agent_link._EndpointHandler.timeout == 0.3
    monkeypatch.undo()
    assert agent_link._EndpointHandler.timeout == agent_link.CONNECTION_IDLE_SECONDS


@needs_sdk
def test_open_says_nothing_on_stdout_but_the_protocol(tmp_path: Path) -> None:
    """``cellucid mcp --open``: the viewer starts while the app is answered,
    and stdout carries only MCP messages (the review's race probe)."""
    data = Path(__file__).parent / "fixtures" / "older_anndata" / "legacy-0.6-twin-current.h5ad"
    for _ in range(3):
        env = {**_server_env(), "CELLUCID_AGENTS_DIR": str(tmp_path / "agents")}
        client = StdioClient(env, "--open", str(data), "--no-browser", "--port", "0", "--web-source-url", "http://127.0.0.1:9", "--web-cache-dir", str(tmp_path / "web"))
        try:
            assert client.initialize()["result"]["serverInfo"]["name"] == "cellucid"
            client.send({"id": 2, "method": "tools/call", "params": {"name": "cellucid_status", "arguments": {}}})
            client.reply(2)
            # The viewer's start, its imports included: until it says it started or could not.
            client.wait_stderr(lambda lines: any("could not open" in line or "Cellucid viewer:" in line for line in lines))
        finally:
            assert client.close() == 0
        assert client.garbage == []


def test_names_written_by_others_are_marked_as_the_tab_marks_them() -> None:
    """One rule in three languages: the package's status result and its own
    sentences mark names (a viewer's, a dataset's) as the web app's tool
    results and refusals do, byte for byte (the examples ``mcp_tools.json``
    carries, made by the web app's own functions, with a name meant to close
    the tags)."""
    rule = core.mcp_tools_document()["resultData"]
    assert rule["tags"] == {"open": "<cellucid-data>", "close": "</cellucid-data>"}
    data = rule["examples"]["data"]
    assert core.result_data_text(data["json"]) == data["text"]
    assert core.result_data_text(core._js_json(json.loads(data["json"]))) == data["text"]
    name = rule["examples"]["name"]
    assert core.quote_name(name["name"]) == name["text"]
    sentence = rule["examples"]["sentence"]
    assert sentence["built"] == core.NOT_FOUND.format(choice=core.quote_name(name["name"]))
    assert core.refusal_text(sentence["built"]) == sentence["text"]
    # Nothing in a name can close a mark: one closing tag per marked name, the real one.
    assert sentence["text"].count(rule["tags"]["close"]) == 1
    assert core.LocalAgentSession.error_result(sentence["built"]) == {"content": [{"type": "text", "text": sentence["text"]}], "isError": True}
    # A quoted name nobody marked gains no mark: its tags are escaped.
    unmarked = rule["examples"]["unmarked"]
    assert core.refusal_text(unmarked["built"]) == unmarked["text"]
    assert rule["tags"]["close"] not in unmarked["text"]
