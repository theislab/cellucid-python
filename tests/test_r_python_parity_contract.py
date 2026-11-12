"""cellucid-python and cellucid-r are peers: one public surface, two languages.

Every public name of this package has its twin among cellucid-r's exports, or
the stated reason it has none; every hand-written export of cellucid-r is the
twin of a name here, or stated as R's own. The actions' functions are made
from the one action catalog on both sides and must agree name for name.
``prepare()`` and ``cellucid_prepare()`` take the same arguments, and a file
served from either states the same launch defaults. A capability added to one
package fails this test until the other has it or the reason is written here.
"""

from __future__ import annotations

import inspect
import json
import re
from pathlib import Path

from _web_root import sibling_file

import cellucid
from cellucid.anndata_server import AnnDataServer, serve_anndata
from cellucid.anndata_session import ApplyPlan
from cellucid.history import Tracker
from cellucid.jupyter import (
    ActionRequest,
    AnnDataViewer,
    LiveHistory,
    LoadSessionRequest,
    SessionRequest,
    ViewerActions,
    show_anndata,
)
from cellucid.prepare_data import prepare
from cellucid.session import CellucidSession

#: Each public name of this package and its cellucid-r twin; a value starting
#: with "-" is the stated reason R has no function of its own for it.
TWINS: dict[str, str] = {
    # cellucid
    "prepare": "cellucid_prepare",
    "restore_original_order": "cellucid_restore_original_order",
    "show_anndata": "cellucid_view",
    "Widget": "- the viewer cellucid_view() returns in a Jupyter notebook (pane = \"notebook\")",
    "AnnDataViewer": "- the viewer object cellucid_view() returns",
    "serve_anndata": "cellucid_serve",
    "AnnDataServer": "- the viewer object cellucid_serve() returns",
    "read_session": "cellucid_read_session",
    "CellucidSession": "- the cellucid_session object cellucid_read_session() returns",
    "apply_cellucid_session_to_anndata": "cellucid_apply_session",
    "ApplyPlan": "- the cellucid_apply_plan object cellucid_apply_session() returns",
    "ApplyReport": "- the object's attribute \"cellucid_apply_report\"",
    "history_applied": "cellucid_history_applied",
    "pull": "cellucid_pull",
    "write_gmt": "cellucid_write_gmt",
    "read_gmt": "cellucid_read_gmt",
    "PullSummary": "- cellucid_pull(return_summary = TRUE)",
    "Selection": "- the list cellucid_selection() returns",
    "ViewerResult": "- the records cellucid_results() returns",
    "clear_web_cache": "cellucid_clear_web_cache",
    "get_web_cache_dir": "- tools::R_user_dir(\"cellucid\", \"cache\")/web, stated on ?cellucid_clear_web_cache",
    "subsample_cells": "cellucid_subsample_cells",
    "compute_transition_drift": "cellucid_transition_drift",
    "add_transition_drift_to_obsm": "- an R object has no obsm: cellucid_prepare(vector_fields = ) stores the drift",
    "__version__": "- utils::packageVersion(\"cellucid\")",
    # A viewer
    "AnnDataViewer.act": "cellucid_act",
    "AnnDataViewer.apply_session_to_anndata": "cellucid_apply_session",
    "AnnDataViewer.clear_hooks": "cellucid_off",
    "AnnDataViewer.clear_web_cache": "cellucid_clear_web_cache",
    "AnnDataViewer.data_mapping": "- viewer$data_mapping",
    "AnnDataViewer.data_mapping_origin": "- viewer$data_mapping$origin",
    "AnnDataViewer.debug_connection": "cellucid_debug_connection",
    "AnnDataViewer.display": "- cellucid_view(pane = ) shows the viewer",
    "AnnDataViewer.ensure_web_ui_cached": "- cellucid_clear_web_cache(), then the next viewer downloads it",
    "AnnDataViewer.gene_sets": "cellucid_gene_sets",
    "AnnDataViewer.history": "cellucid_last_history_change",
    "AnnDataViewer.load_session": "cellucid_load_session",
    "AnnDataViewer.obs_names": "- the object's own cell names (colnames)",
    "AnnDataViewer.on_click": "cellucid_on",
    "AnnDataViewer.on_command_error": "cellucid_on",
    "AnnDataViewer.on_data_mapping": "cellucid_on",
    "AnnDataViewer.on_history_change": "cellucid_on",
    "AnnDataViewer.on_hover": "cellucid_on",
    "AnnDataViewer.on_message": "cellucid_on",
    "AnnDataViewer.on_ready": "cellucid_on",
    "AnnDataViewer.on_selection": "cellucid_on",
    "AnnDataViewer.pull": "cellucid_pull",
    "AnnDataViewer.push": "cellucid_push",
    "AnnDataViewer.register_hook": "cellucid_on",
    "AnnDataViewer.request_session": "cellucid_request_session",
    "AnnDataViewer.results": "cellucid_results",
    "AnnDataViewer.selection": "cellucid_selection",
    "AnnDataViewer.send_message": "- the channel's raw commands are the package's own (every one has its function)",
    "AnnDataViewer.server_url": "- viewer$url",
    "AnnDataViewer.stop": "cellucid_stop",
    "AnnDataViewer.subset": "cellucid_subset",
    "AnnDataViewer.turn": "cellucid_turn",
    "AnnDataViewer.unregister_hook": "cellucid_off",
    "AnnDataViewer.viewer_url": "- viewer$viewer_url",
    "AnnDataViewer.wait_for_event": "cellucid_wait_event",
    "AnnDataViewer.wait_for_ready": "cellucid_wait_ready",
    # Requests
    "ActionRequest.result": "cellucid_action_result",
    "ActionRequest.done": "- printing the request says waiting, answered or failed",
    "ActionRequest.step": "- request$step",
    "ActionRequest.unrecorded": "- request$unrecorded (cellucid_action_result() warns it, class cellucid_action_not_recorded)",
    "ActionRequest.cancel": "- R waits only inside cellucid_action_result(): nothing waits to be canceled",
    "SessionRequest.result": "cellucid_session_result",
    "SessionRequest.done": "- printing the request says waiting, answered or failed",
    "SessionRequest.cancel": "- R waits only inside cellucid_session_result(): nothing waits to be canceled",
    "LoadSessionRequest.result": "cellucid_load_session_result",
    "LoadSessionRequest.done": "- printing the request says waiting, answered or failed",
    "LoadSessionRequest.cancel": "- R waits only inside cellucid_load_session_result(): nothing waits to be canceled",
    "LiveHistory.last_change": "cellucid_last_history_change",
    "LiveHistory.tracker": "cellucid_last_history_change",
    "LiveHistory.head": "cellucid_last_history_change",
    "LiveHistory.branch": "cellucid_last_history_change",
    "LiveHistory.on_change": "cellucid_on",
    "LiveHistory.read": "cellucid_request_session",
    # A session and its history
    "CellucidSession.apply": "cellucid_apply_session",
    "CellucidSession.conversations": "cellucid_conversations",
    "CellucidSession.tracker": "cellucid_tracker",
    "CellucidSession.save": "- cellucid_session_result(file = ) keeps the file; one read from a path is that file",
    "CellucidSession.cells": "- session$cells",
    "CellucidSession.chrome": "- session$chrome",
    "CellucidSession.dataset_fingerprint": "- session$dataset_fingerprint",
    "CellucidSession.genes": "- session$genes",
    "CellucidSession.manifest": "- session$manifest",
    "CellucidSession.requires": "- session$requires",
    "CellucidSession.sha256": "- session$sha256",
    "CellucidSession.source": "- session$source",
    "Tracker.branch": "cellucid_branch",
    "Tracker.group": "cellucid_group",
    "Tracker.methods": "cellucid_methods",
    "Tracker.notes_table": "cellucid_notes",
    "Tracker.path_to": "- cellucid_steps(tracker, branch) is a branch's path from the root",
    "Tracker.pseudocode": "cellucid_pseudocode",
    "Tracker.stale": "cellucid_stale",
    "Tracker.step": "cellucid_step",
    "Tracker.steps": "cellucid_steps",
    "Tracker.turn": "- one row of cellucid_turns(tracker)",
    "ApplyPlan.apply": "cellucid_apply_plan",
    "ApplyPlan.table": "- the cellucid_apply_plan object is the table",
    "ApplyPlan.writes_anything": "- any(plan$status == \"write\")",
}

#: cellucid-r's own exports that are no twin above, and why Python has none.
R_ONLY: dict[str, str] = {
    "cellucid_stop_all": "Python stops each viewer (viewer.stop()); interpreter exit stops them all",
    "cellucid_viewers": "Python keeps the viewer objects it returned",
    "cellucid_mcp": "Python's is the command line's `cellucid mcp`",
    "cellucid_trackers": "CellucidSession.tracker() and the session's manifest list them",
    "cellucid_branches": "Tracker.branches",
    "cellucid_groups": "Tracker.groups",
    "cellucid_journal": "Tracker.journal",
    "cellucid_turns": "Tracker.turns",
    "cellucid_step_state": "Step.state",
    "cellucid_step_diff": "Step.diff()",
}


def _public(cls: type, base: type | None = None) -> set[str]:
    names = {name for name in dir(cls) if not name.startswith("_")}
    if base is not None:
        names -= {name for name in dir(base) if not name.startswith("_")}
    return names


def _python_names() -> set[str]:
    names = set(cellucid.__all__)
    for cls, base in (
        (AnnDataViewer, ViewerActions),
        (ActionRequest, None),
        (SessionRequest, None),
        (LoadSessionRequest, None),
        (LiveHistory, None),
        (CellucidSession, None),
        (ApplyPlan, None),
    ):
        names |= {f"{cls.__name__}.{name}" for name in _public(cls, base)}
    # The history's methods (its fields are the records R's tables hold).
    fields = {field for field in getattr(Tracker, "__dataclass_fields__", {})}
    names |= {f"Tracker.{name}" for name in _public(Tracker) - fields}
    return names


def _r_exports() -> tuple[set[str], set[str]]:
    """cellucid-r's hand-written exports, and the generated block's."""
    namespace = sibling_file("cellucid-r", "NAMESPACE").read_text(encoding="utf-8")
    block = re.search(r"# BEGIN GENERATED ACTION EXPORTS.*?\n(.*?)# END GENERATED ACTION EXPORTS", namespace, re.S)
    assert block is not None
    generated = set(re.findall(r"^export\((\w+)\)$", block.group(1), re.M))
    every = set(re.findall(r"^export\((\w+)\)$", namespace, re.M))
    return every - generated, generated


def test_every_public_name_has_its_r_twin_or_a_stated_reason() -> None:
    hand_written, _generated = _r_exports()
    python_names = _python_names()
    assert set(TWINS) == python_names, (
        f"Python names without a twin decision: {sorted(python_names - set(TWINS))}; "
        f"twins of names Python no longer has: {sorted(set(TWINS) - python_names)}"
    )
    twins = {twin for twin in TWINS.values() if not twin.startswith("- ")}
    assert twins <= hand_written, f"twins cellucid-r does not export: {sorted(twins - hand_written)}"
    assert hand_written == twins | set(R_ONLY), (
        f"cellucid-r exports without a Python twin or a stated reason: {sorted(hand_written - twins - set(R_ONLY))}; "
        f"stated R-only exports it no longer has: {sorted(set(R_ONLY) - hand_written)}"
    )
    assert not twins & set(R_ONLY)


def test_the_actions_functions_are_the_same_in_both_packages() -> None:
    _hand_written, generated = _r_exports()
    methods = {name for name in vars(ViewerActions) if not name.startswith("_")}
    assert {f"cellucid_{name}" for name in methods} == generated
    # Every action of the catalog has its function in each package; the count
    # is said only by the generated action reference (a hand-kept count on the
    # R reference page went stale whenever an action was added or replaced).
    catalog = json.loads((Path(cellucid.__file__).parent / "action_catalog.json").read_text(encoding="utf-8"))
    assert {f"cellucid_{action['tool']}" for action in catalog["actions"]} == generated
    reference = (
        Path(__file__).resolve().parents[1]
        / "docs/user_guide/r_package/g_api_reference_coverage/01_public_functions_and_classes.md"
    ).read_text(encoding="utf-8")
    assert re.search(r"action of the viewer, \d+ of them", reference) is None


def _r_formals(source: str, name: str) -> list[str]:
    body = re.search(rf"^{name} <- function\((.*?)\n\) \{{", source, re.S | re.M)
    assert body is not None, name
    return re.findall(r"^\s+(\w+)(?: = [^\n]*?)?,?$", body.group(1), re.M)


def test_prepare_takes_the_same_arguments_in_both_packages() -> None:
    r_source = sibling_file("cellucid-r", "R/cellucid_prepare.R").read_text(encoding="utf-8")
    python = list(inspect.signature(prepare).parameters)
    r = _r_formals(r_source, "cellucid_prepare")
    assert python == r
    # The R reference page's signature block is the function's.
    reference = (
        Path(__file__).resolve().parents[1]
        / "docs/user_guide/r_package/g_api_reference_coverage/01_public_functions_and_classes.md"
    ).read_text(encoding="utf-8")
    block = re.search(r"```r\ncellucid_prepare\(\n(.*?)\n\)\n```", reference, re.S)
    assert block is not None
    assert re.findall(r"^\s+(\w+)", block.group(1), re.M) == r


def test_a_served_file_states_the_same_launch_defaults_from_both_packages() -> None:
    # Both packages name a served file's launch choices as prepare() and
    # cellucid_prepare() do (var_gene_id_column, embedding_defaults), with no
    # second name for any of them.
    serving = {"data", "port", "host", "open_browser", "quiet", "serve_web_ui", "web_source_url", "web_cache_dir",
               "allowed_hosts", "cache_bytes"}
    python = set(inspect.signature(serve_anndata).parameters) - serving
    # show_anndata(), AnnDataServer and the viewer's show take the same names.
    for function in (show_anndata, AnnDataServer.__init__):
        assert set(inspect.signature(function).parameters) - serving - {"self", "height"} == python, function
    r_source = sibling_file("cellucid-r", "R/viewer-api.R").read_text(encoding="utf-8")
    listed = re.search(r"^\.view_path_arguments <- c\((.*?)\n\)", r_source, re.S | re.M)
    assert listed is not None
    assert python == set(re.findall(r'"(\w+)"', listed.group(1)))


def test_both_packages_say_the_same_words_when_the_web_app_source_is_unreachable() -> None:
    """Offline start parity (R9X-15): the cached build served, or why not and
    how to get the app, in one text, filled the same way in both packages."""
    import re

    from cellucid import web_cache

    source = sibling_file("cellucid-r", "R/web-cache.R").read_text(encoding="utf-8")

    def r_template(name: str) -> str:
        match = re.search(rf"^{re.escape(name)} <- (paste0\((?P<parts>.*?)\n\)|\"(?P<one>[^\"]*)\")", source, re.M | re.S)
        assert match is not None, name
        if match.group("one") is not None:
            return match.group("one")
        return "".join(re.findall(r'"([^"]*)"', match.group("parts")))

    for python_name, r_name in (
        ("WEB_APP_OFFLINE_NOTICE", ".web_app_offline_notice"),
        ("WEB_APP_UNAVAILABLE", ".web_app_unavailable"),
        ("WEB_APP_NOT_CACHED", ".web_app_not_cached"),
        ("WEB_APP_CACHE_UNVERIFIED", ".web_app_cache_unverified"),
    ):
        assert r_template(r_name) == getattr(web_cache, python_name), python_name
