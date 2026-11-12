"""
The viewer's history, live: what changed, notebook turns, sessions pushed in.

The viewer tells this package every change of its history (a step recorded,
the step shown moved, an edit of the history's structure, a tracker started,
used, renamed or removed, a session loaded) as a ``history_changed`` event
(``viewer-protocol.json``): :attr:`LiveHistory.last_change` is the latest,
:meth:`LiveHistory.on_change` runs code on each. The history itself is read
the way a session file is (:meth:`LiveHistory.read`: the viewer's history as a
:class:`~cellucid.CellucidSession`, its trackers, branches, steps, notes and
turns), and changed through the history's actions, the same ones the State
panel and agents use: the viewer's generated ``history_*`` methods
(``viewer.history_undo()``, ``viewer.history_note_add(...)``).

A notebook **turn** (``with viewer.turn("QC"):``) makes every action sent
inside it one turn of the viewer's history: the State panel shows them as
one group and undoes them together.

A session file is pushed into the viewer with :meth:`AnnDataViewer.load_session`:
its trackers are added beside the person's own.
"""

from __future__ import annotations

import copy
from collections.abc import Callable
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ._viewer import AnnDataViewer, SessionRequest, _EventWaiter


class LoadSessionRequest:
    """A session file pushed into the viewer (:meth:`AnnDataViewer.load_session`).

    ``done`` says whether the viewer answered; :meth:`result` returns what it
    shows then: ``{"tracker": {"id", "name"}, "head"}``.
    """

    def __init__(self, viewer: AnnDataViewer, request_id: str, name: str, waiter: _EventWaiter) -> None:
        self._viewer = viewer
        self.request_id = request_id
        #: The file's name, as the loaded tracker records it.
        self.name = name
        self._waiter: _EventWaiter | None = waiter

    @property
    def done(self) -> bool:
        """Whether the viewer answered (with the tracker shown, or with why not)."""
        return self._waiter is not None and self._waiter.payload is not None

    def result(self, timeout: float | None = 60.0) -> dict[str, Any]:
        """
        What the viewer shows once the file is loaded, waiting up to ``timeout``
        seconds (None: forever): ``{"tracker": {"id", "name"}, "head"}``.

        Raises ``RuntimeError`` with the viewer's reason when it did not load
        the file (another dataset's session, a damaged file), ``TimeoutError``
        when no answer arrived in time -- or, in a notebook widget, when the
        answer has not arrived yet: the viewer answers between cells.
        """
        waiter = self._waiter
        if waiter is None:
            raise RuntimeError("This session load was canceled.")
        if waiter.payload is None:
            self._viewer._await(waiter, timeout, f"Timed out waiting for the viewer to load {self.name}")
        assert waiter.payload is not None
        if waiter.event == "command_error":
            raise RuntimeError(f"The viewer did not load {self.name}: {waiter.payload['reason']}")
        return {"tracker": copy.deepcopy(waiter.payload["tracker"]), "head": waiter.payload["head"]}

    def cancel(self) -> None:
        """Stop waiting for the answer (the file may be loaded already)."""
        if self._waiter is None or self._waiter.payload is not None:
            return
        viewer = self._viewer
        with viewer._event_cv:
            viewer._action_waiters.pop(self.request_id, None)
            viewer._action_types.pop(self.request_id, None)
        self._waiter = None

    def __repr__(self) -> str:
        state = "canceled" if self._waiter is None else "answered" if self.done else "waiting"
        return f"LoadSessionRequest({self.name!r}, {state})"


class LiveHistory:
    """The history of a viewer, as it changes (``viewer.history``).

    ``last_change`` is the latest ``history_changed`` the viewer sent
    (``{"kind", "tracker", "head", "step", "branch"}``), with ``tracker``,
    ``head`` and ``branch`` read from it; :meth:`on_change` runs a callback
    on each change; :meth:`read` asks the viewer for its history as a session
    (the same reading model as a file: :class:`~cellucid.history.Tracker`,
    branches, steps as a table). The history is changed through the
    viewer's generated ``history_*`` methods, each returning its
    :class:`~cellucid.jupyter.ActionRequest`; an action the viewer cannot
    carry out yet is refused with the viewer's reason (``request.result()``
    raises it).
    """

    def __init__(self, viewer: AnnDataViewer) -> None:
        self._viewer = viewer

    # -- what the viewer said -----------------------------------------------

    @property
    def last_change(self) -> dict[str, Any] | None:
        """The latest change the viewer told (a copy), or None before the first."""
        change = self._viewer.state.history
        return None if change is None else dict(change)

    @property
    def tracker(self) -> str | None:
        """The id of the tracker in use, as of the latest change."""
        change = self._viewer.state.history
        return None if change is None else str(change["tracker"])

    @property
    def head(self) -> str | None:
        """The step shown, as of the latest change."""
        change = self._viewer.state.history
        return None if change is None else str(change["head"])

    @property
    def branch(self) -> str | None:
        """The id of the active branch, as of the latest change (None when the viewer names none)."""
        change = self._viewer.state.history
        return None if change is None else change["branch"]

    def on_change(self, callback: Callable[[dict[str, Any]], Any]) -> Callable[[dict[str, Any]], Any]:
        """Run ``callback(change)`` on every change of the viewer's history (a decorator too)."""
        return self._viewer.register_hook("history_changed", callback)

    def read(self) -> SessionRequest:
        """Ask the viewer for its history now, as a session: ``.result()`` is
        the :class:`~cellucid.CellucidSession` (in a notebook, collect it in a
        later cell)."""
        return self._viewer.request_session()

    def __repr__(self) -> str:
        change = self._viewer.state.history
        if change is None:
            return "LiveHistory(no change told yet)"
        return f"LiveHistory(tracker {change['tracker']}, head {change['head']}, branch {change['branch']})"
