"""One stand-in for cellucid's AnnData server, for tests that open a viewer
without serving anything.

It holds every attribute and method the viewers read from their server
(``self._server.<name>`` in ``cellucid/jupyter/_viewer.py`` and
``_widget.py``), in one place: a test subclasses it for what it records.
``test_the_stand_in_server_has_what_the_viewers_read`` fails when a viewer
starts reading something it lacks. (A stand-in written by hand in each test
file falls behind the server the day it gains an attribute.)
"""

from __future__ import annotations

from typing import Any


class StandInServer:
    """Records how it was made, serves nothing. ``port_offset`` moves the
    port it reports bound from the one asked (a server bound elsewhere)."""

    port_offset = 0

    def __init__(self, **kwargs: Any) -> None:
        self.host = kwargs.get("host") or "127.0.0.1"
        self.port = int(kwargs.get("port") or 8765) + self.port_offset
        self.dataset_id = "demo"
        self.dataset_name = "Demo"
        # The id the data declared (None: it declares none), as the real server keeps it.
        self.declared_dataset_id: str | None = None
        self.data_format = "anndata"
        self.server_info: dict[str, Any] = {"n_cells": 3}
        self.runtime = "python"
        self.viewer_channel = False
        # With a port not given, the next free one when it is taken (the viewer sets it).
        self.next_free_port = False
        self.stop_calls = 0

    @property
    def url(self) -> str:
        return f"http://{self.host}:{self.port}"

    def start(self) -> None:
        return None

    def start_background(self) -> None:
        return None

    def is_running(self) -> bool:
        return True

    def stop(self) -> None:
        self.stop_calls += 1
