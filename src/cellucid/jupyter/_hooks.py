"""
Event hook registration and the latest-event snapshot.

:class:`HookRegistry` owns the callbacks a notebook registers with
``@viewer.on_selection`` and friends, and :class:`ViewerState` keeps the most
recent payload per event type. Neither touches a viewer, a server, or a socket,
so hook dispatch semantics can be exercised on their own.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

# Type alias for hook callbacks
HookCallback = Callable[[dict[str, Any]], None]


@dataclass
class HookRegistry:
    """
    Manages event hooks for viewer interactions.

    This class provides a clean interface for registering, unregistering,
    and triggering callbacks for various viewer events.

    Supported events:
        - selection: User selected cells (via lasso, click, etc.)
        - hover: User hovering over a cell
        - click: User clicked on a cell
        - ready: Viewer finished initial load
        - command_error: A command sent from the notebook failed in the viewer
        - message: Raw message (catches all events)

    Example:
        hooks = HookRegistry()

        @hooks.on('selection')
        def my_handler(event):
            print(event['cells'])

        # Or without decorator:
        hooks.register('selection', my_handler)

        # Trigger (called internally when frontend sends message):
        hooks.trigger('selection', {'cells': [1, 2, 3]})
    """

    _callbacks: dict[str, list[HookCallback]] = field(default_factory=dict)

    def register(self, event: str, callback: HookCallback) -> HookCallback:
        """
        Register a callback for an event.

        Args:
            event: Event name ('selection', 'hover', 'click', 'ready', 'message')
            callback: Function to call when event fires. Receives event dict.

        Returns:
            The callback (for decorator usage)
        """
        if event not in self._callbacks:
            self._callbacks[event] = []
        self._callbacks[event].append(callback)
        return callback

    def unregister(self, event: str, callback: HookCallback) -> bool:
        """
        Remove a callback from an event.

        Args:
            event: Event name
            callback: The callback function to remove

        Returns:
            True if callback was found and removed, False otherwise
        """
        if event in self._callbacks:
            try:
                self._callbacks[event].remove(callback)
                return True
            except ValueError:
                pass
        return False

    def has_callbacks(self, event: str) -> bool:
        """
        Report whether at least one callback is registered for an event.

        A caller that reports an event itself when nothing is listening -- as
        :class:`~cellucid.jupyter.BaseViewer` does for ``command_error`` -- uses
        this to stand aside once the notebook has taken the event over.

        Args:
            event: Event name

        Returns:
            True if at least one callback is registered for that event
        """
        return bool(self._callbacks.get(event))

    def clear(self, event: str | None = None):
        """
        Clear callbacks for an event, or all events if event is None.

        Args:
            event: Event name to clear, or None to clear all
        """
        if event is None:
            self._callbacks.clear()
        elif event in self._callbacks:
            self._callbacks[event].clear()

    def trigger(self, event: str, data: dict[str, Any]):
        """
        Trigger all callbacks for an event.

        Args:
            event: Event name to trigger
            data: Event data dict passed to callbacks
        """
        message_callbacks = list(self._callbacks.get("message", ()))
        for callback in message_callbacks:
            callback({"event": event, **data})

        if event == "message":
            return

        event_callbacks = list(self._callbacks.get(event, ()))
        for callback in event_callbacks:
            callback(data)

    def on(self, event: str) -> Callable[[HookCallback], HookCallback]:
        """
        Decorator for registering event callbacks.

        Example:
            @hooks.on('selection')
            def handle_selection(event):
                print(event['cells'])
        """

        def decorator(callback: HookCallback) -> HookCallback:
            self.register(event, callback)
            return callback

        return decorator


@dataclass
class ViewerState:
    """
    Thread-safe (read-mostly) snapshot of the latest viewer → Python events.

    This is intentionally small: it stores *the latest* payload per event type.
    For awaiting new events, use `viewer.wait_for_event(...)`.
    """

    ready: dict[str, Any] | None = None
    selection: dict[str, Any] | None = None
    hover: dict[str, Any] | None = None
    click: dict[str, Any] | None = None
    command_error: dict[str, Any] | None = None

    last_event_type: str | None = None
    last_event: dict[str, Any] | None = None
    last_updated_at: float | None = None
