# Cellucid Hooks System - Development Guide

This document explains the communication system between the Python
viewer and the JavaScript frontend. Use this guide when modifying or extending
the hooks system.

## Environment Support

| Direction | Colab | Jupyter | JupyterLab | VSCode |
|-----------|-------|---------|------------|--------|
| **Python → Frontend** | ✅ | ✅ | ✅ | ✅ |
| `highlight_cells()` | ✅ | ✅ | ✅ | ✅ |
| `set_color_by()` | ✅ | ✅ | ✅ | ✅ |
| `send_message()` | ✅ | ✅ | ✅ | ✅ |
| **Frontend → Python** | ✅ | ✅ | ✅ | ✅ |
| `on_selection` hook | ✅ | ✅ | ✅ | ✅ |
| `on_hover` hook | ✅ | ✅ | ✅ | ✅ |
| `on_click` hook | ✅ | ✅ | ✅ | ✅ |

**Both directions work in all environments!**

Frontend → Python works via HTTP POST to the local data server. The viewer iframe
on `cellucid.com` POSTs events to `/_cellucid/events` on the Python server, which
has CORS headers allowing cross-origin requests.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              JUPYTER NOTEBOOK                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────────────┐         ┌──────────────────────────────────────┐ │
│  │   Python Viewer      │         │   Embedded iframe (cellucid.com)     │ │
│  │   + Data Server      │         │                                      │ │
│  │                      │         │  ┌─────────────────────────────────┐│ │
│  │  ┌────────────────┐  │         │  │  JupyterBridgeDataSource        ││ │
│  │  │  HookRegistry  │  │  HTTP   │  │                                 ││ │
│  │  │                │◄─┼─────────┼──│  _postEventToPython(event)      ││ │
│  │  │ - selection    │  │  POST   │  │                                 ││ │
│  │  │ - hover        │  │  to     │  │  Sends (via POST):              ││ │
│  │  │ - click        │  │/_cellucid│  │  - selection events             ││ │
│  │  │ - ready        │  │/events  │  │  - hover events                 ││ │
│  │  │ - message      │  │         │  │  - click events                 ││ │
│  │  └────────────────┘  │         │  │  - ready event                  ││ │
│  │                      │         │  │                                 ││ │
│  │  send_message()──────┼─────────┼─►│  Receives (via postMessage):    ││ │
│  │                      │  POST   │  │  - highlight commands           ││ │
│  │                      │  MSG    │  │  - setColorBy commands          ││ │
│  │                      │         │  │  - setVisibility commands       ││ │
│  └──────────────────────┘         │  └─────────────────────────────────┘│ │
│                                    └──────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Safety Guarantees

The hooks system is designed to be **non-breaking** when not in Jupyter context:

### JavaScript Side
- All `notify*` methods check `if (!this._connected) return;` first
- Using optional chaining (`jupyterSource?.notifyClick(...)`) is safe
- `_postToParent` checks if we're in an iframe before sending
- The `isJupyterContext()` helper returns `false` outside Jupyter

### Python Side
- `HookRegistry.trigger()` catches exceptions in callbacks
- `_handle_frontend_message()` gracefully handles unknown message types
- Hooks with no registered callbacks simply do nothing

### Example: Safe Event Wiring
```javascript
// This code runs safely in both Jupyter and standalone contexts:
function handleCellClick(cellIndex, event) {
    // Business logic always runs
    highlightCell(cellIndex);
    showCellInfo(cellIndex);

    // Jupyter notification - silently no-ops if not in Jupyter
    jupyterSource?.notifyClick(cellIndex, {
        button: event.button,
        shift: event.shiftKey
    });
}
```

## File Locations

### Python Side (cellucid-python)

| File | Purpose |
|------|---------|
| `src/cellucid/jupyter.py` | Main hooks implementation, viewer classes |
| `src/cellucid/_server_base.py` | HTTP event routing infrastructure |
| `src/cellucid/server.py` | CellucidServer with events endpoint |
| `src/cellucid/anndata_server.py` | AnnDataServer with events endpoint |

Key classes:
- `HookRegistry` - Manages callback registration and triggering
- `BaseViewer` - Contains hook decorators (`on_selection`, `on_hover`, etc.)
- `_event_callbacks` - Global dict mapping viewer IDs to event callbacks
- `route_event()` - Routes HTTP POSTed events to viewer callbacks

### JavaScript Side (cellucid)

| File | Purpose |
|------|---------|
| `assets/js/data/jupyter-source.js` | `JupyterBridgeDataSource` class |
| `assets/js/app/state.js` | `DataState` class (highlight state) |
| `assets/js/app/main.js` | Wires up event listeners |

## Adding a New Hook

### Step 1: Define the Hook in Python

In `jupyter.py`, add a new property decorator to `BaseViewer`:

```python
@property
def on_my_new_event(self) -> Callable[[HookCallback], HookCallback]:
    """
    Decorator to register a handler for my_new_event.

    Called when [describe when this fires].

    Event data:
        - field1: type - description
        - field2: type - description

    Example:
        @viewer.on_my_new_event
        def handle(event):
            print(event['field1'])
    """
    return self._hooks.on('my_new_event')
```

### Step 2: Add to Event Map (if needed)

If the frontend sends a different message type name, add mapping in
`_handle_frontend_message`:

```python
def _handle_frontend_message(self, message: dict):
    msg_type = message.get('type', '')

    event_map = {
        'selection': 'selection',
        'hover': 'hover',
        'click': 'click',
        'ready': 'ready',
        'frontendMsgType': 'my_new_event',  # Add mapping here
    }
    # ...
```

### Step 3: Send Event from Frontend

In `jupyter-source.js`, the notification methods are already implemented:

```javascript
// Available methods in JupyterBridgeDataSource:
notifySelection(cellIndices, source)  // User selected cells
notifyHover(cellIndex, position)      // User hovering (debounced 50ms)
notifyClick(cellIndex, options)       // User clicked cell
notifyReady(info)                     // Dataset loaded
notifyCustomEvent(eventType, data)    // Custom events
```

To add a new notification method:

```javascript
// In JupyterBridgeDataSource class
notifyMyNewEvent(data) {
    if (!this._connected) return;  // IMPORTANT: Guard against no-Jupyter context
    const event = {
        type: 'my_new_event',
        viewerId: this._config?.viewerId,
        field1: data.field1,
        field2: data.field2
    };
    this._postEventToPython(event);  // HTTP POST to Python server
}
```

### Step 4: Wire Up in Frontend

In `main.js` or `ui.js`, call the notify methods. Always check if jupyterSource exists:

```javascript
// Get reference to jupyter source (set up during initialization)
let jupyterSource = null;

// During app initialization:
if (isJupyterContext()) {
    jupyterSource = createJupyterBridgeDataSource();
    await jupyterSource.initialize();
}

// Wire up events (safe even if jupyterSource is null):
canvas.addEventListener('click', (e) => {
    const cellIndex = pickCellAtPosition(e.clientX, e.clientY);
    if (cellIndex !== null) {
        jupyterSource?.notifyClick(cellIndex, {
            button: e.button,
            shift: e.shiftKey,
            ctrl: e.ctrlKey || e.metaKey
        });
    }
});

// For hover events:
canvas.addEventListener('mousemove', (e) => {
    const cellIndex = pickCellAtPosition(e.clientX, e.clientY);
    const position = getWorldPosition(e.clientX, e.clientY);
    jupyterSource?.notifyHover(cellIndex, position);
});

// When selection changes (e.g., after lasso tool):
function onSelectionComplete(selectedCells, method) {
    jupyterSource?.notifySelection(selectedCells, method);
}

// When dataset finishes loading:
function onDatasetLoaded(metadata) {
    jupyterSource?.notifyReady({
        nCells: metadata.stats?.n_cells || 0,
        dimensions: metadata.dimensions || 3
    });
}
```

## Event Data Schemas

### selection
```python
{
    'cells': list[int],      # Indices of selected cells
    'source': str,           # 'lasso' | 'click' | 'range'
}
```

### hover
```python
{
    'cell': int | None,      # Cell index, None if not hovering over cell
    'position': {            # World coordinates
        'x': float,
        'y': float,
        'z': float
    }
}
```

### click
```python
{
    'cell': int,             # Clicked cell index
    'button': int,           # 0=left, 1=middle, 2=right
    'shift': bool,           # Shift key held
    'ctrl': bool,            # Ctrl/Cmd key held
}
```

### ready
```python
{
    'n_cells': int,          # Number of cells in dataset
    'dimensions': int,       # Embedding dimensionality (2 or 3)
}
```

## Message Routing

### Python → Frontend (send_message)

```python
viewer.send_message({'type': 'highlight', 'cells': [1,2,3], 'color': '#ff0000'})
```

Flow:
1. `BaseViewer.send_message()` called
2. Executes JavaScript via `IPython.display.Javascript()`
3. JS calls `window.cellucidViewers[viewerId].sendMessage()`
4. `postMessage()` to iframe with origin `https://www.cellucid.com`
5. `JupyterBridgeDataSource._handleMessage()` receives it

### Frontend → Python (hooks) - All Environments

Flow (HTTP POST):
1. Frontend calls `jupyterSource._postEventToPython(event)`
2. HTTP POST to `${serverUrl}/_cellucid/events`
3. Data server receives the POST request
4. `handle_event_post()` parses JSON and extracts viewerId
5. `route_event()` finds the viewer callback by ID
6. `viewer._handle_frontend_message()` called
7. `HookRegistry.trigger()` fires callbacks

This works because:
- The iframe on cellucid.com can make HTTP requests to the Python server
- The Python server has CORS headers allowing cross-origin requests
- No special environment APIs are needed - just standard HTTP

## Testing Hooks

### Manual Testing (All Environments)

```python
# In any Jupyter notebook (Jupyter, JupyterLab, Colab, VSCode)
viewer = show_anndata(adata)

@viewer.on_message
def debug_all(event):
    print(f"[DEBUG] {event}")

# Now interact with viewer - all events will print
```

Works in all environments via HTTP POST to the data server.

### Simulating Events

```python
# Directly trigger hooks for testing
viewer._handle_frontend_message({
    'type': 'selection',
    'cells': [1, 2, 3],
    'source': 'test'
})
```

## Common Modifications

### Adding a New Send Command

To add a new command from Python to frontend:

1. **Python**: Just use `send_message()`
   ```python
   viewer.send_message({'type': 'myCommand', 'param': value})
   ```

2. **Frontend**: Handle in `JupyterBridgeDataSource._handleMessage()`
   ```javascript
   _handleMessage(event) {
       // ...
       if (event.data.type === 'myCommand') {
           this._handleMyCommand(event.data);
       }
   }

   _handleMyCommand(data) {
       // Do something with data.param
   }
   ```

### Changing Event Data Schema

1. Update the frontend to send new fields
2. Update the Python docstring to document new fields
3. Update this documentation

### Debouncing High-Frequency Events

For events like `hover` that fire rapidly:

```javascript
// In jupyter-source.js
constructor() {
    // ...
    this._hoverDebounceTimer = null;
}

notifyHover(cell, position) {
    clearTimeout(this._hoverDebounceTimer);
    this._hoverDebounceTimer = setTimeout(() => {
        const event = {
            type: 'hover',
            viewerId: this._config?.viewerId,
            cell: cell,
            position: position
        };
        this._postEventToPython(event);
    }, 50);  // 50ms debounce
}
```

## Troubleshooting

### Hook not firing

1. **Is the data server running?** Check that the viewer displayed without errors
2. Check browser console for HTTP POST errors to `/_cellucid/events`
3. Verify `viewerId` matches between Python and frontend
4. Check CORS errors in browser console (should not occur with default setup)
5. Try `viewer._handle_frontend_message({'type': 'selection', 'cells': [1,2,3]})` to test hooks directly

### Messages not reaching frontend

1. Verify viewer is displayed (`viewer._displayed == True`)
2. Check iframe is loaded (network tab)
3. Verify postMessage origin is correct

### Multiple viewers interfering

Each viewer has unique `_viewer_id` (hex token). Messages are routed by ID.
Check that frontend includes `viewerId` in all messages.

## HookRegistry API Reference

```python
class HookRegistry:
    def register(self, event: str, callback: Callable[[dict], None]) -> Callable
        """Register a callback for an event."""

    def unregister(self, event: str, callback: Callable) -> bool
        """Remove a callback. Returns True if found."""

    def clear(self, event: str | None = None)
        """Clear callbacks for event, or all if event is None."""

    def trigger(self, event: str, data: dict)
        """Fire all callbacks for event with data."""

    def on(self, event: str) -> Callable
        """Decorator factory for registering callbacks."""
```

## BaseViewer API Reference

```python
class BaseViewer:
    # ========================================
    # HOOKS (Frontend → Python) - All environments
    # ========================================

    # Decorator properties
    @property on_selection -> decorator
    @property on_hover -> decorator
    @property on_click -> decorator
    @property on_ready -> decorator
    @property on_message -> decorator  # Catches all events

    # Programmatic registration
    def register_hook(event: str, callback) -> callback
    def unregister_hook(event: str, callback) -> bool
    def clear_hooks(event: str | None = None)

    # ========================================
    # COMMANDS (Python → Frontend) - All environments
    # ========================================

    # Convenience methods
    def highlight_cells(cell_indices: list[int], color: str = "#ff0000")
    def clear_highlights()
    def set_color_by(field: str)
    def set_visibility(cell_indices: list[int] | None, visible: bool)
    def reset_view()

    # Low-level
    def send_message(message: dict)  # Python → Frontend
    def _handle_frontend_message(message: dict)  # Frontend → Python (internal)
```
