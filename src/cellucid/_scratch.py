"""
Scratch folders: where the package writes the arrays it lays out or decodes;
and the machine's memory, which sizes what is held beside them.

A served store (:mod:`cellucid._served_store`) writes its snapshots, its
gene-major layouts and the decoded copies of compressed source arrays here;
``prepare()`` decodes a compressed matrix it reads from a file
(``backed="r"``) here (:mod:`cellucid.prepare_data._gene_major`). Each folder
is made when first needed, in the temporary folder, and removed when its
owner closes it -- after every memory-mapped array kept in it is let go,
since Windows removes no file that is mapped. It names its process (a file
``owner``), so the next folder made removes those of processes that died
before removing their own (a killed kernel): :func:`remove_stale_scratch`.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import threading
from pathlib import Path
from typing import Protocol

__all__ = ["SCRATCH_PREFIX", "Scratch", "machine_memory_bytes", "remove_stale_scratch"]

#: The prefix of a scratch folder in the temporary folder.
SCRATCH_PREFIX = "cellucid-scratch-"


def machine_memory_bytes() -> int:
    """The physical memory of this machine, in bytes: ``sysconf`` on Linux
    and macOS, ``GlobalMemoryStatusEx`` on Windows; 8 GiB where neither
    answers (the device the web app takes a browser that reports nothing to
    be)."""
    try:
        if sys.platform == "win32":
            import ctypes

            class _Status(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            status = _Status()
            status.dwLength = ctypes.sizeof(_Status)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):  # type: ignore[attr-defined]
                return int(status.ullTotalPhys)
        else:
            pages = os.sysconf("SC_PHYS_PAGES")
            size = os.sysconf("SC_PAGE_SIZE")
            if pages > 0 and size > 0:
                return int(pages * size)
    except (OSError, ValueError, AttributeError):
        pass
    return 8 * 2**30


class Mapped(Protocol):
    """Something holding a memory-mapped file of a scratch folder, let go on :meth:`release`."""

    def release(self) -> None: ...


def _process_alive(pid: int) -> bool:
    """Whether process ``pid`` runs (on this machine), never signalling it
    (``os.kill`` on Windows ends a process, whatever the signal). A process
    that runs but may not be asked (another user's) counts as running."""
    if pid <= 0:
        return False
    if sys.platform == "win32":  # pragma: no cover - Windows only
        import ctypes
        from ctypes import wintypes

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
        kernel32.OpenProcess.restype = wintypes.HANDLE
        kernel32.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
        kernel32.GetExitCodeProcess.restype = wintypes.BOOL
        kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
        kernel32.CloseHandle.restype = wintypes.BOOL
        handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return ctypes.get_last_error() == 5  # ERROR_ACCESS_DENIED: it runs; no such process otherwise
        try:
            code = wintypes.DWORD()
            return bool(kernel32.GetExitCodeProcess(handle, ctypes.byref(code))) and code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def remove_stale_scratch(temporary: Path | None = None) -> list[Path]:
    """Remove the scratch folders whose process is gone (a kernel killed
    before it closed them), and return them. A folder whose owner runs, or
    names none it can read, is left alone."""
    removed = []
    for folder in Path(temporary or tempfile.gettempdir()).glob(f"{SCRATCH_PREFIX}*"):
        try:
            pid = int((folder / "owner").read_text().strip())
        except (OSError, ValueError):
            continue
        if pid != os.getpid() and not _process_alive(pid):
            shutil.rmtree(folder, ignore_errors=True)
            removed.append(folder)
    return removed


class Scratch:
    """One scratch folder, made when first needed and removed by :meth:`close`.

    ``stop``, when set, refuses a new folder (its owner is closing). What maps
    a file of the folder is kept (:meth:`keep`) and let go first on close.
    """

    def __init__(self, stop: threading.Event | None = None) -> None:
        self._stop = stop
        self._lock = threading.Lock()
        self._kept: list[Mapped] = []
        self.folder: Path | None = None

    def path(self) -> Path:
        """The folder, made now if it was not."""
        with self._lock:
            if self._stop is not None and self._stop.is_set():
                raise RuntimeError("The served store is closed.")
            if self.folder is None:
                remove_stale_scratch()
                folder = Path(tempfile.mkdtemp(prefix=SCRATCH_PREFIX))
                (folder / "owner").write_text(str(os.getpid()))
                self.folder = folder
            return self.folder

    def keep(self, mapped: Mapped) -> None:
        """Let ``mapped``'s memory-mapped file go when the folder is closed."""
        with self._lock:
            self._kept.append(mapped)

    def release(self) -> None:
        """Let every memory-mapped array kept here go (each one's own lock taken
        outside this one, which it may take inside its own)."""
        with self._lock:
            kept = list(self._kept)
        for mapped in kept:
            mapped.release()

    def close(self) -> None:
        """Let every mapped array go, then remove the folder."""
        self.release()
        with self._lock:
            folder, self.folder = self.folder, None
        if folder is not None:
            shutil.rmtree(folder, ignore_errors=True)
