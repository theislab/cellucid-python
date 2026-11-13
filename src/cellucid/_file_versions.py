"""Strong validators of a file's version, read from the open file.

One rule for every place Cellucid tells whether a file changed: the server's
``ETag`` and dataset list (:mod:`cellucid.anndata_server`) and the web asset
cache's re-check (:mod:`cellucid.web_cache`). R's ``.file_version()`` builds
the same tag.

Served files are opened with :func:`open_shared`: a file Cellucid reads or
streams can be replaced or removed meanwhile, on Windows too, where Python's
``open`` refuses another process the deletion (a person's ``write_h5ad`` to a
new file moved over the served one would fail while a response is being
sent).
"""

from __future__ import annotations

import io
import os
import sys
from pathlib import Path

if sys.platform == "win32":
    import ctypes
    import msvcrt
    from ctypes import wintypes

    class _FileBasicInfo(ctypes.Structure):
        """``FILE_BASIC_INFO``: a file's times, in 100 ns since 1601."""

        _fields_ = [
            ("CreationTime", ctypes.c_int64),
            ("LastAccessTime", ctypes.c_int64),
            ("LastWriteTime", ctypes.c_int64),
            ("ChangeTime", ctypes.c_int64),
            ("FileAttributes", wintypes.DWORD),
        ]

    _kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    _get_file_information = _kernel32.GetFileInformationByHandleEx
    _get_file_information.argtypes = (wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD)
    _get_file_information.restype = wintypes.BOOL
    _create_file = _kernel32.CreateFileW
    _create_file.argtypes = (
        wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD, ctypes.c_void_p,
        wintypes.DWORD, wintypes.DWORD, wintypes.HANDLE,
    )
    _create_file.restype = wintypes.HANDLE
    _close_handle = _kernel32.CloseHandle
    _close_handle.argtypes = (wintypes.HANDLE,)
    _close_handle.restype = wintypes.BOOL
    _GENERIC_READ = 0x80000000
    _FILE_SHARE_ALL = 0x1 | 0x2 | 0x4  # read, write, delete
    _OPEN_EXISTING = 3
    _FILE_ATTRIBUTE_NORMAL = 0x80
    _INVALID_HANDLE_VALUE = wintypes.HANDLE(-1).value

    def open_shared(path: Path) -> io.FileIO:
        """``path`` opened for reading, unbuffered, sharing read, write and
        delete: as on POSIX, the file may be rewritten, replaced or removed
        while it is open (the open file keeps reading the version it opened)."""
        handle = _create_file(
            os.fspath(path), _GENERIC_READ, _FILE_SHARE_ALL, None,
            _OPEN_EXISTING, _FILE_ATTRIBUTE_NORMAL, None,
        )
        if handle is None or handle == _INVALID_HANDLE_VALUE:
            error = ctypes.get_last_error()
            raise OSError(None, ctypes.FormatError(error).strip(), os.fspath(path), error)
        try:
            fd = msvcrt.open_osfhandle(handle, os.O_RDONLY | os.O_BINARY)
        except BaseException:
            _close_handle(handle)
            raise
        return os.fdopen(fd, "rb", buffering=0)

    def _change_time(fd: int, status: os.stat_result) -> int:
        """The status-change time of the file open as ``fd``: NTFS's change
        time, read from the handle. ``st_ctime`` is no such time on Windows:
        the creation time from ``os.stat`` (and from ``os.fstat`` before
        Python 3.12), which a rewrite in place does not change."""
        info = _FileBasicInfo()
        if not _get_file_information(msvcrt.get_osfhandle(fd), 0, ctypes.byref(info), ctypes.sizeof(info)):  # 0: FileBasicInfo
            raise ctypes.WinError(ctypes.get_last_error())
        return int(info.ChangeTime)

else:

    def _change_time(fd: int, status: os.stat_result) -> int:
        """The status-change time of the file open as ``fd``."""
        return status.st_ctime_ns

    def open_shared(path: Path) -> io.FileIO:
        """``path`` opened for reading, unbuffered (on POSIX a file open for
        reading may be rewritten, replaced or removed meanwhile)."""
        return open(path, "rb", buffering=0)  # noqa: SIM115 - the caller closes it


def file_version(fd: int) -> tuple[int, str]:
    """The size of the file open as ``fd`` and a strong validator of this
    version of it: its size, modification time, status-change time and
    inode, all read from the open file.

    A file rewritten in place gets a new modification and status-change
    time; one replaced by another (a new file moved over it, even one whose
    modification time was copied with it, ``cp -p``, ``rsync -t``) a new
    inode and status-change time. So a browser that kept bytes of the old
    version -- its byte cache, a kept gene index, a converted store, each
    filed under this tag, which ``/_cellucid/datasets`` states beside the
    file's size -- revalidates them and is told they changed, and the
    viewer, reading with ``If-Range``-free range requests, sees the new tag
    and refuses to mix versions. Two writes within the file system's
    timestamp resolution (nanoseconds on APFS, ext4, NTFS; two seconds on
    FAT) that keep the size cannot be told apart by any validator short of
    hashing the whole file.

    Every tag is read from an open file, never from a path: the bytes a
    response sends are the version its tag names, and the tag the dataset
    list states is the one the file's responses carry (on Windows a path's
    ``os.stat`` and a handle's ``os.fstat`` disagree on ``st_ctime``).
    """
    status = os.fstat(fd)
    return status.st_size, f'"{status.st_size:x}-{status.st_mtime_ns:x}-{_change_time(fd, status):x}-{status.st_ino:x}"'


def path_version(path: Path) -> tuple[int, str]:
    """The size and validator (:func:`file_version`) of the file at ``path``, read from it opened."""
    with open_shared(path) as file:
        return file_version(file.fileno())
