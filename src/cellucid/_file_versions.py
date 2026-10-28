"""Strong validators of a file's version, read from the open file.

One rule for every place Cellucid tells whether a file changed: the server's
``ETag`` and dataset list (:mod:`cellucid.anndata_server`) and the web asset
cache's re-check (:mod:`cellucid.web_cache`). R's ``.file_version()`` builds
the same tag.
"""

from __future__ import annotations

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

    _get_file_information = ctypes.WinDLL("kernel32", use_last_error=True).GetFileInformationByHandleEx
    _get_file_information.argtypes = (wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD)
    _get_file_information.restype = wintypes.BOOL

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
    with open(path, "rb", buffering=0) as file:
        return file_version(file.fileno())
