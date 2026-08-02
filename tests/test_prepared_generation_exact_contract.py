from __future__ import annotations

import inspect
import json
import os
import queue
import subprocess
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.cli import _detect_data_format
from cellucid.prepare_data import prepare
from cellucid.server import CORSRequestHandler


def _prepare_kwargs(out_dir: Path, *, dimensions: int) -> dict[str, object]:
    if dimensions == 2:
        embedding = np.array(
            [
                [-3.0, 1.0],
                [0.5, 5.0],
                [8.0, -2.0],
            ],
            dtype=np.float32,
        )
    elif dimensions == 3:
        embedding = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 2.0],
                [0.0, 1.0, 4.0],
            ],
            dtype=np.float32,
        )
    else:
        raise AssertionError("test fixture supports only 2D and 3D")

    return {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        f"X_umap_{dimensions}d": embedding,
        "out_dir": out_dir,
        "dataset_name": "Atomic generation",
        "dataset_id": "atomic-generation",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }


def _snapshot(directory: Path) -> dict[str, bytes]:
    return {
        path.relative_to(directory).as_posix(): path.read_bytes()
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def _assert_only_persistent_lock(parent: Path, target_name: str) -> None:
    lock_path = parent / f".{target_name}.cellucid.lock"
    assert sorted(parent.glob(f".{target_name}.cellucid*")) == [lock_path]
    assert lock_path.is_file()
    assert not lock_path.is_symlink()
    assert lock_path.stat().st_size == 0


@contextmanager
def _live_export_lock(target: Path) -> Iterator[None]:
    process = subprocess.Popen(
        [
            sys.executable,
            "-c",
            (
                "import sys\n"
                "from pathlib import Path\n"
                "from cellucid.prepare_data._locking import _exclusive_export_generation\n"
                "with _exclusive_export_generation(Path(sys.argv[1])):\n"
                "    print('READY', flush=True)\n"
                "    sys.stdin.read(1)\n"
            ),
            str(target),
        ],
        env=os.environ.copy(),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    ready_lines: queue.Queue[str] = queue.Queue()
    reader = threading.Thread(
        target=lambda: ready_lines.put(process.stdout.readline()),
        daemon=True,
    )
    reader.start()
    try:
        # Importing the scientific stack in a fresh interpreter can exceed
        # five seconds on a loaded CI runner. This deadline covers process
        # startup only; the lock-contender assertion below remains bounded by
        # its independent ten-second timeout.
        ready = ready_lines.get(timeout=15).strip()
    except queue.Empty:
        process.kill()
        _stdout, stderr = process.communicate(timeout=5)
        raise AssertionError(f"lock owner did not become ready: {stderr}") from None
    if ready != "READY":
        _stdout, stderr = process.communicate(timeout=5)
        raise AssertionError(f"lock owner did not start: {stderr}")
    try:
        yield
    finally:
        assert process.stdin is not None
        process.stdin.write("\n")
        process.stdin.flush()
        _stdout, stderr = process.communicate(timeout=5)
        reader.join(timeout=5)
        assert process.returncode == 0, stderr


def _attempt_prepare_in_subprocess(target: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys\n"
                "from pathlib import Path\n"
                "import numpy as np\n"
                "import pandas as pd\n"
                "from cellucid.prepare_data import prepare\n"
                "embedding = np.array([[-3.0, 1.0], [0.5, 5.0], [8.0, -2.0]], "
                "dtype=np.float32)\n"
                "try:\n"
                "    prepare(\n"
                "        latent_space=embedding.copy(),\n"
                "        obs=pd.DataFrame({'score': [0.25, 0.5, 0.75]}),\n"
                "        X_umap_2d=embedding,\n"
                "        out_dir=Path(sys.argv[1]),\n"
                "        dataset_name='Atomic generation',\n"
                "        dataset_id='atomic-generation',\n"
                "        obs_categorical_dtype='uint16',\n"
                "        centroid_min_points=1,\n"
                "        force=True,\n"
                "    )\n"
                "except RuntimeError as error:\n"
                "    if 'generation is already active' in str(error):\n"
                "        raise SystemExit(73) from error\n"
                "    raise\n"
            ),
            str(target),
        ],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )


def _attempt_lock_in_subprocess(target: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys\n"
                "from pathlib import Path\n"
                "from cellucid.prepare_data._locking import _exclusive_export_generation\n"
                "try:\n"
                "    with _exclusive_export_generation(Path(sys.argv[1])):\n"
                "        raise SystemExit(3)\n"
                "except RuntimeError as error:\n"
                "    if 'generation is already active' in str(error):\n"
                "        raise SystemExit(0) from error\n"
                "    raise\n"
            ),
            str(target),
        ],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )


def _leave_dead_export_lock(target: Path) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import os, sys\n"
                "from pathlib import Path\n"
                "from cellucid.prepare_data._locking import _exclusive_export_generation\n"
                "owner = _exclusive_export_generation(Path(sys.argv[1]))\n"
                "owner.__enter__()\n"
                "os._exit(73)\n"
            ),
            str(target),
        ],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 73, completed.stderr


def _crash_export_at_rename(target: Path, phase: str) -> None:
    environment = os.environ.copy()
    environment["PYTHONIOENCODING"] = "cp1252"
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import importlib, os, sys\n"
                "from pathlib import Path\n"
                "import numpy as np\n"
                "import pandas as pd\n"
                "module = importlib.import_module('cellucid.prepare_data')\n"
                "transaction = importlib.import_module('cellucid.prepare_data._transaction')\n"
                "target = Path(sys.argv[1])\n"
                "phase = sys.argv[2]\n"
                "real_rename = transaction._rename_export_path\n"
                "def crashing_rename(source, destination):\n"
                "    source = Path(source)\n"
                "    destination = Path(destination)\n"
                "    is_journal = destination.name.endswith('.cellucid-transaction.json')\n"
                "    is_stage_publish = (\n"
                "        source.name.startswith(f'.{target.name}.cellucid-stage-')\n"
                "        and destination == target\n"
                "    )\n"
                "    if phase == 'before-stage-publish' and is_stage_publish:\n"
                "        os._exit(73)\n"
                "    real_rename(source, destination)\n"
                "    is_prior_move = (\n"
                "        source == target\n"
                "        and destination.name.startswith(f'.{target.name}.cellucid-backup-')\n"
                "    )\n"
                "    if phase == 'journal-published' and is_journal:\n"
                "        os._exit(73)\n"
                "    if phase == 'prior-moved' and is_prior_move:\n"
                "        os._exit(73)\n"
                "    if phase == 'stage-published' and is_stage_publish:\n"
                "        os._exit(73)\n"
                "transaction._rename_export_path = crashing_rename\n"
                "embedding = np.array(\n"
                "    [[-3.0, 1.0], [0.5, 5.0], [8.0, -2.0]],\n"
                "    dtype=np.float32,\n"
                ")\n"
                "module.prepare(\n"
                "    latent_space=embedding.copy(),\n"
                "    obs=pd.DataFrame({'score': [0.25, 0.5, 0.75]}),\n"
                "    X_umap_2d=embedding,\n"
                "    out_dir=target,\n"
                "    dataset_name='Atomic generation',\n"
                "    dataset_id='atomic-generation',\n"
                "    created_at='2026-01-02T03:04:05Z',\n"
                "    obs_categorical_dtype='uint16',\n"
                "    centroid_min_points=1,\n"
                "    force=True,\n"
                ")\n"
                "raise SystemExit(4)\n"
            ),
            str(target),
            phase,
        ],
        env=environment,
        capture_output=True,
        encoding="cp1252",
        timeout=20,
        check=False,
    )
    assert completed.returncode == 73, completed.stderr


def _write_exported_dataset(directory: Path, dataset_id: str) -> None:
    directory.mkdir()
    (directory / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": dataset_id,
                "name": dataset_id,
            }
        ),
        encoding="utf-8",
    )
    (directory / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (directory / "points_2d.bin").write_bytes(b"\x00" * 8)


def _write_zarr_v2_root(directory: Path) -> None:
    directory.mkdir()
    (directory / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (directory / ".zattrs").write_text("{}", encoding="utf-8")


def test_prepare_exposes_only_canonical_artifact_names() -> None:
    parameters = inspect.signature(prepare).parameters
    assert "obs_manifest_filename" not in parameters
    assert "obs_binary_dirname" not in parameters
    assert "var_manifest_filename" not in parameters
    assert "var_binary_dirname" not in parameters


def test_force_replacement_removes_every_stale_capability(tmp_path: Path) -> None:
    output = tmp_path / "generation"
    initial = _prepare_kwargs(output, dimensions=3)
    initial["var"] = pd.DataFrame(index=["GeneA", "GeneB"])
    initial["gene_expression"] = np.array(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        dtype=np.float32,
    )
    initial["connectivities"] = np.array(
        [[0.0, 0.5, 0.0], [0.5, 0.0, 1.5], [0.0, 1.5, 0.0]],
        dtype=np.float64,
    )
    prepare(**initial)

    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    prepare(**replacement)

    assert (output / "points_2d.bin").is_file()
    assert not (output / "points_3d.bin").exists()
    assert not (output / "var_manifest.json").exists()
    assert not (output / "var").exists()
    assert not (output / "connectivity_manifest.json").exists()
    assert not (output / "connectivity").exists()
    identity = json.loads((output / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["embeddings"]["available_dimensions"] == [2]
    assert identity["stats"]["n_genes"] == 0
    assert identity["stats"]["has_connectivity"] is False


def test_failed_force_generation_preserves_the_prior_generation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _generation as prepare_module

    output = tmp_path / "generation"
    prepare(**_prepare_kwargs(output, dimensions=3))
    original = _snapshot(output)
    real_write_binary = prepare_module._write_binary
    calls = 0

    def fail_during_generation(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("synthetic staged write failure")
        return real_write_binary(*args, **kwargs)

    monkeypatch.setattr(prepare_module, "_write_binary", fail_during_generation)
    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    with pytest.raises(OSError, match="synthetic staged write failure"):
        prepare(**replacement)

    assert _snapshot(output) == original
    _assert_only_persistent_lock(tmp_path, output.name)


def test_failed_initial_generation_leaves_no_partial_target(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _generation as prepare_module

    output = tmp_path / "generation"

    def fail_write(*_args, **_kwargs):
        raise OSError("synthetic initial write failure")

    monkeypatch.setattr(prepare_module, "_write_binary", fail_write)
    with pytest.raises(OSError, match="synthetic initial write failure"):
        prepare(**_prepare_kwargs(output, dimensions=2))

    assert not output.exists()
    _assert_only_persistent_lock(tmp_path, output.name)


def test_concurrent_generation_for_the_same_target_is_rejected(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    prepare(**_prepare_kwargs(output, dimensions=3))
    original = _snapshot(output)

    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    with _live_export_lock(output):
        contender = _attempt_prepare_in_subprocess(output)
        assert contender.returncode == 73, contender.stderr
        assert _snapshot(output) == original

    prepare(**replacement)
    assert (output / "points_2d.bin").is_file()
    assert not (output / "points_3d.bin").exists()


def test_dead_export_owner_does_not_block_the_next_generation(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    _leave_dead_export_lock(output)

    prepare(**_prepare_kwargs(output, dimensions=2))

    assert (output / "points_2d.bin").is_file()
    assert (
        json.loads((output / "dataset_identity.json").read_text(encoding="utf-8"))["id"]
        == "atomic-generation"
    )


def test_legacy_pid_lock_file_does_not_claim_live_ownership(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    lock_path = tmp_path / ".generation.cellucid.lock"
    legacy_contents = b"12345\n"
    lock_path.write_bytes(legacy_contents)

    prepare(**_prepare_kwargs(output, dimensions=2))

    assert (output / "points_2d.bin").is_file()
    assert lock_path.read_bytes() == legacy_contents


def test_same_process_lock_is_not_reentrant_and_distinct_targets_are_independent(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data._locking import _exclusive_export_generation

    first = tmp_path / "first"
    second = tmp_path / "second"
    with _exclusive_export_generation(first):
        with (
            pytest.raises(RuntimeError, match="generation.*already active"),
            _exclusive_export_generation(first),
        ):
            raise AssertionError("same-target lock unexpectedly entered")
        with _exclusive_export_generation(second):
            pass

    _assert_only_persistent_lock(tmp_path, first.name)
    _assert_only_persistent_lock(tmp_path, second.name)


def test_case_alias_cannot_release_a_same_process_lock(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data._locking import _exclusive_export_generation

    canonical_target = tmp_path / "CaseSensitiveTarget"
    alias_target = tmp_path / "casesensitivetarget"
    alias_lock_path = tmp_path / ".casesensitivetarget.cellucid.lock"

    with _exclusive_export_generation(canonical_target):
        if not alias_lock_path.exists():
            pytest.skip("the test filesystem is case-sensitive")
        with (
            pytest.raises(RuntimeError, match="generation.*already active"),
            _exclusive_export_generation(alias_target),
        ):
            raise AssertionError("case alias lock unexpectedly entered")
        contender = _attempt_lock_in_subprocess(alias_target)
        assert contender.returncode == 0, contender.stderr

    with _exclusive_export_generation(alias_target):
        pass


def test_unsafe_lock_nodes_are_retained_and_do_not_mutate_their_targets(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data._locking import _exclusive_export_generation

    directory_target = tmp_path / "directory-target"
    directory_lock = tmp_path / ".directory-target.cellucid.lock"
    directory_lock.mkdir()
    with (
        pytest.raises(RuntimeError, match="regular non-symbolic file"),
        _exclusive_export_generation(directory_target),
    ):
        raise AssertionError("directory lock unexpectedly entered")
    assert directory_lock.is_dir()
    directory_lock.rmdir()
    with _exclusive_export_generation(directory_target):
        pass

    symlink_target = tmp_path / "symlink-target"
    symlink_lock = tmp_path / ".symlink-target.cellucid.lock"
    victim = tmp_path / "lock-victim"
    victim_contents = b"must remain unchanged"
    victim.write_bytes(victim_contents)
    try:
        symlink_lock.symlink_to(victim)
    except OSError as error:
        pytest.skip(f"symbolic links are unavailable: {error}")
    with (
        pytest.raises(RuntimeError, match="symbolic link"),
        _exclusive_export_generation(symlink_target),
    ):
        raise AssertionError("symbolic lock unexpectedly entered")
    assert symlink_lock.is_symlink()
    assert victim.read_bytes() == victim_contents
    symlink_lock.unlink()
    with _exclusive_export_generation(symlink_target):
        pass

    linked_target = tmp_path / "linked-target"
    linked_alias_target = tmp_path / "linked-alias-target"
    linked_lock = tmp_path / ".linked-target.cellucid.lock"
    linked_alias_lock = tmp_path / ".linked-alias-target.cellucid.lock"
    with _exclusive_export_generation(linked_target):
        try:
            os.link(linked_lock, linked_alias_lock)
        except OSError as error:
            pytest.skip(f"hard links are unavailable: {error}")
        with (
            pytest.raises(RuntimeError, match="non-linked regular"),
            _exclusive_export_generation(linked_alias_target),
        ):
            raise AssertionError("hard-linked lock unexpectedly entered")
        linked_alias_lock.unlink()
        contender = _attempt_lock_in_subprocess(linked_target)
        assert contender.returncode == 0, contender.stderr


def test_existing_lock_replaced_during_open_is_not_adopted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _locking as prepare_module

    target = tmp_path / "generation"
    lock_path = tmp_path / ".generation.cellucid.lock"
    lock_path.write_bytes(b"")
    real_open = prepare_module._open_export_lock_descriptor
    create_modes: list[bool] = []
    attempts = 16

    for attempt in range(attempts):
        replacement = tmp_path / f".generation.cellucid.lock.replacement-{attempt}"
        replacement.write_bytes(b"")
        replaced = False

        def replace_before_open(
            path: Path,
            *,
            create: bool,
            expected_stat: os.stat_result | None,
            replacement_path: Path = replacement,
        ):
            nonlocal replaced
            create_modes.append(create)
            if path == lock_path and not create and not replaced:
                os.replace(replacement_path, path)
                replaced = True
            return real_open(
                path,
                create=create,
                expected_stat=expected_stat,
            )

        with monkeypatch.context() as race_patch:
            race_patch.setattr(
                prepare_module,
                "_open_export_lock_descriptor",
                replace_before_open,
            )
            with (
                pytest.raises(RuntimeError, match="changed while establishing ownership"),
                prepare_module._exclusive_export_generation(target),
            ):
                raise AssertionError("recreated lock unexpectedly entered")
        assert replaced
        assert not replacement.exists()

    assert create_modes == [False] * attempts
    assert lock_path.is_file()
    assert not target.exists()
    with prepare_module._exclusive_export_generation(target):
        pass
    _assert_only_persistent_lock(tmp_path, target.name)


def test_same_inode_lock_mutation_during_open_is_not_adopted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _locking as prepare_module

    target = tmp_path / "generation"
    lock_path = tmp_path / ".generation.cellucid.lock"
    lock_path.write_bytes(b"")
    inspected_stat = lock_path.lstat()
    inspected_generation = prepare_module._export_lock_generation(inspected_stat)
    real_open = prepare_module._open_export_lock_descriptor
    mutated = False

    def mutate_before_open(
        path: Path,
        *,
        create: bool,
        expected_stat: os.stat_result | None,
    ):
        nonlocal mutated
        if path == lock_path and not create and not mutated:
            path.write_bytes(b"same inode, different generation")
            mutated = True
            assert path.stat().st_dev == inspected_stat.st_dev
            assert path.stat().st_ino == inspected_stat.st_ino
            assert prepare_module._export_lock_generation(path.stat()) != inspected_generation
        return real_open(
            path,
            create=create,
            expected_stat=expected_stat,
        )

    with monkeypatch.context() as race_patch:
        race_patch.setattr(
            prepare_module,
            "_open_export_lock_descriptor",
            mutate_before_open,
        )
        with (
            pytest.raises(RuntimeError, match="changed while establishing ownership"),
            prepare_module._exclusive_export_generation(target),
        ):
            raise AssertionError("mutated lock unexpectedly entered")

    assert mutated
    assert lock_path.read_bytes() == b"same inode, different generation"
    lock_path.write_bytes(b"")
    with prepare_module._exclusive_export_generation(target):
        pass
    _assert_only_persistent_lock(tmp_path, target.name)


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows denies replacing the opened lock path at this boundary",
)
def test_lock_replaced_after_os_acquisition_is_rejected_before_entry(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _locking as prepare_module

    target = tmp_path / "generation"
    lock_path = tmp_path / ".generation.cellucid.lock"
    replacement = tmp_path / ".generation.cellucid.lock.replacement"
    lock_path.write_bytes(b"")
    replacement.write_bytes(b"")
    real_acquire = prepare_module._acquire_export_lock_descriptor
    replaced = False

    def replace_after_acquire(descriptor: int) -> None:
        nonlocal replaced
        real_acquire(descriptor)
        os.replace(replacement, lock_path)
        replaced = True

    with monkeypatch.context() as race_patch:
        race_patch.setattr(
            prepare_module,
            "_acquire_export_lock_descriptor",
            replace_after_acquire,
        )
        with (
            pytest.raises(RuntimeError, match="changed while establishing ownership"),
            prepare_module._exclusive_export_generation(target),
        ):
            raise AssertionError("replaced lock unexpectedly entered")

    assert replaced
    assert lock_path.is_file()
    assert not replacement.exists()
    with prepare_module._exclusive_export_generation(target):
        pass
    _assert_only_persistent_lock(tmp_path, target.name)


def test_lock_cleanup_preserves_cancellation_and_body_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.prepare_data import _locking as prepare_module

    cancellation_target = tmp_path / "cancellation-target"

    def interrupt_unlock(_descriptor: int) -> None:
        raise KeyboardInterrupt

    with monkeypatch.context() as cleanup_patch:
        cleanup_patch.setattr(
            prepare_module,
            "_release_export_lock_descriptor",
            interrupt_unlock,
        )
        with (
            pytest.raises(KeyboardInterrupt),
            prepare_module._exclusive_export_generation(cancellation_target),
        ):
            pass
    with prepare_module._exclusive_export_generation(cancellation_target):
        pass

    body_error_target = tmp_path / "body-error-target"

    def fail_unlock(_descriptor: int) -> None:
        raise OSError("synthetic unlock failure")

    with monkeypatch.context() as cleanup_patch:
        cleanup_patch.setattr(
            prepare_module,
            "_release_export_lock_descriptor",
            fail_unlock,
        )
        with (
            pytest.raises(ValueError, match="body failure") as captured,
            prepare_module._exclusive_export_generation(body_error_target),
        ):
            raise ValueError("body failure")
    assert any(
        "synthetic unlock failure" in note for note in getattr(captured.value, "__notes__", ())
    )
    with prepare_module._exclusive_export_generation(body_error_target):
        pass


@pytest.mark.skipif(not hasattr(os, "fork"), reason="POSIX fork is unavailable")
def test_forked_child_closes_inherited_descriptors_and_owns_its_new_lock(
    tmp_path: Path,
) -> None:
    target = tmp_path / "generation"
    environment = os.environ.copy()
    environment["PYTHONWARNINGS"] = "error::DeprecationWarning"
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import errno, os, sys\n"
                "from pathlib import Path\n"
                "from cellucid.prepare_data import _locking as prepare_module\n"
                "target = Path(sys.argv[1])\n"
                "trigger_read, trigger_write = os.pipe()\n"
                "result_read, result_write = os.pipe()\n"
                "release_read, release_write = os.pipe()\n"
                "role = 'parent'\n"
                "child_owner = None\n"
                "with prepare_module._exclusive_export_generation(target):\n"
                "    inherited_descriptor = next(\n"
                "        iter(prepare_module._EXPORT_LOCK_REGISTRY.values())\n"
                "    )[1]\n"
                "    child = os.fork()\n"
                "    if child == 0:\n"
                "        role = 'child'\n"
                "        os.close(trigger_write)\n"
                "        os.close(result_read)\n"
                "        os.close(release_write)\n"
                "        try:\n"
                "            os.fstat(inherited_descriptor)\n"
                "        except OSError as error:\n"
                "            if error.errno != errno.EBADF:\n"
                "                os.write(result_write, b'F')\n"
                "                os._exit(4)\n"
                "        else:\n"
                "            os.write(result_write, b'D')\n"
                "            os._exit(5)\n"
                "        os.read(trigger_read, 1)\n"
                "        child_owner = prepare_module._exclusive_export_generation(target)\n"
                "        child_owner.__enter__()\n"
                "    else:\n"
                "        os.close(trigger_read)\n"
                "        os.close(result_write)\n"
                "        os.close(release_read)\n"
                "if role == 'child':\n"
                "    os.write(result_write, b'A')\n"
                "    os.read(release_read, 1)\n"
                "    child_owner.__exit__(None, None, None)\n"
                "    os._exit(0)\n"
                "os.write(trigger_write, b'G')\n"
                "os.close(trigger_write)\n"
                "result = os.read(result_read, 1)\n"
                "try:\n"
                "    with prepare_module._exclusive_export_generation(target):\n"
                "        parent_contended = False\n"
                "except RuntimeError as error:\n"
                "    parent_contended = 'generation is already active' in str(error)\n"
                "os.write(release_write, b'R')\n"
                "os.close(release_write)\n"
                "os.close(result_read)\n"
                "waited, status = os.waitpid(child, 0)\n"
                "child_ok = waited == child and os.waitstatus_to_exitcode(status) == 0\n"
                "raise SystemExit(0 if result == b'A' and parent_contended and child_ok else 3)\n"
            ),
            str(target),
        ],
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_publish_failure_restores_the_prior_generation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "generation"
    prepare(**_prepare_kwargs(output, dimensions=3))
    original = _snapshot(output)
    real_rename = Path.rename

    def fail_staging_publish(path: Path, target: Path):
        if path.name.startswith(".generation.cellucid-stage-") and Path(target) == output:
            raise OSError("synthetic publish failure")
        return real_rename(path, target)

    monkeypatch.setattr(Path, "rename", fail_staging_publish)
    replacement = _prepare_kwargs(output, dimensions=2)
    replacement["force"] = True
    with pytest.raises(OSError, match="synthetic publish failure"):
        prepare(**replacement)

    assert _snapshot(output) == original
    _assert_only_persistent_lock(tmp_path, output.name)


@pytest.mark.parametrize(
    (
        "had_target",
        "target_exists",
        "stage_exists",
        "backup_exists",
        "expected_owner",
    ),
    [
        (True, True, False, False, "target"),
        (True, True, True, False, "target"),
        (True, False, True, True, "backup"),
        (True, True, False, True, "target"),
        (False, False, False, False, None),
        (False, False, True, False, None),
        (False, True, False, False, "target"),
    ],
    ids=[
        "prior-target-before-stage",
        "prior-target-and-stage-roll-back",
        "prior-moved-restore",
        "new-target-and-backup-commit",
        "initial-before-stage",
        "initial-stage-roll-back",
        "initial-target-commit",
    ],
)
def test_export_transaction_recovers_every_unambiguous_state(
    tmp_path: Path,
    had_target: bool,
    target_exists: bool,
    stage_exists: bool,
    backup_exists: bool,
    expected_owner: str | None,
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    transaction_id = "0123456789abcdef0123456789abcdef"
    journal, journal_temp, stage, backup = prepare_module._export_transaction_paths(
        target,
        transaction_id,
    )
    expected_journal = (
        '{"format":"cellucid-export-transaction","version":1,'
        f'"transaction_id":"{transaction_id}",'
        f'"had_target":{"true" if had_target else "false"}}}\n'
    ).encode("ascii")
    journal.write_bytes(expected_journal)

    for path, exists, owner in (
        (target, target_exists, "target"),
        (stage, stage_exists, "stage"),
        (backup, backup_exists, "backup"),
    ):
        if exists:
            path.mkdir()
            (path / "owner").write_text(owner, encoding="utf-8")

    prepare_module._recover_export_transaction(target)

    assert not journal.exists()
    assert not journal_temp.exists()
    assert not stage.exists()
    assert not backup.exists()
    if expected_owner is None:
        assert not target.exists()
    else:
        assert (target / "owner").read_text(encoding="utf-8") == expected_owner


@pytest.mark.parametrize(
    ("had_target", "state"),
    [
        (True, (False, False, False)),
        (True, (False, False, True)),
        (True, (False, True, False)),
        (True, (True, True, True)),
        (False, (False, False, True)),
        (False, (False, True, True)),
        (False, (True, False, True)),
        (False, (True, True, False)),
        (False, (True, True, True)),
    ],
)
def test_export_transaction_fails_closed_for_ambiguous_states(
    tmp_path: Path,
    had_target: bool,
    state: tuple[bool, bool, bool],
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    transaction_id = "abcdef0123456789abcdef0123456789"
    journal, _journal_temp, stage, backup = prepare_module._export_transaction_paths(
        target,
        transaction_id,
    )
    journal.write_bytes(
        prepare_module._serialize_export_transaction(
            transaction_id,
            had_target=had_target,
        )
    )
    for path, exists in zip((target, stage, backup), state, strict=True):
        if exists:
            path.mkdir()
            (path / "sentinel").write_text(path.name, encoding="utf-8")
    before = {path: _snapshot(path) if path.is_dir() else None for path in (target, stage, backup)}

    with pytest.raises(RuntimeError, match="cannot be recovered without guessing"):
        prepare_module._recover_export_transaction(target)

    assert journal.exists()
    for path in (target, stage, backup):
        if before[path] is None:
            assert not path.exists()
        else:
            assert _snapshot(path) == before[path]


@pytest.mark.parametrize(
    "mutation",
    [
        "malformed-journal",
        "different-transaction",
        "different-target-history",
        "temporary-control-reappeared",
    ],
)
def test_export_transaction_revalidates_journal_before_publication_mutation(
    tmp_path: Path,
    mutation: str,
) -> None:
    from cellucid.prepare_data import _locking as lock_module
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    target.mkdir()
    (target / "owner").write_text("prior", encoding="utf-8")
    prior_snapshot = _snapshot(target)

    with lock_module._exclusive_export_generation(target):
        transaction_id, had_target, stage, backup = prepare_module._begin_export_transaction(target)
        assert had_target is True
        (stage / "owner").write_text("candidate", encoding="utf-8")
        journal, journal_temp, _expected_stage, _expected_backup = (
            prepare_module._export_transaction_paths(target, transaction_id)
        )
        if mutation == "malformed-journal":
            journal.write_bytes(b"{}")
        elif mutation == "different-transaction":
            journal.write_bytes(
                prepare_module._serialize_export_transaction(
                    "abcdef0123456789abcdef0123456789",
                    had_target=True,
                )
            )
        elif mutation == "different-target-history":
            journal.write_bytes(
                prepare_module._serialize_export_transaction(
                    transaction_id,
                    had_target=False,
                )
            )
        else:
            journal_temp.write_bytes(b"unowned")

        with pytest.raises(
            RuntimeError,
            match="journal|temporary control",
        ):
            prepare_module._publish_export_generation(
                stage,
                target,
                transaction_id=transaction_id,
                had_target=had_target,
            )

    assert _snapshot(target) == prior_snapshot
    assert stage.is_dir()
    assert not backup.exists()


@pytest.mark.parametrize(
    "phase",
    [
        "journal-published",
        "before-stage-publish",
        "prior-moved",
        "stage-published",
    ],
)
def test_killed_export_recovers_at_every_rename_boundary(
    tmp_path: Path,
    phase: str,
) -> None:
    target = tmp_path / "generation"
    expected_replacement = tmp_path / "expected-replacement"
    old_kwargs = _prepare_kwargs(target, dimensions=3)
    old_kwargs["created_at"] = "2026-01-02T03:04:05Z"
    prepare(**old_kwargs)
    old_snapshot = _snapshot(target)
    new_kwargs = _prepare_kwargs(expected_replacement, dimensions=2)
    new_kwargs["created_at"] = "2026-01-02T03:04:05Z"
    prepare(**new_kwargs)
    new_snapshot = _snapshot(expected_replacement)

    _crash_export_at_rename(target, phase)
    retry = _prepare_kwargs(target, dimensions=2)
    with pytest.raises(FileExistsError, match="Refusing to replace non-empty"):
        prepare(**retry)

    expected = new_snapshot if phase == "stage-published" else old_snapshot
    assert _snapshot(target) == expected
    _assert_only_persistent_lock(tmp_path, target.name)


def test_transaction_journal_rejects_noncanonical_bytes_without_mutation(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    target.mkdir()
    sentinel = target / "sentinel"
    sentinel.write_bytes(b"prior generation")
    journal, _journal_temp = prepare_module._export_transaction_control_paths(target)
    journal.write_text(
        '{"format":"cellucid-export-transaction","version":1,'
        '"transaction_id":"0123456789abcdef0123456789abcdef",'
        '"had_target":true} ',
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="not canonical"):
        prepare_module._recover_export_transaction(target)

    assert sentinel.read_bytes() == b"prior generation"
    assert journal.exists()


def test_partial_transaction_journal_write_is_recovered_before_target_policy(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    target.mkdir()
    sentinel = target / "sentinel"
    sentinel.write_bytes(b"prior generation")
    _journal, journal_temp = prepare_module._export_transaction_control_paths(target)
    journal_temp.write_bytes(b'{"format":"cellucid-export')

    with pytest.raises(FileExistsError, match="Refusing to replace non-empty"):
        prepare(**_prepare_kwargs(target, dimensions=2))

    assert sentinel.read_bytes() == b"prior generation"
    assert not journal_temp.exists()
    _assert_only_persistent_lock(tmp_path, target.name)


def test_transaction_journal_hard_link_fails_closed(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    target.mkdir()
    transaction_id = "0123456789abcdef0123456789abcdef"
    journal, _journal_temp, _stage, _backup = prepare_module._export_transaction_paths(
        target,
        transaction_id,
    )
    source = tmp_path / "journal-source"
    source.write_bytes(
        prepare_module._serialize_export_transaction(
            transaction_id,
            had_target=True,
        )
    )
    os.link(source, journal)

    with pytest.raises(RuntimeError, match="non-linked"):
        prepare_module._recover_export_transaction(target)

    assert source.exists()
    assert journal.exists()
    assert os.path.samefile(source, journal)


def test_transaction_stage_symlink_fails_closed(
    tmp_path: Path,
) -> None:
    from cellucid.prepare_data import _transaction as prepare_module

    target = tmp_path / "generation"
    target.mkdir()
    transaction_id = "0123456789abcdef0123456789abcdef"
    journal, _journal_temp, stage, _backup = prepare_module._export_transaction_paths(
        target,
        transaction_id,
    )
    journal.write_bytes(
        prepare_module._serialize_export_transaction(
            transaction_id,
            had_target=True,
        )
    )
    outside = tmp_path / "outside"
    outside.mkdir()
    try:
        stage.symlink_to(outside, target_is_directory=True)
    except (NotImplementedError, OSError) as error:
        pytest.skip(f"Directory symlinks are unavailable: {error}")

    with pytest.raises(RuntimeError, match="non-symbolic directory"):
        prepare_module._recover_export_transaction(target)

    assert stage.is_symlink()
    assert outside.is_dir()
    assert journal.exists()


@pytest.mark.parametrize("marker", ["suffix-only", "partial-v2"])
def test_cli_rejects_incomplete_zarr_declarations(
    tmp_path: Path,
    marker: str,
) -> None:
    store = tmp_path / "fixture.zarr"
    store.mkdir()
    if marker == "partial-v2":
        (store / ".zattrs").write_text("{}", encoding="utf-8")
    assert _detect_data_format(store) == "unknown"


def test_cli_accepts_one_complete_zarr_declaration(tmp_path: Path) -> None:
    store = tmp_path / "fixture"
    _write_zarr_v2_root(store)
    assert _detect_data_format(store) == "zarr"


def test_exported_root_rejects_an_orphan_instead_of_omitting_it(
    tmp_path: Path,
) -> None:
    root = tmp_path / "exports"
    root.mkdir()
    _write_exported_dataset(root / "valid", "valid")
    (root / "orphan").mkdir()

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = root
    with pytest.raises(ValueError, match="orphan.*complete exported dataset"):
        handler._list_datasets()
    with pytest.raises(ValueError, match="orphan.*complete exported dataset"):
        _detect_data_format(root)


def test_exported_root_rejects_duplicate_dataset_ids(tmp_path: Path) -> None:
    root = tmp_path / "exports"
    root.mkdir()
    _write_exported_dataset(root / "first", "duplicate")
    _write_exported_dataset(root / "second", "duplicate")

    handler = CORSRequestHandler.__new__(CORSRequestHandler)
    handler.data_dir = root
    with pytest.raises(ValueError, match="duplicate dataset id"):
        handler._list_datasets()
