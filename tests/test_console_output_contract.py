from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _run_with_cp1252(script: str, *arguments: str) -> subprocess.CompletedProcess[bytes]:
    environment = os.environ.copy()
    environment["PYTHONIOENCODING"] = "cp1252:strict"
    return subprocess.run(
        [sys.executable, "-c", script, *arguments],
        env=environment,
        capture_output=True,
        timeout=30,
        check=False,
    )


def test_console_output_escapes_only_unencodable_cp1252_text() -> None:
    completed = _run_with_cp1252(
        "import sys\n"
        "from cellucid._console import console_print\n"
        "from cellucid._server_base import print_detail, print_success\n"
        "console_print('ASCII | → | ✓ | 细胞 | 🧬 | ×')\n"
        "console_print('错误', file=sys.stderr)\n"
        "print_detail('标签', '值🧬')\n"
        "print_success('完成🧬')\n"
    )

    assert completed.returncode == 0, completed.stderr.decode("cp1252")
    stdout = completed.stdout.decode("cp1252")
    stderr = completed.stderr.decode("cp1252")
    assert stdout == os.linesep.join(
        (
            "ASCII | \\u2192 | \\u2713 | \\u7ec6\\u80de | \\U0001f9ec | ×",
            "      \\u6807\\u7b7e: \\u503c\\U0001f9ec",
            "      \\u2713 \\u5b8c\\u6210\\U0001f9ec",
            "",
        )
    )
    assert stderr == f"\\u9519\\u8bef{os.linesep}"


def test_prepare_accepts_unicode_values_with_strict_cp1252_console(
    tmp_path: Path,
) -> None:
    target = tmp_path / "输出-🧬"
    completed = _run_with_cp1252(
        "import sys\n"
        "from pathlib import Path\n"
        "import numpy as np\n"
        "import pandas as pd\n"
        "from cellucid.prepare_data import prepare\n"
        "embedding = np.array([[-3.0, 1.0], [0.5, 5.0], [8.0, -2.0]], "
        "dtype=np.float32)\n"
        "prepare(\n"
        "    latent_space=embedding.copy(),\n"
        "    obs=pd.DataFrame({'score': [0.25, 0.5, 0.75]}),\n"
        "    X_umap_2d=embedding,\n"
        "    out_dir=Path(sys.argv[1]),\n"
        "    dataset_name='细胞 atlas 🧬',\n"
        "    dataset_id='unicode-console',\n"
        "    created_at='2026-01-02T03:04:05Z',\n"
        "    obs_categorical_dtype='uint16',\n"
        "    centroid_min_points=1,\n"
        ")\n",
        str(target),
    )

    stderr = completed.stderr.decode("cp1252")
    assert completed.returncode == 0, stderr
    assert "UnicodeEncodeError" not in stderr
    stdout = completed.stdout.decode("cp1252")
    assert "\\u8f93\\u51fa-\\U0001f9ec" in stdout
    assert "\\u2192" in stdout
    assert "\\u2713" in stdout
    assert (target / "dataset_identity.json").is_file()
