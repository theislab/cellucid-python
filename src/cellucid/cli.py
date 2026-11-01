"""
Cellucid CLI - Unified Command Line Interface

This module provides a single entry point for all cellucid command-line operations:
- serve: Serve a dataset file or folder the viewer reads (AnnData and MuData
  files and Zarr stores, Loom, 10x Cell Ranger ``.h5`` files and ``outs/``
  folders, MatrixMarket folders, R objects; ``input_formats.json``) to the
  viewer, which reads it itself

Usage:
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example
    cellucid serve /path/to/my_dataset.zarr  # written by cellucid.prepare()
    cellucid --version
"""

from __future__ import annotations

import argparse
import contextlib
import errno
import logging
import platform
import socket
import sys
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Final

# Import shared configuration from _server_base
from ._console import _text_for_console_stream, console_print
from ._server_base import CELLUCID_WEB_URL, DEFAULT_HOST, DEFAULT_PORT

if TYPE_CHECKING:
    from _typeshed import SupportsWrite

ISSUE_TRACKER_URL = "https://github.com/theislab/cellucid-python/issues"

# Every way ``cellucid serve`` can end early falls on one of two sides, and the
# side decides whether a traceback is printed.
#
# These types are the ones an *operator condition* arrives as: a flag that is
# missing or does not apply, a path that is not a dataset, a file that does
# not parse, a port another program already holds, an optional package that is
# not installed. Cellucid raises the first group deliberately through the
# exception taxonomy its documentation defines -- ``ValueError``, ``TypeError``,
# ``KeyError``, ``RuntimeError`` -- and the operating system reports its own as
# ``OSError``, of which ``socket.gaierror`` and ``urllib``'s ``URLError`` are
# subclasses. None of them is a defect, and the person reading the terminal is
# usually a biologist, so each is reported as one named line saying what failed
# and what to do about it.
#
# Everything else -- ``AttributeError``, ``IndexError``, ``AssertionError``,
# ``NameError`` and the rest -- means Cellucid did something it did not intend.
# That traceback is the bug report, so it is printed in full, always, even under
# ``--quiet``, and the message says where to send it.
_OPERATOR_ERROR_TYPES: Final = (
    ValueError,
    TypeError,
    KeyError,
    RuntimeError,
    OSError,
    ImportError,
)


class _ConsoleArgumentParser(argparse.ArgumentParser):
    """Argument parser whose help, usage, version, and errors are console-escaped.

    argparse writes to the stream directly, so on a CP1252 Windows console any
    character outside that code page -- in help text, or in a path echoed back
    in an error -- raised ``UnicodeEncodeError`` before a command ran. This
    applies the escaping ``console_print`` gives every other console line.
    Subparsers inherit this class, so every parser in the tree is covered.
    """

    def _print_message(self, message: str, file: SupportsWrite[str] | None = None) -> None:
        if message:
            stream = sys.stderr if file is None else file
            # argparse's own tolerance of a missing stream, e.g. under pythonw.
            with contextlib.suppress(AttributeError, OSError):
                stream.write(_text_for_console_stream(message, stream))


def _create_common_server_parser() -> argparse.ArgumentParser:
    """
    Create a parent parser with common server arguments.

    This parser is shared by server subcommands to ensure consistent behavior
    and reduce code duplication.
    """
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--port",
        "-p",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to serve on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--host",
        "-H",
        type=str,
        default=DEFAULT_HOST,
        help=(
            f"Host to bind to (default: {DEFAULT_HOST}, which no other machine "
            "can reach). Use 0.0.0.0 to accept connections on every interface, "
            "which is what a compute node needs when a tunnel terminates on a "
            "different machine such as a cluster login node. A non-loopback bind "
            "has no authentication: every user who can reach the port reads the "
            "dataset"
        ),
    )
    parser.add_argument(
        "--allowed-host",
        action="append",
        type=str,
        default=None,
        dest="allowed_host",
        metavar="HOST",
        help=(
            "Extra exact Host name this server answers to; repeat once per name. "
            "Needed only behind a reverse proxy such as jupyter-server-proxy. "
            "One bare DNS name or IP address, with no port, scheme or wildcard"
        ),
    )
    parser.add_argument(
        "--cache-bytes",
        type=int,
        default=None,
        dest="cache_bytes",
        metavar="BYTES",
        help=(
            "Bytes a store served in priority order keeps of what it makes "
            "(renumbered genes), least recently used first out; 0 keeps nothing "
            "(default: a 64th of the machine's memory)"
        ),
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically",
    )
    output_group = parser.add_mutually_exclusive_group()
    output_group.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress info messages",
    )
    output_group.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose/debug logging",
    )
    parser.add_argument(
        "--no-web-ui",
        action="store_true",
        help="Serve scientific data endpoints without establishing the web application",
    )
    parser.add_argument(
        "--web-source-url",
        default=CELLUCID_WEB_URL,
        metavar="URL",
        help="Exact origin publishing cellucid-web-assets.json",
    )
    parser.add_argument(
        "--web-cache-dir",
        type=Path,
        default=None,
        metavar="PATH",
        help="Directory for the active verified web build",
    )

    return parser


def _detect_data_format(path: Path) -> str:
    """
    Detect the format of the data at the given path, by the input format
    registry the web app reads too (``classify_anndata_path``).

    Returns:
        The input format id of the registry (``'h5ad'``, ``'zarr'``,
        ``'h5mu'``, ``'mudata-zarr'``, ``'loom'``, ``'tenx-h5'``,
        ``'tenx-outs'``, ``'rds'``, ``'rda'``, ``'mtx'``); ``'unknown'``
        when the path is none of them, or a format Cellucid recognises but
        never opens (an h5Seurat file)
    """
    from ._anndata_files import classify_anndata_path

    if not path.exists():
        return "unknown"
    try:
        return classify_anndata_path(path).kind
    except ValueError:
        return "unknown"


def _add_serve_subparser(subparsers, common_parser: argparse.ArgumentParser) -> None:
    """Add the 'serve' subcommand with auto-detection."""
    serve_parser = subparsers.add_parser(
        "serve",
        parents=[common_parser],
        help="Serve a dataset file or folder to the viewer (AnnData, MuData, Loom, 10x, MatrixMarket, R objects)",
        description="Serve a dataset file or folder to the viewer, which reads it itself.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Formats (detected, by the input format registry the web app reads):
    - .h5ad files; AnnData Zarr store folders (v2 or v3) and zipped stores (.zip)
    - .h5mu files and MuData Zarr stores
    - .loom files; 10x Cell Ranger .h5 files and outs/ folders
    - MatrixMarket folders (a .mtx[.gz] file of any name, such as GEO's GSM..._matrix.mtx.gz)
    - R objects: .rds, .rda, .RData, gzipped or not (.rds.gz from GEO)
    The browser reads every one itself; the server serves the bytes. An
    h5Seurat file is not opened: in R, read it with SeuratDisk::LoadH5Seurat()
    and save it with saveRDS(), then serve the .rds. A store
    written by cellucid.prepare() names itself (uns['cellucid']); any other is
    named after its file (its id the file's name without the extension), unless
    --dataset-name and --dataset-id name it.

Examples:
    # Serve an h5ad file (auto-detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example

    # Serve a zarr store (auto-detected)
    cellucid serve /path/to/data.zarr --dataset-name Example --dataset-id example

    # Serve a Cellucid-ready AnnData written by cellucid.prepare()
    cellucid serve /path/to/my_dataset.zarr

    # Serve only named obs columns, leaving out one Cellucid cannot serve
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \\
        --obs-key cell_type --obs-key total_counts

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Behind jupyter-server-proxy, name the proxy's host explicitly
    cellucid serve /path/to/data --allowed-host hub.example.org

    # Open on a chosen neighbor graph and velocity field
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \\
        --knn connectivities --vector-field-default velocity

    # Choose what is plotted: t-SNE in 2D, PCA components 1 and 2 in a second
    # run, and expression from a layer named by the gene-symbol column
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \\
        --embedding-2d X_tsne --expression layer:log1p --gene-id-column gene_symbols
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \\
        --embedding-2d X_pca:0,1

    # For SSH tunnel access from a remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -N -L 8765:127.0.0.1:8765 you@server
    # Then open the exact Viewer URL printed by Cellucid.

    # On an HPC compute node, where a tunnel can only terminate on the login node:
    # On the compute node: cellucid serve /path/to/data --host 0.0.0.0 --no-browser
    # On local machine:    ssh -N -L 8765:compute-node:8765 you@login-node
    # The login node opens the second connection, so the server has to accept it
    # on every interface. Then open the exact Viewer URL printed by Cellucid.
""",
    )

    serve_parser.add_argument(
        "data_path",
        type=str,
        help=(
            "Path to a file or folder Cellucid reads: .h5ad, AnnData or MuData Zarr "
            "(folder or .zip), .h5mu, .loom, a 10x .h5 or outs/ folder, a MatrixMarket "
            "folder, or an R object (.rds, .rda, .RData, gzipped or not)"
        ),
    )
    serve_parser.add_argument(
        "--latent-key",
        type=str,
        default=None,
        metavar="KEY",
        help="obsm latent space (scVI, PCA) the outlier filter measures in",
    )
    serve_parser.add_argument(
        "--centroid-min-points",
        type=_positive_int,
        default=None,
        metavar="N",
        help=(
            "Fewest cells a category needs for the outlier filter to measure it "
            "(default: the store's own, else 10)"
        ),
    )
    serve_parser.add_argument(
        "--dataset-name",
        type=str,
        default=None,
        metavar="NAME",
        help="Dataset name (default: the store's uns['cellucid']['name'], else the file's name)",
    )
    serve_parser.add_argument(
        "--dataset-id",
        type=str,
        default=None,
        metavar="ID",
        help="Stable dataset identifier (default: the store's uns['cellucid']['id'], else the file's name without its extension)",
    )
    serve_parser.add_argument(
        "--obs-key",
        action="append",
        type=str,
        default=None,
        dest="obs_key",
        metavar="KEY",
        help=(
            "Exact AnnData obs column the viewer shows; repeat once per column, "
            "in order. Every column is shown when this is omitted"
        ),
    )
    serve_parser.add_argument(
        "--modality",
        type=str,
        default=None,
        metavar="NAME",
        help="Modality a multimodal dataset (a MuData) opens on; the default one when omitted",
    )
    serve_parser.add_argument(
        "--vector-field-default",
        type=str,
        default=None,
        metavar="FIELD_ID",
        help="Vector field the overlay opens on: <field> of an obsm '<field>_<basis>' array drawn over 'X_<basis>'",
    )

    for dimension in (1, 2, 3):
        serve_parser.add_argument(
            f"--embedding-{dimension}d",
            type=str,
            default=None,
            dest=f"embedding_{dimension}d",
            metavar="KEY[:COLUMNS]",
            help=(
                f"AnnData obsm array to plot in {dimension}D, optionally with its "
                f"0-based columns, e.g. X_pca:{','.join(str(c) for c in range(dimension))}. "
                "Overrides the default (X_umap_<n>d or X_umap) for this dimension"
            ),
        )

    serve_parser.add_argument(
        "--normalize-embeddings",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Scale every embedding into the viewer's coordinate range; the no- "
            "form shows the coordinates as stored (default: the store's own, "
            "else scaled)"
        ),
    )

    serve_parser.add_argument(
        "--knn",
        type=str,
        default=None,
        metavar="KEY",
        help=(
            "AnnData obsp matrix the viewer draws as the neighbor graph, usually "
            "'connectivities'; read only when a view draws it"
        ),
    )

    serve_parser.add_argument(
        "--gene-id-column",
        type=str,
        default=None,
        metavar="COLUMN",
        help=(
            "AnnData var column (raw.var with --expression raw) whose values "
            "name the genes; the var index is used when omitted"
        ),
    )

    serve_parser.add_argument(
        "--expression",
        type=str,
        default=None,
        metavar="X|raw|layer:NAME",
        help=(
            "Where AnnData gene expression is read: X, raw, or layer:NAME "
            "(default: raw when it holds more genes than X or X is scaled, as Scanpy "
            "reads genes, else X; without X, raw or a layer, reported as guessed)"
        ),
    )

    serve_parser.set_defaults(func=_run_serve)


def _positive_int(value: str) -> int:
    """Parse a positive integer flag value, refusing anything else by name."""
    try:
        number = int(value, 10)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a positive integer; got {value!r}") from None
    if number < 1:
        raise argparse.ArgumentTypeError(f"expected a positive integer; got {value!r}")
    return number


def _parse_embedding_flag(value: str, *, flag: str) -> str | tuple[str, list[int]]:
    """Parse ``KEY`` or ``KEY:c0,c1`` into an ``embedding_keys`` value.

    The columns are 0-based. A key that itself contains ``:`` is read whole
    unless what follows its last ``:`` is a comma-separated list of integers.
    """
    if not value:
        raise ValueError(f"{flag} needs an obsm key, for example {flag} X_tsne.")
    key, separator, columns = value.rpartition(":")
    if separator and key and columns and all(part.isdigit() for part in columns.split(",")):
        return key, [int(part) for part in columns.split(",")]
    return value


def _run_serve(args: argparse.Namespace) -> None:
    """Execute the 'serve' subcommand with auto-detection."""
    # Configure logging based on args
    _configure_logging(args)

    if not args.quiet:
        from . import __version__

        console_print(f"cellucid v{__version__}")

    # Detect format
    data_path = Path(args.data_path)
    data_format = _detect_data_format(data_path)

    if data_format == "unknown":
        if not data_path.exists():
            raise FileNotFoundError(f"Path not found: {data_path}")
        from ._anndata_files import RefusedFormatError, classify_anndata_path

        # The registry's own explanation of why the path is none of its
        # formats, or of how to convert one it recognises but never opens.
        try:
            classify_anndata_path(data_path)
        except RefusedFormatError:
            raise
        except ValueError as error:
            raise ValueError(f"Unable to detect format for: {data_path}. {error}") from error

    # The viewer reads the store; the server names it, as given here or as
    # its uns['cellucid'] declares.
    from .anndata_server import serve_anndata

    embedding_keys = {
        f"{dimension}d": _parse_embedding_flag(
            getattr(args, f"embedding_{dimension}d"),
            flag=f"--embedding-{dimension}d",
        )
        for dimension in (1, 2, 3)
        if getattr(args, f"embedding_{dimension}d") is not None
    }
    serve_anndata(
        data=str(data_path),
        port=args.port,
        host=args.host,
        open_browser=not args.no_browser,
        quiet=args.quiet,
        serve_web_ui=not args.no_web_ui,
        web_source_url=args.web_source_url,
        web_cache_dir=args.web_cache_dir,
        allowed_hosts=args.allowed_host,
        cache_bytes=args.cache_bytes,
        latent_key=args.latent_key,
        gene_id_column=args.gene_id_column,
        normalize_embeddings=args.normalize_embeddings,
        centroid_min_points=args.centroid_min_points,
        dataset_name=args.dataset_name,
        dataset_id=args.dataset_id,
        obs_keys=args.obs_key,
        embedding_keys=embedding_keys or None,
        knn_key=args.knn,
        expression_source=args.expression,
        vector_field_default=args.vector_field_default,
        modality=args.modality,
    )


def _report_error(message: str) -> None:
    """Write one final ``Error:`` line to stderr, after everything already said.

    Progress output goes to stdout and is block-buffered whenever the command is
    piped into a file or a pager, so without this flush the failure line lands
    above the banner and steps it followed. The documentation tells readers to
    read the last ``Error:`` line; that has to hold when the transcript is saved,
    not only when it is watched live.
    """
    sys.stdout.flush()
    console_print(f"Error: {message}", file=sys.stderr)


def _bind_target(args: argparse.Namespace) -> str:
    """Describe the address the operator asked the server to bind."""
    return f"{getattr(args, 'host', DEFAULT_HOST)}:{getattr(args, 'port', DEFAULT_PORT)}"


def _operating_system_message(error: OSError, args: argparse.Namespace) -> str:
    """Name the operating-system condition and the flag that resolves it."""
    host = getattr(args, "host", DEFAULT_HOST)
    port = getattr(args, "port", DEFAULT_PORT)

    if isinstance(error, socket.gaierror) or error.errno == errno.EADDRNOTAVAIL:
        return (
            f"Host {host!r} is not an address this machine can serve from. Use "
            "--host 127.0.0.1 to serve only this computer, which is the default, or "
            "--host 0.0.0.0 to accept connections from other machines."
        )
    if error.errno == errno.EADDRINUSE:
        return (
            f"Port {port} is already in use on {host}, so Cellucid could not start its "
            "server there. Another program is listening on that port, often an earlier "
            "Cellucid that was never stopped. Serve on a different port with "
            "--port 9000, or let the operating system pick a free one with --port 0 and "
            "use the Viewer URL it prints."
        )
    if error.errno in {errno.EACCES, errno.EPERM} and error.filename is None:
        return (
            f"The operating system refused to let Cellucid bind {_bind_target(args)}. "
            "Ports below 1024 need administrator rights on most systems, so choose a "
            "higher port such as --port 8765, or --port 0 for any free port."
        )
    if error.filename is not None:
        detail = error.strerror or "could not be read"
        return (
            f"{detail}: {error.filename}. Check the path and its permissions, then run "
            "the command again."
        )
    return str(error)


def _required_python_range() -> str:
    """Return the Python range this installation declares, as it declares it.

    The range is read from the installed metadata rather than written here, so
    one declaration in pyproject.toml governs both the install and this message.
    """
    from importlib import metadata as _metadata

    try:
        declared = _metadata.metadata("cellucid")["Requires-Python"]
    except (_metadata.PackageNotFoundError, KeyError):
        declared = None
    return f"Python {declared}" if declared else "a supported Python"


def _import_error_message(error: ImportError) -> str:
    """Name the exact import condition, and only advise pip when pip can fix it.

    An ``ImportError`` is not one condition. A module that is absent from the
    environment is an install away. A module that is present but too old to hold
    the name asked of it is a version conflict, and installing it again changes
    nothing. A module the standard library owns is neither: it belongs to the
    running interpreter, and no distribution supplies it.
    """
    name = error.name or ""
    top_level = name.partition(".")[0]

    if not name:
        return (
            f"Cellucid could not import a module it needs for this input: {error}. "
            "Install it into the same environment as cellucid and run the command "
            "again."
        )

    if top_level in sys.stdlib_module_names or top_level in sys.builtin_module_names:
        return (
            f"Cellucid needs the standard-library module {name!r}, and this Python "
            f"build cannot import it. Python {platform.python_version()} at "
            f"{sys.executable} is running Cellucid, which requires "
            f"{_required_python_range()}. Run Cellucid on a Python in that range. "
            "No package installs this module: the interpreter itself provides it."
        )

    if type(error) is not ModuleNotFoundError:
        # The module imported; a name inside it did not. That is a version this
        # Cellucid cannot use, and installing the same version again repeats it.
        location = f" at {error.path}" if error.path else ""
        return (
            f"Cellucid could not read a name it needs from the installed {name!r}"
            f"{location}: {error.msg}. The installed version of that package is not "
            "one this Cellucid can use. Upgrade it in this environment and run the "
            f"command again: pip install --upgrade {top_level}"
        )

    if name != top_level:
        return (
            f"Cellucid needs {name!r}, and the installed {top_level!r} does not "
            f"provide it. Upgrade it in this environment and run the command again: "
            f"pip install --upgrade {top_level}"
        )

    return (
        f"Cellucid needs the {name!r} package to read this input and it is not "
        "installed in this environment. Install it and run the command again: "
        f"pip install {name}"
    )


def _operator_error_message(error: Exception, args: argparse.Namespace) -> str:
    """Return one exact line naming a condition the operator can correct."""
    if isinstance(error, ImportError):
        message = _import_error_message(error)
    elif isinstance(error, OSError):
        message = _operating_system_message(error, args)
    elif type(error) is KeyError and len(error.args) == 1 and type(error.args[0]) is str:
        # ``str(KeyError("no such column"))`` is ``"'no such column'"``. The
        # quotes make a written instruction read like a dumped value.
        message = error.args[0]
    else:
        message = str(error)
    # The documented CLI contract is one actionable line, so that the advice
    # "read the last Error: line" is literally true. Nothing raised on the serve
    # path wraps its message today; this keeps that guarantee if one ever does.
    return " ".join(message.split())


def _configure_logging(args: argparse.Namespace) -> None:
    """Configure logging based on command line arguments."""
    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        )
    elif not args.quiet:
        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
        )


def create_parser() -> argparse.ArgumentParser:
    """
    Create the main argument parser with all subcommands.

    Returns
    -------
    argparse.ArgumentParser
        The configured argument parser with subcommands.
    """
    # Import version here to avoid circular imports.
    from . import __version__

    parser = _ConsoleArgumentParser(
        prog="cellucid",
        description="Cellucid - Interactive Single-Cell Data Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
For more information, see:
    https://docs.cellucid.com
    https://github.com/theislab/cellucid-python

Quick start:
    # Serve any file or folder the web app reads (format detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example
    cellucid serve /path/to/my_dataset.zarr  # written by cellucid.prepare()
    # Formats only the browser reads (R objects, outs/ and MatrixMarket folders,
    # Loom, 10x .h5) carry no name Python can read, so name them:
    cellucid serve /path/to/object.rds --dataset-name PBMC --dataset-id pbmc
    cellucid serve /path/to/outs/ --dataset-name Sample --dataset-id sample

Use 'cellucid serve --help' for more options.
""",
    )

    # Add version flag
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version=f"cellucid {__version__}",
    )

    # Create subparsers for commands
    subparsers = parser.add_subparsers(
        title="commands",
        dest="command",
        description="Available commands. Use 'cellucid <command> --help' for details.",
        metavar="<command>",
    )

    # Create common server parser
    common_parser = _create_common_server_parser()

    # Add subcommands
    _add_serve_subparser(subparsers, common_parser)

    return parser


def main(args: list[str] | None = None) -> int:
    """
    Main entry point for the cellucid CLI.

    Parameters
    ----------
    args : list of str, optional
        Command line arguments. If None, uses sys.argv[1:].

    Returns
    -------
    int
        Exit code (0 for success, non-zero for failure).
    """
    parser = create_parser()

    # Handle case where no command is provided
    if args is None:
        args = sys.argv[1:]

    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        return 2

    try:
        parsed_args = parser.parse_args(args)
    except SystemExit as error:
        return error.code if type(error.code) is int else 1

    # Check if a command was provided
    if not hasattr(parsed_args, "func") or parsed_args.func is None:
        parser.print_help(file=sys.stderr)
        return 2

    try:
        # Execute the command
        parsed_args.func(parsed_args)
        return 0
    except KeyboardInterrupt:
        console_print("\nInterrupted by user", file=sys.stderr)
        return 130
    except _OPERATOR_ERROR_TYPES as error:
        # The stack is written to stderr rather than to a logger so that it does
        # not depend on whether ``--quiet`` skipped ``logging.basicConfig`` or on
        # how an embedding process configured logging.
        if getattr(parsed_args, "verbose", False):
            sys.stdout.flush()
            traceback.print_exception(error, file=sys.stderr)
        _report_error(_operator_error_message(error, parsed_args))
        return 1
    except Exception as error:
        sys.stdout.flush()
        traceback.print_exception(error, file=sys.stderr)
        _report_error(
            f"Cellucid failed unexpectedly with {type(error).__name__}: "
            f"{' '.join(str(error).split())}. This is a defect in Cellucid rather than "
            "a problem with your data or your command, so please report it at "
            f"{ISSUE_TRACKER_URL} and include the traceback printed above."
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
