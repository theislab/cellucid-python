"""
Cellucid CLI - Unified Command Line Interface

This module provides a single entry point for all cellucid command-line operations:
- serve: Serve data (auto-detects h5ad, zarr, or pre-exported)

Usage:
    cellucid serve /path/to/exported_data
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example
    cellucid --version
"""

from __future__ import annotations

import argparse
import errno
import logging
import socket
import sys
import traceback
from pathlib import Path
from typing import Final

# Import shared configuration from _server_base
from ._console import console_print
from ._server_base import CELLUCID_WEB_URL, DEFAULT_HOST, DEFAULT_PORT

ISSUE_TRACKER_URL = "https://github.com/theislab/cellucid-python/issues"

# Every way ``cellucid serve`` can end early falls on one of two sides, and the
# side decides whether a traceback is printed.
#
# These types are the ones an *operator condition* arrives as: a flag that is
# missing or does not apply, a path that is not a dataset, a manifest that does
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
        help=f"Host to bind to (default: {DEFAULT_HOST}). Use 0.0.0.0 for remote access.",
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
    Detect the format of the data at the given path.

    Returns:
        'h5ad' - HDF5-based AnnData file
        'zarr' - Zarr-based AnnData store
        'exported' - Pre-exported cellucid dataset
        'unknown' - Unable to detect format
    """

    if not path.exists():
        return "unknown"

    from .anndata_adapter import _classify_anndata_path

    if path.is_file():
        try:
            return _classify_anndata_path(path)
        except ValueError:
            return "unknown"

    if not path.is_dir():
        return "unknown"

    declares_zarr = (
        path.suffix == ".zarr" or (path / ".zgroup").exists() or (path / ".zattrs").exists()
    )
    if declares_zarr:
        try:
            return _classify_anndata_path(path)
        except ValueError:
            return "unknown"

    from .server._datasets import _list_exported_datasets

    if _list_exported_datasets(path):
        return "exported"

    return "unknown"


def _add_serve_subparser(subparsers, common_parser: argparse.ArgumentParser) -> None:
    """Add the 'serve' subcommand with auto-detection."""
    serve_parser = subparsers.add_parser(
        "serve",
        parents=[common_parser],
        help="Serve data (auto-detects h5ad, zarr, or pre-exported)",
        description="Serve data for visualization. Automatically detects format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Auto-detection:
    - .h5ad files → served directly via AnnData
    - .zarr directories → served directly via AnnData
    - Export directories (or exports roots with multiple datasets) → served as pre-exported data

Examples:
    # Serve an h5ad file (auto-detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example

    # Serve a zarr store (auto-detected)
    cellucid serve /path/to/data.zarr --dataset-name Example --dataset-id example

    # Serve pre-exported data (auto-detected)
    cellucid serve /path/to/exported_dataset

    # Serve only named obs columns, leaving out one Cellucid cannot serve
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example \\
        --obs-key cell_type --obs-key total_counts

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Behind jupyter-server-proxy, name the proxy's host explicitly
    cellucid serve /path/to/data --allowed-host hub.example.org

    # For SSH tunnel access from remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -L 8765:localhost:8765 user@server
    # Then open the exact Viewer URL printed by Cellucid.
""",
    )

    serve_parser.add_argument(
        "data_path",
        type=str,
        help="Path to h5ad file, zarr directory, or pre-exported dataset",
    )
    # AnnData-specific options. Supplying one for another format is an error.
    serve_parser.add_argument(
        "--latent-key",
        type=str,
        default=None,
        metavar="KEY",
        help="Explicit obsm latent-space key for AnnData centroid outliers",
    )
    serve_parser.add_argument(
        "--dataset-name",
        type=str,
        default=None,
        metavar="NAME",
        help="Explicit dataset name (required for h5ad and zarr)",
    )
    serve_parser.add_argument(
        "--dataset-id",
        type=str,
        default=None,
        metavar="ID",
        help="Explicit stable dataset identifier (required for h5ad and zarr)",
    )
    serve_parser.add_argument(
        "--obs-key",
        action="append",
        type=str,
        default=None,
        dest="obs_key",
        metavar="KEY",
        help=(
            "Exact AnnData obs column to serve; repeat once per column. "
            "Every column is served when this is omitted. Naming the columns "
            "leaves out one Cellucid cannot serve, such as a datetime column"
        ),
    )
    serve_parser.add_argument(
        "--vector-field-default",
        type=str,
        default=None,
        metavar="FIELD_ID",
        help="Exact default field id when AnnData contains multiple UMAP vector fields",
    )

    serve_parser.set_defaults(func=_run_serve)


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
        raise ValueError(
            f"Unable to detect format for: {data_path}. Expected an exact .h5ad "
            "file, complete Zarr store, or complete prepared Cellucid generation."
        )

    if data_format == "exported":
        inapplicable = [
            flag
            for flag, supplied in (
                ("--latent-key", args.latent_key is not None),
                ("--dataset-name", args.dataset_name is not None),
                ("--dataset-id", args.dataset_id is not None),
                ("--obs-key", args.obs_key is not None),
                (
                    "--vector-field-default",
                    args.vector_field_default is not None,
                ),
            )
            if supplied
        ]
        if inapplicable:
            raise ValueError(
                f"{', '.join(inapplicable)} may only be used with direct AnnData input. "
                "A prepared export already declares its identity and its columns in "
                "dataset_identity.json, so remove the flag and run the command again."
            )
        # Pre-exported dataset - use standard server
        from .server import serve

        serve(
            data_dir=str(data_path),
            port=args.port,
            host=args.host,
            open_browser=not args.no_browser,
            quiet=args.quiet,
            serve_web_ui=not args.no_web_ui,
            web_source_url=args.web_source_url,
            web_cache_dir=args.web_cache_dir,
            allowed_hosts=args.allowed_host,
        )
    else:
        # AnnData (h5ad or zarr) - use AnnData server.
        #
        # Both flags are reported together. Naming only the first sends someone
        # who supplied neither back to the terminal twice for one mistake.
        missing_identity = [
            flag
            for flag, value in (
                ("--dataset-name", args.dataset_name),
                ("--dataset-id", args.dataset_id),
            )
            if not isinstance(value, str) or not value
        ]
        if missing_identity:
            verb = "are" if len(missing_identity) > 1 else "is"
            raise ValueError(
                f"{' and '.join(missing_identity)} {verb} required when serving h5ad or "
                "zarr data. An .h5ad file and a Zarr store carry no Cellucid identity of "
                "their own, so name the dataset on the command line, for example: "
                '--dataset-name "My study" --dataset-id my-study-v1'
            )
        if not args.quiet:
            console_print(
                "\nImporting dependencies (anndata, numpy, scipy)...",
                end=" ",
                flush=True,
            )

        from .anndata_server import serve_anndata

        if not args.quiet:
            console_print("done")

        adapter_options = {
            "latent_key": args.latent_key,
            "dataset_name": args.dataset_name,
            "dataset_id": args.dataset_id,
            "obs_keys": args.obs_key,
            "vector_field_default": args.vector_field_default,
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
            **adapter_options,
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


def _operator_error_message(error: Exception, args: argparse.Namespace) -> str:
    """Return one exact line naming a condition the operator can correct."""
    if isinstance(error, ImportError):
        if error.name:
            return (
                f"Cellucid needs the {error.name!r} package to read this input and it "
                "is not installed in this environment. Install it and run the command "
                f"again: pip install {error.name}"
            )
        return (
            f"Cellucid could not import a package it needs for this input: {error}. "
            "Install it into the same environment as cellucid and run the command "
            "again."
        )
    if isinstance(error, OSError):
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

    parser = argparse.ArgumentParser(
        prog="cellucid",
        description="Cellucid - Interactive Single-Cell Data Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
For more information, see:
    https://cellucid.readthedocs.io
    https://github.com/theislab/cellucid-python

Quick start:
    # Serve any data (format auto-detected)
    cellucid serve /path/to/data.h5ad --dataset-name Example --dataset-id example
    cellucid serve /path/to/data.zarr --dataset-name Example --dataset-id example
    cellucid serve /path/to/exported_data

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
