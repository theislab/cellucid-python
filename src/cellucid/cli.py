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
import logging
import sys
from pathlib import Path

# Import shared configuration from _server_base
from ._server_base import CELLUCID_WEB_URL, DEFAULT_HOST, DEFAULT_PORT

logger = logging.getLogger("cellucid.cli")


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

    from .server import _list_exported_datasets

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

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

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

        print(f"cellucid v{__version__}")

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
                (
                    "--vector-field-default",
                    args.vector_field_default is not None,
                ),
            )
            if supplied
        ]
        if inapplicable:
            raise ValueError(
                f"{', '.join(inapplicable)} may only be used with direct AnnData input."
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
        )
    else:
        # AnnData (h5ad or zarr) - use AnnData server
        if not isinstance(args.dataset_name, str) or not args.dataset_name:
            raise ValueError("--dataset-name is required when serving h5ad or zarr data")
        if not isinstance(args.dataset_id, str) or not args.dataset_id:
            raise ValueError("--dataset-id is required when serving h5ad or zarr data")
        if not args.quiet:
            print("\nImporting dependencies (anndata, numpy, scipy)...", end=" ", flush=True)

        from .anndata_server import serve_anndata

        if not args.quiet:
            print("done")

        adapter_options = {
            "latent_key": args.latent_key,
            "dataset_name": args.dataset_name,
            "dataset_id": args.dataset_id,
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
            **adapter_options,
        )


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
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 130
    except Exception as e:
        logger.exception("Unexpected error")
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
