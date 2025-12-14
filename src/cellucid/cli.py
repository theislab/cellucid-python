"""
Cellucid CLI - Unified Command Line Interface

This module provides a single entry point for all cellucid command-line operations:
- serve: Serve pre-exported cellucid datasets
- serve-anndata: Serve AnnData (h5ad/zarr) directly without pre-export

Usage:
    cellucid serve /path/to/exported_data
    cellucid serve-anndata /path/to/data.h5ad
    cellucid --version

The CLI consolidates the previously separate `cellucid` and `cellucid-anndata`
commands into a single unified interface with subcommands.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Import shared configuration from _server_base
from ._server_base import DEFAULT_PORT, DEFAULT_HOST

logger = logging.getLogger("cellucid.cli")


def _create_common_server_parser() -> argparse.ArgumentParser:
    """
    Create a parent parser with common server arguments.

    This parser is shared by both 'serve' and 'serve-anndata' subcommands
    to ensure consistent behavior and reduce code duplication.
    """
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--port", "-p",
        type=int,
        default=DEFAULT_PORT,
        help=f"Port to serve on (default: {DEFAULT_PORT})",
    )
    parser.add_argument(
        "--host", "-H",
        type=str,
        default=DEFAULT_HOST,
        help=f"Host to bind to (default: {DEFAULT_HOST}). Use 0.0.0.0 for remote access.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress info messages",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose/debug logging",
    )

    return parser


def _add_serve_subparser(subparsers, common_parser: argparse.ArgumentParser) -> None:
    """Add the 'serve' subcommand for pre-exported datasets."""
    serve_parser = subparsers.add_parser(
        "serve",
        parents=[common_parser],
        help="Serve a pre-exported cellucid dataset",
        description="Serve a pre-exported cellucid dataset directory over HTTP.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Serve a local dataset
    cellucid serve /path/to/my_dataset

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Serve on all interfaces (for remote access)
    cellucid serve /path/to/data --host 0.0.0.0

    # For SSH tunnel access from remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -L 8765:localhost:8765 user@server
    # Then open: https://www.cellucid.com?remote=http://localhost:8765
""",
    )

    serve_parser.add_argument(
        "data_dir",
        type=str,
        help="Path to the pre-exported dataset directory",
    )

    serve_parser.set_defaults(func=_run_serve)


def _add_serve_anndata_subparser(subparsers, common_parser: argparse.ArgumentParser) -> None:
    """Add the 'serve-anndata' subcommand for direct AnnData serving."""
    anndata_parser = subparsers.add_parser(
        "serve-anndata",
        parents=[common_parser],
        help="Serve AnnData (h5ad or zarr) directly",
        description="Serve AnnData (h5ad file or zarr directory) directly for visualization.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Serve an h5ad file
    cellucid serve-anndata /path/to/data.h5ad

    # Serve a zarr store
    cellucid serve-anndata /path/to/data.zarr

    # Serve on a different port
    cellucid serve-anndata /path/to/data.h5ad --port 9000

    # Load entire file into memory (disable lazy loading)
    cellucid serve-anndata /path/to/data.h5ad --no-backed

Supported formats:
    - .h5ad files: HDF5-based AnnData files
    - .zarr directories: Zarr-based AnnData stores

Note:
    Serving data directly from AnnData is convenient but slower than using
    pre-exported data. For production use, consider exporting first:
        # Export to optimized format (coming soon)
        # cellucid prepare /path/to/data.h5ad -o /path/to/export/
        # Then serve the exported directory
        cellucid serve /path/to/export/
""",
    )

    anndata_parser.add_argument(
        "anndata_path",
        type=str,
        help="Path to h5ad file or zarr directory",
    )
    anndata_parser.add_argument(
        "--no-backed",
        action="store_true",
        help="Load entire file into memory (default: lazy loading for h5ad)",
    )
    anndata_parser.add_argument(
        "--latent-key",
        type=str,
        default=None,
        metavar="KEY",
        help="Key in obsm for latent space (auto-detected if not specified)",
    )

    anndata_parser.set_defaults(func=_run_serve_anndata)


def _run_serve(args: argparse.Namespace) -> None:
    """Execute the 'serve' subcommand."""
    # Configure logging based on args
    _configure_logging(args)

    # Import here to avoid circular imports and speed up CLI startup
    from .server import serve

    serve(
        data_dir=args.data_dir,
        port=args.port,
        host=args.host,
        open_browser=not args.no_browser,
        quiet=args.quiet,
    )


def _run_serve_anndata(args: argparse.Namespace) -> None:
    """Execute the 'serve-anndata' subcommand."""
    # Configure logging based on args
    _configure_logging(args)

    # Import here to avoid circular imports and speed up CLI startup
    from .anndata_server import serve_anndata

    serve_anndata(
        data=args.anndata_path,
        port=args.port,
        host=args.host,
        open_browser=not args.no_browser,
        quiet=args.quiet,
        backed=not args.no_backed,
        latent_key=args.latent_key,
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
    # Import version here to avoid circular imports
    try:
        from . import __version__
    except ImportError:
        __version__ = "0.0.0"

    parser = argparse.ArgumentParser(
        prog="cellucid",
        description="Cellucid - Interactive Single-Cell Data Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
For more information, see:
    https://cellucid.readthedocs.io
    https://github.com/theislab/cellucid-python

Quick start:
    # Serve pre-exported data
    cellucid serve /path/to/exported_data

    # Serve AnnData directly (slower, but no export needed)
    cellucid serve-anndata /path/to/data.h5ad

Use 'cellucid <command> --help' for more information on a specific command.
""",
    )

    # Add version flag
    parser.add_argument(
        "--version", "-V",
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
    _add_serve_anndata_subparser(subparsers, common_parser)

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
        parser.print_help()
        return 0

    parsed_args = parser.parse_args(args)

    # Check if a command was provided
    if not hasattr(parsed_args, 'func') or parsed_args.func is None:
        parser.print_help()
        return 0

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


# =============================================================================
# DEPRECATED: Legacy entry points for backwards compatibility
# =============================================================================
# These functions maintain backwards compatibility with the old CLI commands.
# They will be removed in a future version.

def main_serve_legacy() -> int:
    """
    Legacy entry point for 'cellucid' command.

    This is a backwards-compatible wrapper that redirects to the new CLI.
    Users should migrate to 'cellucid serve <path>' instead.

    .. deprecated::
        Use 'cellucid serve <path>' instead.
    """
    import warnings
    import sys

    # If we have arguments, inject 'serve' as the command
    if len(sys.argv) > 1:
        # Check if user is trying to use --help without a path
        if sys.argv[1] in ('-h', '--help'):
            # Show the serve subcommand help
            return main(['serve', '--help'])
        elif sys.argv[1].startswith('-'):
            # Flags before path - need to handle properly
            # Find where the path is and inject 'serve' before everything
            return main(['serve'] + sys.argv[1:])
        else:
            # Path is first argument - inject 'serve' before it
            return main(['serve'] + sys.argv[1:])
    else:
        # No arguments - show general help
        parser = create_parser()
        parser.print_help()
        print("\n" + "=" * 60)
        print("NOTE: The 'cellucid' command now requires a subcommand.")
        print("Use 'cellucid serve <path>' to serve pre-exported data.")
        print("=" * 60)
        return 0


def main_anndata_legacy() -> int:
    """
    Legacy entry point for 'cellucid-anndata' command.

    This is a backwards-compatible wrapper that redirects to the new CLI.
    Users should migrate to 'cellucid serve-anndata <path>' instead.

    .. deprecated::
        Use 'cellucid serve-anndata <path>' instead.
    """
    import sys

    # Inject 'serve-anndata' as the command
    if len(sys.argv) > 1:
        if sys.argv[1] in ('-h', '--help'):
            return main(['serve-anndata', '--help'])
        else:
            return main(['serve-anndata'] + sys.argv[1:])
    else:
        # No arguments - show serve-anndata help
        return main(['serve-anndata', '--help'])


if __name__ == "__main__":
    sys.exit(main())
