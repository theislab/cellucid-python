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


def _detect_data_format(path: Path) -> str:
    """
    Detect the format of the data at the given path.

    Returns:
        'h5ad' - HDF5-based AnnData file
        'zarr' - Zarr-based AnnData store
        'exported' - Pre-exported cellucid dataset
        'unknown' - Unable to detect format
    """
    # Check file extension first (before exists check for better error messages)
    path_str = str(path).lower()
    if path_str.endswith('.h5ad'):
        return 'h5ad' if path.exists() else 'unknown'
    if path_str.endswith('.zarr'):
        return 'zarr' if path.exists() else 'unknown'

    if not path.exists():
        return 'unknown'

    # For directories, check contents
    if path.is_dir():
        # Check for pre-exported dataset (has dataset_identity.json - always created by prepare())
        if (path / 'dataset_identity.json').exists():
            return 'exported'
        # Check for zarr structure (.zattrs or .zgroup at root)
        if (path / '.zattrs').exists() or (path / '.zgroup').exists():
            return 'zarr'

    return 'unknown'


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
    - Directories with dataset_identity.json → served as pre-exported data

Examples:
    # Serve an h5ad file (auto-detected)
    cellucid serve /path/to/data.h5ad

    # Serve a zarr store (auto-detected)
    cellucid serve /path/to/data.zarr

    # Serve pre-exported data (auto-detected)
    cellucid serve /path/to/exported_dataset

    # Serve on a different port
    cellucid serve /path/to/data --port 9000

    # Load entire AnnData into memory (disable lazy loading)
    cellucid serve /path/to/data.h5ad --no-backed

    # For SSH tunnel access from remote server:
    # On the server: cellucid serve /path/to/data
    # On local machine: ssh -L 8765:localhost:8765 user@server
    # Then open: https://www.cellucid.com?remote=http://localhost:8765
""",
    )

    serve_parser.add_argument(
        "data_path",
        type=str,
        help="Path to h5ad file, zarr directory, or pre-exported dataset",
    )
    # AnnData-specific options (ignored for pre-exported data)
    serve_parser.add_argument(
        "--no-backed",
        action="store_true",
        help="Load entire AnnData into memory (default: lazy loading for h5ad)",
    )
    serve_parser.add_argument(
        "--latent-key",
        type=str,
        default=None,
        metavar="KEY",
        help="Key in obsm for latent space (AnnData only, auto-detected if not specified)",
    )

    serve_parser.set_defaults(func=_run_serve)


def _add_serve_anndata_subparser(subparsers, common_parser: argparse.ArgumentParser) -> None:
    """Add the 'serve-anndata' subcommand (deprecated, kept for backwards compatibility)."""
    anndata_parser = subparsers.add_parser(
        "serve-anndata",
        parents=[common_parser],
        help="[Deprecated] Use 'cellucid serve' instead (auto-detects format)",
        description="Serve AnnData directly. DEPRECATED: Use 'cellucid serve' instead.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
DEPRECATED: This command is deprecated. Use 'cellucid serve' instead.
The 'serve' command now auto-detects h5ad, zarr, and pre-exported formats.

Examples:
    # Old way (deprecated):
    cellucid serve-anndata /path/to/data.h5ad

    # New way (recommended):
    cellucid serve /path/to/data.h5ad
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
    """Execute the 'serve' subcommand with auto-detection."""
    # Configure logging based on args
    _configure_logging(args)

    if not args.quiet:
        try:
            from . import __version__
        except ImportError:
            __version__ = "0.0.0"
        print(f"cellucid v{__version__}")

    # Detect format
    data_path = Path(args.data_path)
    data_format = _detect_data_format(data_path)

    if data_format == 'unknown':
        if not data_path.exists():
            print(f"Error: Path not found: {data_path}", file=sys.stderr)
        else:
            print(f"Error: Unable to detect format for: {data_path}", file=sys.stderr)
            print("Expected: .h5ad file, .zarr directory, or pre-exported directory (with dataset_identity.json)", file=sys.stderr)
        sys.exit(1)

    if data_format == 'exported':
        # Pre-exported dataset - use standard server
        from .server import serve
        serve(
            data_dir=str(data_path),
            port=args.port,
            host=args.host,
            open_browser=not args.no_browser,
            quiet=args.quiet,
        )
    else:
        # AnnData (h5ad or zarr) - use AnnData server
        if not args.quiet:
            print("\nImporting dependencies (anndata, numpy, scipy)...", end=" ", flush=True)

        from .anndata_server import serve_anndata

        if not args.quiet:
            print("done")

        serve_anndata(
            data=str(data_path),
            port=args.port,
            host=args.host,
            open_browser=not args.no_browser,
            quiet=args.quiet,
            backed=not args.no_backed,
            latent_key=args.latent_key,
        )


def _run_serve_anndata(args: argparse.Namespace) -> None:
    """Execute the 'serve-anndata' subcommand."""
    # Configure logging based on args
    _configure_logging(args)

    if not args.quiet:
        try:
            from . import __version__
        except ImportError:
            __version__ = "0.0.0"
        print(f"cellucid v{__version__}")
        print("\nImporting dependencies (anndata, numpy, scipy)...", end=" ", flush=True)

    # Import here to avoid circular imports and speed up CLI startup
    # This import is slow because it loads anndata/numpy/scipy
    from .anndata_server import serve_anndata

    if not args.quiet:
        print("done")

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
    # Serve any data (format auto-detected)
    cellucid serve /path/to/data.h5ad
    cellucid serve /path/to/data.zarr
    cellucid serve /path/to/exported_data

Use 'cellucid serve --help' for more options.
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
