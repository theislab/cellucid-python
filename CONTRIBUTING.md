# Contributing

We welcome contributions to Cellucid! This guide will help you get started.

## Development Setup

### Prerequisites

- Python 3.10+
- Git

### Installation

1. Clone the repository:

```bash
git clone https://github.com/theislab/cellucid-python.git
cd cellucid-python
```

2. Create a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install in development mode:

```bash
pip install -e ".[dev,docs]"
```

4. Set up pre-commit hooks:

```bash
pre-commit install
```

## Development Workflow

### Code Style

We use [Ruff](https://github.com/astral-sh/ruff) for linting and formatting:

```bash
ruff check .        # Lint
ruff format .       # Format
```

### Type Checking

```bash
mypy src/cellucid
```

### Testing

```bash
pytest                    # Run all tests
pytest -v --tb=short      # Verbose with short traceback
pytest --cov=cellucid     # With coverage
```

## Building Documentation

```bash
cd docs
make html
```

The built docs will be in `docs/_build/html/`.

## Pull Request Guidelines

1. Fork the repository and create a feature branch
2. Make your changes with clear commit messages
3. Add tests for new functionality
4. Ensure all tests pass and code style checks succeed
5. Submit a pull request with a clear description

## Reporting Issues

Please report bugs and feature requests on the [GitHub issue tracker](https://github.com/theislab/cellucid-python/issues).

When reporting bugs, include:
- Python version
- Operating system
- Steps to reproduce
- Expected vs actual behavior
