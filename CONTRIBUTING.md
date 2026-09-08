# Contributing to EvoFlow

EvoFlow welcomes bug reports, feature requests, documentation improvements, tests, and scientific-method contributions.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
ruff check .
```

Please keep modules focused, document external-tool assumptions, and add tests for new behavior.
