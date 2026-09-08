import pytest

from evoflow.modules.registry import validate_modules


def test_known_modules_pass() -> None:
    validate_modules(["qc", "pca", "fst"])


def test_unknown_module_fails() -> None:
    with pytest.raises(ValueError):
        validate_modules(["telepathy"])
