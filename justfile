default:
    @just --list

clean:
    @find . | grep -E "(__pycache__|\.pyc|\.pyo$)" | xargs rm -rf
    @rm -rf src/*.egg-info/ build/ dist/ .pytest_cache/ .ruff_cache/

test:
    @echo Testing…
    pytest

upload: clean
    @if [ -z "$(git describe --exact-match 2>/dev/null || true)" ]; then echo "Current commit is not a tag; abort task"; exit 1; fi
    python3 -m build --sdist --wheel
    twine upload dist/*
