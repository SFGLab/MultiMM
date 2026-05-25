#!/bin/bash
set -e

uv run isort simulation/ setup.py

uv run black simulation/ setup.py

uv run docformatter --in-place simulation/*.py setup.py

uv run ruff check .
