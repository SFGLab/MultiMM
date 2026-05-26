#!/bin/bash
set -e

uv run isort src/ --float-to-top

uv run black src/ --target-version py311

uv run docformatter --in-place src/multimm/*.py

uv run ruff check .

