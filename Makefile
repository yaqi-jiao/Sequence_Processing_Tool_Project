# -------------------------------------------------------
# Variables
# -------------------------------------------------------
# Prefer venv python if exists
BASE_PYTHON ?= python3.11

ifeq ($(OS),Windows_NT)
    PYTHON := .venv/Scripts/python.exe
else
    PYTHON := .venv/bin/python
endif

OUT ?= dist
PKG_DIR ?= .

.DEFAULT_GOAL := help
SHELL := bash
.PHONY: help venv install clean test lint build run
.SUFFIXES:
.DELETE_ON_ERROR:


# -------------------------------------------------------
# Help
# -------------------------------------------------------
help:
	@echo "Makefile commands:"
	@echo "  make venv     - Create virtual environment"
	@echo "  make install  - Install package in editable mode"
	@echo "  make clean    - Remove build artifacts"
	@echo "  make test     - Run tests"
	@echo "  make lint     - Check code format"
	@echo "  make build    - Build package (wheel)"
	@echo "  make run      - Example of running main script"


# -------------------------------------------------------
# Create virtual environment
# -------------------------------------------------------
venv:
	@echo "[venv] creating .venv with $(PYTHON)"
	@$(BASE_PYTHON) -m venv .venv
	@$(PYTHON) -m pip install --upgrade pip
	@echo "[venv] ready"

# -------------------------------------------------------
# Install package in editable mode
# -------------------------------------------------------
install:
	@$(PYTHON) -m pip install -e .
	@echo "[install] done (python=$$( $(PYTHON) -V ))"

# -------------------------------------------------------
# Clean: remove build artifacts
# -------------------------------------------------------
clean:
	@find . -name '*.pyc' -delete
	@find . -name '__pycache__' -type d -exec rm -rf {} +
	@rm -rf $(OUT) build *.egg-info
	@echo "[clean] removed build/test artifacts"

# -------------------------------------------------------
# Run unittest
# -------------------------------------------------------
test:
	@$(PYTHON) -m unittest discover tests
	@echo "[test] all tests passed"

# -------------------------------------------------------
# Lint: Static checks (ruff)
# -------------------------------------------------------
lint:
	@$(PYTHON) -m pip show ruff >/dev/null 2>&1 || $(PYTHON) -m pip install --upgrade ruff
	@$(PYTHON) -m ruff check . || (echo '[lint] ruff failed' >&2; exit 1)
	@echo "[lint] code formatting checked"

# -------------------------------------------------------
# Build wheel for PyPI
# -------------------------------------------------------
build:
	@$(PYTHON) -m pip show build >/dev/null 2>&1 || $(PYTHON) -m pip install --upgrade build
	@$(PYTHON) -m build -o $(OUT) $(PKG_DIR)
	@echo "[build] wheel/sdist built in ./$(OUT)"

# -------------------------------------------------------
# Run your tool manually (optional)
# -------------------------------------------------------
run:
	@$(PYTHON) -c "from SequenceTool.SeqProcessingTool import SeqProcessingTool; print('Tool loaded OK')"
	@echo "[run] SeqProcessingTool loaded successfully"