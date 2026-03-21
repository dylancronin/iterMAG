#!/usr/bin/env bash
set -euo pipefail

PROJECT="iterMAG"

echo "[info] Creating project: ${PROJECT}"

# Top-level
mkdir -p "${PROJECT}"
cd "${PROJECT}"

touch README.md LICENSE .gitignore pyproject.toml environment.yml

# Python package
mkdir -p src/itermag
touch src/itermag/__init__.py
touch src/itermag/cli.py
touch src/itermag/workflow.py

# Package submodules
mkdir -p src/itermag/utils
touch src/itermag/utils/paths.py
touch src/itermag/utils/config.py
touch src/itermag/utils/logging.py

mkdir -p src/itermag/config
touch src/itermag/config/config.yaml

mkdir -p src/itermag/templates
touch src/itermag/templates/config.yaml.j2

# Snakemake workflow
mkdir -p src/itermag/workflow
touch src/itermag/workflow/Snakefile

mkdir -p src/itermag/workflow/rules
touch src/itermag/workflow/rules/preprocess.smk
touch src/itermag/workflow/rules/assembly.smk
touch src/itermag/workflow/rules/binning.smk
touch src/itermag/workflow/rules/refine.smk
touch src/itermag/workflow/rules/qc.smk

mkdir -p src/itermag/workflow/envs
touch src/itermag/workflow/envs/base.yaml
touch src/itermag/workflow/envs/assembly.yaml
touch src/itermag/workflow/envs/binning.yaml

mkdir -p src/itermag/workflow/scripts
touch src/itermag/workflow/scripts/summarize_bins.py

mkdir -p src/itermag/workflow/schemas
touch src/itermag/workflow/schemas/config.schema.yaml

mkdir -p src/itermag/workflow/resources

# Tests
mkdir -p tests/data/small_example
touch tests/test_cli.py

# Docs
mkdir -p docs
touch docs/usage.md

# Optional profile (for HPC later)
mkdir -p profile
touch profile/config.yaml
touch profile/cluster.yaml

echo "[done] Structure created"
