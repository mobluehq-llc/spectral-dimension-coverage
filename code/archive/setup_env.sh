#!/usr/bin/env bash
# Builds the pinned environment (Python 3.12, exact versions of requirements.lock) in archive/.venv. Requires `uv`.
set -euo pipefail
cd "$(dirname "$0")"
uv venv --python 3.12 --clear .venv
uv pip install --python .venv/bin/python -r requirements.lock
.venv/bin/python -I -c "import numpy,scipy,pandas,flint,highspy;print('env ok',numpy.__version__,scipy.__version__,pandas.__version__,flint.__version__)"
