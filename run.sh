#!/bin/bash
# Wrapper for src/launcher.py. Paths inside the tool are resolved from the
# source file location. This can be invoked from any working directory.

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python3 "$SCRIPT_DIR/src/launcher.py" "$@"
