#!/bin/bash
# Script to run sphinx-autobuild for automatic documentation rebuilding

cd "$(dirname "$0")"
sphinx-autobuild source build/html --open-browser --host 127.0.0.1 --port 8000

