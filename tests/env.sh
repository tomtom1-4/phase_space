#!/bin/bash

# Get the parent directory of the current working directory
PARENT_DIR="$(dirname "$(pwd)")"

# Add the parent directory to LD_LIBRARY_PATH
export LD_LIBRARY_PATH="$PARENT_DIR:$LD_LIBRARY_PATH"

# Optionally, print out the new LD_LIBRARY_PATH for verification
echo "Updated LD_LIBRARY_PATH: $LD_LIBRARY_PATH"