#!/bin/bash

# Initialize and update submodules
git submodule update --init --recursive

# If --latest flag is passed, update to latest versions
if [ "$1" = "--latest" ]; then
    echo "Updating submodules to latest versions..."
    git submodule update --remote --merge
fi

# Ensure nested submodules are also updated
git submodule foreach --recursive git submodule update --init --recursive

# Show current status
git submodule status --recursive