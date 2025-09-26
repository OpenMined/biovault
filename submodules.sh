#!/bin/bash

git submodule update --init --recursive
git submodule foreach --recursive git submodule update --init --recursive
git submodule foreach --recursive git pull --rebase