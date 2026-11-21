#!/bin/bash
# Copy new or changed code from scratch into this repo, commit it and push it to GitHub.
#   ./push_code.sh "what changed"      (message optional; default lists the top-level folders that changed)
# Commits only if something changed. Never deletes files.
set -euo pipefail
cd "$(dirname "$0")"
./sync_from_scratch.sh
git add -A
if git diff --cached --quiet; then echo "nothing new to commit"; exit 0; fi
dirs=$(git diff --cached --name-only | awk -F/ '{print (NF>1 ? $1 : "top level")}' | sort -u | paste -sd, - | sed 's/,/, /g')
n=$(git diff --cached --name-only | wc -l)
git commit -q -m "${1:-Code sync from scratch: $n files ($dirs)}"
git push -q
git log -1 --format='pushed %h: %s'
