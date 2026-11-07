#!/usr/bin/env bash
# Creates the next numbered architecture decision record in docs/adr/.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 \"decision title\"" >&2
  exit 2
fi

dir="docs/adr"
mkdir -p "$dir"

existing=$(find "$dir" -maxdepth 1 -name '[0-9][0-9][0-9][0-9]-*.md' | wc -l)
next=$(printf "%04d" $((existing + 1)))
slug=$(echo "$1" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/-/g; s/^-|-$//g')
file="$dir/$next-$slug.md"

cat > "$file" <<EOF
# $next. $1

Date: $(date +%F)
Status: Proposed

## Context

## Options considered

## Decision

## Consequences
EOF

echo "Created $file"
