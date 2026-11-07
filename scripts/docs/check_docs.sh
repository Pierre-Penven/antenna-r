#!/usr/bin/env bash
# Checks that the required documentation files, arc42 sections and C4 diagrams exist.
set -euo pipefail

fail=0

required_files=(
  "README.md"
  "docs/technical-design.md"
  "docs/arc42"
  "docs/c4"
  "docs/adr"
)

for f in "${required_files[@]}"; do
  if [[ ! -e "$f" ]]; then
    echo "MISSING: $f" >&2
    fail=1
  fi
done

arc42_sections=(
  "Introduction" "Constraints" "Context" "Solution Strategy"
  "Building Block" "Runtime" "Deployment" "Crosscutting"
  "Architecture Decisions" "Quality" "Risks" "Glossary"
)

for s in "${arc42_sections[@]}"; do
  if ! grep -rqi "$s" docs/arc42/ 2>/dev/null; then
    echo "MISSING arc42 section: $s" >&2
    fail=1
  fi
done

if ! grep -rlE "C4Context|C4Container" docs/c4/ >/dev/null 2>&1; then
  echo "MISSING: no C4 context or container diagram in docs/c4/" >&2
  fail=1
fi

if grep -rn "TODO: verify" README.md docs/ 2>/dev/null; then
  echo "WARNING: unverified items remain (see above)" >&2
fi

exit "$fail"
