#!/usr/bin/env python3
"""Report relative Markdown links whose target file does not exist."""

import re
import sys
from pathlib import Path

LINK = re.compile(r"\[[^\]]*\]\(([^)#\s]+)(?:#[^)]*)?\)")
SKIPPED_DIRS = {".git", "target"}


def check(root: Path) -> int:
    broken = 0
    for md in root.rglob("*.md"):
        if SKIPPED_DIRS.intersection(md.parts):
            continue
        for match in LINK.finditer(md.read_text(encoding="utf-8")):
            target = match.group(1)
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            if not (md.parent / target).exists():
                print(f"{md}: broken link -> {target}", file=sys.stderr)
                broken += 1
    return broken


if __name__ == "__main__":
    sys.exit(1 if check(Path(sys.argv[1] if len(sys.argv) > 1 else ".")) else 0)
