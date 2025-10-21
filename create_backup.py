#!/usr/bin/env python3
"""Create backup of files before fixes."""

import shutil
from pathlib import Path

# Define backup directory
backup_dir = Path("before_Claude_151025")
backup_dir.mkdir(exist_ok=True)

# Files and directories to backup
items = [
    "webapp",
    "aladin_scripts",
    "cross_match_ngc2264.py"
]

print(f"Creating backup in {backup_dir.absolute()}")

for item in items:
    src = Path(item)
    if not src.exists():
        print(f"Warning: {item} not found, skipping")
        continue

    dst = backup_dir / item

    if src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        print(f"Copied directory: {item}")
    else:
        shutil.copy2(src, dst)
        print(f"Copied file: {item}")

print("Backup complete!")
