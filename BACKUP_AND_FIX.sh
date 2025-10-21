#!/bin/bash
# Backup script - Run this before applying fixes
# Usage: bash BACKUP_AND_FIX.sh

set -e

BACKUP_DIR="before_Claude_151025"

echo "Creating backup directory: $BACKUP_DIR"
mkdir -p "$BACKUP_DIR"

echo "Backing up files..."
cp -r webapp "$BACKUP_DIR/"
cp -r aladin_scripts "$BACKUP_DIR/"
cp cross_match_ngc2264.py "$BACKUP_DIR/"

echo "Backup complete! Files saved in $BACKUP_DIR"
echo ""
echo "You can now proceed with applying the fixes."
echo "To restore from backup if needed:"
echo "  cp -r $BACKUP_DIR/* ."
