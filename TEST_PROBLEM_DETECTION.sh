#!/bin/bash
# Test script to verify problem detection logging
#
# This script runs the dynamic server with debug logging enabled
# to trace why problem_pages might be empty.
#
# Usage:
#   ./TEST_PROBLEM_DETECTION.sh

set -e

echo "=========================================="
echo "Problem Detection Debug Test"
echo "=========================================="
echo ""
echo "This will start the dynamic server with debug logging enabled."
echo "The logs will show:"
echo "  - Which pages are being generated"
echo "  - Problem detection for each panel"
echo "  - Page-level problem aggregation"
echo "  - Final index JSON contents"
echo ""
echo "Press Ctrl+C to stop the server when done."
echo ""

# Enable debug logging
export PROBLEM_DETECTION_DEBUG=1

# Optional: Limit to first page only for faster testing
# export ALADIN_PAGE_ONLY=1

# Start the server
# The server will generate pages on-demand when you visit them
echo "Starting server at http://127.0.0.1:5050/"
echo ""
echo "To test:"
echo "  1. Open http://127.0.0.1:5050/page/2MASS/1 in your browser"
echo "  2. Watch the console output for debug logs"
echo "  3. Check aladin_scripts/html/2MASS_index.json after page generation"
echo ""
echo "Starting in 3 seconds..."
sleep 3

python webapp/dyn_server.py --host 127.0.0.1 --port 5050
