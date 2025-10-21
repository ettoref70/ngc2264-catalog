# Summary: Debug Logging for Problem Detection

## What Was Added

Comprehensive debug logging has been added to [cross_match_ngc2264.py](cross_match_ngc2264.py) to trace why `problem_pages` is empty.

### Code Changes

1. **Import logging** (line 31)
2. **Configure logger** (lines 89-99) with environment variable control
3. **Log page generation start** (line 2450-2451)
4. **Log page skipping** (line 2457) when `ALADIN_PAGE_ONLY` is set
5. **Log problem detection** (lines 4487-4505):
   - Case A detection (no IDs + nearby master)
   - Case B detection (multiple nearby masters)
   - Case C detection (D² > 1)
   - All panel details
6. **Log page aggregation** (lines 5130-5141):
   - Panel counts per page
   - Problematic panel names
   - Page added to `problem_pages`
7. **Log index writing** (lines 5750-5752):
   - Final index path
   - Total and problem page counts
   - Full list of problem pages

### Files Created

1. **[DEBUG_LOGGING_GUIDE.md](DEBUG_LOGGING_GUIDE.md)** - Comprehensive guide with:
   - Quick start instructions
   - Log level explanations
   - What to look for in output
   - Common scenarios and troubleshooting
   - Example debug session

2. **[TEST_PROBLEM_DETECTION.sh](TEST_PROBLEM_DETECTION.sh)** - Test script that:
   - Enables debug logging
   - Starts the server
   - Provides instructions for testing

3. **[PROBLEM_PAGES_ANALYSIS.md](PROBLEM_PAGES_ANALYSIS.md)** - Updated with logging info

---

## How to Use

### Quick Test

```bash
# Enable debug and run server
bash TEST_PROBLEM_DETECTION.sh

# In another terminal or browser, visit:
# http://127.0.0.1:5050/page/2MASS/1

# Watch the console for detailed logs
```

### Manual Test

```bash
# Enable debug logging
export PROBLEM_DETECTION_DEBUG=1

# Optional: Test just first page
export ALADIN_PAGE_ONLY=1

# Start server
python webapp/dyn_server.py --host 127.0.0.1 --port 5050

# Visit a page to trigger generation
```

---

## What You'll Learn

The logs will reveal:

1. ✅ **Is the problem detection code running?**
   - Look for: `Starting page generation for catalog '2MASS'`
   - Confirms: `aladin_dir=provided`

2. ✅ **Are any panels being marked problematic?**
   - Look for: `Problem detected: Case A/B/C`
   - Shows: Which criteria triggered

3. ✅ **What are the actual values?**
   - See: `num_valid_master`, `near_sep`, `near_thr`, `close_count`, `best_d2` for every panel

4. ✅ **Is the index being written correctly?**
   - Look for: `Wrote index aladin_scripts/html/2MASS_index.json`
   - Shows: Final `problem_pages` list

---

## Expected Outcomes

### If problems ARE being detected:
- You'll see `[INFO] Page X added to problem_pages`
- Index will have non-empty `problem_pages` array
- **Next step**: Test frontend navigation with those pages

### If NO problems are detected:
- All panels will show `is_problem=False`
- Index will have `problem_pages: []`
- **Next steps**:
  1. Review the actual values (are all sources very clean?)
  2. Consider adjusting criteria if too strict
  3. Check if test data exists with known problems

---

## Log Output Locations

All logs go to **stdout** (console where server runs).

To save to file:
```bash
export PROBLEM_DETECTION_DEBUG=1
python webapp/dyn_server.py --host 127.0.0.1 --port 5050 2>&1 | tee problem_detection.log
```

---

## Key Logging Points

| Location | What It Logs | Level |
|----------|--------------|-------|
| [cross_match_ngc2264.py:2450](cross_match_ngc2264.py#L2450) | Page generation start | INFO |
| [cross_match_ngc2264.py:2457](cross_match_ngc2264.py#L2457) | Page skipping (if any) | DEBUG |
| [cross_match_ngc2264.py:4487-4497](cross_match_ngc2264.py#L4487-L4497) | Problem detected (Case A/B/C) | DEBUG |
| [cross_match_ngc2264.py:4502-4505](cross_match_ngc2264.py#L4502-L4505) | Panel details | DEBUG |
| [cross_match_ngc2264.py:5131-5136](cross_match_ngc2264.py#L5131-L5136) | Page aggregation | DEBUG |
| [cross_match_ngc2264.py:5141](cross_match_ngc2264.py#L5141) | Page added to list | INFO |
| [cross_match_ngc2264.py:5750-5752](cross_match_ngc2264.py#L5750-L5752) | Index written | INFO |

---

## Environment Variables Reference

```bash
# Enable debug logging (required for detailed output)
export PROBLEM_DETECTION_DEBUG=1

# Limit to specific page for faster testing (optional)
export ALADIN_PAGE_ONLY=1

# Existing variables that affect behavior:
# - DYN_RADIUS_DEG: Search radius
# - ALADIN_HTML_MODE: all|problems|none
# - DYN_NO_IMAGES: Skip background images
```

---

## Troubleshooting

### No output at all
- Check environment variable is set: `echo $PROBLEM_DETECTION_DEBUG`
- Try unbuffered Python: `python -u webapp/dyn_server.py ...`

### Loop doesn't run
- Check if page is being requested in browser
- Look for errors in server output
- Verify `aladin_dir` is provided (should see `aladin_dir=provided`)

### Problems not detected
- Review DEBUG logs for actual values
- Compare against criteria in [PROBLEM_PAGES_ANALYSIS.md](PROBLEM_PAGES_ANALYSIS.md)
- Consider if data is genuinely clean

---

## Next Steps

After running with debug logging:

1. **Analyze the output** - What do the actual values show?
2. **Compare to criteria** - Are any being triggered?
3. **Adjust if needed** - Criteria, thresholds, or detection logic
4. **Test navigation** - Once problems are detected, verify frontend works

---

**Created**: 2025-10-15
**Status**: Ready for testing
**See Also**: [DEBUG_LOGGING_GUIDE.md](DEBUG_LOGGING_GUIDE.md), [PROBLEM_PAGES_ANALYSIS.md](PROBLEM_PAGES_ANALYSIS.md)
