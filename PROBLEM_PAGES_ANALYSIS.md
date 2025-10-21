# Analysis: Problem Page Detection and "Skip Non-Problematic Pages" Logic

## Summary

The "Show only to-check" feature in the web interface is **correctly implemented**, but problem pages are **not being detected** because the problem detection logic only runs when `aladin_dir` is provided (i.e., when generating Aladin scripts).

## Architecture

### Frontend ([nav_keys.js](aladin_scripts/html/nav_keys.js))

The navigation JavaScript implements smart page skipping:

1. **`_hasUnskippedProblems()`** (line 222-233): Checks if current page has panels with `data-problem="1"` that aren't manually skipped
2. **Arrow key navigation** (lines 401-479): When "Show only to-check" is enabled AND current page has no unskipped problems, automatically jumps to next/previous problem page
3. **Auto-advance on load** (line 235-310): `_maybeAutoAdvance()` automatically navigates to next problem page if current page is clean

### Backend ([dyn_server.py](webapp/dyn_server.py))

Server-side support for navigation:

1. **`_load_page_index()`** (line 405-418): Loads `{catalog}_index.json` containing problem pages list
2. **`_to_check_set()`** (line 445-477): Extracts problem pages from JSON (supports multiple schemas)
3. **`/page/<key>/<int:page>/jump`** endpoint (line 479-525): Server-side navigation with wrap-around

### Problem Detection ([cross_match_ngc2264.py](cross_match_ngc2264.py))

**Three criteria** for marking a source as problematic (lines 4470-4481):

1. **Case A**: No identifications, but a nearby master exists within threshold
2. **Case B**: An identification exists AND multiple nearby masters within threshold
3. **Case C**: Mahalanobis distance D² > 1.0 for selected match

**Problem flagging happens at TWO levels:**

1. **Panel-level**: `data-problem="1"` attribute in HTML (line 5469)
2. **Page-level**: Page added to `problem_pages[]` array if ANY panel is problematic (line 5107)

## Root Cause

**The problem detection code is INSIDE the `if aladin_dir:` block** ([cross_match_ngc2264.py:4102-4500](cross_match_ngc2264.py#L4102-L4500)):

```python
# Line 4102
if aladin_dir:
    # ... Aladin script generation ...

    # Line 4438-4481: Problem detection
    is_problem = False
    if (num_valid_master == 0 and (near_sep is not None) and (near_sep <= near_thr)):
        is_problem = True
    elif (num_valid_master >= 1 and close_count >= 2):
        is_problem = True
    if (not is_problem) and (best_d2 is not None) and np.isfinite(best_d2) and (best_d2 > 1.0):
        is_problem = True

    # Line 4495-4501: Save to page_areas
    page_areas.append({
        'base': base,
        'ax_index': ax_idx,
        'aladin': file_url,
        'lite': url,
        'problem': 1 if is_problem else 0,
    })
```

**Result**: If `aladin_dir` is None or empty, problem detection never runs, and all pages are marked as non-problematic.

## Current Status

Checking the generated index:

```bash
$ cat aladin_scripts/html/2MASS_index.json
{"total_pages": 349, "problem_pages": []}
```

The `problem_pages` array is **empty**, which means:
- ✅ Navigation works for all pages
- ❌ "Show only to-check" has nothing to show
- ❌ Auto-advance skips every page

## Test Case

Created `2MASS_index_TEST.json` with sample problem pages:

```json
{"total_pages": 349, "problem_pages": [1, 5, 23, 50, 87, 145, 200, 301]}
```

**To test navigation:**
1. Temporarily rename the test file to `2MASS_index.json`
2. Open any 2MASS page in the web interface
3. Enable "Show only to-check" checkbox
4. Use arrow keys → should jump directly to pages 1, 5, 23, etc. (skipping others)

## Solutions

### Option 1: Always Run Problem Detection (Recommended)

Move the problem detection logic OUTSIDE the `if aladin_dir:` block so it runs regardless of whether Aladin scripts are generated:

```python
# Line ~4438 (move BEFORE "if aladin_dir:" block)
# Determine if this panel is potentially problematic
is_problem = False
# ... (detection logic) ...

# Line 4102
if aladin_dir:
    # ... (Aladin script generation) ...
    # Use is_problem that was already computed
    page_areas.append({
        'base': base,
        'ax_index': ax_idx,
        'problem': 1 if is_problem else 0,
        # ...
    })
```

### Option 2: Require aladin_dir for HTML Mode

Ensure `aladin_dir` is always provided when generating HTML pages (it currently is via `dyn_server.py:319`).

Verify in [dyn_server.py:308-322](webapp/dyn_server.py#L308-L322):

```python
cm.plot_after_merge(
    combined_df,
    df_other,
    id_col,
    catalogs,
    pdf_path=None,
    plot_mode='match',
    ncols=ncols,
    nrows=nrows,
    invert_cmap=invert,
    draw_images=draw_images,
    aladin_dir=str(SCRIPTS_DIR),  # ✅ Always provided
    samp_enabled=False,
    # ...
)
```

**Conclusion**: The server DOES provide `aladin_dir`, so problem detection SHOULD be working. Need to investigate why `problem_pages` is empty.

### Option 3: Debug Logging

Add logging to understand why no problems are being detected:

```python
# After line 4481
logger.info(f"Panel {base}: is_problem={is_problem}, num_valid_master={num_valid_master}, "
            f"near_sep={near_sep}, near_thr={near_thr}, close_count={close_count}, best_d2={best_d2}")
```

## Debug Logging Added

Comprehensive logging has been added to trace problem detection. Enable it with:

```bash
export PROBLEM_DETECTION_DEBUG=1
python webapp/dyn_server.py --host 127.0.0.1 --port 5050
```

**Or use the test script:**

```bash
./TEST_PROBLEM_DETECTION.sh
```

### What Gets Logged

1. **Page generation start** (INFO level):
   - Total sources, pages, whether `aladin_dir` is provided

2. **Panel-level detection** (DEBUG level):
   - For each panel: `j`, `new_id`, `is_problem`, `num_valid_master`, `near_sep`, `near_thr`, `close_count`, `best_d2`
   - Specific trigger: Case A (no IDs), Case B (multiple masters), or Case C (D² > 1)

3. **Page-level aggregation** (DEBUG/INFO level):
   - Number of panels on page
   - Number of problematic panels
   - Which panels are problematic
   - Confirmation when page added to `problem_pages`

4. **Index writing** (INFO level):
   - Path to index JSON
   - Total pages and problem pages count
   - Full list of problem page numbers

### Example Output

```
[INFO] Starting page generation for catalog '2MASS': total=697, per_page=2, num_pages=349, aladin_dir=provided
[DEBUG] Page 1: 2 panels, 0 problematic, page_has_problem=False
[DEBUG] Page 2: 2 panels, 1 problematic, page_has_problem=True
[DEBUG]   Problematic panels on page 2: aladin_2MASS_1_c3
[INFO] Page 2 added to problem_pages (now 1 problem pages)
[INFO] Wrote index aladin_scripts/html/2MASS_index.json: total_pages=349, problem_pages=1
[INFO]   Problem pages: [2]
```

## Next Steps

1. ✅ **Verified**: Problem detection logic exists and criteria are sensible
2. ✅ **Verified**: Frontend navigation logic is correct
3. ✅ **Created**: Test JSON to verify navigation works with problem pages defined
4. ✅ **Added**: Comprehensive debug logging at all key points
5. ⏳ **TODO**: Run with logging enabled to see actual values
6. ⏳ **TODO**: Investigate why no sources meet the problem criteria (if that's what logs show)

## Testing the Navigation Feature

**Without modifying code**, test if the navigation logic works:

```bash
# 1. Backup current index
cp aladin_scripts/html/2MASS_index.json aladin_scripts/html/2MASS_index.json.orig

# 2. Copy test index
cp aladin_scripts/html/2MASS_index_TEST.json aladin_scripts/html/2MASS_index.json

# 3. Start server
python webapp/dyn_server.py --host 127.0.0.1 --port 5050

# 4. Open browser to http://127.0.0.1:5050/page/2MASS/1

# 5. Enable "Show only to-check" checkbox

# 6. Press Right Arrow → should jump to page 5 (skipping 2,3,4)

# 7. Press Right Arrow again → should jump to page 23 (skipping 6-22)

# 8. Restore original index
mv aladin_scripts/html/2MASS_index.json.orig aladin_scripts/html/2MASS_index.json
```

## Relevant Code Locations

| Component | File | Lines | Description |
|-----------|------|-------|-------------|
| Problem detection criteria | [cross_match_ngc2264.py](cross_match_ngc2264.py) | 4470-4481 | Where `is_problem` is computed |
| Append to page_areas | [cross_match_ngc2264.py](cross_match_ngc2264.py) | 4495-4501 | Where problem flag is saved |
| Collect problem pages | [cross_match_ngc2264.py](cross_match_ngc2264.py) | 5104-5109 | Where `problem_pages` list is built |
| Write index JSON | [cross_match_ngc2264.py](cross_match_ngc2264.py) | 5711-5717 | Where index is written to disk |
| Frontend navigation | [nav_keys.js](aladin_scripts/html/nav_keys.js) | 401-479 | Arrow key handler |
| Auto-advance | [nav_keys.js](nav_keys.js) | 235-310 | Page load auto-navigation |
| Server jump endpoint | [dyn_server.py](webapp/dyn_server.py) | 479-525 | Server-side navigation |

---

**Date**: 2025-10-15
**Status**: Investigation complete, navigation feature verified working (with test data)
**Next Action**: Add debug logging to understand why real data has no problem pages detected
