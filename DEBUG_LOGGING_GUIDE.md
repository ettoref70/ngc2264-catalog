# Debug Logging Guide for Problem Detection

## Quick Start

Enable debug logging to trace why `problem_pages` might be empty:

```bash
# Option 1: Run test script (recommended)
bash TEST_PROBLEM_DETECTION.sh

# Option 2: Manual
export PROBLEM_DETECTION_DEBUG=1
python webapp/dyn_server.py --host 127.0.0.1 --port 5050
```

Then visit `http://127.0.0.1:5050/page/2MASS/1` in your browser and watch the console output.

---

## Environment Variables

| Variable | Values | Effect |
|----------|--------|--------|
| `PROBLEM_DETECTION_DEBUG` | `1`, `true`, `yes`, `on` | Enable DEBUG-level logging for problem detection |
| `ALADIN_PAGE_ONLY` | Page number (e.g., `1`) | Only generate one page (faster testing) |

---

## Log Levels

### INFO (Always shown)

- **Page generation start**: Confirms loop is running
  ```
  [INFO] Starting page generation for catalog '2MASS': total=697, per_page=2,
         num_pages=349, aladin_dir=provided
  ```

- **Problem page added**: Confirms when a page is marked problematic
  ```
  [INFO] Page 23 added to problem_pages (now 5 problem pages)
  ```

- **Index written**: Shows final output
  ```
  [INFO] Wrote index aladin_scripts/html/2MASS_index.json: total_pages=349, problem_pages=5
  [INFO]   Problem pages: [2, 15, 23, 87, 145]
  ```

### DEBUG (Only with `PROBLEM_DETECTION_DEBUG=1`)

- **Panel-level detection**: For every single panel
  ```
  [DEBUG] [aladin_2MASS_1_c23] Panel j=22, new_id=06411481+0944287, is_problem=True,
          num_valid_master=0, near_sep=2.345, near_thr=3.000, close_count=0, best_d2=None
  [DEBUG] [aladin_2MASS_1_c23] Problem detected: Case A (no identifications, nearby master) -
          num_valid_master=0, near_sep=2.345, near_thr=3.000
  ```

- **Page-level aggregation**: Summary per page
  ```
  [DEBUG] Page 12: 2 panels, 1 problematic, page_has_problem=True
  [DEBUG]   Problematic panels on page 12: aladin_2MASS_1_c23
  ```

- **Page skipping** (if `ALADIN_PAGE_ONLY` is set):
  ```
  [DEBUG] Skipping page 2 (ALADIN_PAGE_ONLY=1)
  ```

---

## What to Look For

### 1. Is the loop running?

Look for:
```
[INFO] Starting page generation for catalog '2MASS': total=697, per_page=2, num_pages=349, aladin_dir=provided
```

**If missing**: The function isn't being called or `aladin_dir` is None.

### 2. Are pages being generated?

Look for DEBUG messages like:
```
[DEBUG] Page 1: 2 panels, 0 problematic, page_has_problem=False
```

**If missing**: Check if `ALADIN_PAGE_ONLY` is limiting generation, or if loop exits early.

### 3. Are any panels problematic?

Look for:
```
[DEBUG] [aladin_2MASS_1_c23] Problem detected: Case A/B/C ...
```

**If missing**: No sources meet the three problem criteria:
- **Case A**: `num_valid_master=0` AND `near_sep <= near_thr`
- **Case B**: `num_valid_master >= 1` AND `close_count >= 2`
- **Case C**: `best_d2 > 1.0`

### 4. Is the index being written?

Look for:
```
[INFO] Wrote index aladin_scripts/html/2MASS_index.json: total_pages=349, problem_pages=0
```

**If `problem_pages=0`**: Check panel-level logs to see why no problems detected.

---

## Common Scenarios

### Scenario 1: No log output at all

**Possible causes:**
- Server not running with debug enabled
- Python buffering output (try `python -u`)
- Logs going to a different handler

**Fix:**
```bash
# Ensure unbuffered output
python -u webapp/dyn_server.py --host 127.0.0.1 --port 5050
```

### Scenario 2: Loop runs but no panels logged

**Possible causes:**
- `aladin_dir` not provided (problem detection inside that block)
- Exception caught silently

**Check:**
- Look for `aladin_dir=provided` in startup log
- Add more try/except logging

### Scenario 3: Panels logged but all `is_problem=False`

**Possible causes:**
- Catalog data is genuinely clean (unlikely for 349 pages)
- Problem detection variables not computed correctly
- Criteria too strict

**Debug:**
- Check actual values: `num_valid_master`, `near_sep`, `near_thr`, `close_count`, `best_d2`
- Verify against criteria (see above)

### Scenario 4: Problems detected but not in index

**Possible causes:**
- Index writing failing silently
- Wrong path being read by frontend
- Multiple indices (with radius suffix)

**Check:**
```bash
find aladin_scripts/html -name "*_index.json"
cat aladin_scripts/html/2MASS_index.json
cat aladin_scripts/html/2MASS_r0.1_index.json  # if using radius
```

---

## Example Debug Session

```bash
$ export PROBLEM_DETECTION_DEBUG=1
$ export ALADIN_PAGE_ONLY=1  # Test first page only
$ python webapp/dyn_server.py --host 127.0.0.1 --port 5050

# In browser: visit http://127.0.0.1:5050/page/2MASS/1

# Console output:
[INFO] Starting page generation for catalog '2MASS': total=697, per_page=2, num_pages=349, aladin_dir=provided
[DEBUG] [aladin_2MASS_1_c1] Panel j=0, new_id=06405754+0945563, is_problem=False, num_valid_master=1, near_sep=0.123, near_thr=3.000, close_count=1, best_d2=0.45
[DEBUG] [aladin_2MASS_1_c2] Panel j=1, new_id=06410213+0939049, is_problem=True, num_valid_master=0, near_sep=2.876, near_thr=3.000, close_count=0, best_d2=None
[DEBUG] [aladin_2MASS_1_c2] Problem detected: Case A (no identifications, nearby master) - num_valid_master=0, near_sep=2.876, near_thr=3.000
[DEBUG] Page 1: 2 panels, 1 problematic, page_has_problem=True
[DEBUG]   Problematic panels on page 1: aladin_2MASS_1_c2
[INFO] Page 1 added to problem_pages (now 1 problem pages)
[INFO] Wrote index aladin_scripts/html/2MASS_index.json: total_pages=349, problem_pages=1
[INFO]   Problem pages: [1]
```

---

## Disabling Debug Logging

```bash
# Remove or set to 0
unset PROBLEM_DETECTION_DEBUG
# or
export PROBLEM_DETECTION_DEBUG=0
```

---

## Next Steps After Debugging

1. **If problems ARE detected**: Check frontend JSON loading, verify navigation works
2. **If NO problems detected**: Review criteria, possibly relax thresholds or add new cases
3. **If loop doesn't run**: Check `dyn_server.py` ensures `aladin_dir` is always passed

---

**Date**: 2025-10-15
**Status**: Logging added, ready for testing
