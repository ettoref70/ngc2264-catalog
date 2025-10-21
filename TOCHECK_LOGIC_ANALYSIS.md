# "Show Only To-Check" Logic Analysis

## Overview

The "show only to-check" feature allows users to filter catalog pages to show only problematic sources that need manual review. This analysis covers both client-side (JavaScript) and server-side (Python) implementations.

## How It Works

### 1. **Data Flow**

```
Page Generation (Python)
    ↓
Creates index JSON with to-check flags
    ↓
HTML pages include PROBLEM_PAGES array
    ↓
JavaScript reads checkbox state
    ↓
Filters visible panels OR navigates to problem pages
```

### 2. **Client-Side Logic** (nav_keys.js)

#### Storage Management
```javascript
// Read checkbox state (per-catalog persistence)
function _readOnly() {
    const k = "ONLY_PROB_" + PAGE_KEY;  // e.g., "ONLY_PROB_2MASS"
    return localStorage.getItem(k) === "1";
}

// Write checkbox state
function _writeOnly(v) {
    const k = "ONLY_PROB_" + PAGE_KEY;
    localStorage.setItem(k, v ? "1" : "0");
}
```

#### Panel Visibility (updatePanels)
```javascript
function updatePanels() {
    const only = document.getElementById("only_prob");
    const onlyOn = !!(only && only.checked);
    const pb = document.querySelectorAll(".panelbox");

    pb.forEach(function(p) {
        const prob = p.getAttribute("data-problem") === "1";  // Is this problematic?
        const skip = p.getAttribute("data-skip") === "1";      // Is this skipped?
        const mask = p.querySelector(".mask");

        let show = true;
        if (onlyOn) {
            show = prob && !skip;  // Show only if problematic AND not skipped
        } else {
            show = !skip;          // Show all except skipped
        }

        if (mask) {
            mask.style.display = show ? "none" : "block";
        }
    });
}
```

**Logic**:
- When checkbox is ON: Show panels where `data-problem="1"` AND `data-skip="1"` is NOT set
- When checkbox is OFF: Show all panels except those with `data-skip="1"`

#### Keyboard Navigation (Arrow Keys)

**ArrowLeft** (lines 374-393):
```javascript
if (k === 'ArrowLeft') {
    if (onlyOn) {
        var pb = document.querySelectorAll('.panelbox');
        if (!pb.length) {
            // No panels on page -> go to previous page
            window.location.replace(PAGE_KEY + '_page' + Math.max(1, PAGE_NUM-1) + '.html#seek=prev');
            return;
        }
        // Check if ANY non-skipped problem panels exist
        for (var i=0; i<pb.length; i++) {
            if (isProblemPanel(pb[i]) && pb[i].getAttribute('data-skip') !== '1') {
                break;  // Found one, use normal navigation
            }
            if (i === pb.length-1) {
                // NO problem panels on this page -> go to previous page
                window.location.replace(PAGE_KEY + '_page' + Math.max(1, PAGE_NUM-1) + '.html#seek=prev');
                return;
            }
        }
    }
    // Normal prev navigation
    var prev = findNav('prev');
    if (prev && prev.href) {
        window.location.href = prev.href;
    }
}
```

**ArrowRight** (lines 394-413):
Similar logic for next page navigation.

**Logic**:
- If "only to-check" is ON and current page has NO problem panels → auto-skip to next/prev page
- Otherwise → use normal navigation links

#### Auto-Advance (_maybeAutoAdvance)

Called on page load (line 450) to automatically skip non-problem pages:

```javascript
function _maybeAutoAdvance() {
    var onlyOn = localStorage.getItem('ONLY_PROB_' + PAGE_KEY) === '1';
    if (!onlyOn) return;  // Feature disabled

    // If current page IS a problem page, stop
    if (_isProblemPage(window.PROBLEM_PAGES, PAGE_NUM)) {
        sessionStorage.removeItem('AUTO_ALADIN_SUPPRESS');
        return;
    }

    // Current page is NOT a problem page -> find next problem page
    var dir = extractDirectionFromHash();  // 'next' or 'prev'
    var pages = window.PROBLEM_PAGES.slice().map(Number).filter(Number.isFinite).sort();

    if (dir === 'prev') {
        // Find previous problem page < current, or wrap to last
        target = pages.find(p => p < current) || pages[pages.length-1];
    } else {
        // Find next problem page > current, or wrap to first
        target = pages.find(p => p > current) || pages[0];
    }

    if (target !== current) {
        sessionStorage.setItem('AUTO_ALADIN_SUPPRESS', '1');
        window.location.replace(PAGE_KEY + '_page' + target + '.html#seek=' + dir);
    }
}
```

**Logic**: When "only to-check" is enabled, automatically skip non-problem pages on load.

### 3. **Server-Side Logic** (dyn_server.py)

#### Index JSON Structure

The server generates `<catalog>_index.json` files with:

```json
{
  "total_pages": 25,
  "to_check_pages": [3, 7, 12, 18],
  "pages": [
    {"page": 1, "to_check": false},
    {"page": 2, "to_check": false},
    {"page": 3, "to_check": true},
    ...
  ]
}
```

#### Parsing To-Check Pages (_to_check_set)

```python
def _to_check_set(index: dict) -> set[int]:
    """Extract page numbers marked as to-check from index JSON."""
    pages = set()

    # Method 1: Direct lists
    for k in ('to_check_pages', 'tocheck_pages', 'problematic_pages'):
        lst = index.get(k)
        if isinstance(lst, list):
            pages.update(int(p) for p in lst)

    # Method 2: Structured list with per-page flags
    arr = index.get('pages')
    if isinstance(arr, list):
        for it in arr:
            p = int(it.get('page'))
            flag = bool(it.get('to_check') or it.get('tocheck') or it.get('problematic'))
            if flag:
                pages.add(p)

    return pages
```

**Supports multiple JSON schemas** for flexibility.

#### Page Jump Endpoint (/page/<key>/<int:page>/jump)

Handles programmatic navigation with to-check filtering:

```python
@app.get('/page/<key>/<int:page>/jump')
def page_jump(key: str, page: int):
    """Navigate to next/prev page with optional to-check filtering.

    Query params:
      dir=next|prev           (default: next)
      only=to_check|all       (default: all)
      img=1|0                 (passthrough)
      debug=1|0               (passthrough)
    """
    direction = request.args.get('dir', 'next').lower()
    only = request.args.get('only', 'all').lower()

    index = _load_page_index(key)
    total = _index_total_pages(index)

    if only == 'to_check':
        tset = sorted(_to_check_set(index))
        if not tset:
            target = page  # No problem pages, stay
        else:
            if direction == 'prev':
                prevs = [p for p in tset if p < page]
                target = prevs[-1] if prevs else tset[-1]  # Wrap around
            else:
                nexts = [p for p in tset if p > page]
                target = nexts[0] if nexts else tset[0]    # Wrap around
    else:
        # Normal navigation with wrap-around
        if direction == 'prev':
            target = (page - 1) if page > 1 else total
        else:
            target = (page + 1) if page < total else 1

    return redirect(url_for('page', key=key, page=target, img=want_img, debug=debug))
```

**Logic**:
- Reads index JSON to find problem pages
- Navigates cyclically through problem pages only
- Falls back to current page if no problem pages exist

## Identified Issues & Fixes

### ✅ Issue 1: Arrow Key Logic Bug (FIXED)

**Problem**: Lines 382-389 in nav_keys.js had inverted logic

```javascript
// BEFORE (BUGGY)
for(var i=0;i<pb.length;i++){
    if(isProblemPanel(pb[i]) && pb[i].getAttribute('data-skip') !== '1'){
        break; // Found problem panel
    }
    if(i === pb.length-1){
        // This executes even if we found a problem panel!
        window.location.replace(...); // BUG: Navigates away
        return;
    }
}
```

**Status**: Actually, reviewing more carefully, the logic is **CORRECT**:
- Loop looks for problem panels
- `break` stops the loop if found
- The `if (i === pb.length-1)` check only executes if we reached the end WITHOUT breaking
- This correctly identifies "no problem panels on page"

### ✅ Issue 2: Race Between updatePanels and Navigation

**Problem**: `updatePanels()` is called, but navigation might happen before DOM updates

**Status**: **NOT A BUG** - The navigation checks happen synchronously in the event handler, and `updatePanels()` is only for visual filtering, not navigation decisions.

### ❌ Issue 3: Inconsistent State on Page Load

**Problem**: Checkbox state persisted in localStorage but panels visible before JS loads

**Potential Fix**: Add CSS rule to initially hide all panels, then show via JS:

```css
.panelbox .mask {
    display: block;  /* Hide by default */
}
.panelbox[data-show="1"] .mask {
    display: none;   /* Show when marked */
}
```

Then JS sets `data-show="1"` attribute instead of manipulating display directly.

**Severity**: Low - causes brief flicker on page load

### ❌ Issue 4: TOTAL_PAGES Undefined Reference

**Problem**: Line 399 and 407 reference `TOTAL_PAGES` which might not be defined

```javascript
window.location.replace(PAGE_KEY + '_page' + Math.min(TOTAL_PAGES, PAGE_NUM+1) + '.html#seek=next');
```

**Issue**: If `TOTAL_PAGES` is undefined, `Math.min(undefined, ...)` returns `NaN`

**Fix**:
```javascript
var total = _totalPages() || 999;  // Safe fallback
window.location.replace(PAGE_KEY + '_page' + Math.min(total, PAGE_NUM+1) + '.html#seek=next');
```

**Severity**: Medium - Can cause navigation to break

### ❌ Issue 5: Missing Error Handling in _maybeAutoAdvance

**Problem**: If `PROBLEM_PAGES` is malformed or `PAGE_NUM` is invalid, auto-advance fails silently

**Current**: Wrapped in try-catch, but errors are swallowed

**Improvement**: Log errors for debugging:

```javascript
function _maybeAutoAdvance() {
    try {
        // ... existing logic ...
    } catch(e) {
        console.error('[nav_keys] Auto-advance failed:', e);
    }
}
```

**Severity**: Low - Already handled, just hard to debug

## Recommendations

### 🟢 High Priority

1. **Fix TOTAL_PAGES undefined reference** (Issue #4)
   - Use `_totalPages()` function instead of direct variable access
   - Add safe fallback value

### 🟡 Medium Priority

2. **Add error logging** (Issue #5)
   - Help debugging in production

3. **Document data-problem attribute**
   - Clarify how pages are marked as problematic
   - Document in code generation

### 🔵 Low Priority

4. **Add CSS-based initial state** (Issue #3)
   - Prevent flash of wrong content
   - Better UX on slow connections

5. **Add visual indicator**
   - Show "X of Y problem pages" counter
   - Make it clear which mode is active

## Testing Checklist

- [ ] Checkbox persists across page reloads
- [ ] Panels filter correctly when checkbox toggled
- [ ] Arrow keys skip to next/prev problem page when enabled
- [ ] Auto-advance works on page load
- [ ] Wrap-around navigation works (last → first, first → last)
- [ ] Works when PROBLEM_PAGES is empty
- [ ] Works when PROBLEM_PAGES is undefined
- [ ] Works when TOTAL_PAGES is undefined
- [ ] Multiple catalogs maintain separate checkbox states
- [ ] Skip feature works in combination with only-to-check

## Summary

The "show only to-check" logic is **mostly correct** with a few edge cases:

✅ **Working Correctly:**
- Panel filtering based on data-problem attribute
- Checkbox state persistence per catalog
- Auto-advance to skip non-problem pages
- Keyboard navigation with problem-page awareness
- Server-side index generation and parsing
- Wrap-around navigation

❌ **Needs Fixing:**
- TOTAL_PAGES undefined reference (Medium severity)
- Initial page load flicker (Low severity)
- Error logging for debugging (Low severity)

The core logic is sound, but the implementation could be more robust in edge cases.
