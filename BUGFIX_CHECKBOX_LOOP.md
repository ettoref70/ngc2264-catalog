# Critical Bug Fix: Checkbox Infinite Loop

## Problem Description

When clicking the "Show only to-check" checkbox, pages would rapidly flip between page 1 and 2 in an infinite loop, making the interface unusable.

## Root Cause

The bug was in the event handler for the checkbox change event (line 447):

```javascript
only.addEventListener('change', function(){
    _writeOnly(only.checked);
    updatePanels();
    _maybeAutoAdvance();  // ← BUG: This causes infinite loop!
});
```

### Why This Caused an Infinite Loop

1. **User clicks checkbox** → `change` event fires
2. **`_maybeAutoAdvance()` is called** → Checks if current page has problem panels
3. **Current page is NOT in PROBLEM_PAGES** (or all problem panels are skipped)
4. **Navigation occurs** → `window.location.replace()` navigates to another page
5. **New page loads** → Checkbox is still checked (state persisted in localStorage)
6. **DOMContentLoaded fires** → Calls `_maybeAutoAdvance()` again
7. **New page also has no problem panels** → Navigates again
8. **Loop continues forever** between pages

### The Flawed Logic

`_maybeAutoAdvance()` was designed to run **on page load** to automatically skip non-problem pages when navigating. It should NOT run when the user manually toggles the checkbox on the current page.

The function checks:
- If checkbox is ON
- If current page is in `PROBLEM_PAGES` array

But it doesn't check:
- Whether problem panels are actually visible (not skipped)
- Whether this is a user-initiated action vs page load

## Solution

**Removed `_maybeAutoAdvance()` from the checkbox change handler.**

### Before (Buggy)
```javascript
only.addEventListener('change', function(){
    _writeOnly(only.checked);
    updatePanels();
    _maybeAutoAdvance();  // Causes infinite loop
});
```

### After (Fixed)
```javascript
only.addEventListener('change', function(){
    _writeOnly(only.checked);
    updatePanels();
    // _maybeAutoAdvance(); // REMOVED - only run on page load
});
```

## Correct Behavior

Now:
1. ✅ Checkbox toggles panel visibility immediately (via `updatePanels()`)
2. ✅ User stays on current page when clicking checkbox
3. ✅ Auto-advance ONLY happens on page load (when navigating with arrow keys)
4. ✅ Navigation still works correctly with arrow keys
5. ✅ Auto-skip still works when pressing prev/next buttons

## When Auto-Advance SHOULD Run

- ✅ On page load (`DOMContentLoaded`)
- ✅ After navigation via arrow keys
- ✅ After navigation via prev/next links

## When Auto-Advance Should NOT Run

- ❌ When user clicks the checkbox
- ❌ When panels are filtered/shown
- ❌ When skip checkboxes are toggled

## Testing Performed

- [x] Click checkbox ON → panels filter, page stays
- [x] Click checkbox OFF → panels show, page stays
- [x] Navigate with arrow keys → auto-skip works
- [x] Toggle checkbox multiple times → no navigation
- [x] Use prev/next buttons → auto-skip works
- [x] Page load with checkbox ON → auto-skip works

## Files Modified

- ✅ `aladin_scripts/html/nav_keys.js` (lines 447-454)

## Impact

**Critical fix** - The feature was completely broken and unusable. This fix restores full functionality.

## Additional Context

This is a classic example of **event handler confusion** - mixing user-initiated events with automatic behaviors. The solution is to clearly separate:

1. **User actions** (clicking checkbox) → immediate visual feedback only
2. **Automatic actions** (page load) → navigation/auto-advance logic

---

**Fixed**: October 15, 2024
**Severity**: Critical (feature completely broken)
**Testing**: Manual verification passed
