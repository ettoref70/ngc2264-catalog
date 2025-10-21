# "Show Only To-Check" Logic - Fixes Applied

## Summary

Fixed **3 critical bugs** in the "show only to-check" navigation feature in [nav_keys.js](aladin_scripts/html/nav_keys.js).

## Issues Fixed

### 0. ✅ **CRITICAL: Infinite Loop When Clicking Checkbox**

**Problem**: Clicking the "show only to-check" checkbox caused pages to rapidly flip between page 1 and 2 in an infinite loop, making the feature **completely unusable**.

**Root Cause**: The checkbox change event handler was calling `_maybeAutoAdvance()`, which would navigate to another page, which would load with the checkbox still checked, triggering auto-advance again, creating an infinite navigation loop.

**Location**: `aladin_scripts/html/nav_keys.js:447-454`

**Before**:
```javascript
only.addEventListener('change', function(){
    _writeOnly(only.checked);
    updatePanels();
    _maybeAutoAdvance();  // BUG: Causes infinite navigation loop!
});
```

**After**:
```javascript
only.addEventListener('change', function(){
    _writeOnly(only.checked);
    updatePanels();
    // _maybeAutoAdvance(); // REMOVED - only run on page load, not on checkbox toggle
});
```

**Impact**:
- **CRITICAL FIX** - Feature was completely broken and unusable
- User can now safely toggle checkbox without triggering navigation
- Auto-advance still works correctly on page load and arrow key navigation
- Checkbox now only filters panels visually, doesn't cause navigation

**Severity**: **CRITICAL** (feature completely broken)

---

### 1. ✅ TOTAL_PAGES Undefined Reference (Medium Severity)

**Problem**: Lines 399 and 407 directly referenced `TOTAL_PAGES` variable which may not be defined, causing `Math.min(undefined, ...)` to return `NaN` and break navigation.

**Location**: `aladin_scripts/html/nav_keys.js:399-407`

**Before**:
```javascript
else if(k === 'ArrowRight'){
  if(onlyOn){
    var pb2 = document.querySelectorAll('.panelbox');
    if(!pb2.length){
      e.preventDefault();
      // BUG: TOTAL_PAGES might be undefined
      window.location.replace(PAGE_KEY + '_page' + Math.min(TOTAL_PAGES, PAGE_NUM+1) + '.html#seek=next');
      return;
    }
```

**After**:
```javascript
else if(k === 'ArrowRight'){
  if(onlyOn){
    var pb2 = document.querySelectorAll('.panelbox');
    var totalPages = _totalPages() || 999;  // Safe fallback for unknown page count
    if(!pb2.length){
      e.preventDefault();
      window.location.replace(PAGE_KEY + '_page' + Math.min(totalPages, PAGE_NUM+1) + '.html#seek=next');
      return;
    }
```

**Impact**:
- Prevents navigation breaking when `TOTAL_PAGES` is undefined
- Uses proper `_totalPages()` function which handles edge cases
- Falls back to 999 if page count is unknown (safe upper limit)

---

### 2. ✅ Missing Error Logging (Low Severity)

**Problem**: Errors in `_maybeAutoAdvance` and localStorage access were silently swallowed, making debugging difficult.

**Location**: `aladin_scripts/html/nav_keys.js:228, 286`

**Before**:
```javascript
function _maybeAutoAdvance(){
  try{
    var onlyOn = (function(){
      try{
        var key = 'ONLY_PROB_' + (PAGE_KEY || '');
        return localStorage.getItem(key) === '1';
      }catch(e){ return false; }  // Silent failure
    })();
    // ...
  }catch(e){}  // Silent failure
}
```

**After**:
```javascript
function _maybeAutoAdvance(){
  try{
    var onlyOn = (function(){
      try{
        var key = 'ONLY_PROB_' + (PAGE_KEY || '');
        return localStorage.getItem(key) === '1';
      }catch(e){
        console.error('[nav_keys] Error reading only-prob state:', e);
        return false;
      }
    })();
    // ...
  }catch(e){
    console.error('[nav_keys] Auto-advance error:', e);
  }
}
```

**Impact**:
- Errors now logged to console for debugging
- Doesn't change functionality, just improves observability
- Helps diagnose issues in production

---

## Logic Verified as Correct

During analysis, the following logic was verified and confirmed to be working correctly:

### ✅ Arrow Key Navigation
The loop logic for checking if a page has problem panels is **correct**:

```javascript
for(var i=0; i<pb.length; i++){
  if(isProblemPanel(pb[i]) && pb[i].getAttribute('data-skip') !== '1'){
    break;  // Found a problem panel, stop loop
  }
  if(i === pb.length-1){
    // Only executes if we reached the end WITHOUT breaking
    window.location.replace(...);  // Navigate away (correct!)
    return;
  }
}
```

This correctly identifies "no problem panels on current page" and auto-skips.

### ✅ Panel Visibility Filter
The `updatePanels()` function correctly shows/hides panels:

```javascript
if(onlyOn){
  show = prob && !skip;  // Show only problematic AND not skipped
}else{
  show = !skip;          // Show all except skipped
}
```

### ✅ Server-Side Navigation
The `/page/<key>/<int:page>/jump` endpoint correctly handles wrap-around navigation through problem pages only.

---

## Testing Performed

- [x] Navigation with arrow keys when "only to-check" is enabled
- [x] Navigation when `TOTAL_PAGES` is undefined
- [x] Navigation when `PROBLEM_PAGES` is empty
- [x] Auto-advance on page load
- [x] Error logging in browser console
- [x] Checkbox state persistence across reloads

---

## Additional Recommendations (Not Implemented)

These are low-priority enhancements for future consideration:

### 1. CSS-Based Initial State (Low Priority)
Add CSS to prevent flash of wrong content on slow page loads:

```css
/* Hide all panels initially */
.panelbox .mask {
    display: block;
}

/* Show when marked by JavaScript */
.panelbox[data-show="1"] .mask {
    display: none;
}
```

Then update JavaScript to set `data-show` attribute instead of directly manipulating `display` style.

### 2. Visual Indicator (Low Priority)
Add a counter showing "Page X of Y problem pages" when filter is active.

### 3. URL State Sync (Low Priority)
Persist "only to-check" state in URL query parameter for shareable links.

---

## Files Modified

- ✅ `aladin_scripts/html/nav_keys.js` - Fixed TOTAL_PAGES reference and added error logging

## Documentation Created

- ✅ `TOCHECK_LOGIC_ANALYSIS.md` - Comprehensive analysis of the feature
- ✅ `TOCHECK_FIXES.md` - This file (summary of fixes)

---

## Conclusion

The "show only to-check" feature is now more robust with:
- ✅ Fixed navigation edge cases
- ✅ Better error visibility for debugging
- ✅ Verified correct logic throughout

All critical and medium-severity issues have been resolved.
