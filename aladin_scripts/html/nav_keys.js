// Keyboard navigation: ←/→ pages; 'a' to trigger Aladin. Bounds-aware.
(function(){
  var _mx, _my; document.addEventListener('mousemove', function(ev){ _mx=ev.clientX; _my=ev.clientY; }, {passive:true});
  function isEditable(el){ if(!el) return false; var tag=(el.tagName||'').toUpperCase(); return el.isContentEditable||tag==='INPUT'||tag==='TEXTAREA'||tag==='SELECT'; }
  function closest(el, sel){ while(el && el.nodeType===1){ if(el.matches && el.matches(sel)) return el; el = el.parentElement; } return null; }
  function _totalPages(){ try{ return (typeof TOTAL_PAGES==='number' && TOTAL_PAGES>0) ? TOTAL_PAGES : 1; }catch(e){ return 1; } }
  function findNav(which){
    var navs = document.querySelectorAll('.nav');
    for(var i=0;i<navs.length;i++){
      var as = navs[i].querySelectorAll('a[href]');
      for(var j=0;j<as.length;j++){
        var a = as[j]; var t = (a.textContent||'').trim().toLowerCase();
        if(which==='prev' && t.indexOf('prev page') !== -1) return a;
        if(which==='next' && t.indexOf('next page') !== -1) return a;
      }
    }
    var m = (location.pathname||'').match(/([^\/]+)_page(\d+)\.html$/);
    if(m){ var base = m[1]; var n = parseInt(m[2],10)||1; var total=_totalPages(); var target = which==='prev' ? (n-1) : (n+1); if(target<1 || target>total) return null; var a = document.createElement('a'); a.setAttribute('href', base + '_page' + target + '.html'); return a; }
    return null;
  }
  function triggerAladin(){
    var cand=null;
    if(typeof _mx==='number' && typeof _my==='number'){ var el=document.elementFromPoint(_mx,_my); cand = closest(el, 'a.btn.al[data-ajs]'); }
    if(!cand) cand = document.querySelector('a.btn.al[data-ajs]');
    if(cand){ try{ cand.click(); return true; }catch(e){} }
    try{ if(typeof ALADIN_ITEMS!=='undefined' && ALADIN_ITEMS.length){ var idx=(window._srcIdx===undefined?0:window._srcIdx); var ajs = ALADIN_ITEMS[idx] || ALADIN_ITEMS[0]; var u=(typeof SAMP_URL!=='undefined'?SAMP_URL:'http://127.0.0.1:8765'); var base=(u && u.endsWith && u.endsWith('/'))?u.slice(0,-1):u; var url = base + '/run_samp?file=' + encodeURIComponent(ajs) + '&_t=' + Date.now(); var sink=document.getElementsByName('aladin_sink')[0]; if(sink){ try{ sink.src = url; return true; }catch(e){} } try{ window.open(url, '_blank'); return true; }catch(e){} } }catch(e){}
    return false;
  }
  function onKey(e){
    if(isEditable(e.target)) return;
    var k = (e.key||'');
    var onlyOn = false; try{ onlyOn = (typeof _readOnly==='function') ? _readOnly() : false; }catch(_e){}
    if(k === 'ArrowLeft'){
      if(onlyOn && Array.isArray(window.PROBLEM_PAGES) && window.PROBLEM_PAGES.length){
        try{ var arr=window.PROBLEM_PAGES.slice().sort(function(a,b){return a-b}); var target=null; for(var i=arr.length-1;i>=0;i--){ if(arr[i] < PAGE_NUM){ target=arr[i]; break; } } if(target===null) target=arr[arr.length-1]; if(target && target!==PAGE_NUM){ e.preventDefault(); location.href = PAGE_KEY + '_page' + target + '.html#seek=prev'; return; } }catch(_e){}
      }
      var a = findNav('prev'); if(a && a.getAttribute('href')){ e.preventDefault(); location.href = a.getAttribute('href'); }
    }
    else if(k === 'ArrowRight'){
      if(onlyOn && Array.isArray(window.PROBLEM_PAGES) && window.PROBLEM_PAGES.length){
        try{ var arr=window.PROBLEM_PAGES.slice().sort(function(a,b){return a-b}); var target=null; for(var i=0;i<arr.length;i++){ if(arr[i] > PAGE_NUM){ target=arr[i]; break; } } if(target===null) target=arr[0]; if(target && target!==PAGE_NUM){ e.preventDefault(); location.href = PAGE_KEY + '_page' + target + '.html#seek=next'; return; } }catch(_e){}
      }
      var b = findNav('next'); if(b && b.getAttribute('href')){ e.preventDefault(); location.href = b.getAttribute('href'); }
    }
    else if(k && k.toLowerCase() === 'a' && !e.ctrlKey && !e.metaKey && !e.altKey){ if(triggerAladin()){ e.preventDefault(); } }
    else if(k && k.toLowerCase() === 'g' && !e.ctrlKey && !e.metaKey && !e.altKey){
      // 'g' for goto: focus the page input if present; otherwise prompt and navigate
      try{
        var inp = document.querySelector('.nav form.jump input[type=number]');
        if(inp){ inp.focus(); inp.select(); e.preventDefault(); return; }
      }catch(_e){}
      try{ var total=_totalPages(); var defv=String((typeof PAGE_NUM==='number')?PAGE_NUM:1); var v = prompt('Go to page (1..'+ total +'):', defv); if(v!==null){ e.preventDefault(); gotoPage(v); return; } }catch(_e){}
    }
  }
  document.addEventListener('keydown', onKey, true);
})();
