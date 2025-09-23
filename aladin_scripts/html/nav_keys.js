// Keyboard navigation: ←/→ pages; 'a' to trigger Aladin.
(function(){
  var _mx, _my; document.addEventListener('mousemove', function(ev){ _mx=ev.clientX; _my=ev.clientY; }, {passive:true});
  function isEditable(el){ if(!el) return false; var tag=(el.tagName||'').toUpperCase(); return el.isContentEditable||tag==='INPUT'||tag==='TEXTAREA'||tag==='SELECT'; }
  function closest(el, sel){ while(el && el.nodeType===1){ if(el.matches && el.matches(sel)) return el; el = el.parentElement; } return null; }
  function findNav(which){
    var navs = document.querySelectorAll('.nav');
    for(var i=0;i<navs.length;i++){
      var as = navs[i].querySelectorAll('a[href]');
      for(var j=0;j<as.length;j++){
        var a = as[j]; var t = (a.textContent||'').trim();
        if(which==='prev' && t.indexOf('Prev Page') !== -1) return a;
        if(which==='next' && t.indexOf('Next Page') !== -1) return a;
      }
    }
    var m = (location.pathname||'').match(/([^\/]+)_page(\d+)\.html$/);
    if(m){ var base = m[1]; var n = parseInt(m[2],10); var target = which==='prev' ? (n-1) : (n+1); if(target>=1){ var a = document.createElement('a'); a.setAttribute('href', base + '_page' + target + '.html'); return a; } }
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
    if(k === 'ArrowLeft'){ var a = findNav('prev'); if(a && a.getAttribute('href')){ e.preventDefault(); location.href = a.getAttribute('href'); } }
    else if(k === 'ArrowRight'){ var b = findNav('next'); if(b && b.getAttribute('href')){ e.preventDefault(); location.href = b.getAttribute('href'); } }
    else if(k && k.toLowerCase() === 'a' && !e.ctrlKey && !e.metaKey && !e.altKey){ if(triggerAladin()){ e.preventDefault(); } }
  }
  document.addEventListener('keydown', onKey, true);
})();

// Position the per-panel Skip checkbox along the upper edge,
// vertically aligned with the Aladin/Aladin Lite buttons, and
// horizontally centered over the panel.
(function(){
  function alignSkipBoxes(){
    try{
      var wrap = document.querySelector('.wrap');
      if(!wrap) return;
      var wrect = wrap.getBoundingClientRect();
      var yTarget = null;
      var al = document.querySelector('a.btn.al');
      if(al){ try{ yTarget = al.getBoundingClientRect().top - wrect.top; }catch(e){} }
      if(yTarget === null){
        var lite = document.querySelector('a.btn.lite');
        if(lite){ try{ yTarget = lite.getBoundingClientRect().top - wrect.top; }catch(e){} }
      }
      if(yTarget === null){
        // Fallback to ~4.8% of image height if buttons are missing
        yTarget = wrect.height * 0.048;
      }
      var panels = document.querySelectorAll('.panelbox');
      for(var i=0;i<panels.length;i++){
        var p = panels[i];
        var prect = p.getBoundingClientRect();
        var topWithin = yTarget - (prect.top - wrect.top);
        if(!isFinite(topWithin)) topWithin = 8;
        // Clamp within panel bounds
        if(topWithin < 4) topWithin = 4;
        var maxTop = prect.height - 4;
        if(topWithin > maxTop) topWithin = Math.max(4, maxTop);
        var cb = p.querySelector('input.skipbox');
        if(!cb) continue;
        cb.style.position = 'absolute';
        cb.style.left = '50%';
        cb.style.top = String(topWithin) + 'px';
        cb.style.transform = 'translateX(-50%)';
        cb.style.margin = '0';
        cb.style.zIndex = '100';
        // Ensure the overlay mask does not block clicks on the checkbox
        var mask = p.querySelector('.mask');
        if(mask){ mask.style.pointerEvents = 'none'; }
      }
    }catch(e){}
  }
  if(document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', alignSkipBoxes);
  } else {
    alignSkipBoxes();
  }
  window.addEventListener('resize', function(){ try{ alignSkipBoxes(); }catch(e){} }, {passive:true});
})();
