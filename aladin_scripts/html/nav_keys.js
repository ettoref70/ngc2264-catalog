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
