
// Keyboard navigation: ←/→ pages; 'a' to trigger Aladin. Bounds-aware.
(function(){
  var KEY_DEBUG_URL = '/api/key_debug';
  var _mx = null;
  var _my = null;

  document.addEventListener('mousemove', function(ev){
    _mx = ev.clientX;
    _my = ev.clientY;
  }, {passive:true});

  function isEditable(el){
    if(!el) return false;
    var tag = (el.tagName || '').toUpperCase();
    return el.isContentEditable || tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT';
  }

  function closest(el, sel){
    while(el && el.nodeType === 1){
      if(el.matches && el.matches(sel)) return el;
      el = el.parentElement;
    }
    return null;
  }

  function _autoKey(){
    try{
      var base = (typeof PAGE_KEY === 'string' && PAGE_KEY) ? PAGE_KEY : 'global';
      return 'AUTO_ALADIN_' + base;
    }catch(e){
      return 'AUTO_ALADIN_global';
    }
  }

  function _consumeAutoSuppress(){
    var suppressed = false;
    try{
      suppressed = sessionStorage.getItem('AUTO_ALADIN_SUPPRESS') === '1';
      if(suppressed){
        sessionStorage.removeItem('AUTO_ALADIN_SUPPRESS');
      }
    }catch(e){
      suppressed = false;
    }
    return suppressed;
  }

  function readAuto(){
    try{
      var key = _autoKey();
      var val = null;
      if(window.localStorage){
        val = localStorage.getItem(key);
        if(val === null){
          val = localStorage.getItem('AUTO_ALADIN');
        }
      }
      return val === '1';
    }catch(e){
      return false;
    }
  }

  function writeAuto(on){
    try{
      var key = _autoKey();
      if(window.localStorage){
        localStorage.setItem(key, on ? '1' : '0');
        localStorage.setItem('AUTO_ALADIN', on ? '1' : '0');
      }
    }catch(e){}
    if(!on){
      try{ sessionStorage.removeItem('AUTO_ALADIN_SUPPRESS'); }catch(_e){}
    }
  }

  function triggerAladin(){
    var cand = null;
    try{
      if(typeof _mx === 'number' && typeof _my === 'number'){
        var el = document.elementFromPoint(_mx, _my);
        cand = closest(el, 'a.btn.al[data-ajs]');
      }
    }catch(e){}
    if(!cand){
      try{ cand = document.querySelector('a.btn.al[data-ajs]'); }catch(e){}
    }
    if(!cand){
      try{ cand = document.querySelector('[data-auto-aladin="1"]'); }catch(e){}
    }
    if(cand){
      try{
        cand.click();
        return true;
      }catch(e){}
    }
    try{
      if(Array.isArray(window.ALADIN_ITEMS) && window.ALADIN_ITEMS.length){
        var idx = (typeof window._srcIdx === 'number') ? window._srcIdx : 0;
        if(idx < 0 || idx >= window.ALADIN_ITEMS.length){
          idx = 0;
        }
        var ajs = window.ALADIN_ITEMS[idx] || window.ALADIN_ITEMS[0];
        var base = (typeof window.SAMP_URL === 'string' && window.SAMP_URL) ? window.SAMP_URL : 'http://127.0.0.1:8765';
        if(base && base.charAt(base.length - 1) === '/'){
          base = base.slice(0, -1);
        }
        var url = base + '/run_samp?file=' + encodeURIComponent(ajs) + '&_t=' + Date.now();
        var sink = document.getElementsByName('aladin_sink')[0];
        if(sink){
          try{ sink.src = url; return true; }catch(e){}
        }
        window.open(url, '_blank');
        return true;
      }
    }catch(e){}
    return false;
  }

  function findNav(dir){
    dir = (dir === 'prev') ? 'prev' : 'next';
    try{
      var anchors = document.querySelectorAll('.nav a[href]');
      for(var i=0;i<anchors.length;i++){
        var a = anchors[i];
        if(!a) continue;
        var label = (a.textContent || '').toLowerCase();
        if(a.dataset && a.dataset.nav){
          var ds = String(a.dataset.nav || '').toLowerCase();
          if(ds === dir) return a;
        }
        if(dir === 'prev'){
          if(label.indexOf('prev page') !== -1 || label.indexOf('previous page') !== -1) return a;
          if(label.indexOf('« prev') !== -1 || label.indexOf('prev «') !== -1) return a;
        }else if(dir === 'next'){
          if(label.indexOf('next page') !== -1) return a;
          if(label.indexOf('» next') !== -1 || label.indexOf('next »') !== -1) return a;
        }
      }
    }catch(e){}
    try{
      var btn = document.querySelector('.nav [data-nav="' + dir + '"]');
      if(btn){
        var href = btn.getAttribute('href') || btn.getAttribute('data-href');
        if(href){
          var auto = document.createElement('a');
          auto.setAttribute('href', href);
          return auto;
        }
      }
    }catch(e){}
    try{
      if(typeof PAGE_KEY === 'string' && PAGE_KEY){
        var current = Number(PAGE_NUM || 0);
        if(!Number.isFinite(current) || current < 1) current = 1;
        var total = _totalPages();
        var target = dir === 'prev' ? (current - 1) : (current + 1);
        if(Number.isFinite(total) && total > 0){
          if(target < 1 || target > total) return null;
        }else if(target < 1){
          return null;
        }
        if(target === current) return null;
        var a = document.createElement('a');
        var hash = '#seek=' + dir;
        a.setAttribute('href', PAGE_KEY + '_page' + target + '.html' + hash);
        return a;
      }
    }catch(e){}
    try{
      var path = String(location.pathname || '');
      var m = path.match(/([^\/]+)_page(\d+)\.html$/);
      if(m){
        var base = m[1];
        var num = parseInt(m[2], 10) || 0;
        var target2 = dir === 'prev' ? (num - 1) : (num + 1);
        if(target2 >= 1){
          var stub = document.createElement('a');
          stub.setAttribute('href', base + '_page' + target2 + '.html#seek=' + dir);
          return stub;
        }
      }
    }catch(e){}
    return null;
  }

  function _totalPages(){
    try{
      return (typeof TOTAL_PAGES === 'number' && TOTAL_PAGES > 0) ? TOTAL_PAGES : 1;
    }catch(e){
      return 1;
    }
  }

  function logKey(evt){
    try{
      var params = [
        'key=' + encodeURIComponent(evt.key || ''),
        'code=' + encodeURIComponent(evt.code || ''),
        'ctrl=' + (evt.ctrlKey ? '1' : '0'),
        'alt=' + (evt.altKey ? '1' : '0'),
        'shift=' + (evt.shiftKey ? '1' : '0'),
        'meta=' + (evt.metaKey ? '1' : '0'),
        'ts=' + Date.now()
      ].join('&');
      var img = new Image();
      img.src = KEY_DEBUG_URL + '?' + params;
    }catch(e){}
  }

  function _isProblemPage(list, page){
    if(!Array.isArray(list)) return false;
    var target = Number(page);
    if(!Number.isFinite(target)) return false;
    for(var i=0;i<list.length;i++){
      if(Number(list[i]) === target) return true;
    }
    return false;
  }

  function _hasUnskippedProblems(){
    try{
      var pb = document.querySelectorAll('.panelbox');
      for(var i=0;i<pb.length;i++){
        var p = pb[i];
        if(p.getAttribute('data-problem') === '1' && p.getAttribute('data-skip') !== '1'){
          return true;
        }
      }
    }catch(e){}
    return false;
  }

  function _maybeAutoAdvance(){
    try{
      var onlyOn = (function(){
        try{
          var key = 'ONLY_PROB_' + (PAGE_KEY || '');
          return localStorage.getItem(key) === '1';
        }catch(e){ return false; }
      })();
      if(!onlyOn) return false;
      if(_hasUnskippedProblems()) return false;
      var dir = (function(){
        try{
          var h = String(location.hash || '');
          var m = h.match(/seek=(next|prev)/);
          return m ? m[1] : 'next';
        }catch(e){ return 'next'; }
      })();
      var total = _totalPages();
      var current = Number(PAGE_NUM || 0);
      if(!Number.isFinite(current) || current < 1){ current = 1; }
      var target = (dir === 'prev') ? (current - 1) : (current + 1);
      if(target < 1 || target > total){
        return false;
      }
      if(target !== current){
        window.location.replace(PAGE_KEY + '_page' + target + '.html#seek=' + dir);
        return true;
      }
    }catch(e){}
    return false;
  }

  function maybeTriggerAuto(){
    try{
      var cb = document.getElementById('auto_aladin');
      if(cb && cb.checked){
        setTimeout(function(){ triggerAladin(); }, 200);
      }
    }catch(e){}
  }

  function collectToggles(){
    var buttons = Array.prototype.slice.call(document.querySelectorAll('.btn.tgl'));
    if(!buttons.length) return buttons;
    var panel = null;
    if(typeof _mx === 'number' && typeof _my === 'number'){
      var el = document.elementFromPoint(_mx, _my);
      panel = closest(el, '.panelbox');
      if(!panel){
        var btn = closest(el, '.btn.tgl');
        if(btn){ panel = closest(btn, '.panelbox'); }
      }
    }
    if(panel){
      var base = panel.getAttribute('data-base');
      var filtered = buttons.filter(function(btn){
        if(base){
          var id = btn.getAttribute('data-panel');
          if(id) return id === base;
        }
        return panel.contains(btn);
      });
      if(filtered.length) buttons = filtered;
    }
    buttons.sort(function(a, b){
      var ra = a.getBoundingClientRect();
      var rb = b.getBoundingClientRect();
      if(Math.abs(ra.top - rb.top) > 4) return ra.top - rb.top;
      return ra.left - rb.left;
    });
    return buttons;
  }

  function triggerToggle(rank){
    var buttons = collectToggles();
    if(!buttons.length) return false;
    var idx = Math.min(Math.max(rank - 1, 0), buttons.length - 1);
    var btn = buttons[idx];
    if(!btn) return false;
    var href = btn.getAttribute('data-href-force') || btn.getAttribute('href');
    if(!href){
      try{ btn.click(); return true; }catch(e){ return false; }
    }
    window.location.href = href;
    return true;
  }

  function setupAuto(skipImmediate){
    try{
      var cb = document.getElementById('auto_aladin');
      if(!cb) return;
      var on = readAuto();
      var suppressed = _consumeAutoSuppress();
      cb.checked = on;
      cb.addEventListener('change', function(){
        writeAuto(cb.checked);
        if(cb.checked){
          try{ sessionStorage.removeItem('AUTO_ALADIN_SUPPRESS'); }catch(_e){}
          var hopped = _maybeAutoAdvance();
          if(!hopped){
            maybeTriggerAuto();
          }
        }
      });
      if(on && !suppressed && !skipImmediate){
        maybeTriggerAuto();
      }
    }catch(e){}
  }

  function isProblemPanel(el){
    if(!el) return false;
    try{ return el.getAttribute('data-problem') === '1'; }catch(e){ return false; }
  }

  if(!Array.isArray(window.PROBLEM_PAGES)){
    window.PROBLEM_PAGES = [];
  }

  document.addEventListener('keydown', function(e){
    if(isEditable(e.target)) return;
    logKey(e);
    var k = (e.key || '');
    if(!e.ctrlKey && !e.metaKey && !e.altKey){
      if(k === '1' || k === '2' || k === '3'){
        if(triggerToggle(parseInt(k, 10))){
          e.preventDefault();
          return;
        }
      }
    }
    var only = document.getElementById('only_prob');
    var onlyOn = !!(only && only.checked);
    if(k === 'ArrowLeft'){
      if(onlyOn){
        e.preventDefault();
        var currentL = Number(PAGE_NUM || 0);
        if(!Number.isFinite(currentL) || currentL < 1){ currentL = 1; }
        if(currentL > 1){
          window.location.href = PAGE_KEY + '_page' + (currentL - 1) + '.html#seek=prev';
        }
        return;
      }
      var prev = findNav('prev');
      if(prev && prev.href){
        e.preventDefault();
        window.location.href = prev.href;
      }
    }
    else if(k === 'ArrowRight'){
      if(onlyOn){
        e.preventDefault();
        var currentR = Number(PAGE_NUM || 0);
        if(!Number.isFinite(currentR) || currentR < 1){ currentR = 1; }
        var totalPages = _totalPages();
        if(currentR < totalPages){
          window.location.href = PAGE_KEY + '_page' + (currentR + 1) + '.html#seek=next';
        }
        return;
      }
      var next = findNav('next');
      if(next && next.href){
        e.preventDefault();
        window.location.href = next.href;
      }
    }
    else if(k.toLowerCase() === 'a' && !e.ctrlKey && !e.metaKey && !e.altKey){
      if(triggerAladin()){ e.preventDefault(); }
    }
    else if(k.toLowerCase() === 'g' && !e.ctrlKey && !e.metaKey && !e.altKey){
      try{
        var inp = document.querySelector('.nav form.jump input[type=number]');
        if(inp){ inp.focus(); inp.select(); e.preventDefault(); return; }
      }catch(_e){}
      try{
        var total = _totalPages();
        var defv = String((typeof PAGE_NUM === 'number') ? PAGE_NUM : 1);
        var v = prompt('Go to page (1..' + total + '):', defv);
        if(v !== null){ e.preventDefault(); gotoPage(v); return; }
      }catch(_e){}
    }
  }, true);

  document.addEventListener('DOMContentLoaded', function(){
    try{
      const sk=_readSkip();
      document.querySelectorAll('.panelbox').forEach(function(p){
        const b=p.getAttribute('data-base');
        if(sk[b]){
          p.setAttribute('data-skip','1');
          const cb=p.querySelector('.skipbox');
          if(cb) cb.checked=true;
        }
      });
    }catch(e){}
    const only=document.getElementById('only_prob');
    if(only){
      only.checked = _readOnly();
      only.addEventListener('change', function(){
        _writeOnly(only.checked);
        updatePanels();
        if(only.checked){
          var hopped = _maybeAutoAdvance();
          if(!hopped){
            maybeTriggerAuto();
          }
        }
      });
    }
    updatePanels();
    var autoHopped = _maybeAutoAdvance();
    setupAuto(autoHopped);
    if(!autoHopped){
      maybeTriggerAuto();
    }
  }, true);
})();
