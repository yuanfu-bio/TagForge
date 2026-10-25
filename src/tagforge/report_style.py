"""Shared appearance and theme controls for TagForge HTML reports."""

REPORT_CSS = 'body{margin:0;background:#f5f7fb;color:#172033;font:15px Arial,sans-serif}main{max-width:1180px;margin:auto;padding:36px}h1{font-size:36px;margin:0}.section-title{margin:36px 0 16px;padding-top:28px;border-top:1px solid #dce3ee}.sub{color:#64748b}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:14px;margin:28px 0}.mini-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:14px;margin:14px 0}.pi-column{min-width:0}.pi-metrics{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px;margin:0 0 14px}.pi-top-panel h3{margin:0 0 8px}.pi-top-panel h4{margin:8px 0 0}.report-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:14px;align-items:stretch;margin:16px 0}.report-grid .panel{margin:0;display:flex;flex-direction:column;min-width:0;overflow:hidden}#sat,#features{width:100%;min-width:0}.scroll-window{height:520px;overflow:auto;border:1px solid #e6ebf2;border-radius:8px}.scroll-window th{position:sticky;top:0;z-index:1}.icon-button{display:inline-flex;align-items:center;justify-content:center;width:36px;height:36px;padding:0;border:1px solid #94a3b8;border-radius:8px;background:white;color:#172033;cursor:pointer;font-size:19px}.card,.panel{background:white;border:1px solid #dce3ee;border-radius:14px;padding:18px}.value{font-size:28px;font-weight:750;color:#155e75}.panel{margin:16px 0;overflow:auto}table{border-collapse:collapse;width:100%}td,th{padding:9px;border-bottom:1px solid #e6ebf2;text-align:left}th{color:#475569;background:#f8fafc}'

DARK_CSS = '.theme{float:right}body.dark{background:#07111f;color:#e5edf8}body.dark .card,body.dark .panel{background:#0f1d30;border-color:#24344d}body.dark .sub{color:#94a3b8}body.dark th{background:#17263b;color:#93c5fd}body.dark td,body.dark th,body.dark .scroll-window{border-color:#24344d}body.dark .value{color:#5eead4}body.dark .section-title{border-color:#24344d}body.dark .icon-button{background:#0f1d30;color:#e5edf8;border-color:#64748b}'

THEME_JS = (
    "function toggleTheme(button,ids){const dark=document.body.classList.toggle('dark');"
    "button.textContent=dark?'☀':'☾';"
    "button.setAttribute('aria-label',dark?'Switch to light theme':'Switch to dark theme');button.title=button.getAttribute('aria-label');"
    "const grid=dark?'#24344d':'#e2e8f0';const surface=dark?'#0f1d30':'#ffffff';const font=dark?'#cbd5e1':'#172033';"
    "ids.forEach(id=>{if(document.getElementById(id))Plotly.relayout(id,{'paper_bgcolor':surface,'plot_bgcolor':surface,"
    "'font.color':font,'font.family':'Arial','xaxis.gridcolor':grid,'yaxis.gridcolor':grid,'yaxis2.gridcolor':grid,'yaxis3.gridcolor':grid})})}"
)
