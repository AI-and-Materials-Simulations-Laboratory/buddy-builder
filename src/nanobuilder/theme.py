"""Visual theme for the GUI: one scoped stylesheet, light and dark palettes.

Everything is scoped to ``.nb-app`` so the tool never restyles the rest of a
notebook. The palette re-points Jupyter's own widget CSS variables, which both
ipywidgets 7 (classic Notebook / nanoHUB AppMode) and 8 (JupyterLab, Notebook
7, Voila) read, plus a few direct selectors for the parts those variables do
not reach (tabs, buttons, sliders).

Themes: ``nb-buddy`` (pastel yellow, the default), ``nb-buddy2`` (mint
green), ``nb-light`` and ``nb-dark``. The accent colour (green in the Buddy
themes, gold in Light/Dark) is reserved for the action to take next
(Generate, Download ZIP), the current tab and checked controls.
"""

from __future__ import annotations

FONT_URL = "https://fonts.googleapis.com/css2?family=Fredoka:wght@400;500;600;700&display=swap"
FONT_FILE = __import__("pathlib").Path(__file__).resolve().parent / "assets" / "fonts" / "Fredoka-latin.woff2"


def font_face_css() -> str:
    """@font-face for the bundled Fredoka (SIL OFL), so titles use it even offline."""
    try:
        import base64
        data = base64.b64encode(FONT_FILE.read_bytes()).decode("ascii")
    except OSError:
        return ""
    return ("@font-face { font-family: 'Fredoka'; font-style: normal; font-weight: 300 700; "
            "font-stretch: 75% 125%; font-display: swap; "
            f"src: url(data:font/woff2;base64,{data}) format('woff2'); }}")
# Small text (labels, captions, inputs): Arial. Titles, headings, tabs and buttons: Fredoka.
FONT_STACK = 'Arial, "Helvetica Neue", Helvetica, sans-serif'
DISPLAY_STACK = '"Fredoka", "Baloo 2", "Nunito", "Arial Rounded MT Bold", Arial, sans-serif'

LIGHT = {
    "bg": "#EDF0F4",        # page: cool lab grey
    "panel": "#FFFFFF",     # cards, tab body
    "input": "#FFFFFF",
    "ink": "#1A2230",
    "muted": "#5C6778",
    "rule": "#D5DBE3",
    "au": "#8A5A00",        # deep gold: primary action / selection
    "au_ink": "#FFFFFF",    # text on gold
    "au_soft": "#F6EBD3",
    "stage": "#FFFFFF",     # viewer background
    "info": "#2C5AA0",
    "ok": "#1F6B3A",
    "warn": "#94520A",
    "err": "#B3261E",
}

DARK = {
    "bg": "#0E141C",
    "panel": "#151D28",
    "input": "#0F1620",
    "ink": "#E4E9F0",
    "muted": "#97A3B4",
    "rule": "#2A3646",
    "au": "#E3A83B",
    "au_ink": "#1B1405",
    "au_soft": "#3A2C10",
    "stage": "#0B1017",
    "info": "#8DB2EE",
    "ok": "#79D29E",
    "warn": "#F2B76E",
    "err": "#FF8F86",
}


# Buddy: the logo's pastel yellow, pink and green; deep green carries the primary action.
BUDDY = {
    "bg": "#FBF6DF",        # page: pale butter yellow
    "panel": "#FFFDF6",     # cards, tab body: warm white
    "input": "#FFFFFF",
    "ink": "#2B3A2F",       # deep green-grey text
    "muted": "#66705F",
    "rule": "#E8DFBE",
    "au": "#2E7D55",        # primary action / selection (green, 5:1 on white)
    "au_ink": "#FFFFFF",
    "au_soft": "#DDF1E2",   # pastel green
    "stage": "#FFFDF6",
    "info": "#3B6EA5",
    "ok": "#2E7D55",
    "warn": "#A0661A",
    "err": "#B3261E",
    "pink": "#F7D9E1",      # pastel pink accents
    "yellow": "#FFF1B8",
}


# Buddy 2: the same family on a light mint-green page.
BUDDY2 = dict(BUDDY, bg="#E6F4EA", panel="#FAFEFB", rule="#CBE4D3", ink="#22362A",
              muted="#5A6D60", au_soft="#D2EEDA", stage="#FAFEFB")


def _vars(p: dict) -> str:
    own = "".join(f"--nb-{k.replace('_', '-')}:{v};" for k, v in p.items())
    jupyter = (
        "--jp-widgets-color:var(--nb-ink);--jp-widgets-label-color:var(--nb-ink);"
        "--jp-widgets-readout-color:var(--nb-ink);--jp-widgets-font-size:13px;"
        "--jp-widgets-input-color:var(--nb-ink);--jp-widgets-input-background-color:var(--nb-input);"
        "--jp-widgets-input-border-color:var(--nb-rule);--jp-widgets-input-focus-border-color:var(--nb-au);"
        "--jp-widgets-slider-handle-border-color:var(--nb-muted);"
        "--jp-widgets-slider-handle-background-color:var(--nb-panel);"
        "--jp-widgets-slider-active-handle-color:var(--nb-au);"
        "--jp-layout-color0:var(--nb-panel);--jp-layout-color1:var(--nb-panel);"
        "--jp-layout-color2:var(--nb-bg);--jp-layout-color3:var(--nb-rule);"
        "--jp-border-color0:var(--nb-rule);--jp-border-color1:var(--nb-rule);--jp-border-color2:var(--nb-rule);"
        "--jp-ui-font-color0:var(--nb-ink);--jp-ui-font-color1:var(--nb-ink);--jp-ui-font-color2:var(--nb-muted);"
        "--jp-content-font-color1:var(--nb-ink);--jp-brand-color1:var(--nb-au);"
        f"--jp-ui-font-family:{FONT_STACK};--jp-widgets-font-family:{FONT_STACK};"
    )
    return own + jupyter


def _page_dark(scope: str) -> str:
    """Dark page background around the tool in Voila, JupyterLab and classic AppMode."""
    containers = ("#main, .jp-LabShell, #jp-main-content-panel, #rendered_cells, .jp-Notebook, .jp-Cell, .jp-OutputArea, .jp-OutputArea-output, "
                  "#notebook, #notebook-container, .container, .output_area, .output_subarea, "
                  ".jp-RenderedHTMLCommon, .rendered_html")
    sel = ", ".join(f"{scope} {c.strip()}" for c in containers.split(","))
    return (f"{scope} {{ background: {DARK['bg']} !important; color: {DARK['ink']}; }} "
            f"{sel} {{ background: transparent !important; box-shadow: none !important; "
            f"color: {DARK['ink']}; }}")


def stylesheet() -> str:
    light, dark, buddy, buddy2 = _vars(LIGHT), _vars(DARK), _vars(BUDDY), _vars(BUDDY2)
    return f"""<style>
{font_face_css()}
.nb-app {{ {light} }}
.nb-app.nb-dark {{ {dark} }}
.nb-app.nb-buddy {{ {buddy} }}
.nb-app.nb-buddy2 {{ {buddy2} }}
.nb-app.nb-buddy2 .nb-head {{ background: linear-gradient(100deg, #E3F6E8 0%, {BUDDY['yellow']} 55%, {BUDDY['pink']} 100%);
  border-radius: 16px; padding: 10px 16px; border: 1px solid {BUDDY2['rule']}; }}
.nb-app.nb-buddy2 .jupyter-button {{ border-radius: 10px; }}
/* Buddy theme: pastel header band and soft cards */
.nb-app.nb-buddy .nb-head {{ background: linear-gradient(100deg, {BUDDY['yellow']} 0%, {BUDDY['pink']} 55%, {BUDDY['au_soft']} 100%);
  border-radius: 16px; padding: 10px 16px; border: 1px solid {BUDDY['rule']}; }}
.nb-app.nb-buddy .widget-tab > .widget-tab-contents, .nb-app.nb-buddy .jupyter-widget-tab > .widget-tab-contents {{
  border-radius: 0 0 14px 14px; }}
.nb-app.nb-buddy .jupyter-button {{ border-radius: 10px; }}
.nb-app.nb-buddy .jupyter-button.mod-active {{ background: {BUDDY['au_soft']}; border-color: {BUDDY['au']}; }}
.nb-app.nb-buddy .nb-sec {{ background: {BUDDY['pink']}33; }}
@media (prefers-color-scheme: dark) {{ .nb-app.nb-auto {{ {dark} }} }}
/* Page around the tool follows the theme (Voila / AppMode show only the tool). */
{_page_dark("body:has(.nb-app.nb-dark)")}
@media (prefers-color-scheme: dark) {{ {_page_dark("body:has(.nb-app.nb-auto)")} }}

.nb-app {{ background: var(--nb-bg); color: var(--nb-ink); font-family: {FONT_STACK};
  font-size: 13px; line-height: 1.45; padding: 18px 22px 22px; border-radius: 14px; box-sizing: border-box; }}
.nb-app, .nb-app .widget-html-content, .nb-app .widget-label, .nb-app label {{ color: var(--nb-ink); }}
.nb-app .widget-html-content {{ font-family: inherit; font-size: 13px; line-height: 1.45; }}
.nb-app small {{ color: var(--nb-muted); font-size: 12px; }}
.nb-app a {{ color: var(--nb-info); }}

/* Header */
.nb-app .nb-head {{ align-items: center; justify-content: space-between; margin-bottom: 14px; gap: 12px; flex-wrap: wrap; }}
.nb-app .nb-brand {{ display: flex; align-items: center; gap: 14px; }}
.nb-app .nb-logo {{ width: 72px; height: 72px; border-radius: 50%; flex: 0 0 auto;
  box-shadow: 0 1px 4px rgba(0,0,0,0.12); }}
.nb-app .nb-title {{ font-size: 24px; font-weight: 600; letter-spacing: -0.01em; margin: 0; }}
.nb-app .nb-sub {{ color: var(--nb-muted); margin: 2px 0 0; max-width: 72ch; }}
.nb-app .nb-switches {{ display: flex; flex-direction: column; align-items: flex-end; gap: 6px; }}
.nb-app .nb-switch-label {{ color: var(--nb-muted); font-size: 11px; margin-right: 6px; }}

/* Sections inside a tab; Advanced expander */
.nb-app .nb-sec-h {{ font-size: 12px; font-weight: 600; color: var(--nb-muted); text-transform: uppercase;
  letter-spacing: 0.04em; margin: 12px 0 4px; }}
.nb-app .nb-adv-h {{ font-size: 12px; font-weight: 600; color: var(--nb-au); margin: 14px 0 2px;
  border-top: 1px dashed var(--nb-rule); padding-top: 8px; }}
.nb-app .nb-row {{ flex-wrap: wrap; gap: 4px 12px; }}
.nb-app .widget-toggle-buttons {{ flex-wrap: wrap; height: auto; }}
.nb-app .widget-toggle-buttons .widget-toggle-button {{ margin: 2px 6px 2px 0; }}
.nb-app .nb-choices .jupyter-button {{ height: 34px; padding: 0 12px; display: inline-flex;
  align-items: center; justify-content: center; gap: 7px; }}
.nb-app .nb-choices .jupyter-button i {{ order: -1; margin: 0; }}

/* Console */
.nb-app .bb-log {{ height: 170px; overflow: auto; display: flex; flex-direction: column-reverse;
  background: var(--nb-input); border: 1px solid var(--nb-rule); border-radius: 8px; padding: 6px 8px;
  font: 12px/1.45 ui-monospace, Menlo, Consolas, monospace; box-sizing: border-box; }}
.nb-app .bb-log > pre {{ margin: 0; font: inherit; white-space: pre-wrap; flex: 0 0 auto; }}
.nb-app .bb-log .bb-in {{ color: var(--nb-au); }}
.nb-app .bb-log .bb-err {{ color: var(--nb-err); }}
.nb-app .bb-log .bb-echo {{ color: var(--nb-muted); }}
.nb-app .bb-console-input input {{ font-family: ui-monospace, Menlo, Consolas, monospace !important; }}

/* Columns */
.nb-app .nb-left {{ padding-right: 20px; box-sizing: border-box; }}
.nb-app .nb-right {{ position: sticky; top: 10px; align-self: flex-start; box-sizing: border-box; }}
.nb-app .nb-h {{ font-size: 13px; font-weight: 600; margin: 16px 0 4px; }}

/* Labels: left-aligned and quiet, values (and checkbox text) in ink */
.nb-app .widget-label {{ color: var(--nb-muted); }}
.nb-app .widget-label-basic, .nb-app .widget-label-basic span {{ color: var(--nb-ink); }}
.nb-app .widget-inline-hbox > .widget-label {{ text-align: left; padding-right: 10px; }}
.nb-app .widget-checkbox {{ width: auto; }}

/* Inputs */
.nb-app input[type=text], .nb-app input[type=number], .nb-app input[type=password],
.nb-app select, .nb-app textarea {{
  background: var(--nb-input); color: var(--nb-ink); border: 1px solid var(--nb-rule);
  border-radius: 6px; font-family: inherit; font-size: 13px; box-sizing: border-box; }}
.nb-app input[type=text], .nb-app input[type=number], .nb-app input[type=password],
.nb-app .widget-dropdown select {{ height: 30px; padding: 0 8px; }}
.nb-app .widget-select select {{ padding: 4px; height: auto; }}
.nb-app .widget-select select option {{ padding: 3px 6px; }}
.nb-app .widget-select select option:checked {{ background: var(--nb-au-soft); color: var(--nb-ink); }}
.nb-app input:focus-visible, .nb-app select:focus-visible, .nb-app textarea:focus-visible,
.nb-app button:focus-visible {{ outline: 2px solid var(--nb-au); outline-offset: 1px; }}
.nb-app input[type=checkbox], .nb-app input[type=radio] {{ accent-color: var(--nb-au); }}
.nb-app .widget-dropdown select {{ appearance: auto; }}

/* Buttons: neutral by default, gold only for the primary action */
.nb-app .widget-button, .nb-app .jupyter-button {{ width: auto; min-width: 0; }}
.nb-app .widget-vbox > .widget-button {{ align-self: flex-start; }}
.nb-app .jupyter-button {{ background: var(--nb-panel); color: var(--nb-ink); border: 1px solid var(--nb-rule);
  border-radius: 6px; height: 32px; padding: 0 14px; font-family: inherit; font-size: 13px;
  font-weight: 500; box-shadow: none; }}
.nb-app .jupyter-button:hover:enabled {{ border-color: var(--nb-muted); box-shadow: none; }}
.nb-app .jupyter-button.mod-primary {{ background: var(--nb-au); border-color: var(--nb-au); color: var(--nb-au-ink); }}
.nb-app .jupyter-button.mod-primary:hover:enabled {{ filter: brightness(1.08); border-color: var(--nb-au); }}
.nb-app .jupyter-button.mod-info, .nb-app .jupyter-button.mod-success, .nb-app .jupyter-button.mod-warning {{
  background: var(--nb-panel); color: var(--nb-ink); }}
.nb-app .jupyter-button.mod-active {{ background: var(--nb-au-soft); border-color: var(--nb-au); color: var(--nb-ink); }}
.nb-app .jupyter-button:disabled {{ opacity: 0.5; }}

/* Theme switch: a compact segmented control */
.nb-app .nb-theme .widget-toggle-buttons {{ gap: 0; }}
.nb-app .nb-theme .jupyter-button {{ border-radius: 0; margin: 0 -1px 0 0; height: 28px; padding: 0 10px; width: auto; }}
.nb-app .nb-theme .jupyter-button:first-child {{ border-radius: 6px 0 0 6px; }}
.nb-app .nb-theme .jupyter-button:last-child {{ border-radius: 0 6px 6px 0; }}

/* Tabs: clear buttons with bold labels; the current step is filled */
.nb-app .lm-TabBar, .nb-app .p-TabBar {{ border-bottom: none; min-height: 0; height: auto !important;
  overflow: visible; }}
.nb-app .lm-TabBar-content, .nb-app .p-TabBar-content {{ display: grid !important;
  grid-template-columns: repeat(auto-fill, minmax(140px, 1fr)); gap: 6px; padding: 2px 0 8px;
  border-bottom: none !important; }}
.nb-app .lm-TabBar-tab, .nb-app .p-TabBar-tab {{ background: var(--nb-panel) !important;
  border: 1.5px solid var(--nb-rule) !important; border-radius: 10px !important;
  color: var(--nb-ink) !important; font-family: {DISPLAY_STACK} !important; font-weight: 600 !important;
  font-size: 14.5px !important; min-width: 0 !important; max-width: none !important;
  padding: 0 8px !important; margin: 0 !important; height: 38px; box-sizing: border-box;
  box-shadow: 0 1px 2px rgba(0,0,0,0.06); cursor: pointer;
  display: flex !important; align-items: center !important; justify-content: center !important; }}
/* Only the label: the (empty) icon and close boxes would push it off centre */
.nb-app .lm-TabBar-tabIcon, .nb-app .p-TabBar-tabIcon,
.nb-app .lm-TabBar-tabCloseIcon, .nb-app .p-TabBar-tabCloseIcon {{ display: none !important; }}
.nb-app .lm-TabBar-tabLabel, .nb-app .p-TabBar-tabLabel {{ flex: 0 1 auto !important; text-align: center;
  line-height: 1.1 !important; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
.nb-app .lm-TabBar-tab:hover, .nb-app .p-TabBar-tab:hover {{ border-color: var(--nb-au) !important; }}
.nb-app .lm-TabBar-tab.lm-mod-current, .nb-app .p-TabBar-tab.p-mod-current {{
  background: var(--nb-au) !important; border-color: var(--nb-au) !important;
  color: var(--nb-au-ink) !important; transform: none !important; }}
.nb-app .lm-TabBar-tab.lm-mod-current::before, .nb-app .p-TabBar-tab.p-mod-current::before {{ display: none; }}
.nb-app .lm-TabBar-tabLabel, .nb-app .p-TabBar-tabLabel {{ font-weight: 600 !important; }}
.nb-app .widget-tab > .widget-tab-contents, .nb-app .jupyter-widget-tab > .widget-tab-contents {{
  background: var(--nb-panel); border: 1px solid var(--nb-rule); border-radius: 12px;
  padding: 16px 16px 14px; overflow: visible !important; }}
.nb-app .jupyter-widget-tab, .nb-app .widget-tab, .nb-app .lm-StackedPanel, .nb-app .p-StackedPanel {{
  overflow: visible !important; }}

/* Browse menu: hover a group, then a material */
.nb-app .bb-menu {{ position: relative; display: inline-block; z-index: 40; }}
.nb-app .bb-menu-btn {{ font-family: {DISPLAY_STACK}; font-size: 14px; font-weight: 600; color: var(--nb-ink);
  background: var(--nb-input);
  border: 1.5px solid var(--nb-au); border-radius: 10px; height: 34px; padding: 0 14px; cursor: pointer; }}
.nb-app .bb-menu ul {{ list-style: none; margin: 0; padding: 4px 0; background: var(--nb-panel);
  border: 1px solid var(--nb-rule); border-radius: 10px; box-shadow: 0 6px 18px rgba(0,0,0,0.14); }}
.nb-app .bb-groups {{ display: none; position: absolute; left: 0; top: 34px; min-width: 250px; }}
.nb-app .bb-menu:hover .bb-groups, .nb-app .bb-menu:focus-within .bb-groups {{ display: block; }}
.nb-app .bb-menu.bb-closed .bb-groups {{ display: none !important; }}
.nb-app .bb-group {{ position: relative; padding: 7px 14px; cursor: default; white-space: nowrap; }}
.nb-app .bb-group > span::after {{ content: "▸"; float: right; margin-left: 14px; color: var(--nb-muted); }}
.nb-app .bb-group:hover {{ background: var(--nb-au-soft); }}
.nb-app .bb-items {{ display: none; position: absolute; left: 100%; top: -5px; min-width: 260px;
  max-height: 420px; overflow-y: auto; }}
.nb-app .bb-group:hover > .bb-items {{ display: block; }}
.nb-app .bb-items a {{ display: block; padding: 6px 14px; color: var(--nb-ink); text-decoration: none;
  cursor: pointer; white-space: nowrap; }}
.nb-app .bb-items a:hover {{ background: var(--nb-au-soft); }}

/* Sad Buddy: shows for a few seconds when a build fails (with why and how to fix), then fades */
.nb-app .bb-sad {{ position: fixed; left: 50%; top: 45%; z-index: 10000; pointer-events: none;
  text-align: center; color: var(--nb-ink); width: min(440px, 90vw);
  animation: bbsad 4.2s ease forwards; }}
.nb-app .bb-sad img {{ width: 200px; height: 200px; border-radius: 50%; display: block; margin: 0 auto 10px;
  box-shadow: 0 10px 30px rgba(0,0,0,0.25); background: #fff; }}
.nb-app .bb-sad-card {{ background: var(--nb-panel); border: 1.5px solid var(--nb-rule); border-radius: 14px;
  padding: 10px 16px 6px; box-shadow: 0 8px 24px rgba(0,0,0,0.18); text-align: left;
  font-family: {FONT_STACK}; font-size: 13px; line-height: 1.4; }}
.nb-app .bb-sad-title {{ font-family: {DISPLAY_STACK}; font-size: 18px; font-weight: 600;
  text-align: center; margin-bottom: 6px; }}
.nb-app .bb-sad-card p {{ margin: 0 0 6px; }}
@keyframes bbsad {{ 0% {{ opacity: 0; transform: translate(-50%, -50%) scale(0.8); }}
  6% {{ opacity: 1; transform: translate(-50%, -50%) scale(1); }}
  88% {{ opacity: 1; transform: translate(-50%, -50%) scale(1); }}
  100% {{ opacity: 0; visibility: hidden; transform: translate(-50%, -50%) scale(0.95); }} }}

/* Display font: titles, headings, tabs and buttons (Fredoka); small text stays Arial */
.nb-app .nb-title, .nb-app .nb-sec-h, .nb-app .nb-h, .nb-app .nb-adv-h,
.nb-app .jupyter-button, .nb-app .widget-button, .nb-app .widget-toggle-button,
.nb-app .nb-dl a, .nb-app .bb-menu-btn {{ font-family: {DISPLAY_STACK} !important; }}
.nb-app .nb-title {{ font-size: 30px; font-weight: 600; letter-spacing: 0; }}
.nb-app .nb-sec-h {{ font-size: 14px; font-weight: 600; letter-spacing: 0.02em; text-transform: none;
  color: var(--nb-ink); }}
.nb-app .nb-h {{ font-size: 16px; font-weight: 600; }}
.nb-app .nb-adv-h {{ font-size: 14px; }}
.nb-app .jupyter-button, .nb-app .widget-toggle-button {{ font-size: 14px !important; font-weight: 500; }}
.nb-app .jupyter-button.mod-primary {{ font-weight: 600; }}
.nb-app .nb-dl a {{ font-size: 14px; }}
.nb-app .nb-dl .nb-size {{ font-family: {FONT_STACK}; }}

/* Unit-cell table */
.nb-app .nb-cell td:first-child {{ width: 34%; }}
.nb-app .bb-console-input textarea {{ font-family: ui-monospace, Menlo, Consolas, monospace !important;
  font-size: 12px; min-height: 40px; resize: vertical; }}

/* Slider, progress */
.nb-app .noUi-connect, .nb-app .ui-slider-range {{ background: var(--nb-au); }}
.nb-app .widget-hprogress .progress {{ background: var(--nb-rule); border-radius: 3px; height: 6px; }}
.nb-app .widget-hprogress .progress-bar {{ background: var(--nb-au); border-radius: 3px; }}
.nb-app .widget-hprogress {{ align-items: center; }}

/* Messages: a coloured rule says what kind; text stays readable */
.nb-app .nb-msg {{ border-left: 3px solid var(--nb-muted); background: var(--nb-panel); padding: 7px 10px;
  border-radius: 0 6px 6px 0; margin: 6px 0; }}
.nb-app .nb-msg ul {{ margin: 0 0 0 16px; padding: 0; }}
.nb-app .nb-msg li + li {{ margin-top: 2px; }}
.nb-app .nb-msg-ok {{ border-left-color: var(--nb-ok); }}
.nb-app .nb-msg-info {{ border-left-color: var(--nb-info); }}
.nb-app .nb-msg-warning {{ border-left-color: var(--nb-warn); }}
.nb-app .nb-msg-error {{ border-left-color: var(--nb-err); color: var(--nb-err); }}
.nb-app details summary {{ cursor: pointer; color: var(--nb-muted); }}
.nb-app pre {{ white-space: pre-wrap; font-size: 11px; }}

/* Results, downloads, viewer stage */
.nb-app .nb-table {{ border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }}
.nb-app .nb-table td {{ padding: 5px 0; border-bottom: 1px solid var(--nb-rule); vertical-align: top; }}
.nb-app .nb-table td:first-child {{ color: var(--nb-muted); width: 44%; padding-right: 12px; }}
.nb-app .nb-dl {{ display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }}
.nb-app .nb-dl a {{ display: inline-flex; align-items: baseline; gap: 6px; padding: 6px 12px;
  border: 1px solid var(--nb-rule); border-radius: 6px; background: var(--nb-panel);
  color: var(--nb-ink); text-decoration: none; font-weight: 500; }}
.nb-app .nb-dl a:hover {{ border-color: var(--nb-muted); }}
.nb-app .nb-dl a.nb-primary {{ background: var(--nb-au); border-color: var(--nb-au); color: var(--nb-au-ink); }}
.nb-app .nb-dl .nb-size {{ font-size: 11px; color: var(--nb-muted); font-weight: 400; font-variant-numeric: tabular-nums; }}
.nb-app .nb-dl a.nb-primary .nb-size {{ color: inherit; opacity: 0.8; }}
.nb-app iframe.nb-stage {{ border: 1px solid var(--nb-rule); border-radius: 10px; background: var(--nb-stage); display: block; }}
.nb-app .nb-placeholder {{ border: 1px dashed var(--nb-rule); border-radius: 10px; color: var(--nb-muted);
  background: var(--nb-panel); display: flex; align-items: center; justify-content: center;
  text-align: center; padding: 12px; box-sizing: border-box; }}
.nb-app .nb-svg svg {{ background: #fff; border-radius: 8px; border: 1px solid var(--nb-rule); max-width: 100%; height: auto; }}

/* Ask Buddy chat (same tokens and fonts as the rest of the app) */
.nb-app .bb-chat {{ background: var(--nb-panel); border: 1.5px solid var(--nb-rule); border-radius: 14px;
  padding: 4px 12px 10px; box-sizing: border-box; }}
.nb-app .bb-chat .bb-chat-h {{ margin: 8px 0 4px; }}
.nb-app .bb-chat-scroll {{ max-height: 300px; min-height: 90px; overflow: auto; display: flex;
  flex-direction: column-reverse; gap: 8px; background: var(--nb-bg); border: 1px solid var(--nb-rule);
  border-radius: 10px; padding: 8px; box-sizing: border-box; }}
.nb-app .bb-turn {{ display: flex; flex-direction: column; max-width: 92%; flex: 0 0 auto; }}
.nb-app .bb-turn-user {{ align-self: flex-end; align-items: flex-end; }}
.nb-app .bb-turn-buddy {{ align-self: flex-start; }}
.nb-app .bb-who {{ font-family: {DISPLAY_STACK}; font-size: 11px; font-weight: 600; color: var(--nb-muted);
  margin: 0 6px 2px; }}
.nb-app .bb-bubble {{ background: var(--nb-panel); border: 1px solid var(--nb-rule); border-radius: 12px;
  padding: 7px 11px; line-height: 1.45; overflow-wrap: anywhere; }}
.nb-app .bb-turn-user .bb-bubble {{ background: var(--nb-au-soft); border-color: var(--nb-au); }}
.nb-app .bb-bubble ul {{ margin: 2px 0 2px 16px; padding: 0; }}
.nb-app .bb-bubble .nb-msg {{ margin: 6px 0 2px; }}
.nb-app .bb-sub {{ font-family: {DISPLAY_STACK}; font-size: 12px; font-weight: 600; color: var(--nb-muted);
  margin: 6px 0 2px; }}
.nb-app pre.bb-cmds {{ margin: 0; padding: 6px 8px; background: var(--nb-input); border: 1px solid var(--nb-rule);
  border-radius: 8px; font: 12px/1.45 ui-monospace, Menlo, Consolas, monospace; color: var(--nb-au); }}
.nb-app .bb-checks {{ list-style: none; margin: 0 !important; }}
.nb-app .bb-checks .bb-ok {{ color: var(--nb-ok); }}
.nb-app .bb-checks .bb-bad {{ color: var(--nb-err); }}
.nb-app .bb-answer {{ white-space: pre-wrap; }}
.nb-app .bb-turn-thinking .bb-bubble {{ color: var(--nb-muted); font-style: italic; border-style: dashed; }}
.nb-app .bb-dots::after {{ content: ""; animation: bbdots 1.2s steps(4, end) infinite; }}
@keyframes bbdots {{ 0% {{ content: ""; }} 25% {{ content: "."; }} 50% {{ content: ".."; }} 75% {{ content: "..."; }} }}
.nb-app .bb-badge {{ display: inline-block; font-size: 11px; padding: 2px 9px; border-radius: 999px;
  border: 1px solid var(--nb-rule); color: var(--nb-muted); background: var(--nb-bg); }}
.nb-app .bb-badge-ai {{ border-color: var(--nb-ok); color: var(--nb-ok); }}
.nb-app .bb-chat-input textarea {{ font-size: 13px; min-height: 40px; resize: vertical; }}

@media (prefers-reduced-motion: reduce) {{ .nb-app * {{ transition: none !important; }} }}
</style>"""
