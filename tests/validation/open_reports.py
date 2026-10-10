"""Open Excel and render HTML in real Chrome; save only PDF visual artifacts."""
import argparse
import base64
import json
import os
from pathlib import Path
import subprocess
import time
import urllib.request
import warnings

import openpyxl
import websocket


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.out.resolve()
    output.mkdir(parents=True, exist_ok=False)
    workbook_checks = []
    for path in root.rglob("*.xlsx"):
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            workbook = openpyxl.load_workbook(path, data_only=True)
        sheets = []
        for sheet in workbook:
            assert sheet.freeze_panes == "A2", (path, sheet.title, sheet.freeze_panes)
            rows = list(sheet.iter_rows(values_only=True))
            assert rows, (path, sheet.title)
            sheets.append({"name": sheet.title, "rows": len(rows), "columns": max(map(len, rows))})
        workbook_checks.append({"path": str(path), "sheets": sheets})
        workbook.close()
    assert workbook_checks
    profile = output / "chrome-profile"
    fonts = output / "fonts.conf"
    fonts.write_text(f'<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd"><fontconfig><dir>/data/project/yanhuixu/yfu/01_source/fonts</dir><dir>/usr/share/fonts</dir><cachedir>{output / "fonts-cache"}</cachedir></fontconfig>')
    environment = dict(os.environ, FONTCONFIG_FILE=str(fonts))
    # websocket-client only reads its no_proxy option with an explicit proxy
    # host; environment bypass also covers the default-proxy path.
    os.environ["no_proxy"] = os.environ["NO_PROXY"] = "127.0.0.1,localhost"
    request_id, errors = 0, []
    results = []
    with (output / "chrome.log").open("w") as log:
        browser = subprocess.Popen(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--remote-debugging-port=0", "--no-first-run", "--no-default-browser-check", "--user-data-dir=" + str(profile), "about:blank"], stdout=log, stderr=log, env=environment)
        connection = None
        try:
            active_port = profile / "DevToolsActivePort"
            deadline = time.monotonic() + 20
            while not active_port.exists():
                if browser.poll() is not None:
                    raise RuntimeError("Chrome exited before opening its debug endpoint")
                if time.monotonic() >= deadline:
                    raise TimeoutError("Chrome debug endpoint did not start")
                time.sleep(0.1)
            port = active_port.read_text().splitlines()[0]
            local_http = urllib.request.build_opener(urllib.request.ProxyHandler({}))
            with local_http.open(f"http://127.0.0.1:{port}/json/list") as response:
                pages = json.load(response)
            page = next(page for page in pages if page["type"] == "page")
            connection = websocket.create_connection(page["webSocketDebuggerUrl"], suppress_origin=True, timeout=10, http_no_proxy=["127.0.0.1", "localhost"])
            def call(method, params=None):
                nonlocal request_id
                request_id += 1
                connection.send(json.dumps({"id": request_id, "method": method, "params": params or {}}))
                while True:
                    message = json.loads(connection.recv())
                    if message.get("method") == "Runtime.exceptionThrown":
                        errors.append(message["params"]["exceptionDetails"])
                    if message.get("id") == request_id:
                        if "error" in message:
                            raise RuntimeError(message["error"])
                        return message.get("result", {})
            call("Runtime.enable")
            call("Page.enable")
            paths = list(root.glob("*/02_output/*/06_report/*.report.html"))
            paths += list(root.glob("*/02_output/00_report/*.html"))
            paths += list(root.glob("*/02_output/*/05_downsample/*.downsample.html"))
            paths += list(root.glob("*/02_output/*/06_report/*.quick_test.html"))
            assert paths
            expression = "JSON.stringify({ready:document.readyState,plotly:typeof Plotly!=='undefined',charts:Array.from(document.querySelectorAll('.js-plotly-plot')).map(e=>({id:e.id,traces:(e.data||[]).length,visible:e.offsetWidth>0,sd:(e.data||[]).map(t=>t.error_y&&t.error_y.array)})),title:document.title,text:document.body.innerText})"
            for index, path in enumerate(paths):
                errors.clear()
                needs_plotly = "Plotly.newPlot" in path.read_text()
                call("Page.navigate", {"url": path.as_uri()})
                deadline = time.monotonic() + 40
                while True:
                    result = call("Runtime.evaluate", {"expression": expression, "returnByValue": True})
                    raw = result.get("result", {}).get("value")
                    state = json.loads(raw) if raw else {}
                    if state.get("ready") == "complete" and (not needs_plotly or state.get("plotly") and state.get("charts")):
                        break
                    if time.monotonic() >= deadline:
                        raise TimeoutError(f"Report did not render: {path}; state={state}; exceptions={errors}")
                    time.sleep(0.25)
                assert not errors, (path, errors)
                call("Runtime.evaluate", {"expression": "document.fonts.ready", "awaitPromise": True})
                pdf = call("Page.printToPDF", {"printBackground": True, "paperWidth": 12, "paperHeight": 16, "marginTop": 0.3, "marginBottom": 0.3})
                pdf_path = output / f"{index:02d}-{path.stem}.pdf"
                pdf_path.write_bytes(base64.b64decode(pdf["data"]))
                assert pdf_path.stat().st_size > 1000
                state.pop("text", None)
                results.append({"html": str(path), "pdf": str(pdf_path), "runtime_exceptions": [], "state": state})
                (output / "partial.json").write_text(json.dumps(results, indent=2) + "\n")
            summary = {"result": "OK", "workbooks": workbook_checks, "html": results, "font": "Arial from user-provided font directory", "visual_format": "PDF"}
            (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            print(json.dumps({"result": "OK", "workbooks": len(workbook_checks), "html": len(results)}), flush=True)
        finally:
            if connection is not None:
                connection.close()
            browser.terminate()
            try:
                browser.wait(timeout=10)
            except subprocess.TimeoutExpired:
                browser.kill()
                browser.wait()


if __name__ == "__main__":
    main()
