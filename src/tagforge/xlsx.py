"""Tiny dependency-free XLSX writer for report tables."""
from __future__ import annotations

import html
import re
import zipfile
from pathlib import Path


def _column(n):
    value = ""
    while n:
        n, rem = divmod(n - 1, 26); value = chr(65 + rem) + value
    return value


def _cell(value, ref, header=False):
    style = ' s="1"' if header else ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f'<c r="{ref}"{style}><v>{value}</v></c>'
    text = html.escape("" if value is None else str(value))
    return f'<c r="{ref}" t="inlineStr"{style}><is><t>{text}</t></is></c>'


def write_xlsx(path: Path, sheets):
    """Write ``[(sheet_name, rows), ...]`` where rows are sequences."""
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_sheets = []
    used = set()
    for name, rows in sheets:
        name = re.sub(r"[\\/*?:\[\]]", "_", str(name))[:31] or "Sheet"
        original = name; i = 2
        while name in used:
            suffix = f"_{i}"; name = original[:31-len(suffix)] + suffix; i += 1
        used.add(name); safe_sheets.append((name, list(rows)))
    tmp = Path(str(path) + ".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", '<?xml version="1.0" encoding="UTF-8"?>'
          '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
          '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
          '<Default Extension="xml" ContentType="application/xml"/>'
          '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
          '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' +
          ''.join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1,len(safe_sheets)+1)) + '</Types>')
        z.writestr("_rels/.rels", '<?xml version="1.0" encoding="UTF-8"?>'
          '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
          '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        sheet_defs = ''.join(f'<sheet name="{html.escape(name)}" sheetId="{i}" r:id="rId{i}"/>' for i,(name,_) in enumerate(safe_sheets,1))
        z.writestr("xl/workbook.xml", '<?xml version="1.0" encoding="UTF-8"?>'
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>' + sheet_defs + '</sheets></workbook>')
        rels = ''.join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1,len(safe_sheets)+1))
        rels += f'<Relationship Id="rId{len(safe_sheets)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
        z.writestr("xl/_rels/workbook.xml.rels", '<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' + rels + '</Relationships>')
        z.writestr("xl/styles.xml", '<?xml version="1.0" encoding="UTF-8"?><styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF334155"/><bgColor indexed="64"/></patternFill></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs></styleSheet>')
        for index, (_, rows) in enumerate(safe_sheets, 1):
            data = []
            max_cols = max((len(row) for row in rows), default=1)
            for r, row in enumerate(rows, 1):
                cells = ''.join(_cell(value, f"{_column(c)}{r}", r == 1) for c, value in enumerate(row, 1))
                data.append(f'<row r="{r}">{cells}</row>')
            widths = ''.join(f'<col min="{c}" max="{c}" width="18" customWidth="1"/>' for c in range(1,max_cols+1))
            xml = '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><cols>' + widths + '</cols><sheetData>' + ''.join(data) + '</sheetData><autoFilter ref="A1:' + _column(max_cols) + '1"/><freezePane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></worksheet>'
            z.writestr(f"xl/worksheets/sheet{index}.xml", xml)
    tmp.replace(path)

