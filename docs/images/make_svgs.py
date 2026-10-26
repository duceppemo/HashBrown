"""The HashBrown logo (reads pressed into one assembly, like a hash brown) and the pipeline diagram.

Run from docs/images/: python make_svgs.py
"""

LIGHT = {"text": "#1d2433", "muted": "#5d6678", "box": "#f6f7f9", "edge": "#d7dbe2", "arrow": "#8a93a3",
         "reads": "#8a93a3", "ref": "#5d6678"}
DARK = {"text": "#e6e8ec", "muted": "#9aa3b2", "box": "#1d2129", "edge": "#363c47", "arrow": "#737d8e",
        "reads": "#737d8e", "ref": "#9aa3b2"}
# One adaptive file: colors are CSS variables, flipped by a media query inside the SVG. Works wherever
# the SVG is shown in an <img>, GitHub wikis included.
AUTO = {k: f"var(--{k})" for k in LIGHT}

GOLDEN = "#c87f3a"  # the crispy accent: the assembled genome
FONT = "system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"


def style():
    light = "; ".join(f"--{k}: {v}" for k, v in LIGHT.items())
    dark = "; ".join(f"--{k}: {v}" for k, v in DARK.items())
    return (f"\n  <style>:root {{ {light} }} "
            f"@media (prefers-color-scheme: dark) {{ :root {{ {dark} }} }}</style>")


# The icon: shredded reads (long strips and short dashes) pressed into one golden patty.
# Each row is (y, [(x_start, x_end, kind)]) with kind 'long' (nanopore) or 'short' (illumina).
ROWS = [
    (34, [(38, 92, 'long'), (98, 116, 'short')]),
    (52, [(28, 44, 'short'), (50, 118, 'long')]),
    (70, [(24, 86, 'long'), (92, 108, 'short'), (114, 130, 'short')]),
    (88, [(26, 42, 'short'), (48, 64, 'short'), (70, 126, 'long')]),
    (106, [(32, 98, 'long'), (104, 122, 'short')]),
    (124, [(42, 58, 'short'), (64, 114, 'long')]),
]


def logo():
    t = AUTO
    parts = [
        # The patty: one piece out of many shreds
        f'<rect x="12" y="20" width="130" height="118" rx="26" fill="none" '
        f'stroke="{GOLDEN}" stroke-width="5"/>'
    ]
    for y, strips in ROWS:
        for x0, x1, kind in strips:
            color = t["reads"] if kind == 'long' else t["ref"]
            width = 7 if kind == 'long' else 6
            parts.append(f'<path d="M{x0},{y} H{x1}" fill="none" stroke="{color}" '
                         f'stroke-width="{width}" stroke-linecap="round"/>')
    body = "\n  ".join(parts)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 158" width="640" height="158"
  font-family="{FONT}">
  <title>HashBrown: long-read, short-read and hybrid bacterial genome assembly</title>{style()}
  {body}
  <text x="172" y="88" font-size="54" font-weight="800" letter-spacing="2" fill="var(--text)">Hash<tspan fill="{GOLDEN}">Brown</tspan></text>
  <text x="174" y="120" font-size="14.5" letter-spacing="0.6" fill="var(--muted)">Long-read, short-read and hybrid bacterial genome assembly</text>
</svg>
'''


def report_icon():
    """Icon-only logo variant shipped inside the package (hashbrown/logo.svg) for the reports.

    Inlined into report HTML, so: no <style>, no :root variables (they would fight the page's own
    CSS variables) and fixed theme-neutral mid-tone colors, legible on the light and dark page alike.
    """
    reads, ref = "#8a93a3", "#737d8e"
    parts = [f'<rect x="12" y="20" width="130" height="118" rx="26" fill="none" '
             f'stroke="{GOLDEN}" stroke-width="7"/>']
    for y, strips in ROWS:
        for x0, x1, kind in strips:
            color = reads if kind == 'long' else ref
            width = 7 if kind == 'long' else 6
            parts.append(f'<path d="M{x0},{y} H{x1}" fill="none" stroke="{color}" '
                         f'stroke-width="{width}" stroke-linecap="round"/>')
    body = "\n  ".join(parts)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 154 158" width="27" height="28"
  role="img" aria-label="HashBrown">
  {body}
</svg>
'''


def box(x, y, w, h, title, sub, t, dashed=False):
    dash = ' stroke-dasharray="5 4"' if dashed else ''
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="9" fill="{t["box"]}" '
            f'stroke="{t["edge"]}" stroke-width="1.2"{dash}/>'
            f'<text x="{x + w / 2}" y="{y + 22}" font-size="13" font-weight="700" text-anchor="middle" '
            f'fill="{t["text"]}">{title}</text>'
            f'<text x="{x + w / 2}" y="{y + 40}" font-size="11" text-anchor="middle" '
            f'fill="{t["muted"]}">{sub}</text>')


def arrow(d, t, marker='ah', dashed=False):
    dash = ' stroke-dasharray="5 4"' if dashed else ''
    stroke = GOLDEN if marker == 'ahg' else t["arrow"]
    return (f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="1.8"{dash} '
            f'marker-end="url(#{marker})"/>')


def pipeline():
    # Two columns flowing top to bottom, a routing gutter between them, and two full-width boxes
    # (QC, then the report) across the bottom. Long reads go down the left column; short reads down
    # the right. In hybrid runs the trimmed short reads cross the gutter (golden dashed connector)
    # to polish the long-read assembly, and both columns end with an optional reorientation.
    t = AUTO
    bw, bh = 250, 52           # box size
    lx, rx = 40, 400           # column left edges
    lc, rc = lx + bw / 2, rx + bw / 2  # column centers (165, 525)
    gx = (lx + bw + rx) / 2    # gutter channel x (345)
    rows = [34, 120, 206, 292, 378]  # box tops; centers are top + 26
    qy, qh = 464, 56           # QC box
    ry, rh = 554, 56           # report box

    parts = [
        f'<defs>'
        f'<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{t["arrow"]}"/></marker>'
        f'<marker id="ahg" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{GOLDEN}"/></marker>'
        f'</defs>'
    ]

    # Column headers
    parts.append(f'<text x="{lc}" y="20" font-size="12" font-weight="700" text-anchor="middle" '
                 f'fill="{t["muted"]}">Long reads (Nanopore)</text>')
    parts.append(f'<text x="{rc}" y="20" font-size="12" font-weight="700" text-anchor="middle" '
                 f'fill="{t["muted"]}">Short reads (Illumina · Ion Torrent)</text>')

    # Left column: the long-read path
    parts.append(box(lx, rows[0], bw, bh, 'Prepare', 'Porechop · Filtlong', t))
    parts.append(box(lx, rows[1], bw, bh, 'Assemble', 'Flye or Shasta', t))
    parts.append(box(lx, rows[2], bw, bh, 'Polish (long reads)', 'Medaka · optional', t))
    parts.append(box(lx, rows[3], bw, bh, 'Polish (short reads)', 'Polypolish → pyPOLCA · hybrid only', t,
                     dashed=True))
    parts.append(box(lx, rows[4], bw, bh, 'Reorient', 'dnaapler · optional', t, dashed=True))
    for y0, y1 in zip(rows, rows[1:], strict=False):
        parts.append(arrow(f'M{lc},{y0 + bh + 2} L{lc},{y1 - 2}', t))
    parts.append(arrow(f'M{lc},{rows[4] + bh + 2} L{lc},{qy - 2}', t))

    # Right column: the short-read path, with Reorient aligned to the left column's
    parts.append(box(rx, rows[0], bw, bh, 'Trim', 'fastp', t))
    parts.append(box(rx, rows[1], bw, bh, 'Assemble', 'SKESA or SPAdes', t))
    parts.append(box(rx, rows[4], bw, bh, 'Reorient', 'dnaapler · optional', t, dashed=True))
    parts.append(arrow(f'M{rc},{rows[0] + bh + 2} L{rc},{rows[1] - 2}', t))
    parts.append(arrow(f'M{rc},{rows[1] + bh + 2} L{rc},{rows[4] - 2}', t))
    parts.append(arrow(f'M{rc},{rows[4] + bh + 2} L{rc},{qy - 2}', t))
    parts.append(f'<text x="{rc + 8}" y="284" font-size="11" font-style="italic" '
                 f'fill="{t["muted"]}">short reads only</text>')

    # Hybrid: trimmed short reads cross the gutter to polish the long-read assembly
    y_trim, y_pol = rows[0] + 26, rows[3] + 26
    parts.append(arrow(f'M{rx - 2},{y_trim} H{gx + 6} Q{gx},{y_trim} {gx},{y_trim + 6} '
                       f'V{y_pol - 6} Q{gx},{y_pol} {gx - 6},{y_pol} H{lx + bw + 2}',
                       t, marker='ahg', dashed=True))
    parts.append(f'<text x="{gx + 8}" y="190" font-size="11" font-style="italic" '
                 f'fill="{GOLDEN}">hybrid</text>')

    # QC and report boxes spanning both columns: every run ends here
    parts.append(box(lx, qy, rx + bw - lx, qh, 'QC — every run',
                     'QUAST · flagstat · mosdepth · dot plots · assembly graphs', t))
    parts.append(arrow(f'M{gx},{qy + qh + 2} L{gx},{ry - 2}', t))
    parts.append(box(lx, ry, rx + bw - lx, rh, 'Report — &lt;sample&gt;.html + .pdf',
                     'read &amp; assembly funnels · QV · species ID · completeness · typing · AMR · '
                     'annotation · contamination', t))

    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 690 626" width="690" height="626"
  font-family="{FONT}">
  <title>HashBrown pipeline: long reads are prepared, assembled and polished; short reads polish them (hybrid) or are assembled on their own; assemblies are optionally reoriented, and everything ends in QC and an HTML + PDF report per sample</title>{style()}
  {"".join(parts)}
</svg>
'''


if __name__ == '__main__':
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'hashbrown_logo.svg'), 'w') as f:
        f.write(logo())
    with open(os.path.join(here, 'pipeline.svg'), 'w') as f:
        f.write(pipeline())
    package_icon = os.path.join(here, '..', '..', 'hashbrown', 'logo.svg')
    with open(package_icon, 'w') as f:
        f.write(report_icon())
    print('Wrote hashbrown_logo.svg, pipeline.svg and hashbrown/logo.svg')
