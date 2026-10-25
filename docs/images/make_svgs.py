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


def box(x, y, w, h, title, sub, t):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{t["box"]}" '
            f'stroke="{t["edge"]}"/>'
            f'<text x="{x + w / 2}" y="{y + 22}" font-size="13" font-weight="700" text-anchor="middle" '
            f'fill="{t["text"]}">{title}</text>'
            f'<text x="{x + w / 2}" y="{y + 40}" font-size="11" text-anchor="middle" '
            f'fill="{t["muted"]}">{sub}</text>')


def arrow(x0, y0, x1, y1, t):
    return (f'<path d="M{x0},{y0} L{x1},{y1}" stroke="{t["arrow"]}" stroke-width="1.8" '
            f'marker-end="url(#ah)"/>')


def pipeline():
    t = AUTO
    parts = [
        f'<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{t["arrow"]}"/></marker></defs>'
    ]
    # Long-read lane (y=20): trim/filter -> assemble -> medaka -> short polish -> QC
    parts.append(f'<text x="80" y="16" font-size="12" font-weight="700" fill="{t["muted"]}">Long reads (Nanopore)</text>')
    parts.append(box(20, 24, 120, 50, '1 Prepare', 'Porechop · Filtlong', t))
    parts.append(box(170, 24, 120, 50, '2 Assemble', 'Flye or Shasta', t))
    parts.append(box(320, 24, 120, 50, '3 Polish long', 'Medaka (or none)', t))
    parts.append(box(470, 24, 150, 50, '4 Polish short', 'Polypolish → pyPOLCA', t))
    parts.append(arrow(140, 49, 168, 49, t))
    parts.append(arrow(290, 49, 318, 49, t))
    parts.append(arrow(440, 49, 468, 49, t))
    parts.append(arrow(620, 49, 658, 84, t))

    # Short-read lane (y=120)
    parts.append(f'<text x="80" y="112" font-size="12" font-weight="700" fill="{t["muted"]}">Short reads (Illumina, Ion Torrent)</text>')
    parts.append(box(20, 120, 120, 50, 'Trim', 'fastp', t))
    parts.append(box(170, 120, 120, 50, 'Assemble', 'SKESA or SPAdes', t))
    parts.append(arrow(140, 145, 168, 145, t))
    # Short reads joining the hybrid lane
    parts.append(arrow(140, 132, 468, 62, t))
    # Short-only to QC
    parts.append(arrow(290, 145, 658, 116, t))

    # QC box spanning both lanes
    parts.append(box(662, 70, 150, 60, 'QC', 'QUAST · flagstat ·', t))
    parts.append(f'<text x="737" y="124" font-size="11" text-anchor="middle" fill="{t["muted"]}">'
                 f'mosdepth · dot plots · graphs</text>')

    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 830 190" width="830" height="190"
  font-family="{FONT}">
  <title>HashBrown pipeline: long reads are prepared, assembled and polished; short reads polish them or are assembled on their own; everything ends in QC</title>{style()}
  {"".join(parts)}
</svg>
'''


if __name__ == '__main__':
    with open('hashbrown_logo.svg', 'w') as f:
        f.write(logo())
    with open('pipeline.svg', 'w') as f:
        f.write(pipeline())
    print('Wrote hashbrown_logo.svg and pipeline.svg')
