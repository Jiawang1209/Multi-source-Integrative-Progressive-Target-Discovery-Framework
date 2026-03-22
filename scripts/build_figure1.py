#!/usr/bin/env python3

import argparse
import os
import subprocess
import textwrap
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
TMP_DIR = ROOT / "tmp" / "figure1_builder"
OUTPUT_DIR = ROOT / "output" / "figure1"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assemble a Figure 1 style summary from local offline results."
    )
    parser.add_argument("--case-dir", required=True, help="Case directory with PDF outputs.")
    parser.add_argument("--compound-name", required=True, help="Name shown in panel a.")
    parser.add_argument("--compound-subtitle", default="", help="Optional subtitle under the name.")
    parser.add_argument(
        "--methods",
        default="SEA, PPB2, SwissTargetPrediction, ChEMBL",
        help="Methods shown in the target-prediction funnel.",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Output PNG path. Defaults to output/figure1/<case-name>_figure1.png",
    )
    return parser.parse_args()


def load_font(size, bold=False):
    candidates = []
    if bold:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica.ttc",
                "/System/Library/Fonts/Supplemental/PingFang.ttc",
            ]
        )
    else:
        candidates.extend(
            [
                "/System/Library/Fonts/Supplemental/Arial.ttf",
                "/System/Library/Fonts/Supplemental/Helvetica.ttc",
                "/System/Library/Fonts/Supplemental/PingFang.ttc",
            ]
        )
    for path in candidates:
        if os.path.exists(path):
            try:
                return ImageFont.truetype(path, size=size)
            except OSError:
                continue
    return ImageFont.load_default()


def trim_white(image):
    bg = Image.new(image.mode, image.size, "white")
    diff = ImageChops.difference(image, bg)
    bbox = diff.getbbox()
    if not bbox:
        return image
    return image.crop(bbox)


def render_pdf_first_page(pdf_path, stem):
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    prefix = TMP_DIR / stem
    env = os.environ.copy()
    env["XDG_CACHE_HOME"] = str(TMP_DIR / "cache")
    env["FONTCONFIG_PATH"] = env.get("FONTCONFIG_PATH", "/opt/homebrew/etc/fonts")
    env["FONTCONFIG_FILE"] = env.get("FONTCONFIG_FILE", "/opt/homebrew/etc/fonts/fonts.conf")
    (TMP_DIR / "cache").mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["pdftoppm", "-f", "1", "-l", "1", "-png", str(pdf_path), str(prefix)],
        check=True,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    image_path = prefix.parent / f"{prefix.name}-1.png"
    return trim_white(Image.open(image_path).convert("RGB"))


def fetch_top_targets(case_dir, limit=5):
    venn_path = Path(case_dir) / "venn.rds"
    if not venn_path.exists():
        return []

    r_code = """
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(args[1])
all_targets <- sort(unique(unlist(x)))
votes <- sapply(all_targets, function(g) sum(sapply(x, function(v) g %in% v)))
df <- data.frame(target = all_targets, votes = votes)
df <- df[order(-df$votes, df$target), ]
df <- df[df$votes >= 2, ]
if (nrow(df) == 0) {
  df <- head(data.frame(target = all_targets, votes = votes), 5)
} else {
  df <- head(df, 5)
}
apply(df, 1, function(row) {
  cat(row[["target"]], row[["votes"]], sep = "\\t")
  cat("\\n")
})
"""
    result = subprocess.run(
        ["Rscript", "-e", r_code, str(venn_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    targets = []
    for line in result.stdout.splitlines():
        if not line.strip() or "\t" not in line:
            continue
        name, votes = line.split("\t")
        targets.append((name.strip(), int(votes)))
    return targets[:limit]


def contain_image(image, size):
    framed = Image.new("RGB", size, "white")
    fitted = image.copy()
    fitted.thumbnail(size, Image.Resampling.LANCZOS)
    x = (size[0] - fitted.size[0]) // 2
    y = (size[1] - fitted.size[1]) // 2
    framed.paste(fitted, (x, y))
    return framed


def draw_panel_label(draw, label, xy):
    font = load_font(42, bold=True)
    draw.text(xy, label, fill="#111111", font=font)


def draw_centered_text(draw, box, text, font, fill="#111111"):
    left, top, right, bottom = box
    bbox = draw.multiline_textbbox((0, 0), text, font=font, align="center", spacing=6)
    w = bbox[2] - bbox[0]
    h = bbox[3] - bbox[1]
    x = left + (right - left - w) / 2
    y = top + (bottom - top - h) / 2
    draw.multiline_text((x, y), text, fill=fill, font=font, align="center", spacing=6)


def draw_arrow(draw, start, end, color="#b9b9b9", width=6):
    draw.line([start, end], fill=color, width=width)
    ex, ey = end
    sx, sy = start
    dx = ex - sx
    dy = ey - sy
    length = max((dx * dx + dy * dy) ** 0.5, 1)
    ux = dx / length
    uy = dy / length
    px = -uy
    py = ux
    arrow = 18
    p1 = (ex, ey)
    p2 = (ex - ux * arrow + px * 10, ey - uy * arrow + py * 10)
    p3 = (ex - ux * arrow - px * 10, ey - uy * arrow - py * 10)
    draw.polygon([p1, p2, p3], fill=color)


def make_panel_a(width, height, compound_name, compound_subtitle, methods):
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_panel_label(draw, "a", (10, 6))

    title_font = load_font(34, bold=True)
    subtitle_font = load_font(22)
    title_text = textwrap.fill(compound_name, width=24)
    draw_centered_text(draw, (40, 34, width - 40, 120), title_text, title_font)
    if compound_subtitle:
        draw_centered_text(draw, (50, 104, width - 50, 150), compound_subtitle, subtitle_font, fill="#555555")

    draw_arrow(draw, (width // 2 - 70, 165), (width // 2 + 10, 195), color="#d6d6d6", width=5)

    funnel = [
        (
            (40, 210),
            (width - 40, 210),
            (width - 95, 305),
            (95, 305),
            "#f6b3ab",
            "Target Prediction\n(" + textwrap.fill(methods, width=24) + ")",
        ),
        ((82, 325), (width - 82, 325), (width - 120, 396), (120, 396), "#f5a05f", "KEGG Analysis"),
        ((116, 414), (width - 116, 414), (width - 145, 478), (145, 478), "#f6da72", "Multi-task QSAR"),
        ((150, 496), (width - 150, 496), (width - 170, 550), (170, 550), "#9fcd7e", "High Confidence\nTargets"),
    ]
    box_font = load_font(22, bold=True)
    for left_top, right_top, right_bottom, left_bottom, color, text in funnel:
        draw.polygon([left_top, right_top, right_bottom, left_bottom], fill=color)
        draw_centered_text(
            draw,
            (left_top[0], left_top[1], right_top[0], left_bottom[1]),
            text,
            box_font,
            fill="#333333",
        )
    return image


def draw_molecule_stub(draw, box, color="#7ec7de"):
    left, top, right, bottom = box
    points = [
        (left + 22, top + 50),
        (left + 70, top + 20),
        (left + 124, top + 42),
        (left + 136, top + 98),
        (left + 86, top + 124),
        (left + 30, top + 104),
    ]
    draw.rounded_rectangle(box, radius=18, outline="#c9dcec", width=3)
    for idx in range(len(points)):
        draw.line([points[idx], points[(idx + 1) % len(points)]], fill=color, width=5)
    extra = [(left + 136, top + 98), (right - 18, top + 128), (left + 70, top + 20), (left + 40, top - 2)]
    draw.line(extra[:2], fill=color, width=5)
    draw.line(extra[2:], fill=color, width=5)
    for x, y in points + [extra[1], extra[3]]:
        draw.ellipse((x - 7, y - 7, x + 7, y + 7), fill="#ffb39f", outline="#5b91ab")


def draw_network(draw, center_x, center_y):
    layers = [
        [(center_x - 90, center_y - 80), (center_x - 90, center_y - 25), (center_x - 90, center_y + 30), (center_x - 90, center_y + 85)],
        [(center_x, center_y - 105), (center_x, center_y - 50), (center_x, center_y + 5), (center_x, center_y + 60), (center_x, center_y + 115)],
        [(center_x + 95, center_y - 70), (center_x + 95, center_y - 15), (center_x + 95, center_y + 40), (center_x + 95, center_y + 95)],
    ]
    layer_colors = ["#7bc8e2", "#f9b4aa", "#f5b56b"]
    for source, target in zip(layers, layers[1:]):
        for x1, y1 in source:
            for x2, y2 in target:
                draw.line((x1, y1, x2, y2), fill="#d9d9d9", width=2)
    for color, layer in zip(layer_colors, layers):
        for x, y in layer:
            draw.ellipse((x - 14, y - 14, x + 14, y + 14), fill=color, outline="#7d7d7d", width=2)


def make_panel_e(width, height, top_targets):
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_panel_label(draw, "e", (10, 6))

    draw_molecule_stub(draw, (40, 94, 210, 258))
    draw_molecule_stub(draw, (40, 300, 210, 464))

    box_font = load_font(20)
    bold_font = load_font(24, bold=True)
    draw.text((246, 196), "Message Passing", fill="#555555", font=box_font)
    draw.text((242, 304), "Hidden Layer", fill="#555555", font=box_font)
    draw.rounded_rectangle((228, 156, 398, 404), radius=18, outline="#96bce0", width=3)
    for i in range(5):
        x0 = 260 + i * 18
        draw.rounded_rectangle((x0, 234, x0 + 12, 322), radius=3, fill="#7bc8e2")
    draw.rounded_rectangle((325, 220, 370, 248), radius=6, fill="#f8c07a")
    draw.rounded_rectangle((325, 254, 370, 282), radius=6, fill="#f9b4aa")
    draw.rounded_rectangle((325, 288, 370, 316), radius=6, fill="#9dd38c")
    draw.text((430, 236), "Readout", fill="#555555", font=box_font)
    draw.line((400, 278, 468, 278), fill="#b7b7b7", width=5)
    draw.polygon([(468, 278), (450, 268), (450, 288)], fill="#b7b7b7")
    draw.rounded_rectangle((476, 246, 518, 310), radius=8, outline="#9bbbe6", width=3)
    draw.text((487, 264), "S", fill="#3c78b9", font=bold_font)

    draw_network(draw, 660, 276)
    draw.text((610, 70), "Multi-task QSAR", fill="#111111", font=load_font(28, bold=True))

    draw.line((820, 192, 880, 192), fill="#d0d0d0", width=4)
    draw.line((820, 360, 880, 360), fill="#d0d0d0", width=4)
    draw.text((886, 172), "Task 1", fill="#555555", font=box_font)
    draw.text((886, 340), "Task n", fill="#555555", font=box_font)

    list_box = (1030, 125, width - 30, 430)
    draw.rounded_rectangle(list_box, radius=16, outline="#7da4d7", width=3)
    list_font = load_font(24)
    header_font = load_font(26, bold=True)
    draw.text((1060, 92), "High-confidence targets", fill="#3465a4", font=header_font)

    y = 160
    for name, votes in top_targets:
        text = f"{name}  (vote={votes})"
        draw.text((1060, y), text, fill="#333333", font=list_font)
        y += 48

    if not top_targets:
        draw.text((1060, 190), "No target summary available", fill="#777777", font=list_font)
    return image


def assemble(case_dir, compound_name, compound_subtitle, methods, output_path):
    case_path = Path(case_dir).resolve()
    venn_pdf = case_path / "1.p_venn.pdf"
    gokegg_pdf = case_path / "2.GO_KEGG.pdf"
    circos_pdf = case_path / "3.KEGG_circos.pdf"

    for pdf in [venn_pdf, gokegg_pdf, circos_pdf]:
        if not pdf.exists():
            raise FileNotFoundError(f"Missing required panel PDF: {pdf}")

    panel_a = make_panel_a(540, 640, compound_name, compound_subtitle, methods)
    panel_b = render_pdf_first_page(venn_pdf, "venn_panel")
    panel_c = render_pdf_first_page(gokegg_pdf, "gokegg_panel")
    panel_d = render_pdf_first_page(circos_pdf, "circos_panel")
    top_targets = fetch_top_targets(case_path)
    panel_e = make_panel_e(1360, 520, top_targets)

    canvas = Image.new("RGB", (1400, 2220), "white")
    draw = ImageDraw.Draw(canvas)

    def paste_panel(image, box):
        left, top, width, height = box
        panel = contain_image(image, (width, height))
        canvas.paste(panel, (left, top))

    paste_panel(panel_a, (20, 20, 520, 620))
    paste_panel(panel_c, (560, 20, 820, 620))
    draw_panel_label(draw, "c", (560, 20))

    draw_panel_label(draw, "b", (20, 650))
    paste_panel(panel_b, (20, 690, 520, 600))

    draw_panel_label(draw, "d", (560, 650))
    paste_panel(panel_d, (560, 690, 820, 760))

    paste_panel(panel_e, (20, 1490, 1360, 520))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_path, format="PNG")


def main():
    args = parse_args()
    case_path = Path(args.case_dir)
    if not case_path.is_absolute():
        case_path = (ROOT / case_path).resolve()

    output_path = Path(args.output) if args.output else OUTPUT_DIR / f"{case_path.name}_figure1.png"
    if not output_path.is_absolute():
        output_path = (ROOT / output_path).resolve()

    assemble(
        case_path,
        args.compound_name,
        args.compound_subtitle,
        args.methods,
        output_path,
    )
    print(output_path)


if __name__ == "__main__":
    main()
