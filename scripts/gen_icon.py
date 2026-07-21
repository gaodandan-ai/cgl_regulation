"""
gen_icon.py — Generate a sky-blue radial-gradient icon with white "Cgl" text.

Gradient: center is the deepest sky-blue (#0284c7), fading outward to a very
light azure (#bae6fd), giving a "glowing core" look.
Shape   : rounded-square (squircle-ish) for .ico multi-resolution, and a
          full-circle PNG for web/manifest use.
Output  : icon.png  (512×512, circle on transparent)
          icon.ico  (256/128/64/48/32/16 px, squircle mask)
"""

import math
import os
import struct
import zlib
from PIL import Image, ImageDraw, ImageFont, ImageFilter

# -- Colour palette (lighter sky-blue family) --------------------------------
CENTER_COLOR  = (56,  189, 248)   # #38bdf8 -- sky-400 (lighter core)
EDGE_COLOR    = (240, 249, 255)   # #f0f9ff -- sky-50  (near-white edge)
MID_COLOR     = (125, 211, 252)   # #7dd3fc -- sky-300 (mid)
TEXT_COLOR    = (255, 255, 255)   # white

# ── Canvas ────────────────────────────────────────────────────────────────────
SIZE = 512          # working canvas
CORNER_R = 112      # rounded-corner radius for .ico squircle (≈22% of size)


def lerp_color(c1, c2, t):
    """Linear interpolate between two RGB tuples; t in [0, 1]."""
    return tuple(int(c1[i] + (c2[i] - c1[i]) * t) for i in range(3))


def radial_gradient(size, center_col, mid_col, edge_col):
    """
    Build an RGBA image with a 3-stop radial gradient:
      0.0 → center_col
      0.5 → mid_col
      1.0 → edge_col
    """
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    pixels = img.load()
    cx = cy = size / 2
    max_r = size / 2  # radius to corner

    for y in range(size):
        for x in range(size):
            dx, dy = x - cx, y - cy
            dist = math.sqrt(dx * dx + dy * dy)
            t = min(dist / max_r, 1.0)

            if t < 0.5:
                col = lerp_color(center_col, mid_col, t / 0.5)
            else:
                col = lerp_color(mid_col, edge_col, (t - 0.5) / 0.5)

            pixels[x, y] = (*col, 255)

    return img


def rounded_rect_mask(size, radius):
    """Alpha mask: white rounded rectangle on black."""
    mask = Image.new("L", (size, size), 0)
    draw = ImageDraw.Draw(mask)
    draw.rounded_rectangle([0, 0, size - 1, size - 1], radius=radius, fill=255)
    return mask


def circle_mask(size):
    """Alpha mask: white circle on black."""
    mask = Image.new("L", (size, size), 0)
    draw = ImageDraw.Draw(mask)
    draw.ellipse([0, 0, size - 1, size - 1], fill=255)
    return mask


def draw_text(img, text, size, color):
    """
    Draw centred text on img.  Falls back to default font if no system font
    is found, scaling the font size to roughly 40% of the icon size.
    """
    draw = ImageDraw.Draw(img)
    target_h = int(size * 0.42)   # desired cap-height ≈ 42% of icon

    font = None
    # Try common Windows/Linux bold sans-serif fonts
    font_candidates = [
        "arialbd.ttf", "Arial Bold.ttf",
        "segoeuib.ttf", "calibrib.ttf",
        "DejaVuSans-Bold.ttf", "LiberationSans-Bold.ttf",
        "NotoSans-Bold.ttf",
    ]
    for fname in font_candidates:
        try:
            font = ImageFont.truetype(fname, target_h)
            break
        except (IOError, OSError):
            continue

    if font is None:
        # Pillow built-in bitmap font — won't scale nicely but always works
        font = ImageFont.load_default()

    # Measure text
    bbox = draw.textbbox((0, 0), text, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    x = (size - tw) // 2 - bbox[0]
    y = (size - th) // 2 - bbox[1]

    # Subtle drop-shadow for depth
    shadow_offset = max(2, size // 80)
    draw.text((x + shadow_offset, y + shadow_offset), text, font=font,
              fill=(0, 80, 140, 100))
    draw.text((x, y), text, font=font, fill=color)


# ── Build icons ───────────────────────────────────────────────────────────────

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def build_png(size=512):
    """Squircle PNG for web/manifest (transparent background)."""
    grad = radial_gradient(size, CENTER_COLOR, MID_COLOR, EDGE_COLOR)
    mask = rounded_rect_mask(size, radius=int(size * 0.22))
    grad.putalpha(mask)
    draw_text(grad, "Cgl", size, TEXT_COLOR)
    return grad


def build_squircle(size):
    """Squircle (rounded-rect) RGBA image for .ico frames."""
    grad = radial_gradient(size, CENTER_COLOR, MID_COLOR, EDGE_COLOR)
    mask = rounded_rect_mask(size, radius=max(4, int(size * 0.22)))
    grad.putalpha(mask)
    draw_text(grad, "Cgl", size, TEXT_COLOR)
    return grad


def save_png(path, size=512):
    img = build_png(size)
    img.save(path, "PNG")
    print(f"  Saved PNG  → {path}  ({size}×{size})")


def save_ico(path):
    """Save a multi-resolution .ico with squircle frames."""
    sizes = [256, 128, 64, 48, 32, 16]
    frames = []
    for s in sizes:
        frame = build_squircle(s).convert("RGBA")
        frames.append(frame)

    # Save using the largest as the primary, embedding all sizes
    frames[0].save(
        path,
        format="ICO",
        sizes=[(s, s) for s in sizes],
        append_images=frames[1:],
    )
    print(f"  Saved ICO  → {path}  (sizes: {sizes})")


if __name__ == "__main__":
    out_png = os.path.join(ROOT, "icon.png")
    out_ico = os.path.join(ROOT, "icon.ico")

    print("Generating sky-blue gradient icon…")
    save_png(out_png)
    save_ico(out_ico)
    print("Done ✓")
