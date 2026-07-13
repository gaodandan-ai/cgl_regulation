"""
scripts/generate_ico.py
=======================
Programmatically creates a beautiful high-resolution systems biology app icon
using Pillow (PIL) and converts it to a standard Windows icon.ico file.

Features:
- Dark indigo background with a radial gradient effect.
- Glowing interconnected nodes (neon teal & neon purple).
- Stylized double helix in the center.
- Outputs multi-size bundled icon.ico.
"""

import os
import math
from PIL import Image, ImageDraw, ImageFilter

def create_gradient_background(size):
    # Create dark indigo radial gradient
    base = Image.new("RGBA", (size, size), (15, 23, 42, 255)) # tailwind slate-900
    draw = ImageDraw.Draw(base)
    
    # Draw radial gradient layers
    center = size // 2
    for r in range(size // 2, 0, -8):
        # Calculate fade
        ratio = r / (size // 2)
        alpha = int(80 * (1.0 - ratio))
        color = (59, 130, 246, alpha) # blue-500 with alpha
        draw.ellipse([center - r, center - r, center + r, center + r], fill=color)
        
    return base

def draw_network_nodes(image, size):
    draw = ImageDraw.Draw(image)
    center = size // 2
    
    # Colors
    teal = (20, 184, 166, 255) # neon teal
    purple = (168, 85, 247, 255) # neon purple
    white = (255, 255, 255, 255)
    
    # 1. Draw network edges (connecting lines)
    # Define some coordinates for node centers
    nodes = [
        (center - 120, center - 120, teal),
        (center + 120, center - 100, purple),
        (center - 140, center + 100, purple),
        (center + 110, center + 120, teal),
        (center - 50, center - 160, teal),
        (center + 50, center + 160, purple),
        (center, center, white) # Central node (Cgl)
    ]
    
    # Connect nodes
    connections = [
        (0, 4), (4, 1), (1, 6), (0, 6), (2, 6), (3, 6), (2, 5), (5, 3), (0, 2), (1, 3)
    ]
    
    for i, j in connections:
        n1 = nodes[i]
        n2 = nodes[j]
        # Draw thick connection line with gradient-like look (using average color)
        avg_color = tuple((n1[2][c] + n2[2][c]) // 2 for c in range(3)) + (120,)
        draw.line([n1[0], n1[1], n2[0], n2[1]], fill=avg_color, width=4)

    # 2. Draw nodes (glowing circles)
    # Create a separate layer for glows to apply gaussian blur
    glow_layer = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    glow_draw = ImageDraw.Draw(glow_layer)
    
    for x, y, col in nodes:
        # Draw glow
        r_glow = 25
        glow_col = col[:3] + (100,)
        glow_draw.ellipse([x - r_glow, y - r_glow, x + r_glow, y + r_glow], fill=glow_col)
        
    # Blur the glows
    glow_layer = glow_layer.filter(ImageFilter.GaussianBlur(10))
    image = Image.alpha_composite(image, glow_layer)
    
    # Draw solid node centers on top
    draw = ImageDraw.Draw(image)
    for x, y, col in nodes:
        # Outer ring
        draw.ellipse([x - 12, y - 12, x + 12, y + 12], fill=col)
        # Inner white core
        draw.ellipse([x - 5, y - 5, x + 5, y + 5], fill=(255, 255, 255, 255))
        
    return image

def draw_dna_helix(image, size):
    # Draw a stylized DNA double helix in the center background/foreground
    draw = ImageDraw.Draw(image)
    center = size // 2
    
    # Generate helix points
    points_a = []
    points_b = []
    
    y_start = center - 80
    y_end = center + 80
    
    for y in range(y_start, y_end, 4):
        # Calculate sine wave offset
        angle = (y - y_start) * 0.05
        x_offset = math.sin(angle) * 45
        
        points_a.append((center + x_offset, y))
        points_b.append((center - x_offset, y))
        
    # Draw DNA rungs (connecting lines between strands)
    for idx in range(0, len(points_a), 6):
        pa = points_a[idx]
        pb = points_b[idx]
        # Alternate rung colors
        col = (20, 184, 166, 180) if idx % 12 == 0 else (168, 85, 247, 180)
        draw.line([pa[0], pa[1], pb[0], pb[1]], fill=col, width=3)
        
    # Draw DNA strands (lines connecting points)
    for i in range(len(points_a) - 1):
        draw.line([points_a[i][0], points_a[i][1], points_a[i+1][0], points_a[i+1][1]], fill=(255, 255, 255, 220), width=6)
        draw.line([points_b[i][0], points_b[i][1], points_b[i+1][0], points_b[i+1][1]], fill=(255, 255, 255, 220), width=6)
        
    # Draw DNA base pair circles
    for idx in range(0, len(points_a), 6):
        pa = points_a[idx]
        pb = points_b[idx]
        draw.ellipse([pa[0] - 4, pa[1] - 4, pa[0] + 4, pa[1] + 4], fill=(255, 255, 255, 255))
        draw.ellipse([pb[0] - 4, pb[1] - 4, pb[0] + 4, pb[1] + 4], fill=(255, 255, 255, 255))
        
    return image

def main():
    size = 512
    print("Generating base vector artwork...")
    img = create_gradient_background(size)
    img = draw_dna_helix(img, size)
    img = draw_network_nodes(img, size)
    
    # Ensure outputs/ directory exists or write to repository root
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    png_path = os.path.join(ROOT, "icon.png")
    ico_path = os.path.join(ROOT, "icon.ico")
    
    # Save high-res PNG
    img.save(png_path, "PNG")
    print(f"Saved high-res icon PNG to: {png_path}")
    
    # Save standard Windows multi-size ICO
    # Recommended standard sizes for Windows
    sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
    img.save(ico_path, format="ICO", sizes=sizes)
    print(f"Saved bundled multi-resolution Windows ICO to: {ico_path}")
    
if __name__ == "__main__":
    main()
