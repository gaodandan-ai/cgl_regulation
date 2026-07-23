"""
scripts/generate_ico.py
=======================
Programmatically creates a beautiful high-resolution systems biology app icon
using Pillow (PIL) and converts it to a standard Windows icon.ico file.

Features:
- Solid/gradient royal blue background inside a modern rounded squircle.
- "Cgl" text in bold white sans-serif font centered on the canvas.
- Subtle inner border and drop shadow for a premium 3D desktop app look.
"""

import os
import sys
from PIL import Image, ImageDraw, ImageFilter, ImageFont

def create_app_icon(size):
    # 1. Create a transparent base image
    base = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(base)
    
    # 2. Draw a modern rounded rectangle (squircle) background with gradient
    # Define bounds (padding of 20px for a clean spacing on the desktop)
    pad = 24
    rect_bounds = [pad, pad, size - pad, size - pad]
    radius = 110 # smooth round corner
    
    # Generate gradient fill
    gradient_img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    g_draw = ImageDraw.Draw(gradient_img)
    
    # Draw rounded rectangle base
    g_draw.rounded_rectangle(rect_bounds, radius=radius, fill=(29, 78, 216, 255)) # Royal Blue base (blue-700)
    
    # Radial royal blue gradient overlay
    center = size // 2
    for r in range(size // 2, 0, -4):
        ratio = r / (size // 2)
        alpha = int(140 * (1.0 - ratio))
        color = (96, 165, 250, alpha) # Bright blue (blue-400)
        # Use a temporary image to mask to the rounded rectangle
        mask_layer = Image.new("L", (size, size), 0)
        m_draw = ImageDraw.Draw(mask_layer)
        m_draw.rounded_rectangle(rect_bounds, radius=radius, fill=255)
        
        # Apply radial glow layers
        g_draw.ellipse([center - r, center - r, center + r, center + r], fill=color)
        
    # Mask the gradient image to only the rounded rect
    mask = Image.new("L", (size, size), 0)
    m_draw = ImageDraw.Draw(mask)
    m_draw.rounded_rectangle(rect_bounds, radius=radius, fill=255)
    
    final_bg = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    final_bg.paste(gradient_img, (0, 0), mask=mask)
    
    # 3. Add a fine, semi-transparent white inner border for a premium glassmorphism outline
    draw_bg = ImageDraw.Draw(final_bg)
    draw_bg.rounded_rectangle(rect_bounds, radius=radius, outline=(255, 255, 255, 40), width=4)
    
    # 4. Load bold clean sans-serif font
    font_size = 180
    font = None
    
    # Try various clean bold sans-serif system fonts
    font_paths = [
        "C:\\Windows\\Fonts\\segoeuib.ttf",  # Segoe UI Bold (Standard Win 10/11)
        "C:\\Windows\\Fonts\\arialbd.ttf",   # Arial Bold
        "C:\\Windows\\Fonts\\trebucbd.ttf",  # Trebuchet MS Bold
        "C:\\Windows\\Fonts\\calibrib.ttf",   # Calibri Bold
    ]
    
    for path in font_paths:
        if os.path.exists(path):
            try:
                font = ImageFont.truetype(path, font_size)
                print(f"Loaded font: {os.path.basename(path)}")
                break
            except Exception:
                pass
                
    if font is None:
        font = ImageFont.load_default()
        print("Using default fallback font.")

    # 5. Draw text "Cgl" in the center with a soft drop shadow
    # Text content
    text_str = "Cgl"
    
    # Calculate text bounding box to align perfectly in center
    draw_txt = ImageDraw.Draw(final_bg)
    try:
        bbox = draw_txt.textbbox((0, 0), text_str, font=font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]
        tx = (size - text_w) // 2 - bbox[0]
        ty = (size - text_h) // 2 - bbox[1] - 10 # slightly offset up to visually center
    except AttributeError:
        # Older PIL versions fallback
        text_w, text_h = draw_txt.textsize(text_str, font=font)
        tx = (size - text_w) // 2
        ty = (size - text_h) // 2
        
    # Draw soft shadow layer for text depth
    shadow_layer = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    s_draw = ImageDraw.Draw(shadow_layer)
    s_draw.text((tx + 4, ty + 6), text_str, fill=(15, 23, 42, 100), font=font) # Dark drop shadow
    shadow_layer = shadow_layer.filter(ImageFilter.GaussianBlur(8))
    
    # Composite background, shadow, and final white text
    final_bg = Image.alpha_composite(final_bg, shadow_layer)
    draw_final = ImageDraw.Draw(final_bg)
    draw_final.text((tx, ty), text_str, fill=(255, 255, 255, 255), font=font)
    
    return final_bg

def main():
    size = 512
    print("Generating Cgl royal blue app icon...")
    img = create_app_icon(size)
    
    # Save target paths
    ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    png_path = os.path.join(ROOT, "icon.png")
    ico_path = os.path.join(ROOT, "icon.ico")
    
    # Save high-res PNG
    img.save(png_path, "PNG")
    print(f"Saved high-res icon PNG to: {png_path}")
    
    # Save standard Windows multi-size ICO
    sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
    img.save(ico_path, format="ICO", sizes=sizes)
    print(f"Saved bundled multi-resolution Windows ICO to: {ico_path}")

if __name__ == "__main__":
    main()
