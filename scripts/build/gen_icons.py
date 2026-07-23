#!/usr/bin/env python3
"""
Regenerate crisp multi-size ICO from icon.png (512x512 source).
PIL's ICO append_images is unreliable — we write ICO binary manually.
"""
from PIL import Image
import struct, io, os

src = Image.open('icon.png').convert('RGBA')

# ── 1. Generate individual PNGs for web/ ──────────────────────────────────────
for s in [16, 32, 48, 64, 96, 128, 192, 256, 512]:
    img = src.resize((s, s), Image.LANCZOS)
    img.save('web/icon-{}.png'.format(s), 'PNG', optimize=True)

src.resize((16, 16), Image.LANCZOS).save('web/favicon-16.png', 'PNG', optimize=True)
src.resize((32, 32), Image.LANCZOS).save('web/favicon-32.png', 'PNG', optimize=True)
src.resize((64, 64), Image.LANCZOS).save('web/favicon.png',   'PNG', optimize=True)
print('PNGs generated.')

# ── 2. Build proper multi-size ICO manually ────────────────────────────────────
ICO_SIZES = [16, 24, 32, 48, 64, 128, 256]

def make_png_bytes(img, size):
    buf = io.BytesIO()
    resized = img.resize((size, size), Image.LANCZOS)
    resized.save(buf, format='PNG')
    return buf.getvalue()

images_bytes = [make_png_bytes(src, s) for s in ICO_SIZES]

# ICO format:
#   Header: 6 bytes
#   Directory entries: n * 16 bytes
#   Image data: concatenated PNG blobs
n = len(ICO_SIZES)
header = struct.pack('<HHH', 0, 1, n)   # reserved, type=1 (ICO), count

dir_offset = 6 + n * 16  # offset of first image data
entries = b''
data_parts = b''
current_offset = dir_offset

for i, (s, blob) in enumerate(zip(ICO_SIZES, images_bytes)):
    w = 0 if s == 256 else s   # ICO uses 0 to mean 256
    h = 0 if s == 256 else s
    # width, height, colorCount, reserved, planes, bitCount, sizeInBytes, offset
    entry = struct.pack('<BBBBHHII', w, h, 0, 0, 1, 32, len(blob), current_offset)
    entries += entry
    data_parts += blob
    current_offset += len(blob)

ico_bytes = header + entries + data_parts

for path in ['web/favicon.ico', 'icon.ico']:
    with open(path, 'wb') as f:
        f.write(ico_bytes)
    print('  Saved {} ({:,} bytes, {} sizes: {})'.format(
        path, os.path.getsize(path), n, ICO_SIZES))

print('Done.')
