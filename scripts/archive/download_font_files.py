import urllib.request, os

font_dir = r"f:\cgl_regulation\web\lib\webfonts"
os.makedirs(font_dir, exist_ok=True)

font_files = [
    "fa-solid-900.woff2", "fa-solid-900.ttf",
    "fa-brands-400.woff2", "fa-brands-400.ttf",
    "fa-regular-400.woff2", "fa-regular-400.ttf",
    "fa-v4compatibility.woff2", "fa-v4compatibility.ttf"
]

base_url = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/webfonts/"

for font in font_files:
    target = os.path.join(font_dir, font)
    if os.path.exists(target) and os.path.getsize(target) > 500:
        continue
    url = base_url + font
    print(f"Downloading {url} -> {target}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
        with open(target, "wb") as f:
            f.write(data)
        print(f"Saved {font} ({len(data)} bytes)")
    except Exception as e:
        print(f"Failed {font}: {e}")
