import urllib.request, os

vendor_dir = r"f:\cgl_regulation\web\lib\vendor"
os.makedirs(vendor_dir, exist_ok=True)

downloads = [
    ("https://cdn.jsdelivr.net/npm/hammerjs@2.0.8/hammer.min.js", os.path.join(vendor_dir, "hammer.min.js")),
    ("https://cdn.jsdelivr.net/npm/chartjs-plugin-zoom@2.0.1/dist/chartjs-plugin-zoom.min.js", os.path.join(vendor_dir, "chartjs-plugin-zoom.min.js")),
]

for url, target in downloads:
    if os.path.exists(target) and os.path.getsize(target) > 500:
        print(f"Already exists: {target}")
        continue
    print(f"Downloading {url} -> {target}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
        with open(target, "wb") as f:
            f.write(data)
        print(f"Saved {target} ({len(data)} bytes)")
    except Exception as e:
        print(f"Failed {url}: {e}")
