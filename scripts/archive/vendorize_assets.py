import urllib.request, os

vendor_dir = r"f:\cgl_regulation\web\lib\vendor"
os.makedirs(vendor_dir, exist_ok=True)

downloads = [
    ("https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js", os.path.join(vendor_dir, "chart.min.js")),
    ("https://cdn.jsdelivr.net/npm/mermaid@10.6.1/dist/mermaid.min.js", os.path.join(vendor_dir, "mermaid.min.js")),
    ("https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css", os.path.join(vendor_dir, "fontawesome.min.css")),
]

for url, target in downloads:
    if os.path.exists(target) and os.path.getsize(target) > 1000:
        print(f"Already exists: {target}")
        continue
    print(f"Downloading {url} -> {target}...")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = resp.read()
        with open(target, "wb") as f:
            f.write(data)
        print(f"Successfully saved {target} ({len(data)} bytes)")
    except Exception as e:
        print(f"Download failed for {url}: {e}")
