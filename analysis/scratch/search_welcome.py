with open("web/app.js", "r", encoding="utf-8") as f:
    for idx, line in enumerate(f):
        if "canvas-overlay" in line or "canvasOverlay" in line:
            print(f"Line {idx+1}: {line.strip()}")
