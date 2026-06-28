with open("web/app.js", "r", encoding="utf-8") as f:
    lines = f.readlines()

for idx in range(10170, 10210):
    if idx < len(lines):
        print(f"{idx+1}: {lines[idx].strip()}")
