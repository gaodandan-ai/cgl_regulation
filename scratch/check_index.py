with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

for i in range(168800, 169000):
    if i < len(content):
        print(f"Index {i}: {repr(content[i])} code: {ord(content[i])}")
