import re

with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

# Find querySingleGene definition
for m in re.finditer(r'function querySingleGene', content):
    start = max(0, content.rfind('\n', 0, m.start()))
    end = content.find('\n', m.end())
    line_num = content.count('\n', 0, m.start()) + 1
    print(f"Line {line_num}: {content[start:end].strip()}")
