with open('web/app.js', 'r', encoding='utf-8') as f:
    lines = f.readlines()

line = lines[7671 - 1]
print(f"Line 7671: {repr(line)}")
for idx, char in enumerate(line):
    print(f"Char {idx}: {repr(char)} code: {ord(char)}")
