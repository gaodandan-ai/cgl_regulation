with open('web/app.js', 'r', encoding='utf-8') as f:
    lines = f.readlines()

line = lines[7508 - 1]
print(f"Line 7508: {repr(line)}")
for idx, char in enumerate(line):
    print(f"Char {idx}: {repr(char)} code: {ord(char)}")
