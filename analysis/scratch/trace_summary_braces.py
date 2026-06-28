with open('web/app.js', 'r', encoding='utf-8') as f:
    lines = f.readlines()

depth = 0
for idx in range(4131, 4940):
    line = lines[idx-1]
    clean_line = line.split('//')[0]
    # Simple character loop ignoring strings
    in_string = False
    string_char = None
    escaped = False
    for char in clean_line:
        if in_string:
            if escaped: escaped = False
            elif char == '\\': escaped = True
            elif char == string_char: in_string = False
        else:
            if char in ('\"', '\'', '`'):
                in_string = True
                string_char = char
            elif char == '{':
                depth += 1
                print(f"Line {idx} depth +1 -> {depth}: {line.strip()}")
            elif char == '}':
                depth -= 1
                print(f"Line {idx} depth -1 -> {depth}: {line.strip()}")
