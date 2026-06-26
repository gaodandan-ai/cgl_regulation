with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

in_string = False
string_char = None
escaped = False
in_line_comment = False
in_block_comment = False

idx = 0
n = len(content)
line_num = 1

while idx < n:
    char = content[idx]
    if char == '\n':
        line_num += 1
        in_line_comment = False
        
    if line_num >= 5940 and line_num <= 5975:
        print(f"L{line_num} Char {idx} {repr(char)}: in_str={in_string}, string_char={string_char}, escaped={escaped}")

    if in_line_comment:
        idx += 1
        continue
    if in_block_comment:
        if char == '*' and idx + 1 < n and content[idx+1] == '/':
            in_block_comment = False
            idx += 2
        else:
            idx += 1
        continue
    if in_string:
        if escaped:
            escaped = False
        elif char == '\\':
            escaped = True
        elif char == string_char:
            in_string = False
            string_char = None
        idx += 1
        continue
        
    if char == '/' and idx + 1 < n and content[idx+1] == '/':
        in_line_comment = True
        idx += 2
        continue
    if char == '/' and idx + 1 < n and content[idx+1] == '*':
        in_block_comment = True
        idx += 2
        continue
    if char in ('\"', '\'', '`'):
        in_string = True
        string_char = char
        idx += 1
        continue
        
    idx += 1
