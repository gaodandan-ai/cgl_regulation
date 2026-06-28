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
            # print(f"Out of string at line {line_num}")
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
        string_start_line = line_num
        idx += 1
        continue
        
    idx += 1

print("At end:")
print("in_string:", in_string, "char:", string_char, "started at line:", string_start_line if in_string else None)
print("in_line_comment:", in_line_comment)
print("in_block_comment:", in_block_comment)
