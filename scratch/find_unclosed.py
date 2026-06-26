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

open_strings = [] # stack of (line_num, char_idx, quote_char)

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
            start_line, start_idx, _ = open_strings.pop()
            if line_num - start_line > 50:
                print(f"Long string closed: opened at line {start_line}, closed at line {line_num}")
                snippet = content[start_idx:start_idx+100].replace('\n', ' ')
                print(f"   Snippet: ... {snippet} ...")
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
        open_strings.append((line_num, idx, char))
        idx += 1
        continue
        
    idx += 1

print(f"Remaining open strings in stack: {len(open_strings)}")
for start_line, start_idx, q in open_strings:
    snippet = content[start_idx:start_idx+100].replace('\n', ' ')
    print(f"Unclosed string of type {q} opened at line {start_line}: ... {snippet} ...")
