with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

depth = 0
stack = []

in_string = False
string_char = None
escaped = False
in_line_comment = False
in_block_comment = False

idx = 0
n = len(content)

line_num = 1
col_num = 1

while idx < n:
    char = content[idx]
    
    if char == '\n':
        line_num += 1
        col_num = 1
        in_line_comment = False
        idx += 1
        continue
    
    if in_line_comment:
        idx += 1
        col_num += 1
        continue
        
    if in_block_comment:
        if char == '*' and idx + 1 < n and content[idx+1] == '/':
            in_block_comment = False
            idx += 2
            col_num += 2
        else:
            idx += 1
            col_num += 1
        continue
        
    if in_string:
        if escaped:
            escaped = False
        elif char == '\\':
            escaped = True
        elif char == string_char:
            in_string = False
        idx += 1
        col_num += 1
        continue
        
    # Check for comments
    if char == '/' and idx + 1 < n and content[idx+1] == '/':
        in_line_comment = True
        idx += 2
        col_num += 2
        continue
        
    if char == '/' and idx + 1 < n and content[idx+1] == '*':
        in_block_comment = True
        idx += 2
        col_num += 2
        continue
        
    # Check for string literals
    if char in ('\"', '\'', '`'):
        in_string = True
        string_char = char
        idx += 1
        col_num += 1
        continue
        
    # Check for braces
    if char == '{':
        depth += 1
        # Extract surrounding context
        start = max(0, idx - 40)
        end = min(n, idx + 40)
        context = content[start:end].replace('\n', ' ')
        stack.append((line_num, col_num, context))
    elif char == '}':
        depth -= 1
        if stack:
            stack.pop()
        else:
            print(f"Excess closing brace at line {line_num}, col {col_num}")
            
    idx += 1
    col_num += 1

print(f"Final brace depth: {depth}")
if stack:
    print("Unclosed braces:")
    for l, c, ctx in stack:
        print(f"Line {l}, Col {c}: ... {ctx.strip()} ...")
