import re

with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

# Let's replace comments first
def strip_comments_and_strings(js_code):
    # Regex to match:
    # 1. block comments: /* ... */
    # 2. line comments: // ...
    # 3. double quoted strings: " ... "
    # 4. single quoted strings: ' ... '
    # 5. template literals: ` ... ` (we will handle ${} by replacing backtick strings)
    # 6. regexes: /.../
    
    # We will parse character by character to be 100% correct, including template literal nesting.
    
    in_string = False
    string_char = None
    escaped = False
    in_line_comment = False
    in_block_comment = False
    
    # We will also keep track of regex. In JS, a slash starts a regex if we are not in a string/comment
    # and the previous non-whitespace character was one of:
    # =, (, ,, [, !, &, |, ?, :, {, ;, +, -, *, ^, ~, %, >, <, return, throw, typeof, yield, delete, void, instanceof, in, case
    # Let's keep a history of non-whitespace characters to determine if a / is division or regex.
    
    processed = []
    idx = 0
    n = len(js_code)
    
    line_num = 1
    
    # To determine if / is regex, we can track the last non-whitespace word/char.
    last_tokens = []
    
    # Let's build a simple helper to check if previous token allows regex
    regex_triggers = {
        '=', '(', ',', '[', '!', '&', '|', '?', ':', '{', ';', '+', '-', '*', '^', '~', '%', '>', '<',
        'return', 'throw', 'typeof', 'yield', 'delete', 'void', 'instanceof', 'in', 'case', '}'
    }
    
    while idx < n:
        char = js_code[idx]
        
        if char == '\n':
            line_num += 1
            in_line_comment = False
            processed.append('\n')
            idx += 1
            continue
            
        if in_line_comment:
            idx += 1
            continue
            
        if in_block_comment:
            if char == '*' and idx + 1 < n and js_code[idx+1] == '/':
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
            
        # Check comments
        if char == '/' and idx + 1 < n and js_code[idx+1] == '/':
            in_line_comment = True
            idx += 2
            continue
            
        if char == '/' and idx + 1 < n and js_code[idx+1] == '*':
            in_block_comment = True
            idx += 2
            continue
            
        # Check string literals
        if char in ('\"', '\'', '`'):
            in_string = True
            string_char = char
            idx += 1
            continue
            
        # Check regex
        if char == '/':
            # Is this division or regex?
            # Let's look backward at last non-whitespace chars/words
            is_regex = False
            # Find last non-whitespace character in processed
            last_char = None
            for c in reversed(processed):
                if c not in (' ', '\n', '\t', '\r'):
                    last_char = c
                    break
            
            # Simple heuristic: if last_char is in regex_triggers or is None, it is a regex!
            if last_char is None or last_char in regex_triggers or last_char.isalnum() is False:
                # Except if last_char is close paren/bracket/identifier, then division
                # wait, in JS, `1 / 2` or `x / y` is division.
                # So if last_char is alphanumeric, or ) or ], it's division.
                if last_char in (')', ']') or (last_char and (last_char.isalnum() or last_char == '_')):
                    is_regex = False
                else:
                    is_regex = True
            
            if is_regex:
                # Scan until end of regex literal
                # regex ends at next unescaped /
                reg_idx = idx + 1
                reg_escaped = False
                in_char_class = False # [ / ] inside character class doesn't close regex
                while reg_idx < n:
                    r_char = js_code[reg_idx]
                    if r_char == '\n':
                        break # regex cannot span lines unless escaped (rare)
                    if reg_escaped:
                        reg_escaped = False
                    elif r_char == '\\':
                        reg_escaped = True
                    elif r_char == '[':
                        in_char_class = True
                    elif r_char == ']':
                        in_char_class = False
                    elif r_char == '/' and not in_char_class:
                        reg_idx += 1
                        break
                    reg_idx += 1
                idx = reg_idx
                continue
                
        # Keep track of last non-whitespace char in processed
        processed.append(char)
        idx += 1
        
    return "".join(processed)

stripped = strip_comments_and_strings(content)

# Now check braces in stripped code
depth = 0
stack = []
line_num = 1
for idx, char in enumerate(stripped):
    if char == '\n':
        line_num += 1
    elif char == '{':
        depth += 1
        start = max(0, idx - 40)
        end = min(len(stripped), idx + 40)
        context = stripped[start:end].replace('\n', ' ')
        stack.append((line_num, context))
    elif char == '}':
        depth -= 1
        if stack:
            stack.pop()
        else:
            print(f"Excess closing brace at stripped line {line_num}")

print(f"Final brace depth: {depth}")
if stack:
    print("Unclosed braces:")
    for l, ctx in stack:
        print(f"Stripped Line {l}: ... {ctx.strip()} ...")
