# inspect_app.py
with open('web/app.js', 'r', encoding='utf-8') as f:
    content = f.read()

# Let's search for function showNodeDetails(locusTag) and its scope
# To find the end of the function, let's search for showNodeDetails in function names
# We know it starts at line 1885. Let's search for the start index
start_idx = content.find("function showNodeDetails(locusTag) {")
print("Starts at index:", start_idx)

# Let's see surrounding text at the end of the function
# From the previous python script output, we saw detail table populating or links loading
# Let's search for some text that is after showNodeDetails, like function showOperonDetails
next_func_idx = content.find("function showOperonDetails")
print("Next function starts at index:", next_func_idx)

# Let's print the code between lines 2500 and 2700
with open('web/app.js', 'r', encoding='utf-8') as f:
    lines = f.readlines()
for idx in range(2500, 2600):
    if idx < len(lines):
        print(f"Line {idx+1}: {lines[idx]}", end="")
