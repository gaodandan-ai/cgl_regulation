"""scripts/add_thermo_conf_div.py — adds thermo confidence div to engineering sim modal in index.html"""

THERMO_DIV = '\n                        <!-- Thermodynamic Confidence -->\n                        <div id="engineering-sim-thermo-confidence" style="font-size:11px;color:var(--text-secondary);"></div>'

with open('web/index.html', 'r', encoding='utf-8') as f:
    content = f.read()

# The unique string right before the Close button div
MARKER = '                        </div>\n                        <div style="display: flex; justify-content: flex-end; margin-top: 4px;">'
REPLACEMENT = '                        </div>' + THERMO_DIV + '\n                        <div style="display: flex; justify-content: flex-end; margin-top: 4px;">'

count = content.count(MARKER)
print(f"Marker occurrences: {count}")

if count == 1:
    content = content.replace(MARKER, REPLACEMENT, 1)
    with open('web/index.html', 'w', encoding='utf-8') as f:
        f.write(content)
    print("SUCCESS: thermo confidence div added to sim modal")
else:
    print("ERROR: marker not found or ambiguous")
    idx = content.find('engineering-sim-summary-table')
    print(f"  Found sim table at char {idx}")
