import os

app_js_path = os.path.join("web", "app.js")

with open(app_js_path, "r", encoding="utf-8") as f:
    content = f.read()

# Let's define the function we want to insert
resizer_code = """

function initSidebarResizer() {
    const resizer = document.getElementById('sidebar-resizer');
    if (!resizer || !rightSidebar) return;

    let isResizing = false;

    resizer.addEventListener('mousedown', (e) => {
        isResizing = true;
        document.body.classList.add('sidebar-no-transition');
        resizer.classList.add('resizing');
        e.preventDefault();
    });

    document.addEventListener('mousemove', (e) => {
        if (!isResizing) return;
        
        const newWidth = window.innerWidth - e.clientX;
        
        if (newWidth >= 280 && newWidth <= 800) {
            document.documentElement.style.setProperty('--right-sidebar-width', `${newWidth}px`);
        }
    });

    document.addEventListener('mouseup', () => {
        if (isResizing) {
            isResizing = false;
            document.body.classList.remove('sidebar-no-transition');
            resizer.classList.remove('resizing');
        }
    });
}
"""

# Let's search for the end of initEventListeners function
# It has:
#         forwardBtn.addEventListener('click', () => {
#             navigateHistory('forward');
#         });
#     }
# }

# Let's try both \r\n and \n line endings for locating
search_str_crlf = "forwardBtn.addEventListener('click', () => {\r\n\r\n            navigateHistory('forward');\r\n\r\n        });\r\n\r\n    }\r\n\r\n}"
search_str_lf = "forwardBtn.addEventListener('click', () => {\n\n            navigateHistory('forward');\n\n        });\n\n    }\n\n}"

idx = content.find(search_str_crlf)
if idx != -1:
    print("Found with CRLF at:", idx)
    insert_pos = idx + len(search_str_crlf)
    new_content = content[:insert_pos] + resizer_code.replace("\n", "\r\n") + content[insert_pos:]
else:
    idx = content.find(search_str_lf)
    if idx != -1:
        print("Found with LF at:", idx)
        insert_pos = idx + len(search_str_lf)
        new_content = content[:insert_pos] + resizer_code + content[insert_pos:]
    else:
        print("Could not find the target position using standard templates.")
        # Fallback search
        fallback_crlf = "navigateHistory('forward');\r\n\r\n        });\r\n\r\n    }\r\n\r\n}"
        idx = content.find(fallback_crlf)
        if idx != -1:
            print("Found with fallback CRLF at:", idx)
            insert_pos = idx + len(fallback_crlf)
            new_content = content[:insert_pos] + resizer_code.replace("\n", "\r\n") + content[insert_pos:]
        else:
            new_content = None

if new_content:
    with open(app_js_path, "w", encoding="utf-8") as f:
        f.write(new_content)
    print("Successfully inserted initSidebarResizer function.")
else:
    print("Error: Could not modify file.")
