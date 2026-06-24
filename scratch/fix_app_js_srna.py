import os

def main():
    filepath = os.path.join("web", "app.js")
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    # 1. Update loadNetworkData initial panel display
    old_load = "srnaThresholdPanel.classList.remove('hidden');"
    new_load = """if (filterSrna.checked) {
                srnaThresholdPanel.classList.remove('hidden');
            } else {
                srnaThresholdPanel.classList.add('hidden');
            }"""
            
    if old_load in content:
        content = content.replace(old_load, new_load)
        print("Updated initial load slider logic.")
    else:
        print("Warning: could not find old load target.")

    # 2. Update filterSrna event listener
    old_listener = "filterSrna.addEventListener('change', reRender);"
    new_listener = """filterSrna.addEventListener('change', () => {
        if (filterSrna.checked) {
            srnaThresholdPanel.classList.remove('hidden');
        } else {
            srnaThresholdPanel.classList.add('hidden');
        }
        reRender();
    });"""

    if old_listener in content:
        content = content.replace(old_listener, new_listener)
        print("Updated filterSrna change listener logic.")
    else:
        print("Warning: could not find old listener target.")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)
    print("app.js updated successfully!")

if __name__ == "__main__":
    main()
