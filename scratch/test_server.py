import subprocess
import time
import urllib.request
import urllib.parse
import json
import sys

def test_api_endpoint(url, headers):
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            return data
    except Exception as e:
        return {"error_request": str(e)}

def main():
    print("Starting verification server on port 8005...")
    # Modify port temporarily to avoid conflict with running instances
    # We can launch run_server.py by passing PORT environment variable or modifying the file.
    # Actually, in run_server.py: PORT = 8000.
    # We can write a temporary test_server_run.py that changes PORT to 8005 and runs it.
    with open("run_server.py", "r", encoding="utf-8") as f:
        code = f.read()
    
    test_code = code.replace("PORT = 8000", "PORT = 8005")
    # Disable webbrowser.open in the test run to avoid popping up browser tabs
    test_code = test_code.replace("webbrowser.open(url)", "print('MOCK: webbrowser open', url)")
    
    with open("scratch/test_server_run.py", "w", encoding="utf-8") as f:
        f.write(test_code)
        
    server_process = subprocess.Popen([sys.executable, "scratch/test_server_run.py"])
    time.sleep(2.0) # Wait for server to boot up
    
    success = True
    try:
        # Test 1: Summarize gene with Google Gemini (DummyKey)
        print("Test 1: Querying summarize with Google (DummyKey)...")
        headers = {
            "X-AI-API-Key": "DummyKey",
            "X-AI-Provider": "google"
        }
        res = test_api_endpoint("http://localhost:8005/api/summarize?gene=cg0350&name=whiB4", headers)
        print("Result:", json.dumps(res, ensure_ascii=False))
        if "summary" not in res or "DUMMY_MODE" not in res["summary"]:
            print("FAILED Test 1")
            success = False
            
        # Test 2: Gene Assistant with DeepSeek (DummyKey)
        print("Test 2: Querying gene assistant with DeepSeek (DummyKey)...")
        headers = {
            "X-AI-API-Key": "DummyKey",
            "X-AI-Provider": "deepseek"
        }
        res = test_api_endpoint("http://localhost:8005/api/gene_assistant?query=" + urllib.parse.quote("抗逆"), headers)
        print("Result:", json.dumps(res, ensure_ascii=False))
        if "summary" not in res or "应激反应机制" not in res["summary"]:
            print("FAILED Test 2")
            success = False

        # Test 3: Pathway analysis with Qwen (DummyKey)
        print("Test 3: Querying pathway analysis with Qwen (DummyKey)...")
        headers = {
            "X-AI-API-Key": "DummyKey",
            "X-AI-Provider": "qwen"
        }
        res = test_api_endpoint("http://localhost:8005/api/pathway?pathway=biotin", headers)
        print("Result:", json.dumps(res, ensure_ascii=False))
        if "summary" not in res or "生物素" not in res["summary"]:
            print("FAILED Test 3")
            success = False

        # Test 4: Ollama with empty key (DummyKey fallback in server)
        # Note: In mock/dummy key mode, if we pass "DummyKey" as API Key it returns DUMMY_MODE.
        # But if we pass empty key for Ollama, it will attempt actual HTTP call to http://localhost:11434/v1.
        # Since Ollama is likely not running on the test host, it should fail with "Connection refused" instead of "未提供 API Key".
        # This confirms that the key validation was successfully bypassed for Ollama!
        print("Test 4: Querying pathway with Ollama (Empty Key)...")
        headers = {
            "X-AI-API-Key": "",
            "X-AI-Provider": "ollama"
        }
        res = test_api_endpoint("http://localhost:8005/api/pathway?pathway=biotin", headers)
        print("Result:", json.dumps(res, ensure_ascii=False))
        if "Connection refused" in res.get("error", "") or "connection refused" in res.get("error", "").lower():
            print("Success: Bypassed API Key requirement for Ollama as expected, and attempted connection to local port.")
        elif "未提供 API Key" in res.get("error", ""):
            print("FAILED Test 4: Ollama requires API Key but it should be optional.")
            success = False
        else:
            print("Ollama test warning/info:", res)

    finally:
        print("Stopping verification server...")
        server_process.terminate()
        server_process.wait()
        print("Server stopped.")
        
    if success:
        print("\nAll backend checks PASSED successfully!")
    else:
        print("\nSome backend checks FAILED.")
        sys.exit(1)

if __name__ == "__main__":
    main()
