# test_http_get.py
import urllib.request
import json

headers = {
    'X-AI-API-Key': 'DummyKey',
    'X-AI-Provider': 'google'
}

def test_url(url):
    print(f"Testing GET {url}...")
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            print("Status: Success")
            print("Response summary keys:", list(data.keys()))
            if "summary" in data:
                print("Summary preview (100 chars):", data["summary"][:100] + "...")
            if "error" in data:
                print("Error in response:", data["error"])
    except Exception as e:
        print("Status: Failed")
        print("Error:", e)

test_url("http://localhost:8000/api/protein_domain?gene=cg0350")
print()
test_url("http://localhost:8000/api/binding_site?gene=cg0350")
