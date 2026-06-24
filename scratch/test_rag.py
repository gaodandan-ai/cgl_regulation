import urllib.request
import urllib.parse
import json

url = "http://localhost:8000/api/summarize?gene=sigH&name=sigH"
headers = {
    "X-AI-API-Key": "DummyKey-Test",
    "X-AI-Provider": "google"
}

req = urllib.request.Request(url, headers=headers)
try:
    with urllib.request.urlopen(req) as resp:
        res = json.loads(resp.read().decode('utf-8'))
        print("Status Code: 200")
        print("Response keys:", list(res.keys()))
        print("Summary:", res.get("summary"))
        print("PubMed Papers count:", len(res.get("papers", [])))
        print("RAG Sources:", res.get("rag_sources"))
except Exception as e:
    print("Request failed:", e)
