import urllib.request, json
r = urllib.request.urlopen("http://localhost:8000/api/network?gene=cg0337")
d = json.loads(r.read())
print("keys:", list(d.keys()))
if d.get("edges"):
    print("edge keys:", list(d["edges"][0].keys()))
    print("first edge:", d["edges"][0])
if d.get("nodes"):
    print("node keys:", list(d["nodes"][0].keys()))
    print("first node:", d["nodes"][0])
