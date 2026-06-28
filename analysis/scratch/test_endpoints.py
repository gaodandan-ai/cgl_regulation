# test_endpoints.py
import sys
import os

# Add root directory to python path
sys.path.insert(0, os.getcwd())

from run_server import CustomHTTPRequestHandler

class DummyRequest:
    def __init__(self):
        self.headers = {}

# Instantiate request handler
handler_instance = CustomHTTPRequestHandler.__new__(CustomHTTPRequestHandler)

print("--- Testing perform_protein_domain_analysis ---")
res_domain = handler_instance.perform_protein_domain_analysis("cg0350", "DummyKey")
print("WhiB4 domain result:")
print(res_domain.get("summary")[:200] + "...")

print("\n--- Testing perform_binding_site_analysis ---")
res_binding = handler_instance.perform_binding_site_analysis("cg0350", "DummyKey")
print("WhiB4 binding site result:")
print(res_binding.get("summary")[:200] + "...")
print(f"Total sites: {res_binding.get('total_sites')}")
