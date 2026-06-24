import os
import sys

# Add the parent directory to Python path so we can import run_server
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)
if PARENT_DIR not in sys.path:
    sys.path.insert(0, PARENT_DIR)

from run_server import CustomHTTPRequestHandler

class handler(CustomHTTPRequestHandler):
    """
    Vercel serverless function handler.
    Inherits the GET routes (/api/summarize, /api/test_ai, etc.) from run_server.
    """
    pass
