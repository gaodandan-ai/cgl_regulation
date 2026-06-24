import threading
import time
import urllib.request
import sys
import http.server
import socketserver

# Import run_server's classes
sys.path.append('.')
import run_server

def run_server_thread():
    server_address = ("", 8009)
    socketserver.TCPServer.allow_reuse_address = True
    try:
        server = run_server.ThreadingHTTPServer(server_address, run_server.CustomHTTPRequestHandler)
        print("Server started on port 8009 in thread")
        server.serve_forever()
    except Exception as e:
        print("Server thread exception:", e)

def main():
    t = threading.Thread(target=run_server_thread)
    t.daemon = True
    t.start()
    
    time.sleep(2.0) # wait for server
    
    print("Sending request to http://localhost:8009/index.html ...")
    try:
        with urllib.request.urlopen("http://localhost:8009/index.html", timeout=5) as resp:
            print("Response Code:", resp.getcode())
            print("Response Headers:", resp.info())
            print("Content (first 100 bytes):", resp.read(100))
    except Exception as e:
        print("Request failed with exception:", e)
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
