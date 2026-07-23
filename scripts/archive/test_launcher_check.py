import sys, os, time, urllib.request, threading
sys.path.insert(0, ".")

def poll_health():
    for i in range(20):
        time.sleep(0.5)
        try:
            r = urllib.request.urlopen("http://127.0.0.1:8000/api/check-update")
            print("API RESPONSE:", r.status, r.read().decode('utf-8'))
            os._exit(0)
        except Exception as e:
            pass
    print("Health check timed out")
    os._exit(1)

t = threading.Thread(target=poll_health, daemon=True)
t.start()

import launcher
launcher.main()
