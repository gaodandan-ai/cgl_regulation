#!/usr/bin/env python3
"""Post-deployment smoke checks for the public Cgl Regulation site and API."""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.parse
import urllib.request


def fetch(base_url: str, path: str) -> tuple[int, dict, dict]:
    url = urllib.parse.urljoin(base_url.rstrip("/") + "/", path.lstrip("/"))
    request = urllib.request.Request(url, headers={"User-Agent": "CglDeploymentSmoke/1"})
    with urllib.request.urlopen(request, timeout=20) as response:
        body = response.read()
        payload = json.loads(body.decode("utf-8")) if "json" in response.headers.get_content_type() else {}
        return response.status, dict(response.headers), payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", default="https://cgl-regulation.vercel.app")
    parser.add_argument("--expected-version", default="")
    args = parser.parse_args()

    checks = [
        ("health", "/api/health"),
        ("condition regulation", "/api/condition-regulation/runs"),
        ("target priorities", "/api/intervention-targets?limit=1"),
    ]
    failures = []
    for label, path in checks:
        try:
            failures_before = len(failures)
            status, headers, payload = fetch(args.base_url, path)
            if status != 200:
                failures.append(f"{label}: HTTP {status}")
                continue
            if not payload:
                failures.append(f"{label}: expected JSON but received a non-API response")
                continue
            if label == "health":
                if payload.get("app") != "cgl-regulation" or payload.get("status") != "ok":
                    failures.append("health: invalid application identity")
                if args.expected_version and payload.get("version") != args.expected_version:
                    failures.append(
                        f"health: version {payload.get('version')!r} != {args.expected_version!r}"
                    )
                if headers.get("X-Content-Type-Options", "").lower() != "nosniff":
                    failures.append("health: security headers missing")
            elif label == "condition regulation" and "runs" not in payload:
                failures.append("condition regulation: invalid response schema")
            elif label == "target priorities" and "targets" not in payload:
                failures.append("target priorities: invalid response schema")
            if len(failures) == failures_before:
                print(f"PASS {label}: HTTP {status}")
        except (OSError, ValueError, urllib.error.HTTPError) as exc:
            failures.append(f"{label}: {exc}")

    if failures:
        for failure in failures:
            print(f"FAIL {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
