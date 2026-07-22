import sys
from pathlib import Path

import pytest


BACKEND = Path(__file__).resolve().parents[1] / "backend"
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

import security


def test_public_custom_endpoint_requires_allowlist(monkeypatch):
    monkeypatch.setattr(security, "PUBLIC_DEPLOYMENT", True)
    monkeypatch.delenv("CGL_AI_ALLOWED_HOSTS", raising=False)
    with pytest.raises(ValueError, match="not allowed"):
        security.validate_outbound_url("https://example.com/v1", "custom")


def test_public_endpoint_rejects_loopback(monkeypatch):
    monkeypatch.setattr(security, "PUBLIC_DEPLOYMENT", True)
    monkeypatch.setenv("CGL_AI_ALLOWED_HOSTS", "localhost")
    with pytest.raises(ValueError):
        security.validate_outbound_url("https://localhost/v1", "custom")


def test_public_known_provider_is_allowed(monkeypatch):
    monkeypatch.setattr(security, "PUBLIC_DEPLOYMENT", True)
    monkeypatch.setattr(
        security.socket,
        "getaddrinfo",
        lambda *args, **kwargs: [(None, None, None, None, ("104.18.7.192", 443))],
    )
    assert security.validate_outbound_url("https://api.openai.com/v1", "openai") == "https://api.openai.com/v1"
