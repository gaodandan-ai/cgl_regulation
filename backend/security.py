"""Security helpers shared by the API and outbound AI integrations."""

from __future__ import annotations

import ipaddress
import os
import socket
from urllib.parse import urlparse


PUBLIC_DEPLOYMENT = (
    os.environ.get("VERCEL") == "1"
    or "AWS_LAMBDA_FUNCTION_NAME" in os.environ
    or os.environ.get("CGL_PUBLIC_DEPLOYMENT", "").lower() in {"1", "true", "yes"}
)

_PROVIDER_HOSTS = {
    "openai": {"api.openai.com"},
    "deepseek": {"api.deepseek.com"},
    "qwen": {"dashscope.aliyuncs.com"},
    "kimi": {"api.moonshot.cn"},
    "zhipu": {"open.bigmodel.cn"},
    "google": {"generativelanguage.googleapis.com"},
}


def _is_public_address(address: str) -> bool:
    ip = ipaddress.ip_address(address)
    return not (
        ip.is_private
        or ip.is_loopback
        or ip.is_link_local
        or ip.is_multicast
        or ip.is_reserved
        or ip.is_unspecified
    )


def validate_outbound_url(url: str, provider: str = "custom") -> str:
    """Validate an outbound AI endpoint and reject server-side network pivots."""
    parsed = urlparse((url or "").strip())
    if parsed.scheme not in {"http", "https"} or not parsed.hostname:
        raise ValueError("AI Base URL must be an absolute HTTP(S) URL")
    if parsed.username or parsed.password:
        raise ValueError("Credentials are not allowed in AI Base URLs")

    hostname = parsed.hostname.rstrip(".").lower()
    provider = (provider or "custom").lower()

    if PUBLIC_DEPLOYMENT:
        if parsed.scheme != "https":
            raise ValueError("Public deployments require HTTPS AI endpoints")
        if provider == "ollama":
            raise ValueError("Ollama is available only in the desktop/local application")

        allowed = set(_PROVIDER_HOSTS.get(provider, set()))
        if provider == "custom":
            allowed.update(
                item.strip().lower()
                for item in os.environ.get("CGL_AI_ALLOWED_HOSTS", "").split(",")
                if item.strip()
            )
        if not allowed or hostname not in allowed:
            raise ValueError("This AI endpoint is not allowed by the public deployment")

        try:
            addresses = {item[4][0] for item in socket.getaddrinfo(hostname, parsed.port or 443)}
        except socket.gaierror as exc:
            raise ValueError("AI endpoint hostname could not be resolved") from exc
        if not addresses or any(not _is_public_address(address) for address in addresses):
            raise ValueError("Private or reserved network addresses are not allowed")

    return url.rstrip("/")
