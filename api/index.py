"""Vercel entrypoint for the lightweight public data API."""

import os

os.environ.setdefault("CGL_PUBLIC_DEPLOYMENT", "true")

from fastapi import FastAPI

from api.public_api import configure_public_api


# Keep this explicit assignment: Vercel discovers Python entrypoints statically.
app = FastAPI(
    title="Cgl Regulation Public Data API",
    docs_url=None,
    redoc_url=None,
    openapi_url=None,
)
configure_public_api(app)
