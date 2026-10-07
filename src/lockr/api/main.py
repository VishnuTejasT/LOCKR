from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from starlette.responses import Response

from .errors import install_error_handlers
from .routes import (assembly, dosing, foldchange, graft, inverse, kopen, landscape,
                     scan, sweep, temperature_window)

WEB_DIR = Path(__file__).resolve().parent.parent / "web"


class _RevalidatingStatic(StaticFiles):
    """Static files the browser must revalidate."""
    # no-cache still 304s, but stops a stale landscape.js running against a new index.html

    def file_response(self, *args, **kwargs) -> Response:
        response = super().file_response(*args, **kwargs)
        response.headers["Cache-Control"] = "no-cache"
        return response

app = FastAPI(title="LOCKR API")
install_error_handlers(app)
app.include_router(scan.router)
app.include_router(foldchange.router)
app.include_router(sweep.router)
app.include_router(assembly.router)
app.include_router(kopen.router)
app.include_router(graft.router)
app.include_router(dosing.router)
app.include_router(inverse.router)
app.include_router(landscape.router)
app.include_router(temperature_window.router)

app.mount("/static", _RevalidatingStatic(directory=WEB_DIR / "static"), name="static")


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(WEB_DIR / "index.html")
