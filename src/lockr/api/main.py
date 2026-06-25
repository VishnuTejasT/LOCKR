from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .errors import install_error_handlers
from .routes import assembly, foldchange, scan, sweep

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

app = FastAPI(title="LOCKR API")
install_error_handlers(app)
app.include_router(scan.router)
app.include_router(foldchange.router)
app.include_router(sweep.router)
app.include_router(assembly.router)

app.mount("/static", StaticFiles(directory=WEB_DIR / "static"), name="static")


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(WEB_DIR / "index.html")
