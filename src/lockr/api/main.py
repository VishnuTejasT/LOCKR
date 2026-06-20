from fastapi import FastAPI

from .errors import install_error_handlers
from .routes import assembly, foldchange, scan, sweep

app = FastAPI(title="LOCKR API")
install_error_handlers(app)
app.include_router(scan.router)
app.include_router(foldchange.router)
app.include_router(sweep.router)
app.include_router(assembly.router)
