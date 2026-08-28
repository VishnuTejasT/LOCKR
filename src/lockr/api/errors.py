"""Shared error envelope: {"error": {"code", "field", "message"}}, per spec 6.5."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse


class ApiError(Exception):
    def __init__(self, code: str, message: str, field: str | None = None, status_code: int = 400):
        self.code = code
        self.field = field
        self.message = message
        self.status_code = status_code


def _envelope(code: str, message: str, field: str | None = None) -> dict:
    return {"error": {"code": code, "field": field, "message": message}}


# Pydantic v2 prefixes validator errors with "Value error, ", an implementation detail users shouldn't see.
def _clean_pydantic_message(msg: str) -> str:
    prefix = "Value error, "
    return msg[len(prefix):] if msg.startswith(prefix) else msg


def install_error_handlers(app: FastAPI) -> None:
    @app.exception_handler(ApiError)
    async def _handle_api_error(request: Request, exc: ApiError):
        return JSONResponse(status_code=exc.status_code,
                            content=_envelope(exc.code, exc.message, exc.field))

    @app.exception_handler(RequestValidationError)
    async def _handle_validation_error(request: Request, exc: RequestValidationError):
        first = exc.errors()[0]
        field = ".".join(str(p) for p in first["loc"] if p != "body")
        message = _clean_pydantic_message(first["msg"])
        return JSONResponse(status_code=400,
                            content=_envelope("VALIDATION_ERROR", message, field))

    @app.exception_handler(Exception)
    async def _handle_unexpected(request: Request, exc: Exception):
        return JSONResponse(status_code=500,
                            content=_envelope("INTERNAL_ERROR", f"{type(exc).__name__}: {exc}"))
