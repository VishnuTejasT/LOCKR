"""`lockr serve` -- runs the FastAPI app via uvicorn."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(prog="lockr")
    subparsers = parser.add_subparsers(dest="command", required=True)

    serve = subparsers.add_parser("serve", help="run the LOCKR API locally")
    serve.add_argument("--port", type=int, default=8420)

    args = parser.parse_args()
    if args.command == "serve":
        import uvicorn
        print(f"LOCKR API running at http://127.0.0.1:{args.port}")
        uvicorn.run("lockr.api.main:app", host="127.0.0.1", port=args.port)
