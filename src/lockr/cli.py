"""`lockr serve` / `lockr scan` / `lockr fc` -- CLI entry points."""

from __future__ import annotations

import argparse
import json
import sys


def _err(msg: str) -> None:
    print(f"error: {msg}", file=sys.stderr)
    sys.exit(1)


def _parse_window(s: str) -> tuple[int, int]:
    try:
        left, right = s.split(":")
        return int(left), int(right)
    except (ValueError, AttributeError):
        _err(f"--window must be START:END (e.g. 1:17), got: {s!r}")


def _parse_preserve(s: str) -> list[int]:
    try:
        return [int(p) for p in s.split(",")]
    except ValueError:
        _err(f"--preserve must be comma-separated integers, got: {s!r}")


def _cmd_scan(args) -> None:
    from lockr.api.errors import ApiError
    from lockr.api.routes.scan import _scan_one, _validate_window_and_preserve
    from lockr.api.schemas.common import validate_sequence
    from lockr.api.schemas.scan import Window
    from lockr.engine.parse_batch import parse_batch_input

    if args.sequence and args.file:
        _err("provide either a sequence argument or --file, not both")
    if not args.sequence and not args.file:
        _err("provide a sequence argument or --file")

    if args.file:
        try:
            text = open(args.file).read()
        except OSError as e:
            _err(str(e))
        records, errors = parse_batch_input(text)
        if errors:
            for e in errors:
                print(f"parse error (line {e['line_num']}): {e['message']}", file=sys.stderr)
            sys.exit(1)
        if not records:
            _err("no valid sequences found in file")
    else:
        records = [{"id": "seq_1", "sequence": args.sequence}]

    preserve = _parse_preserve(args.preserve) if args.preserve else []
    results = []

    for rec in records:
        try:
            seq = validate_sequence(rec["sequence"])
        except ValueError as e:
            _err(str(e))

        length = len(seq)
        start, end = _parse_window(args.window) if args.window else (1, length)
        window = Window(start=start, end=end)

        try:
            _validate_window_and_preserve(window, preserve, length,
                                          window_field="sensitive_window",
                                          preserve_field="preserve_positions")
        except ApiError as e:
            _err(e.message)

        start_c, end_c = window.clamped(length)
        result = _scan_one(seq, start_c, end_c, args.ph, args.policy, preserve)
        result.id = rec["id"]
        results.append(result)

    if args.json:
        data = {"results": [r.model_dump(by_alias=True) for r in results]}
        print(json.dumps(data, indent=2))
    else:
        for r in results:
            _print_scan_result(r, suggest=args.suggest)


def _print_scan_result(r, suggest: bool) -> None:
    print(f"Sequence:   {r.sequence}")
    if r.id and r.id != "seq_1":
        print(f"ID:         {r.id}")
    print(f"Length:     {r.length}")
    print(f"Liability:  {r.liability_score:.3f}  [{r.liability_band}]")
    print(f"Net charge: {r.net_charge:.1f}")
    kck_m = r.estimated_kck_nm * 1e-9
    print(f"K_CK:       {kck_m:.3e} M")
    for w in r.warnings:
        print(f"Warning:    {w}")
    if suggest and r.suggested_variants:
        v = r.suggested_variants[0]
        subs = ", ".join(f"{s.from_}{s.position}{s.to}" for s in v.substitutions)
        print(f"\nSuggested:  {v.sequence}")
        if subs:
            print(f"Mutations:  {subs}")
        print(f"Liability:  {v.liability_score:.3f}  [{v.liability_band}]")
        print(f"K_CK:       {v.estimated_kck_nm * 1e-9:.3e} M")
    print()


def _cmd_fc(args) -> None:
    from lockr.api.errors import ApiError
    from lockr.api.routes.foldchange import foldchange
    from lockr.api.schemas.foldchange import FoldChangeRequest
    from pydantic import ValidationError

    try:
        request = FoldChangeRequest(
            k_ck=args.k_ck,
            k_open=args.k_open,
            pull=args.pull,
            luckey=args.luckey,
            k_target=args.k_target,
            target_conc=args.target,
        )
    except ValidationError as e:
        first = e.errors()[0]
        _err(first["msg"])

    try:
        response = foldchange(request)
    except ApiError as e:
        _err(e.message)

    if args.json:
        print(json.dumps(response.model_dump(), indent=2))
    else:
        _print_fc_result(response)


def _print_fc_result(r) -> None:
    print(f"Fold-change:     {r.fold_change:.2f}x")
    print(f"Dominance ratio: {r.dominance_ratio:.2f}")
    print(f"Fraction of DR:  {r.fraction_of_dominance_ratio:.4f}")
    print(f"Regime:          {r.regime.replace('_', '-')}")
    print()
    print(f"Verdict: {r.verdict}")
    print("Recommendations:")
    for rec in r.recommendations:
        print(f"  - {rec}")
    if r.warnings:
        print()
        for w in r.warnings:
            print(f"Warning: {w}")


def main() -> None:
    parser = argparse.ArgumentParser(prog="lockr")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # lockr serve
    serve = subparsers.add_parser("serve", help="run the LOCKR API locally")
    serve.add_argument("--port", type=int, default=8000)

    # lockr scan
    scan = subparsers.add_parser("scan", help="scan a sequence for CK-binding liability")
    scan.add_argument("sequence", nargs="?", help="amino-acid sequence")
    scan.add_argument("--file", metavar="FILE", help="FASTA or raw-sequence file")
    scan.add_argument("--ph", type=float, default=7.4)
    scan.add_argument("--window", metavar="START:END")
    scan.add_argument("--preserve", metavar="POS,POS,...")
    scan.add_argument("--policy", choices=["conservative", "neutralizing"], default="conservative")
    scan.add_argument("--suggest", action=argparse.BooleanOptionalAction, default=True)
    scan.add_argument("--json", action="store_true", help="output raw JSON")

    # lockr fc
    fc = subparsers.add_parser("fc", help="compute fold-change for a LOCKR sensor")
    fc.add_argument("--k-ck", type=float, required=True, dest="k_ck", metavar="FLOAT")
    fc.add_argument("--k-open", type=float, required=True, dest="k_open", metavar="FLOAT")
    fc.add_argument("--pull", type=float, required=True)
    fc.add_argument("--luckey", type=float, required=True)
    fc.add_argument("--k-target", type=float, default=None, dest="k_target", metavar="FLOAT")
    fc.add_argument("--target", type=float, default=None, metavar="FLOAT")
    fc.add_argument("--json", action="store_true", help="output raw JSON")

    args = parser.parse_args()

    if args.command == "serve":
        import uvicorn
        print(f"LOCKR API running at http://127.0.0.1:{args.port}")
        uvicorn.run("lockr.api.main:app", host="127.0.0.1", port=args.port)
    elif args.command == "scan":
        _cmd_scan(args)
    elif args.command == "fc":
        _cmd_fc(args)


if __name__ == "__main__":
    main()
