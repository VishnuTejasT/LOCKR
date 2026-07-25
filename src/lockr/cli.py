"""`lockr serve` / `lockr scan` / `lockr fc`, CLI entry points."""

from __future__ import annotations

import argparse
import json
import sys
import textwrap


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


_SEQ_WRAP_WIDTH = 60


def _print_labeled_seq(label: str, seq: str) -> None:
    # Long sequences print as one giant line otherwise, wrap and indent
    # continuation lines under the label instead.
    indent = " " * len(label)
    lines = textwrap.wrap(seq, _SEQ_WRAP_WIDTH) or [""]
    print(f"{label}{lines[0]}")
    for line in lines[1:]:
        print(f"{indent}{line}")


def _print_scan_result(r, suggest: bool) -> None:
    _print_labeled_seq("Sequence:   ", r.sequence)
    if r.id and r.id != "seq_1":
        print(f"ID:         {r.id}")
    print(f"Length:     {r.length}")
    print(f"Liability:  {r.liability_score:.3f}  [{r.liability_band}]")
    print(f"Net charge: {r.net_charge:.1f}")
    kck_m = r.estimated_kck_nm * 1e-9
    print(f"K_CK:       {kck_m:.3e} M  ({r.estimated_kck_nm:.3g} nM)")
    for w in r.warnings:
        print(f"Warning:    {w}")
    if suggest and r.suggested_variants:
        v = r.suggested_variants[0]
        if v.sequence == r.sequence:
            # nothing was flagged, so the "suggestion" is just the input, skip the noise
            print()
            return
        subs = ", ".join(f"{s.from_}{s.position}{s.to}" for s in v.substitutions)
        print()
        _print_labeled_seq("Suggested:  ", v.sequence)
        if subs:
            print(f"Mutations:  {subs}")
        print(f"Liability:  {v.liability_score:.3f}  [{v.liability_band}]")
        print(f"K_CK:       {v.estimated_kck_nm * 1e-9:.3e} M  ({v.estimated_kck_nm:.3g} nM)")
    print()


def _cmd_fc(args) -> None:
    from lockr.api.errors import ApiError
    from lockr.api.routes.foldchange import foldchange
    from lockr.api.schemas.foldchange import FoldChangeRequest
    from pydantic import ValidationError

    from lockr.api.errors import _clean_pydantic_message

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
        _err(_clean_pydantic_message(first["msg"]))

    try:
        response = foldchange(request)
    except ApiError as e:
        _err(e.message)

    if args.json:
        print(json.dumps(response.model_dump(), indent=2))
    else:
        _print_fc_result(response)


def _print_fc_result(r) -> None:
    print(f"Fold-change:       {r.fold_change:.2f}x  [{r.quality}]")
    print(f"  {r.interpretation}")
    print(f"Best case (max):   {r.max_fold_change:.2f}x  (saturating target, this pull)")
    print(f"Key vs cage ratio: {r.dominance_ratio:.2f}  (lucKey / K_CK, how hard the key competes)")
    print(f"Limited by:        {r.regime.replace('_', '-')}")
    print()
    print(f"What this means: {r.verdict}")
    print("How to improve:")
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
    serve.add_argument("--port", type=int, default=8000,
                       help="local port to serve the web UI and API on (default: 8000)")

    # lockr scan
    scan = subparsers.add_parser("scan", help="scan a sequence for CK-binding liability")
    scan.add_argument("sequence", nargs="?",
                      help="amino-acid sequence to scan, standard residues only (e.g. LISDAELEAIFAEELDC). "
                           "Omit this and use --file instead for a batch of sequences.")
    scan.add_argument("--file", metavar="FILE",
                      help="path to a FASTA file, or a plain-text file with one sequence per line; "
                           "used instead of the sequence argument for scanning many sequences at once")
    scan.add_argument("--ph", type=float, default=7.4,
                      help="pH to compute net charge at (default: 7.4, physiological pH)")
    scan.add_argument("--window", metavar="START:END",
                      help="1-indexed residue range treated as the sensitive/graft region, e.g. 1:17 "
                           "(default: the whole sequence). Acidic residues outside this range are still "
                           "reported but don't count toward the liability score.")
    scan.add_argument("--preserve", metavar="POS,POS,...",
                      help="1-indexed residue positions to never suggest mutating, e.g. 1,2,11,12,15, "
                           "use this for residues that contact your target, or Suggest may mutate them")
    scan.add_argument("--policy", choices=["conservative", "neutralizing"], default="conservative",
                      help="how to fix flagged acidic residues in the suggested variant: "
                           "conservative = D->N/E->Q (keeps shape/H-bonding), "
                           "neutralizing = D->A/E->A (ECLIPSE's original fix) (default: conservative)")
    scan.add_argument("--suggest", action=argparse.BooleanOptionalAction, default=True,
                      help="also print a charge-optimized variant suggestion")
    scan.add_argument("--json", action="store_true", help="output raw JSON instead of a plain-text summary")

    # lockr fc
    fc = subparsers.add_parser("fc", help="compute fold-change for a LOCKR sensor")
    fc.add_argument("--k-ck", type=float, required=True, dest="k_ck", metavar="FLOAT",
                    help="cage-key dissociation constant, in nM, how tightly lucKey binds the open "
                         "cage; lower = tighter = brighter signal (typical 1-100 nM)")
    fc.add_argument("--k-open", type=float, required=True, dest="k_open", metavar="FLOAT",
                    help="basal latch-opening equilibrium, dimensionless, chance the cage pops open "
                         "with no target present (typical ~0.001 = shut 99.9%% of the time)")
    fc.add_argument("--pull", type=float, required=True,
                    help="allosteric coupling strength, dimensionless, how much target binding "
                         "increases the cage's odds of opening (e.g. 10 = ~10x more likely; fit this "
                         "from a wet-lab luminescence titration, typical 10-30)")
    fc.add_argument("--luckey", type=float, required=True,
                    help="reporter key concentration used in the assay, in nM (typical 100-1000 nM)")
    fc.add_argument("--k-target", type=float, default=None, dest="k_target", metavar="FLOAT",
                    help="binder-target affinity, in nM. Must be paired with --target, or omitted "
                         "entirely to assume the target is fully saturating (best case)")
    fc.add_argument("--target", type=float, default=None, metavar="FLOAT",
                    help="target concentration actually present in the sample, in nM. Must be paired "
                         "with --k-target, or omitted entirely to assume saturating target")
    fc.add_argument("--json", action="store_true", help="output raw JSON instead of a plain-text summary")

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
