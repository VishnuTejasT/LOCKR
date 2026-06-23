"""Parse raw/FASTA/mixed sequence input -- Python port of scanner.js parseBatchInput.

Returns (records, errors) where:
  records = [{"id": str, "sequence": str}]
  errors  = [{"line_num": int, "message": str}]
"""

from __future__ import annotations

_STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def parse_batch_input(text: str) -> tuple[list[dict], list[dict]]:
    lines = text.splitlines()
    records: list[dict] = []
    errors: list[dict] = []
    auto_idx = 0
    i = 0

    while i < len(lines):
        trimmed = lines[i].strip()
        if not trimmed:
            i += 1
            continue

        if trimmed.startswith(">"):
            header_line_num = i + 1
            id_ = trimmed[1:].strip()
            i += 1
            seq_parts: list[str] = []
            while i < len(lines) and not lines[i].strip().startswith(">") and lines[i].strip():
                seq_parts.append("".join(lines[i].split()).upper())
                i += 1
            seq = "".join(seq_parts)
            if not seq:
                label = id_ or "(empty header)"
                errors.append({"line_num": header_line_num,
                                "message": f'FASTA record "{label}": no sequence found'})
            else:
                bad = next((c for c in seq if c not in _STANDARD_AA), None)
                if bad:
                    errors.append({"line_num": header_line_num,
                                   "message": f'FASTA record "{id_}": non-standard amino acid {bad!r}'})
                else:
                    auto_idx += 1
                    records.append({"id": id_ or f"seq_{auto_idx}", "sequence": seq})
        else:
            line_num = i + 1
            seq = "".join(trimmed.split()).upper()
            bad = next((c for c in seq if c not in _STANDARD_AA), None)
            if bad:
                preview = trimmed[:40] + "…" if len(trimmed) > 40 else trimmed
                errors.append({"line_num": line_num,
                                "message": f"Line {line_num}: non-standard character {bad!r} in \"{preview}\""})
            else:
                auto_idx += 1
                records.append({"id": f"seq_{auto_idx}", "sequence": seq})
            i += 1

    return records, errors
