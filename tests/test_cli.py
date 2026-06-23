"""CLI integration tests -- run lockr as a subprocess so we catch import and exit-code issues."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import os


def _run(*args: str, input: str | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "lockr.cli", *args],
        capture_output=True, text=True, input=input,
    )


ECLIPSE_HIGH = "LISDAELEAIFAEELDC"
ECLIPSE_LOW  = "LISAAALAAIFAAALAC"

# ── lockr scan ────────────────────────────────────────────────────────────────

def test_scan_eclipse_high_liability():
    r = _run("scan", ECLIPSE_HIGH)
    assert r.returncode == 0
    assert "High" in r.stdout
    # K_CK ~3.32e-05 M -- check order of magnitude and first significant digit
    assert "3.3" in r.stdout or "3.32" in r.stdout


def test_scan_eclipse_high_kck_value():
    r = _run("scan", ECLIPSE_HIGH, "--json")
    assert r.returncode == 0
    data = json.loads(r.stdout)
    kck_nm = data["results"][0]["estimated_kck_nm"]
    kck_m = kck_nm * 1e-9
    assert abs(kck_m - 3.32e-5) / 3.32e-5 < 0.01, f"K_CK = {kck_m:.3e}, expected ~3.32e-5"


def test_scan_eclipse_low_liability():
    r = _run("scan", ECLIPSE_LOW)
    assert r.returncode == 0
    assert "Low" in r.stdout


def test_scan_eclipse_low_kck_value():
    r = _run("scan", ECLIPSE_LOW, "--json")
    assert r.returncode == 0
    data = json.loads(r.stdout)
    kck_nm = data["results"][0]["estimated_kck_nm"]
    kck_m = kck_nm * 1e-9
    assert abs(kck_m - 1.00e-8) / 1.00e-8 < 0.01, f"K_CK = {kck_m:.3e}, expected ~1.00e-8"


def test_scan_json_is_valid_and_matches_human_readable():
    r_text = _run("scan", ECLIPSE_HIGH)
    r_json = _run("scan", ECLIPSE_HIGH, "--json")
    assert r_text.returncode == 0
    assert r_json.returncode == 0
    data = json.loads(r_json.stdout)
    result = data["results"][0]
    # Both outputs reference the same K_CK order of magnitude
    kck_m = result["estimated_kck_nm"] * 1e-9
    assert f"{kck_m:.3e}" in r_text.stdout


def test_scan_bad_sequence_exits_nonzero():
    r = _run("scan", "LISDAELXAIFAEELDC")
    assert r.returncode != 0
    assert "non-standard" in r.stderr.lower() or "error" in r.stderr.lower()
    assert not r.stdout.strip()


def test_scan_file_fasta():
    fasta = ">high\nLISDAELEAIFAEELDC\n>low\nLISAALAAAIFAAALAC\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        f.write(fasta)
        path = f.name
    try:
        r = _run("scan", "--file", path)
        assert r.returncode == 0
        assert "High" in r.stdout
    finally:
        os.unlink(path)


def test_scan_file_raw():
    raw = "LISDAELEAIFAEELDC\nLISAAALAAAIFAAALAC\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write(raw)
        path = f.name
    try:
        r = _run("scan", "--file", path)
        assert r.returncode == 0
    finally:
        os.unlink(path)


def test_scan_no_suggest_flag():
    r = _run("scan", ECLIPSE_HIGH, "--no-suggest")
    assert r.returncode == 0
    assert "Suggested" not in r.stdout


# ── lockr fc ──────────────────────────────────────────────────────────────────

def test_fc_eclipse_foldchange():
    r = _run("fc", "--k-ck", "10", "--k-open", "0.001", "--pull", "10", "--luckey", "500")
    assert r.returncode == 0
    # fold-change ~11x
    assert "11." in r.stdout or "11x" in r.stdout
    assert "key-limited" in r.stdout


def test_fc_json_valid_and_matches_human_readable():
    r_text = _run("fc", "--k-ck", "10", "--k-open", "0.001", "--pull", "10", "--luckey", "500")
    r_json = _run("fc", "--k-ck", "10", "--k-open", "0.001", "--pull", "10", "--luckey", "500", "--json")
    assert r_text.returncode == 0
    assert r_json.returncode == 0
    data = json.loads(r_json.stdout)
    fc = data["fold_change"]
    # human-readable output should contain the same fold-change to 2 dp
    assert f"{fc:.2f}" in r_text.stdout


def test_fc_json_regime_key_limited():
    r = _run("fc", "--k-ck", "10", "--k-open", "0.001", "--pull", "10", "--luckey", "500", "--json")
    data = json.loads(r.stdout)
    assert data["regime"] in ("key_limited", "key-limited")


def test_fc_negative_kck_exits_nonzero():
    r = _run("fc", "--k-ck", "-5", "--k-open", "0.001", "--pull", "10", "--luckey", "500")
    assert r.returncode != 0
    assert r.stderr.strip()
    assert not r.stdout.strip()


def test_fc_zero_kopen_exits_nonzero():
    r = _run("fc", "--k-ck", "10", "--k-open", "0", "--pull", "10", "--luckey", "500")
    assert r.returncode != 0
    assert r.stderr.strip()
