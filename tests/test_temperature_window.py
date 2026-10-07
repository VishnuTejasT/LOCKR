"""Temperature operating window: the fold-change convention, window logic, caps and refusals.

Every number in here is deliberately synthetic and round. The shipped ECLIPSE preset is
placeholder data and is only touched by the one test that checks the preset still runs.
"""

import copy
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient

from lockr.api.main import app
from lockr.api.schemas.temperature_window import ECLIPSE_PRESET
from lockr.engine import temperature_window as tw
from lockr.engine import thermo
from lockr.engine.models import SensorParams
from lockr.engine.temperature import R_KCAL, ThermoSpec, dG_at

client = TestClient(app)

WEB_JS = Path(__file__).resolve().parents[1] / "src/lockr/web/static/js"

# Synthetic construct: round numbers, no relation to any real binder.
SYNTH = {
    "name": "synthetic cage",
    "dG_latch_cage": -4.0,
    "dH_latch_cage": -10.0,
    "dG_key_cage": -11.0,
    "dH_key_cage": -6.0,
    "target_conc": 1e-6,
    "binder": {
        "name": "synthetic binder",
        "dG_target_binder": -12.0,
        "dG_target_binder_sd": 0.0,
        "dH_target_binder": -8.0,
    },
}

NO_SPREAD = {k: {"mean": v["mean"], "sd": 0.0} for k, v in tw.DEFAULT_DH_PRIORS.items()}


def synth(**over):
    system = copy.deepcopy(SYNTH)
    binder = over.pop("binder", {})
    system.update(over)
    system["binder"].update(binder)
    return system


def run(systems, **over):
    kwargs = dict(t_min_c=4.0, t_max_c=45.0, t_step_c=1.0, t_ref_c=25.0,
                  n_samples=200, seed=11, dH_priors=NO_SPREAD)
    kwargs.update(over)
    return tw.scan(systems, **kwargs)


# --- 1. the median at T_ref is the existing fold-change function, not a new model ------

def test_median_at_t_ref_matches_existing_fold_change():
    out = run([synth()], t_min_c=25.0, t_max_c=25.0 + 1.0, t_step_c=1.0)
    result = out["systems"][0]

    T = tw.c_to_k(25.0)
    RT = R_KCAL * T
    params = SensorParams(K_open=np.exp(SYNTH["dG_latch_cage"] / RT),
                          K_CK=np.exp(SYNTH["dG_key_cage"] / RT),
                          lucKey=SensorParams().lucKey, RT=RT)
    expected = thermo.fold_change(SYNTH["target_conc"],
                                  np.exp(SYNTH["binder"]["dG_target_binder"] / RT),
                                  tw.DEFAULT_PULL, params)

    assert result["median_fc"][0] == pytest.approx(expected, rel=1e-12)


def test_vant_hoff_shift_matches_temperature_dG_at_with_zero_dCp():
    """dG_at_T is temperature.dG_at with dCp = 0, so the two must not drift apart."""
    T_ref = tw.c_to_k(25.0)
    spec = ThermoSpec(dG_ref=-4.0, dH_ref=-10.0, dCp=0.0, T_ref=T_ref)
    for t_c in (4.0, 25.0, 45.0):
        T = tw.c_to_k(t_c)
        assert tw.dG_at_T(-4.0, -10.0, T, T_ref) == pytest.approx(dG_at(spec, T), abs=1e-12)


# --- 2 + 11. the boundary convention, asserted so a future change fails loudly ---------

def test_fold_change_is_one_at_zero_target_and_undefined_at_zero_key():
    """Convention this module relies on, checked against the engine rather than assumed.

    Zero target -> exactly 1.0 (baseline, "no change"), NOT 0. Zero Key -> 0/0, which is
    undefined and must stay undefined; it is never silently 0 or 1. p_works counts a nan
    as a failure. If this test fails, the convention moved and this module needs updating.
    """
    assert thermo.fold_change(0.0, 1e-9, 10.0) == 1.0

    no_key = SensorParams(K_open=1e-3, K_CK=1e-8, lucKey=0.0)
    assert thermo._f_open(1e-3, no_key) == 0.0
    with pytest.raises(ZeroDivisionError):
        thermo.fold_change(1e-6, 1e-9, 10.0, no_key)

    # the vectorised path reports that same undefined case as nan, never as 0 or 1
    fc = tw._fold_change(np.array([[1e-3]]), np.array([[1e-8]]), 0.0, 10.0,
                         np.array([[1e-9]]), 1e-6, np.array([[298.15]]))
    assert np.isnan(fc).all()

    # and a nan never counts as working
    window = tw.operating_window(np.array([0.0, 1.0]), np.array([np.nan, np.nan]), 0.8)
    assert window["window"] is None


def test_zero_variance_is_fully_deterministic():
    result = run([synth()], n_samples=300)["systems"][0]
    assert set(result["p_works"]) <= {0.0, 1.0}
    # with no spread every sample is the same curve, so the band has no width
    assert result["p10_fc"] == result["median_fc"] == result["p90_fc"]


# --- 3. determinism -------------------------------------------------------------------

def test_same_seed_reproduces_and_different_seeds_differ():
    spread = {"binder": {"dG_target_binder_sd": 1.5}}
    a = run([synth(**spread)], seed=5, dH_priors=None)["systems"][0]
    b = run([synth(**spread)], seed=5, dH_priors=None)["systems"][0]
    c = run([synth(**spread)], seed=6, dH_priors=None)["systems"][0]
    assert a["median_fc"] == b["median_fc"]
    assert a["median_fc"] != c["median_fc"]


def test_missing_seed_is_generated_and_returned():
    out = tw.scan([synth()], n_samples=200)
    assert isinstance(out["seed"], int)
    # the returned seed reproduces the run it came from
    assert tw.scan([synth()], n_samples=200, seed=out["seed"])["systems"][0]["median_fc"] \
        == out["systems"][0]["median_fc"]


# --- 4. window logic, on hand-built p_works arrays ------------------------------------

TEMPS = np.arange(10.0)


def test_window_contiguous():
    p = np.array([0, 0, .9, .95, .9, 0, 0, 0, 0, 0], dtype=float)
    out = tw.operating_window(TEMPS, p, 0.8)
    assert out["window"] == {"t_low_c": 2.0, "t_high_c": 4.0}
    assert out["other_segments"] == []
    assert out["window_touches_grid_edge"] is False
    assert out["reason"] is None


def test_window_multiple_segments_returns_longest_and_lists_rest():
    p = np.array([.9, 0, .9, .9, .9, 0, .85, .85, 0, .9], dtype=float)
    out = tw.operating_window(TEMPS, p, 0.8)
    assert out["window"] == {"t_low_c": 2.0, "t_high_c": 4.0}
    assert out["other_segments"] == [{"t_low_c": 0.0, "t_high_c": 0.0},
                                     {"t_low_c": 6.0, "t_high_c": 7.0},
                                     {"t_low_c": 9.0, "t_high_c": 9.0}]


def test_window_none_returns_null_and_a_reason():
    out = tw.operating_window(TEMPS, np.full(10, 0.4), 0.8)
    assert out["window"] is None
    assert out["other_segments"] == []
    assert "0.8" in out["reason"] and "0.40" in out["reason"]


def test_window_touching_grid_edge_sets_the_flag():
    assert tw.operating_window(TEMPS, np.ones(10), 0.8)["window_touches_grid_edge"] is True
    low = np.array([.9, .9, 0, 0, 0, 0, 0, 0, 0, 0], dtype=float)
    assert tw.operating_window(TEMPS, low, 0.8)["window_touches_grid_edge"] is True
    high = np.array([0, 0, 0, 0, 0, 0, 0, 0, .9, .9], dtype=float)
    assert tw.operating_window(TEMPS, high, 0.8)["window_touches_grid_edge"] is True


# --- 5. stability cap ------------------------------------------------------------------

def test_binder_tm_below_system_cap_sets_the_cap():
    result = run([synth(t_max_stability_c=40.0, binder={"tm_c": 30.0})])["systems"][0]
    assert result["stability_cap_c"] == 30.0
    assert result["stability_cap_source"] == "binder"


def test_system_cap_below_binder_tm_sets_the_cap():
    result = run([synth(t_max_stability_c=30.0, binder={"tm_c": 40.0})])["systems"][0]
    assert result["stability_cap_c"] == 30.0
    assert result["stability_cap_source"] == "system"


def test_p_works_is_exactly_zero_above_the_cap():
    out = run([synth(binder={"tm_c": 30.0})], fc_threshold=1.0)
    temps = np.array(out["temperatures_c"])
    p_works = np.array(out["systems"][0]["p_works"])
    assert (p_works[temps > 30.0] == 0.0).all()
    # and a threshold of 1.0 is met everywhere below the cap, so the cap is the only zero
    assert (p_works[temps <= 30.0] > 0.0).all()


def test_t_opt_never_lands_above_the_cap():
    uncapped = run([synth()])["systems"][0]
    capped = run([synth(binder={"tm_c": 20.0})])["systems"][0]
    assert capped["t_opt_c"] <= 20.0
    assert capped["stability_cap_source"] == "binder"
    assert uncapped["stability_cap_c"] is None


# --- 6. luciferase activity ------------------------------------------------------------

def test_luciferase_factor_of_one_changes_nothing():
    plain = run([synth()])["systems"][0]
    flat = run([synth()], luciferase_activity=[[0.0, 1.0], [60.0, 1.0]])["systems"][0]
    assert flat["median_fc"] == pytest.approx(plain["median_fc"])
    assert flat["p_works"] == plain["p_works"]


def test_luciferase_factor_of_zero_zeroes_the_signal_there():
    out = run([synth()], fc_threshold=1.0,
              luciferase_activity=[[0.0, 1.0], [30.0, 1.0], [30.001, 0.0], [60.0, 0.0]])
    temps = np.array(out["temperatures_c"])
    result = out["systems"][0]
    median = np.array([v if v is not None else np.nan for v in result["median_fc"]])
    assert (median[temps > 31.0] == 0.0).all()
    assert (median[temps <= 30.0] > 0.0).all()
    assert (np.array(result["p_works"])[temps > 31.0] == 0.0).all()


# --- 7. two binders on the same cage ---------------------------------------------------

def test_different_binder_dG_and_dH_give_different_curves_and_t_opt():
    tight = synth(name="shared cage", binder={"name": "tight",
                                              "dG_target_binder": -13.0,
                                              "dH_target_binder": -20.0})
    loose = synth(name="shared cage", binder={"name": "loose",
                                              "dG_target_binder": -9.0,
                                              "dH_target_binder": 5.0})
    out = run([tight, loose])
    a, b = out["systems"]

    # same cage, same Key, different binder: separate entries that overlay
    assert a["name"] == b["name"] and a["binder_name"] != b["binder_name"]
    assert a["median_fc"] != b["median_fc"]
    assert a["t_opt_c"] != b["t_opt_c"]


def test_binder_is_inert_without_a_target_concentration():
    """theta is pinned at 1 with no target, so Kd cancels. The result must say so."""
    a = run([synth(target_conc=None, binder={"dG_target_binder": -13.0})])
    b = run([synth(target_conc=None, binder={"dG_target_binder": -6.0})])
    assert a["systems"][0]["median_fc"] == b["systems"][0]["median_fc"]
    assert a["systems"][0]["binder_affects_result"] is False
    assert any("does NOT affect this result" in n for n in a["notes"])


# --- 8. graft shift --------------------------------------------------------------------

def test_graft_shift_moves_effective_dG_by_exactly_that_much_and_raises_leak():
    shift = 1.3
    plain = run([synth()])["systems"][0]
    grafted = run([synth(binder={"graft_dG_latch_cage_shift": shift})])["systems"][0]

    assert grafted["dG_latch_cage_raw"] == SYNTH["dG_latch_cage"]
    assert grafted["dG_latch_cage_effective"] == pytest.approx(
        SYNTH["dG_latch_cage"] + shift, abs=1e-12)
    assert grafted["graft_dG_latch_cage_shift"] == shift

    # A less negative latch-cage dG is a looser cage: more background glow (leak) and so
    # less headroom left for the target to buy, i.e. a lower fold change.
    RT = R_KCAL * tw.c_to_k(25.0)
    leak = lambda dG: thermo._f_open(
        np.exp(dG / RT), SensorParams(K_open=np.exp(dG / RT),
                                      K_CK=np.exp(SYNTH["dG_key_cage"] / RT),
                                      lucKey=SensorParams().lucKey, RT=RT))
    assert leak(grafted["dG_latch_cage_effective"]) > leak(plain["dG_latch_cage_effective"])
    assert grafted["peak_median_fc"] < plain["peak_median_fc"]


# --- 9. Kd and dG are the same input ---------------------------------------------------

def test_kd_nM_converts_to_the_same_dG_as_a_direct_dG():
    T_ref = tw.c_to_k(25.0)
    kd_nM = 4.0
    dG = tw.dG_from_kd_nM(kd_nM, T_ref)
    assert dG == pytest.approx(R_KCAL * T_ref * np.log(kd_nM * 1e-9), rel=1e-12)

    body = _request(binder_over={"dG_target_binder": None, "kd_nM": kd_nM})
    via_kd = _post(body)["systems"][0]
    via_dG = _post(_request(binder_over={"dG_target_binder": dG}))["systems"][0]
    assert via_kd["dG_target_binder"] == pytest.approx(via_dG["dG_target_binder"], rel=1e-12)
    assert via_kd["median_fc"] == via_dG["median_fc"]


# --- API plumbing ----------------------------------------------------------------------

def _request(binder_over=None, **over):
    body = {
        "systems": [{
            "name": "synthetic cage",
            "dG_latch_cage": -4.0,
            "dH_latch_cage": -10.0,
            "dG_key_cage": -11.0,
            "dH_key_cage": -6.0,
            "target_conc": 1000.0,      # nM at the API boundary
            "binder": {
                "name": "synthetic binder",
                "dG_target_binder": -12.0,
                "dG_target_binder_sd": 0.0,
                "dH_target_binder": -8.0,
            },
        }],
        "t_min_c": 4.0, "t_max_c": 45.0, "t_step_c": 1.0, "t_ref_c": 25.0,
        "n_samples": 200, "seed": 11,
        "dH_priors": {k: {"mean": v["mean"], "sd": 0.0}
                      for k, v in tw.DEFAULT_DH_PRIORS.items()},
    }
    body.update(over)
    if binder_over:
        body["systems"][0]["binder"].update(binder_over)
    return body


def _post(body):
    response = client.post("/api/temperature/window", json=body)
    assert response.status_code == 200, response.text
    return response.json()


def test_endpoint_returns_a_window_and_echoes_what_it_ran():
    data = _post(_request())
    assert len(data["temperatures_c"]) == 42
    assert data["seed"] == 11
    assert data["echo"]["n_samples"] == 200
    assert data["echo"]["t_ref_c"] == 25.0
    assert data["units"]["target_conc"] == "nM"
    assert data["notes"]

    result = data["systems"][0]
    assert len(result["median_fc"]) == 42 == len(result["p_works"])
    assert result["dG_latch_cage_effective"] == -4.0


def test_nm_to_molar_conversion_happens_at_the_boundary():
    """1000 nM through the API must equal 1e-6 M straight into the engine."""
    api = _post(_request())["systems"][0]
    engine = run([synth(target_conc=1e-6)])["systems"][0]
    assert api["median_fc"] == pytest.approx(engine["median_fc"], rel=1e-12)


def test_performance_ten_systems_full_sweep():
    import time
    body = _request(n_samples=2000)
    body["systems"] = [dict(copy.deepcopy(body["systems"][0]), name=f"s{i}") for i in range(10)]
    body["dH_priors"] = None
    start = time.perf_counter()
    data = _post(body)
    elapsed = time.perf_counter() - start
    assert len(data["systems"]) == 10
    assert len(data["temperatures_c"]) == 42
    assert elapsed < 2.0, f"10 x 42 x 2000 took {elapsed:.2f}s"


# --- 10. every refusal is a 422 --------------------------------------------------------

@pytest.mark.parametrize("label,over", [
    ("t_min above t_max", {"t_min_c": 50.0}),
    ("zero step", {"t_step_c": 0.0}),
    ("negative step", {"t_step_c": -1.0}),
    ("too many grid points", {"t_min_c": 0.0, "t_max_c": 100.0, "t_step_c": 0.01}),
    ("n_samples too low", {"n_samples": 10}),
    ("n_samples too high", {"n_samples": 50_000}),
    ("confidence zero", {"confidence": 0.0}),
    ("confidence one", {"confidence": 1.0}),
    ("confidence above one", {"confidence": 1.4}),
    ("no systems", {"systems": []}),
    ("unknown dH prior", {"dH_priors": {"nope": {"mean": 1.0, "sd": 1.0}}}),
    ("negative prior sd", {"dH_priors": {"latch_cage": {"mean": -10.0, "sd": -1.0}}}),
])
def test_request_level_rejections_are_422(label, over):
    response = client.post("/api/temperature/window", json=_request(**over))
    assert response.status_code == 422, f"{label}: {response.text}"
    assert response.json()["error"]["message"]


def test_more_than_ten_systems_is_422():
    body = _request()
    body["systems"] = body["systems"] * 11
    assert client.post("/api/temperature/window", json=body).status_code == 422


@pytest.mark.parametrize("label,field", [("latch", "dG_latch_cage"), ("key", "dG_key_cage")])
def test_missing_required_dG_is_422(label, field):
    body = _request()
    body["systems"][0].pop(field)
    response = client.post("/api/temperature/window", json=body)
    assert response.status_code == 422
    assert field in response.json()["error"]["message"]


def test_missing_binder_is_422():
    body = _request()
    body["systems"][0].pop("binder")
    response = client.post("/api/temperature/window", json=body)
    assert response.status_code == 422
    assert "binder is required" in response.json()["error"]["message"]


@pytest.mark.parametrize("label,binder_over", [
    ("neither affinity", {"dG_target_binder": None}),
    ("both affinities", {"kd_nM": 5.0}),
    ("negative kd", {"dG_target_binder": None, "kd_nM": -1.0}),
    ("negative sd", {"dG_target_binder_sd": -0.5}),
    ("tm above 100", {"tm_c": 150.0}),
    ("tm below 0", {"tm_c": -5.0}),
])
def test_binder_rejections_are_422(label, binder_over):
    response = client.post("/api/temperature/window",
                           json=_request(binder_over=binder_over))
    assert response.status_code == 422, f"{label}: {response.text}"


@pytest.mark.parametrize("sequence", ["LISDAELX", "LISDAEL1", "lisdaelz", "LISD-AEL", "LISDAEL*"])
def test_non_standard_sequence_character_is_422(sequence):
    response = client.post("/api/temperature/window",
                           json=_request(binder_over={"sequence": sequence}))
    assert response.status_code == 422, response.text
    assert "sequence" in response.json()["error"]["field"]


def test_lowercase_sequence_is_accepted():
    """Case-insensitive: the sequence is provenance, not a prediction input."""
    data = _post(_request(binder_over={"sequence": "lisdaeleaifaeeldc"}))
    assert data["systems"][0]["median_fc"]


# --- 12. the shipped preset still runs -------------------------------------------------

def test_eclipse_preset_runs_and_is_flagged_as_placeholder():
    served = client.get("/api/temperature/window/preset")
    assert served.status_code == 200
    preset = served.json()
    assert preset["placeholder"] is True
    assert preset["placeholder_fields"]
    assert "replace with real data" in preset["badge"]

    data = _post(preset["request"])
    assert len(data["systems"]) == 1
    assert data["systems"][0]["median_fc"]
    assert any("PLACEHOLDER" in note for note in data["notes"]), data["notes"]


def test_preset_is_served_not_duplicated_in_the_frontend():
    """The browser must fetch the preset, so its placeholder numbers live in one place."""
    source = (WEB_JS / "temperature_window.js").read_text()
    assert "/api/temperature/window/preset" in source
    assert "LISDAELEAIFAEELDC" not in source, "preset sequence duplicated into JS"
    assert str(ECLIPSE_PRESET["systems"][0]["dG_latch_cage"]) not in source


# --- 13. the new view's JS at least parses ---------------------------------------------

@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_temperature_window_js_parses():
    for name in ("temperature_window.js", "shared.js"):
        result = subprocess.run(["node", "--check", str(WEB_JS / name)],
                                capture_output=True, text=True)
        assert result.returncode == 0, f"{name}: {result.stderr}"


def test_view_is_wired_into_the_page():
    index = (WEB_JS.parents[1] / "index.html").read_text()
    assert 'data-tab="temperature"' in index
    assert "static/js/temperature_window.js" in index
