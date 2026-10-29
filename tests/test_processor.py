"""
Unit tests on a synthetic fixture (tests/make_fixture.py); no Docker, no real data.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pytest

from processor import channels as C
from processor import fcsio
from processor import inputs as I
from processor import main as M
from processor import panel as P
from processor import pregate as G
from tests import make_fixture

N = 1500


@pytest.fixture(scope="module")
def fixture_dir(tmp_path_factory) -> Path:
    return make_fixture.make(tmp_path_factory.mktemp("fixture"), n_events=N)


def run(monkeypatch, input_dir: Path, output_dir: Path, **env) -> int:
    for k in ("SAMPLE_SHEET", "COMPENSATE", "PREGATE", "COFACTOR", "COFACTORS", "PANEL_OVERRIDES",
              "MAX_EVENTS_PER_SAMPLE", "SEED"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    return M.main()


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_fcs(path: Path):
    import flowio

    fd = flowio.FlowData(str(path))
    t = fd.text
    n = int(t["par"])
    pnn = [t[f"p{i}n"] for i in range(1, n + 1)]
    pns = [t.get(f"p{i}s", "") for i in range(1, n + 1)]
    return fd, pnn, pns, fd.as_array(preprocess=False)


# --- channels -------------------------------------------------------------------

@pytest.mark.parametrize("channel,marker,expected", [
    ("FSC-A", "", "scatter"), ("SSC-B-H", "", "scatter"), ("Time", "", "time"), ("AF-A", "", "other"),
    ("BUV395-A", "CD4", "marker"), ("FITC-A", "", "marker"), ("Yb171Di", "", "marker"), ("171Yb", "", "marker"),
    ("Event_length", "", "other"), ("Ir191Di", "191Ir_DNA1", "marker"), ("Unknown-A", "", "other"),
    ("B2-A", "", "marker"), ("Alexa Fluor 488-A", "HELIOS", "marker"), ("AF-A", "Autofluorescence", "other"),
])
def test_classify(channel, marker, expected):
    assert C.classify(channel, marker) == expected


def test_cytometer_type():
    assert C.cytometer_type({"cyt": "Aurora"}, ["FSC-A", "BV421-A"]) == "spectral"
    assert C.cytometer_type({"cyt": "ID7000"}, []) == "spectral"
    assert C.cytometer_type({"cyt": "Helios"}, []) == "mass"
    assert C.cytometer_type({"cyt": ""}, ["Time", "Yb171Di", "Er168Di", "Nd146Di", "Event_length"]) == "mass"
    assert C.cytometer_type({"cyt": "LSRFortessa"}, ["FSC-A", "FITC-A"]) == "flow"
    assert C.DEFAULT_COFACTOR == {"mass": 5.0, "spectral": 6000.0, "flow": 150.0}


@pytest.mark.parametrize("marker,flags", [
    ("CD4", (True, False)), ("FOXP3", (True, False)), ("HELIOS", (True, False)), ("CD25", (False, True)),
    ("HLA-DR", (False, True)), ("PD-1", (False, True)), ("Ki-67", (False, True)), ("GzmB", (False, True)),
    ("pSTAT5", (False, True)), ("CD45", (False, False)), ("CD45RA", (True, False)), ("191Ir_DNA1", (False, False)),
    ("171Yb_CD3", (True, False)), ("CD69_Nd144Di", (False, True)), ("Zombie NIR", (False, False)),
])
def test_default_flags(marker, flags):
    assert C.default_flags(marker) == flags


# CyTOF marker names (Kowli et al. panel and common spellings)
@pytest.mark.parametrize("marker,channel,flags", [
    # viability: neither clustering nor state
    ("115In_Dead", "In115Di", (False, False)), ("Dead", None, (False, False)), ("LD", None, (False, False)),
    ("L/D", None, (False, False)), ("Live", None, (False, False)), ("LIVE/DEAD Aqua", None, (False, False)),
    ("Viability", None, (False, False)), ("Viab", None, (False, False)), ("195Pt_Cisplatin", "Pt195Di", (False, False)),
    ("Pt198Di", "Pt198Di", (False, False)), ("103Rh_Intercalator", "Rh103Di", (False, False)),
    ("Zombie NIR", None, (False, False)),
    # ... without catching real markers
    ("LDLR", None, (True, False)), ("LDHA", None, (True, False)), ("CD57", None, (True, False)),
    ("194Pt_CD8", "Pt194Di", (True, False)), ("Pt194Di", "Pt194Di", (True, False)),
    # aliases: '_X_', '(X)', '/X'; state when any alias is
    ("168Er_CD154__CD40L_", "Er168Di", (False, True)), ("CD154 (CD40L)", None, (False, True)),
    ("CD40L/CD154", None, (False, True)), ("CD3_(UCHT1)", None, (True, False)),
    # DNA damage and PD-L1
    ("172Yb_pH2A.X", "Yb172Di", (False, True)), ("pH2AX", None, (False, True)), ("gH2AX", None, (False, True)),
    ("γH2AX", None, (False, True)), ("175Lu_PD-L1", "Lu175Di", (False, True)), ("CD274", None, (False, True)),
    # DNA, beads, CD45
    ("191Ir_DNA1", "Ir191Di", (False, False)), ("Ir193Di", "Ir193Di", (False, False)),
    ("Intercalator-Ir", "Ir191Di", (False, False)), ("140Ce_Beads", "Ce140Di", (False, False)),
    ("89Y_CD45", "Y89Di", (False, False)), ("CD45RA", None, (True, False)),
])
def test_default_flags_cytof(marker, channel, flags):
    assert C.default_flags(marker, channel) == flags


def test_gating_roles_and_aliases():
    assert C.gating_role("115In_Dead", "In115Di") == "viability"
    assert C.gating_role("191Ir_DNA1", "Ir191Di") == "dna"
    assert C.gating_role("Intercalator-Rh", "Rh103Di") == "viability"
    assert C.gating_role("103Rh_Intercalator") == "viability" and C.gating_role("191Ir_Intercalator") == "dna"
    assert C.strip_mass_tag("Ki67_Er168Di") == "Ki67" and C.strip_mass_tag("B220_Er168") == "B220"
    assert C.gating_role("140Ce_Beads", "Ce140Di") == "beads"
    assert C.gating_role("CD4") is None
    assert C.marker_aliases("168Er_CD154__CD40L_") == ["CD154CD40L", "CD154", "CD40L"]
    assert C.marker_aliases("CD154 (CD40L)") == ["CD154CD40L", "CD154", "CD40L"]
    assert C.marker_aliases("CD3") == ["CD3"]
    assert C.strip_mass_tag("CD3_Yb171Di") == "CD3" and C.strip_mass_tag("171Yb_CD3") == "CD3"


@pytest.mark.parametrize("channel", ["Event_length", "Center", "Offset", "Width", "Residual"])
def test_helios_parameters_are_other(channel):
    assert C.classify(channel, "") == "other" and C.classify(channel, channel) == "other"


# --- panel ----------------------------------------------------------------------

def test_panel_intersection_conflicts_and_spelling():
    h1 = P.Header("a", ["FSC-A", "BV750-A", "BV785-A", "PerCP-A", "Time"], ["", "HLADR", "CD28", "CD163", ""], {})
    h2 = P.Header("b", ["FSC-A", "BV750-A", "BV785-A", "PerCP-A", "Time"], ["", "HLADR", "CD28", "CD163", ""], {})
    h3 = P.Header("c", ["FSC-A", "BV750-A", "PerCP-A", "Time"], ["", "HLA-DR", "CD28", ""], {})
    panel = P.build([h1, h2, h3])
    assert [c.channel for c in panel.channels] == ["FSC-A", "BV750-A", "Time"]
    assert panel.by_channel()["BV750-A"].marker == "HLADR"
    assert {d["channel"] for d in panel.dropped["a"]} == {"BV785-A", "PerCP-A"}
    assert {d["channel"] for d in panel.dropped["c"]} == {"PerCP-A"}
    assert panel.renamed["c"] == [{"channel": "BV750-A", "from": "HLA-DR", "to": "HLADR"}]
    assert panel.renamed["a"] == []


def test_configure_overrides():
    h = P.Header("a", ["FSC-A", "PE-A", "BV421-A"], ["", "CD25", "CD4"], {})
    panel = P.build([h])
    warnings = P.configure(panel, 6000, {"CD4": 3000}, {"cd25": {"clustering": True}, "CD999": {"state": True}})
    by = panel.by_channel()
    assert by["BV421-A"].cofactor == 3000 and by["PE-A"].cofactor == 6000
    assert (by["PE-A"].clustering, by["PE-A"].state) == (True, True)
    assert by["FSC-A"].cofactor is None and not by["FSC-A"].clustering
    assert warnings == ["PANEL_OVERRIDES: 'CD999' matches no marker channel in the panel"]
    with pytest.raises(ValueError):
        P.configure(panel, 6000, {"CD4": -1}, {})
    with pytest.raises(ValueError):
        P.parse_json_object("[1]", "COFACTORS")


# --- inputs ---------------------------------------------------------------------

def test_sanitise_and_batch():
    assert I.sanitise("B01_tcs3_CD45+CD14-") == "B01_tcs3_CD45_CD14-"
    assert I.sanitise("a b..c") == "a_b_c"
    assert I.batch_from_date("11-Mar-2022") == "2022-03-11"
    assert I.batch_from_date("") == "NA"


def test_compensation_plan():
    ident = (["A", "B"], np.array([[1, 1e-6], [0, 1]]))
    spill = (["A", "B"], np.array([[1, 0.2], [0.05, 1]]))
    assert M.compensation_plan("auto", "spectral", ident)[0] is False
    assert "identity" in M.compensation_plan("auto", "spectral", ident)[1]
    assert M.compensation_plan("auto", "spectral", spill)[0] is False
    assert M.compensation_plan("auto", "flow", spill)[0] is True
    assert M.compensation_plan("true", "spectral", spill)[0] is True
    assert M.compensation_plan("false", "flow", spill)[0] is False
    assert M.compensation_plan("auto", "flow", None)[0] is False


def test_compensation_recovers_true_values(tmp_path):
    rng = np.random.default_rng(3)
    true = np.column_stack([rng.uniform(1e5, 2e5, 500), rng.lognormal(9, 1, 500), rng.lognormal(8, 1, 500),
                            np.arange(500, dtype=float)])
    spill = np.array([[1.0, 0.15], [0.02, 1.0]])
    observed = true.copy()
    observed[:, 1:3] = true[:, 1:3] @ spill
    meta = {"cyt": "LSRFortessa", "spillover": make_fixture.spill_text(["FITC-A", "PE-A"], spill)}
    make_fixture.write_fcs(tmp_path / "c.fcs", observed.astype(np.float32), ["FSC-A", "FITC-A", "PE-A", "Time"],
                           ["", "CD3", "CD4", ""], meta)
    text, *_ = fcsio.read_header(tmp_path / "c.fcs")
    parsed = fcsio.parse_spill(text)
    assert parsed[0] == ["FITC-A", "PE-A"] and np.allclose(parsed[1], spill)
    events, pnn, off_by_one = fcsio.read_events(tmp_path / "c.fcs", parsed)
    assert np.allclose(events[:, 1:3], true[:, 1:3], rtol=1e-3, atol=1.0)


def test_subsample_index_is_seeded_and_sorted():
    a = M.subsample_index(1000, 100, 42, "s1")
    assert len(a) == 100 and np.all(np.diff(a) > 0)
    assert np.array_equal(a, M.subsample_index(1000, 100, 42, "s1"))
    assert not np.array_equal(a, M.subsample_index(1000, 100, 42, "s2"))
    assert M.subsample_index(50, 100, 42, "s1") is None and M.subsample_index(50, 0, 42, "s1") is None


# --- pre-gating -----------------------------------------------------------------

def cytof_panel() -> P.Panel:
    pnn = [c for c, _ in make_fixture.CYTOF_PANEL]
    pns = [m for _, m in make_fixture.CYTOF_PANEL]
    panel = P.build([P.Header("a", pnn, pns, {})])
    P.configure(panel, 5.0, {}, {})
    return panel


def test_find_channels_and_plan():
    gc = G.find_channels(cytof_panel().channels)
    assert gc.dna == ["Ir191Di", "Ir193Di"] and gc.beads == ["Ce140Di"] and gc.viability == ["In115Di"]
    assert G.plan("auto", "mass", gc)[0] is True
    assert G.plan("auto", "spectral", gc)[0] is False
    assert G.plan("false", "mass", gc)[0] is False
    assert G.plan("true", "spectral", gc)[0] is True
    assert G.plan("true", "mass", G.GateChannels())[0] is False
    # spectral demo-like panel: no gating channel (CD45 is not one)
    spectral = P.build([P.Header("a", ["FSC-A", "BUV395-A", "Spark Blue 550-A"], ["", "CD4", "CD45"], {})])
    P.configure(spectral, 6000, {}, {})
    assert not G.find_channels(spectral.channels).any()
    # an unlabelled 140Ce is a bead channel; a labelled 151Eu antibody is not
    p = P.build([P.Header("a", ["Ce140Di", "Eu151Di", "Eu153Di"], ["", "151Eu_CD107a", "153Eu_Beads"], {})])
    P.configure(p, 5.0, {}, {})
    assert G.find_channels(p.channels).beads == ["Ce140Di", "Eu153Di"]


def gate_synthetic(n=40000, seed=3, **kw):
    ev = make_fixture.cytof_events(n, np.random.default_rng(seed), **kw)
    pnn = [c for c, _ in make_fixture.CYTOF_PANEL]
    gc = G.GateChannels(dna=["Ir191Di", "Ir193Di"], beads=["Ce140Di"], viability=["In115Di"])
    cols = {c: ev[:, pnn.index(c)].astype(float) for c in gc.dna + gc.beads + gc.viability}
    return G.gate(cols, gc, {c: 5.0 for c in cols}), gc, cols


def test_gates_on_synthetic_cytof():
    res, gc, _ = gate_synthetic()
    n = res.events_in
    f = make_fixture.CYTOF_FRACTIONS
    assert abs(res.beads_removed / n - f["beads"]) < 0.005
    assert res.beads_free == res.beads_removed  # the synthetic beads have no DNA
    assert abs(res.dna_removed / n - f["debris"]) < 0.01
    assert abs(res.dead_removed / n - f["dead"]) < 0.01
    assert res.events_out == n - res.beads_removed - res.dna_removed - res.dead_removed
    th = res.thresholds
    assert 2.0 <= th["bead"] <= 6.7
    assert abs(th["dna_mode"] - 4.5) < 0.05
    assert th["dna_low"] <= th["dna_mode"] - G.DNA_LOW_MIN + 1e-6
    assert th["dna_high"] >= th["dna_mode"] + G.DNA_HIGH_MIN - 1e-6  # generous on the high side (doublets)
    assert 3.6 < th["dead"] < 5.0 and "valley" in th["dead_method"]
    assert not res.warnings
    # deterministic
    again, *_ = gate_synthetic()
    assert np.array_equal(res.keep, again.keep) and again.thresholds == th


def test_viability_unimodal_removes_little():
    rng = np.random.default_rng(5)
    v = np.where(rng.random(20000) < 0.5, 0, rng.normal(2.8, 0.5, 20000)).clip(0)
    thr, info = G.viability_threshold(v)
    assert (v > thr).mean() < 0.02 and info["dead_method"].startswith("upper fence")
    # a channel that is high for most events never loses more than MAX_DEAD_FRACTION
    v = np.concatenate([rng.normal(1.0, 0.2, 4000), rng.normal(5.0, 0.3, 16000)])
    thr, info = G.viability_threshold(v)
    assert (v > thr).mean() <= G.MAX_DEAD_FRACTION + 1e-9


def test_gate_warns_when_it_removes_most_events():
    n = 20000
    rng = np.random.default_rng(1)
    s = np.where(rng.random(n) < 0.3, rng.normal(4.5, 0.15, n), rng.uniform(0.3, 3.0, n))
    cols = {"Ir191Di": 5 * np.sinh(s), "Ir193Di": 5 * np.sinh(s + 0.3)}
    res = G.gate(cols, G.GateChannels(dna=["Ir191Di", "Ir193Di"]), {})
    assert res.dna_removed / n > 0.6
    assert len(res.warnings) == 1 and "DNA gate removed" in res.warnings[0]
    assert res.events_out > 0.25 * n  # applied anyway


def test_bead_gate_without_beads():
    rng = np.random.default_rng(2)
    thr, why = G.bead_threshold(np.arcsinh(rng.uniform(0, 3, 10000) / 5))
    assert thr is None and "no bead population" in why


def test_offset_off_by_one_is_retried(tmp_path):
    out = make_fixture.make_cytof(tmp_path, n_events=500)
    import flowio
    from flowio.exceptions import FCSParsingError

    with pytest.raises(FCSParsingError):
        flowio.FlowData(str(out / "pt01_pre_tx.fcs"))
    events, pnn, off = fcsio.read_events(out / "pt01_pre_tx.fcs", None)
    assert off is True and events.shape[0] == 500 and pnn[0] == "Time"
    # a file that is fine is read without the option
    make_fixture.write_fcs(tmp_path / "ok.fcs", np.ones((10, 2), np.float32), ["A", "B"], ["", ""], {})
    assert fcsio.read_events(tmp_path / "ok.fcs", None)[2] is False
    # any other offset problem still fails
    raw = bytearray((out / "pt01_month_6.fcs").read_bytes())
    end = int(raw[34:42].decode())
    raw[34:42] = str(end + 2).rjust(8).encode()
    (tmp_path / "bad.fcs").write_bytes(bytes(raw))
    with pytest.raises(Exception):
        fcsio.read_events(tmp_path / "bad.fcs", None)


def test_end_to_end_cytof(tmp_path, monkeypatch):
    fixture = make_fixture.make_cytof(tmp_path / "fx", n_events=6000)
    out = tmp_path / "out"
    assert run(monkeypatch, fixture, out, MAX_EVENTS_PER_SAMPLE="4000") == 0
    meta = out / "metadata"
    rows = {r["channel"]: r for r in read_csv(meta / "panel.csv")}
    assert (rows["In115Di"]["clustering"], rows["In115Di"]["state"]) == ("false", "false")
    assert rows["Er168Di"]["state"] == "true" and rows["Yb172Di"]["state"] == "true"
    assert rows["Lu175Di"]["state"] == "true" and rows["Nd145Di"]["clustering"] == "true"
    for ch in ("Event_length", "Center", "Offset", "Width", "Residual"):
        assert rows[ch]["type"] == "other" and rows[ch]["cofactor"] == ""
    prep = json.loads((meta / "prep.json").read_text())
    assert prep["cytometerType"] == "mass"
    assert prep["pregate"]["applied"] is True and prep["pregate"]["mode"] == "auto"
    assert prep["pregate"]["channels"] == {"dna": ["Ir191Di", "Ir193Di"], "beads": ["Ce140Di"],
                                           "viability": ["In115Di"]}
    gating = {r["sample_id"]: r for r in read_csv(meta / "gating.csv")}
    assert set(gating) == {"pt01_pre_tx", "pt01_month_6"}
    for sid, g in gating.items():
        s = prep["samples"][sid]
        assert s["ignore_offset_error"] is True
        pg = s["pregate"]
        assert int(g["events_in"]) == s["events_in"] and int(g["events_out"]) == pg["events_out"]
        assert int(g["events_in"]) - int(g["beads_removed"]) - int(g["dna_removed"]) - int(g["dead_removed"]) \
            == int(g["events_out"])
        assert 0.10 < 1 - int(g["events_out"]) / int(g["events_in"]) < 0.22
        assert "bead=" in g["thresholds"] and "dead=" in g["thresholds"]
        assert (meta / pg["plot"]).is_file() and (meta / pg["plot"]).stat().st_size < 200_000
        # subsample drawn from the gated events
        assert s["events_out"] == min(4000, pg["events_out"])
    fd, pnn, _, ev = read_fcs(out / "fcs" / "pt01_pre_tx.fcs")
    assert "pregate=beads/dna/live" in fd.text["pennsieve_prep"]
    assert np.all(np.diff(ev[:, pnn.index("Time")]) >= 0)
    # the gated events are gone: no bead events left
    assert ev[:, pnn.index("Ce140Di")].max() < prep["samples"]["pt01_pre_tx"]["pregate"]["thresholds"]["bead"]

    out2 = tmp_path / "out2"
    assert run(monkeypatch, fixture, out2, PREGATE="false") == 0
    assert not (out2 / "metadata" / "gating.csv").exists() and not (out2 / "metadata" / "gates").exists()
    prep2 = json.loads((out2 / "metadata" / "prep.json").read_text())
    assert prep2["pregate"]["applied"] is False
    assert prep2["samples"]["pt01_pre_tx"]["events_out"] == prep2["samples"]["pt01_pre_tx"]["events_in"]
    assert run(monkeypatch, fixture, tmp_path / "out3", PREGATE="sometimes") == 2


# --- end to end -----------------------------------------------------------------

def test_end_to_end(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir.parent, out) == 0
    meta, fcs = out / "metadata", out / "fcs"
    assert sorted(p.name for p in fcs.glob("*.fcs")) == ["batch1_s1.fcs", "batch1_s2.fcs", "batch2_s1.fcs",
                                                          "extra_s1.fcs"]
    samples = read_csv(meta / "samples.csv")
    assert list(samples[0].keys()) == ["sample_id", "file", "batch", "condition", "n_events", "donor"]
    by = {r["sample_id"]: r for r in samples}
    assert by["batch1_s2"] == {"sample_id": "batch1_s2", "file": "b1_s2.fcs", "batch": "batch1",
                               "condition": "stimulated", "n_events": str(N), "donor": "D1"}
    assert by["extra_s1"]["batch"] == "2022-05-05" and by["extra_s1"]["condition"] == "NA"
    assert by["extra_s1"]["donor"] == ""

    panel = read_csv(meta / "panel.csv")
    chans = [r["channel"] for r in panel]
    assert "BV785-A" not in chans and "PerCP-A" not in chans
    rows = {r["channel"]: r for r in panel}
    assert rows["BV750-A"]["marker"] == "HLADR"
    assert rows["AF-A"]["type"] == "other" and rows["AF-A"]["cofactor"] == ""
    assert rows["Time"]["type"] == "time" and rows["FSC-A"]["type"] == "scatter"
    assert rows["PE-A"] == {"channel": "PE-A", "marker": "CD25", "type": "marker", "cofactor": "6000",
                            "clustering": "false", "state": "true"}
    assert rows["Pacific Blue-A"]["clustering"] == "true"

    prep = json.loads((meta / "prep.json").read_text())
    s = prep["samples"]["batch1_s1"]
    assert s["compensation"]["applied"] is False and "identity" in s["compensation"]["reason"]
    assert {d["channel"] for d in s["channels_dropped"]} == {"BV785-A", "PerCP-A"}
    assert prep["samples"]["batch2_s1"]["markers_renamed"][0]["from"] == "HLA-DR"
    assert prep["cytometerType"] == "spectral"
    assert prep["pregate"]["applied"] is False and not (meta / "gating.csv").exists()
    assert "pregate" not in s and s["ignore_offset_error"] is False
    assert {k["file"] for k in prep["skipped"]} == {"empty.fcs", "broken.fcs"} or \
        {Path(k["file"]).name for k in prep["skipped"]} == {"empty.fcs", "broken.fcs"}
    for d in (meta, fcs):
        prov = json.loads((d / "provenance.json").read_text())
        assert prov["app"] == "processor-fcs-prep" and "flowkit" in prov["tools"]
        assert prov["parameters"]["COFACTOR"] == "auto"

    # values: markers arcsinh(x / 6000) of the source, the rest untouched; keywords kept
    fd, pnn, pns, ev = read_fcs(fcs / "batch2_s1.fcs")
    assert fd.version == "3.1" and fd.text["datatype"] == "F"
    assert pnn == chans
    assert pns[pnn.index("BV750-A")] == "HLADR" and pns[pnn.index("PE-A")] == "CD25"
    assert fd.text["pennsieve_prep"].startswith("transform=arcsinh;cofactor=6000;compensation=none")
    assert fd.text["date"] == "05-May-2022" and fd.text["cyt"] == "Aurora"
    _, spnn, _, src = read_fcs(fixture_dir / "b2_s1.fcs")
    for ch in ("FSC-A", "AF-A", "Time"):
        assert np.allclose(ev[:, pnn.index(ch)], src[:, spnn.index(ch)])
    assert np.allclose(ev[:, pnn.index("PE-A")], np.arcsinh(src[:, spnn.index("PE-A")] / 6000), atol=1e-5)
    # $PnR covers every value written (flowCore truncates at $PnR by default)
    for i, ch in enumerate(pnn, start=1):
        assert float(fd.text[f"p{i}r"]) >= ev[:, i - 1].max()


def test_parameters_and_subsample(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir, out, MAX_EVENTS_PER_SAMPLE="300", COFACTOR="150",
               COFACTORS='{"CD3": 500}', PANEL_OVERRIDES='{"HLA-DR": {"clustering": true, "state": false}}') == 0
    rows = {r["marker"]: r for r in read_csv(out / "metadata" / "panel.csv")}
    assert rows["CD3"]["cofactor"] == "500" and rows["CD4"]["cofactor"] == "150"
    assert (rows["HLADR"]["clustering"], rows["HLADR"]["state"]) == ("true", "false")
    assert {r["n_events"] for r in read_csv(out / "metadata" / "samples.csv")} == {"300"}
    _, pnn, _, ev = read_fcs(out / "fcs" / "batch1_s1.fcs")
    assert np.all(np.diff(ev[:, pnn.index("Time")]) >= 0)


def test_no_sheet_defaults(fixture_dir, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(monkeypatch, fixture_dir, out, SAMPLE_SHEET="none.csv") == 0
    samples = read_csv(out / "metadata" / "samples.csv")
    assert list(samples[0].keys()) == ["sample_id", "file", "batch", "condition", "n_events"]
    assert {r["sample_id"] for r in samples} == {"b1_s1", "b1_s2", "b2_s1", "extra_s1"}
    assert {r["condition"] for r in samples} == {"NA"}


def test_symlinked_input_and_duplicates(fixture_dir, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    (inp / "a").mkdir(parents=True)
    (inp / "b").mkdir()
    os.symlink(fixture_dir / "b1_s1.fcs", inp / "a" / "x.fcs")
    os.symlink(fixture_dir / "b1_s2.fcs", inp / "b" / "x.fcs")       # same stem: duplicate sample_id
    os.symlink(fixture_dir, inp / "linked")                         # a symlinked folder
    out = tmp_path / "out"
    assert run(monkeypatch, inp, out, SAMPLE_SHEET="none.csv") == 0
    prep = json.loads((out / "metadata" / "prep.json").read_text())
    assert "x" not in prep["samples"]
    assert sum(1 for s in prep["skipped"] if s["reason"] == "duplicate sample_id") == 2
    assert "b2_s1" in prep["samples"]


def test_failures(tmp_path, monkeypatch):
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "readme.txt").write_text("nothing")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 1
    (tmp_path / "in" / "bad.fcs").write_bytes(b"garbage")
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out") == 1
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out", COMPENSATE="maybe") == 2
    assert run(monkeypatch, tmp_path / "in", tmp_path / "out", COFACTORS="{bad") == 2
