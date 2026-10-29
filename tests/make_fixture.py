"""
Synthetic FCS fixture shaped like the demo data (Cytek Aurora, spectrally
unmixed, two panels across batches), small enough for unit tests.

  python -m tests.make_fixture data/fixture [--events 2000]

writes data/fixture/linked/fcs/ (what a linked `fcs` input looks like) with:
  b1_s1.fcs, b1_s2.fcs     batch 1 panel: BV785-A=CD28, PerCP-A=CD163, BV750-A=HLADR,
                           an identity $SPILLOVER with one 1e-6 entry (as SpectroFlo writes it)
  b2_s1.fcs                batch 2 panel: no BV785-A, PerCP-A=CD28, BV750-A=HLA-DR, no $SPILLOVER
  extra_s1.fcs             not in the sample sheet (processed with defaults)
  empty.fcs                0 events (skipped)
  broken.fcs               not an FCS file (skipped)
  samples.csv              sample sheet: file, sample_id, batch, condition, donor

make_cytof() writes raw (not pre-gated) mass cytometry files shaped like the Kowli
et al. CyTOF data: beads (140Ce), debris, DNA+ cells (Ir191/Ir193), dead cells
bright on 115In, the Helios Gaussian parameters, and $ENDDATA off by one byte as
the CyTOF software and bead normaliser write it.
"""
from __future__ import annotations

import argparse
import csv
from array import array
from pathlib import Path

import numpy as np

SCATTER = ["FSC-A", "FSC-H", "SSC-A", "SSC-H"]
PANEL_B1 = [("AF-A", ""), ("BUV395-A", "CD4"), ("BUV563-A", "CD3"), ("BUV805-A", "CD8"), ("PE-A", "CD25"),
            ("Pacific Blue-A", "FOXP3"), ("BV750-A", "HLADR"), ("BV785-A", "CD28"), ("PerCP-A", "CD163"),
            ("Spark Blue 550-A", "CD45")]
PANEL_B2 = [("AF-A", ""), ("BUV395-A", "CD4"), ("BUV563-A", "CD3"), ("BUV805-A", "CD8"), ("PE-A", "CD25"),
            ("Pacific Blue-A", "FOXP3"), ("BV750-A", "HLA-DR"), ("PerCP-A", "CD28"), ("Spark Blue 550-A", "CD45")]


def events_for(n: int, n_markers: int, rng: np.random.Generator) -> np.ndarray:
    """Scatter, AF, markers (two populations, negative ones spread around 0) and a monotonic Time."""
    scatter = rng.normal(8e5, 1.5e5, size=(n, len(SCATTER))).clip(1e4)
    pos = rng.random((n, n_markers)) < 0.4
    markers = np.where(pos, rng.lognormal(np.log(4e4), 0.5, (n, n_markers)), rng.normal(0, 800, (n, n_markers)))
    time = np.sort(rng.uniform(0, 1500, n))
    return np.column_stack([scatter, markers, time]).astype(np.float32)


def write_fcs(path: Path, events: np.ndarray, pnn: list[str], pns: list[str], meta: dict[str, str]) -> None:
    import flowio

    data = array("f")
    data.frombytes(np.ascontiguousarray(events, dtype="<f4").tobytes())
    with open(path, "wb") as fh:
        flowio.create_fcs(fh, data, pnn, opt_channel_names=pns, metadata_dict=meta)


def spill_text(detectors: list[str], matrix: np.ndarray) -> str:
    return ",".join([str(len(detectors)), *detectors, *(f"{v:g}" for v in matrix.ravel())])


def make(root: Path, n_events: int = 2000, seed: int = 1) -> Path:
    rng = np.random.default_rng(seed)
    out = root / "linked" / "fcs"
    out.mkdir(parents=True, exist_ok=True)
    files = [("b1_s1.fcs", PANEL_B1, "11-Mar-2022", True), ("b1_s2.fcs", PANEL_B1, "11-Mar-2022", True),
             ("b2_s1.fcs", PANEL_B2, "05-May-2022", False), ("extra_s1.fcs", PANEL_B1, "05-May-2022", True)]
    for name, panel, date, spill in files:
        pnn = SCATTER + [c for c, _ in panel] + ["Time"]
        pns = [""] * len(SCATTER) + [m for _, m in panel] + [""]
        meta = {"cyt": "Aurora", "creator": "SpectroFlo 3.0.3", "date": date, "timestep": "0.0001"}
        if spill:
            det = [c for c, _ in panel if c != "AF-A"]
            m = np.eye(len(det))
            m[1, 0] = 1e-6
            meta["spillover"] = spill_text(det, m)
        write_fcs(out / name, events_for(n_events, len(panel), rng), pnn, pns, meta)
    pnn = SCATTER + [c for c, _ in PANEL_B1] + ["Time"]
    pns = [""] * len(SCATTER) + [m for _, m in PANEL_B1] + [""]
    write_fcs(out / "empty.fcs", np.zeros((0, len(pnn)), np.float32), pnn, pns, {"cyt": "Aurora"})
    (out / "broken.fcs").write_bytes(b"this is not an FCS file\n" * 10)
    with open(out / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "sample_id", "batch", "condition", "donor"])
        w.writerow(["b1_s1.fcs", "batch1_s1", "batch1", "unstimulated", "D1"])
        w.writerow(["b1_s2.fcs", "batch1_s2", "batch1", "stimulated", "D1"])
        w.writerow(["b2_s1.fcs", "batch2_s1", "batch2", "unstimulated", "D2"])
        w.writerow(["missing.fcs", "missing", "batch2", "stimulated", "D3"])
    return out


CYTOF_PANEL = [("Time", ""), ("Event_length", ""), ("In115Di", "115In_Dead"), ("Ce140Di", "140Ce_Beads"),
               ("Pr141Di", "141Pr_HLA-DR"), ("Nd145Di", "145Nd_CD4"), ("Sm154Di", "154Sm_CD3"),
               ("Er168Di", "168Er_CD154__CD40L_"), ("Yb172Di", "172Yb_pH2A.X"), ("Lu175Di", "175Lu_PD-L1"),
               ("Ir191Di", "191Ir_DNA1"), ("Ir193Di", "193Ir_DNA2"), ("Center", ""), ("Offset", ""), ("Width", ""),
               ("Residual", "")]
# Population fractions of make_cytof's events.
CYTOF_FRACTIONS = {"beads": 0.02, "debris": 0.08, "dead": 0.06}


def cytof_events(n: int, rng: np.random.Generator, dna_mode: float = 4.5) -> np.ndarray:
    """Raw (untransformed) CyTOF-like events in CYTOF_PANEL order; values are 5 * sinh(arcsinh-level)."""
    raw = lambda a: 5 * np.sinh(np.clip(a, 0, None))  # noqa: E731
    kind = rng.choice(4, size=n, p=[1 - sum(CYTOF_FRACTIONS.values()), CYTOF_FRACTIONS["beads"],
                                    CYTOF_FRACTIONS["debris"], CYTOF_FRACTIONS["dead"]])  # 0 live, 1 bead, 2 debris, 3 dead
    cell = (kind == 0) | (kind == 3)
    s = np.where(cell, rng.normal(dna_mode, 0.15, n), rng.uniform(0.3, 3.0, n))
    s[kind == 1] = rng.uniform(0, 0.8, (kind == 1).sum())
    d1 = s - 0.3 + rng.normal(0, 0.05, n)
    d2 = s + 0.3 + rng.normal(0, 0.05, n)
    bead = np.where(kind == 1, rng.normal(7.2, 0.1, n), np.where(rng.random(n) < 0.8, 0, rng.uniform(0, 1, n)))
    live_v = np.where(rng.random(n) < 0.5, 0, rng.normal(2.8, 0.5, n))
    viab = np.where(kind == 3, rng.normal(5.3, 0.3, n), live_v)
    markers = [np.where(rng.random(n) < 0.4, rng.normal(4, 0.5, n), rng.uniform(0, 1, n)) for _ in range(6)]
    gauss = [rng.normal(500, 100, n), rng.normal(100, 20, n), rng.normal(40, 8, n), rng.normal(80, 20, n)]
    cols = [np.sort(rng.uniform(0, 6e5, n)), rng.integers(10, 30, n).astype(float), raw(viab), raw(bead),
            *[raw(m) for m in markers], raw(d1), raw(d2), *gauss]
    return np.column_stack(cols).astype(np.float32)


def off_by_one(path: Path) -> None:
    """Make $ENDDATA (HEADER and TEXT) point one byte past the DATA segment, as some CyTOF exports do."""
    import re

    raw = bytearray(path.read_bytes())
    end = int(raw[34:42].decode().strip())
    new = str(end + 1)
    if len(new) != len(str(end)):
        raise ValueError("pick another event count: $ENDDATA + 1 changes its width")
    raw[34:42] = new.rjust(8).encode()
    text_end = int(raw[18:26].decode().strip())
    text = bytes(raw[:text_end + 1])
    m = re.search(rb"([^A-Za-z0-9])\$ENDDATA\1(\d+)", text)
    if not m:
        raise ValueError("no $ENDDATA in TEXT")
    raw[m.start(2):m.end(2)] = new.encode()
    path.write_bytes(bytes(raw))


def make_cytof(root: Path, n_events: int = 6000, seed: int = 7) -> Path:
    """Two raw CyTOF files (off-by-one $ENDDATA) and a sample sheet under root/linked/fcs/."""
    rng = np.random.default_rng(seed)
    out = root / "linked" / "fcs"
    out.mkdir(parents=True, exist_ok=True)
    pnn = [c for c, _ in CYTOF_PANEL]
    pns = [m for _, m in CYTOF_PANEL]
    for i, name in enumerate(("pt01_pre_tx.fcs", "pt01_month_6.fcs")):
        n = n_events + 7 * i
        meta = {"cyt": "DVSSCIENCES-CYTOF-6.7.1014", "date": "19-Apr-2018"}
        for attempt in range(5):
            write_fcs(out / name, cytof_events(n + attempt, rng), pnn, pns, meta)
            try:
                off_by_one(out / name)
                break
            except ValueError:
                continue
    with open(out / "samples.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "sample_id", "batch", "condition", "timepoint"])
        w.writerow(["pt01_pre_tx.fcs", "pt01_pre_tx", "batch01", "patient", "pre_tx"])
        w.writerow(["pt01_month_6.fcs", "pt01_month_6", "batch01", "patient", "month_6"])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--events", type=int, default=2000)
    args = ap.parse_args()
    out = make(args.root, args.events)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
