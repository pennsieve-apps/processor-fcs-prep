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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--events", type=int, default=2000)
    args = ap.parse_args()
    out = make(args.root, args.events)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
