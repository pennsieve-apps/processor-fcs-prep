"""
processor-fcs-prep: the first step of the immune cytometry workflows.

Finds the FCS files (and an optional sample sheet) under INPUT_DIR, works out
the panel common to every sample, compensates when that is needed, applies
arcsinh(x / cofactor) to the marker channels and writes

  OUTPUT_DIR/fcs/<sample_id>.fcs        FCS 3.1 float32, panel channels only
  OUTPUT_DIR/fcs/provenance.json
  OUTPUT_DIR/metadata/samples.csv       sample_id, file, batch, condition, n_events, sheet columns
  OUTPUT_DIR/metadata/panel.csv         channel, marker, type, cofactor, clustering, state
  OUTPUT_DIR/metadata/prep.json         per-sample details (cytometer, compensation, dropped channels, ...)
  OUTPUT_DIR/metadata/gating.csv        per-sample pre-gating counts and thresholds (when PREGATE applies)
  OUTPUT_DIR/metadata/gates/<id>.png    per-sample gate plots (when PREGATE applies)
  OUTPUT_DIR/metadata/provenance.json

Files are processed one at a time; only headers are read up front.
"""
from __future__ import annotations

import csv
import json
import logging
import os
import platform
import sys
import time
import traceback
import zlib
from collections import Counter
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

import numpy as np

from processor import channels as C
from processor import fcsio
from processor import inputs as I
from processor import panel as P
from processor import pregate as G

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("processor-fcs-prep")
log.setLevel(logging.INFO)
for noisy in ("flowkit", "flowio", "bokeh", "matplotlib"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "processor-fcs-prep"
OUTPUTS = ("fcs", "metadata")
REPO_ROOT = Path(__file__).resolve().parent.parent
# Keywords copied from the source file (when present) into the output file.
KEEP_KEYWORDS = ("date", "btim", "etim", "cyt", "cytsn", "inst", "op", "src", "smno", "proj", "timestep", "vol",
                 "creator", "fil")


def version_of(root: Path) -> str:
    try:
        return (root / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for pkg in ("FlowKit", "FlowIO", "FlowUtils", "numpy", "pandas", "scipy", "matplotlib"):
        try:
            out[pkg.lower()] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def env_params() -> dict[str, str]:
    return {
        "SAMPLE_SHEET": os.environ.get("SAMPLE_SHEET", "samples.csv"),
        "COMPENSATE": (os.environ.get("COMPENSATE") or "auto").strip().lower(),
        "PREGATE": (os.environ.get("PREGATE") or "auto").strip().lower(),
        "COFACTOR": (os.environ.get("COFACTOR") or "auto").strip(),
        "COFACTORS": os.environ.get("COFACTORS") or "{}",
        "PANEL_OVERRIDES": os.environ.get("PANEL_OVERRIDES") or "{}",
        "MAX_EVENTS_PER_SAMPLE": (os.environ.get("MAX_EVENTS_PER_SAMPLE") or "0").strip(),
        "SEED": (os.environ.get("SEED") or "42").strip(),
    }


def compensation_plan(mode: str, cyt_type: str, spill) -> tuple[bool, str]:
    """(apply?, why) for one file."""
    if mode == "false":
        return False, "COMPENSATE=false"
    if spill is None:
        return False, "no $SPILLOVER in the file"
    if fcsio.is_identity(spill[1]):
        off = float(abs(spill[1] - np.eye(spill[1].shape[0])).max())
        return False, f"$SPILLOVER is the identity matrix (largest deviation {off:.2g})"
    if mode == "true":
        return True, "COMPENSATE=true: applied $SPILLOVER"
    if cyt_type == "spectral":
        return False, "spectral cytometer: data are already unmixed; $SPILLOVER not applied"
    if cyt_type == "mass":
        return False, "mass cytometer: $SPILLOVER not applied (use COMPENSATE=true to force)"
    return True, "applied $SPILLOVER from the file"


def subsample_index(n: int, max_events: int, seed: int, sample_id: str) -> np.ndarray | None:
    """Seeded random subset of event indices, kept in acquisition order (Time stays monotonic)."""
    if max_events <= 0 or n <= max_events:
        return None
    rng = np.random.default_rng([seed, zlib.crc32(sample_id.encode())])
    return np.sort(rng.choice(n, size=max_events, replace=False))


def ranges(text: dict, pnn: list[str], chans: list, out: np.ndarray) -> list[str]:
    """
    $PnR per output channel: the source's range, raised to cover the largest value
    written. Some exports (SpectroFlo's Time: $PnR=129, values up to 1e6) understate
    it, and flowCore's read.FCS truncates at $PnR by default.
    """
    src = {c: text.get(f"p{i + 1}r", "") for i, c in enumerate(pnn)}
    result = []
    for j, ch in enumerate(chans):
        try:
            r = float(src.get(ch.channel) or 0)
        except ValueError:
            r = 0.0
        top = float(np.nanmax(out[:, j])) if out.shape[0] else 0.0
        if not np.isfinite(top):
            top = 0.0
        need = int(np.ceil(top)) + 1
        result.append(str(int(r)) if r >= need and r == int(r) else str(max(int(np.ceil(r)), need)))
    return result


def gating_cofactors(panel: P.Panel, gc: G.GateChannels, main_type: str) -> dict[str, float]:
    """arcsinh cofactor per gating channel: 5 for mass cytometry (the gates are tuned on it), else the panel's."""
    by = panel.by_channel()
    return {c: (5.0 if main_type == "mass" else (by[c].cofactor or 5.0)) for c in gc.dna + gc.beads + gc.viability}


def threshold_text(th: dict) -> str:
    keys = ("bead", "dna_low", "dna_high", "dna_diagonal_center", "dna_diagonal_halfwidth", "dead")
    return ";".join(f"{k}={'NA' if th.get(k) is None else th[k]}" for k in keys if k in th)


def process_sample(rec: I.SampleRecord, header: P.Header, cyt_type: str, panel: P.Panel, params: dict,
                   fcs_dir: Path, app_version: str, gating: dict | None = None) -> dict[str, Any]:
    """
    Compensate, pre-gate (when `gating` is given), select, transform and write one
    sample. Raises on failure.
    """
    t0 = time.monotonic()
    text = header.text
    spill = fcsio.parse_spill(text)
    apply_comp, comp_reason = compensation_plan(params["COMPENSATE"], cyt_type, spill)
    if apply_comp:
        missing = [d for d in spill[0] if d not in header.pnn]
        if missing:
            raise ValueError(f"$SPILLOVER names detectors not in the file: {', '.join(missing)}")
    events, pnn, off_by_one = fcsio.read_events(rec.path, spill if apply_comp else None)
    if off_by_one:
        log.warning("%s: the file's DATA end offset ($ENDDATA) is off by one byte (a known writer bug, e.g. the "
                    "CyTOF software and bead normaliser); read with ignore_offset_error=True", rec.sample_id)
    n_in = events.shape[0]
    if n_in == 0:
        raise ValueError("no events")
    col = {c: i for i, c in enumerate(pnn)}
    gate_summary = None
    rows = None  # event rows to write (None = all)
    if gating is not None:
        gc = gating["channels"]
        res = G.gate({c: events[:, col[c]] for c in gc.dna + gc.beads + gc.viability}, gc, gating["cofactors"],
                     mass=gating["mass"])
        gate_summary = res.summary()
        for w in res.warnings:
            log.warning("%s: %s", rec.sample_id, w)
        th = res.thresholds
        log.info("%s: pregate %d -> beads -%d, DNA -%d, dead -%d -> %d events (%.1f%% kept; %s)", rec.sample_id,
                 n_in, res.beads_removed, res.dna_removed, res.dead_removed, res.events_out,
                 100.0 * res.events_out / n_in, threshold_text(th))
        try:
            G.plot(res, gc, gating["labels"], rec.sample_id, gating["plot_dir"] / f"{rec.sample_id}.png",
                   seed=params["seed"])
            gate_summary["plot"] = f"gates/{rec.sample_id}.png"
        except Exception as exc:  # noqa: BLE001 - a plot must not fail the sample
            log.warning("%s: gate plot failed: %s", rec.sample_id, exc)
        rows = np.flatnonzero(res.keep)
        if len(rows) == 0:
            raise ValueError("no events left after pre-gating")
    n_kept = n_in if rows is None else len(rows)
    idx = subsample_index(n_kept, params["max_events"], params["seed"], rec.sample_id)
    if rows is not None:
        idx = rows if idx is None else rows[idx]
    chans = panel.channels
    out = np.empty((n_in if idx is None else len(idx), len(chans)), dtype=np.float32)
    for j, ch in enumerate(chans):
        x = events[:, col[ch.channel]] if idx is None else events[idx, col[ch.channel]]
        if ch.type == "marker":
            x = np.arcsinh(x / ch.cofactor)
        out[:, j] = x
    del events
    n_out = out.shape[0]
    n_nonfinite = int((~np.isfinite(out)).sum())
    if n_nonfinite:
        log.warning("%s: %d non-finite values in the output", rec.sample_id, n_nonfinite)

    cof_text = ",".join(f"{c.channel}={c.cofactor:g}" for c in chans if c.type == "marker")
    distinct = sorted({c.cofactor for c in chans if c.type == "marker"})
    prep = (f"transform=arcsinh;cofactor={'+'.join(f'{v:g}' for v in distinct) or 'none'};"
            f"compensation={'applied' if apply_comp else 'none'};events={n_out}"
            + (f";pregate=beads+dna+live;pregated_from={n_in}" if gate_summary is not None else "")
            + (f";subsampled_from={n_kept}" if n_out < n_kept else "") + f";app={APP_NAME}@{app_version}")
    keywords = {k: text[k] for k in KEEP_KEYWORDS if text.get(k)}
    keywords["fil"] = f"{rec.sample_id}.fcs"
    keywords.update({
        "PENNSIEVE_PREP": prep,
        "PENNSIEVE_SAMPLE_ID": rec.sample_id,
        "PENNSIEVE_SOURCE_FILE": rec.file,
        "PENNSIEVE_COFACTORS": cof_text or "none",
    })
    fcsio.write_fcs(fcs_dir / f"{rec.sample_id}.fcs", out,
                    [c.channel for c in chans], [c.marker if c.marker != c.channel else "" for c in chans],
                    keywords, ranges(text, header.pnn, chans, out))
    result = {
        "events_in": n_in, "events_out": n_out, "subsampled": n_out < n_kept,
        "ignore_offset_error": off_by_one,
        "compensation": {"applied": apply_comp, "reason": comp_reason,
                         "detectors": spill[0] if (spill and apply_comp) else []},
        "nonfinite_values": n_nonfinite, "seconds": round(time.monotonic() - t0, 2),
    }
    if gate_summary is not None:
        result["pregate"] = gate_summary
    return result


def write_csv(path: Path, header: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    started = datetime.now(timezone.utc)
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2

    params = env_params()
    try:
        if params["COMPENSATE"] not in ("auto", "true", "false"):
            raise ValueError(f"COMPENSATE must be auto, true or false (got {params['COMPENSATE']!r})")
        if params["PREGATE"] not in ("auto", "true", "false"):
            raise ValueError(f"PREGATE must be auto, true or false (got {params['PREGATE']!r})")
        cofactor_param = None
        if params["COFACTOR"].lower() != "auto":
            cofactor_param = float(params["COFACTOR"])
            if cofactor_param <= 0:
                raise ValueError("COFACTOR must be > 0 or 'auto'")
        cofactors = P.parse_json_object(params["COFACTORS"], "COFACTORS")
        overrides = P.parse_json_object(params["PANEL_OVERRIDES"], "PANEL_OVERRIDES")
        params["max_events"] = int(params["MAX_EVENTS_PER_SAMPLE"])
        params["seed"] = int(params["SEED"])
        if params["max_events"] < 0:
            raise ValueError("MAX_EVENTS_PER_SAMPLE must be >= 0")
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2

    fcs_dir, meta_dir = output_dir / "fcs", output_dir / "metadata"
    fcs_dir.mkdir(parents=True, exist_ok=True)
    meta_dir.mkdir(parents=True, exist_ok=True)
    app_version = version_of(REPO_ROOT)

    fcs_files, sheet_path = I.find_inputs(input_dir, params["SAMPLE_SHEET"])
    rel = lambda p: str(p.relative_to(input_dir))  # noqa: E731
    log.info("input: %d FCS file(s), sample sheet: %s", len(fcs_files), rel(sheet_path) if sheet_path else "none")
    if not fcs_files:
        listing = [rel(p) for p in I.walk_files(input_dir)]
        log.error("no *.fcs files under INPUT_DIR; it holds %d file(s): %s", len(listing),
                  ", ".join(listing[:40]) + (" …" if len(listing) > 40 else ""))
        return 1
    try:
        sheet = I.read_sheet(sheet_path)
    except (ValueError, OSError, UnicodeDecodeError) as exc:
        log.error("cannot use the sample sheet %s: %s", sheet_path, exc)
        return 1

    # 1. headers and sample identities
    skipped: list[dict] = []
    candidates: list[tuple[I.SampleRecord, P.Header, str, str, int]] = []
    for path in fcs_files:
        try:
            text, pnn, pns, tot, version = fcsio.read_header(path)
        except Exception as exc:  # noqa: BLE001
            log.error("%s: unreadable FCS header, skipped: %s", rel(path), exc)
            skipped.append({"file": rel(path), "reason": f"unreadable: {type(exc).__name__}: {exc}"})
            continue
        rec = I.make_record(path, text, sheet)
        if tot <= 0:
            log.error("%s: no events ($TOT=%s), skipped", rel(path), text.get("tot"))
            skipped.append({"file": rel(path), "sample_id": rec.sample_id, "reason": "empty: no events"})
            continue
        for note in rec.notes:
            log.warning("%s: %s", rel(path), note)
        cyt_type = C.cytometer_type(text, pnn)
        candidates.append((rec, P.Header(rec.sample_id, pnn, pns, text), cyt_type, version, tot))
    if sheet.path:
        found_names = {c[0].file.lower() for c in candidates} | {Path(s["file"]).name.lower() for s in skipped}
        for name in sheet.rows:
            if name.lower() not in found_names:
                log.warning("sample sheet row %s: no such FCS file under INPUT_DIR", name)

    ids = Counter(c[0].sample_id for c in candidates)
    dups = {s for s, k in ids.items() if k > 1}
    if dups:
        for rec, *_ in candidates:
            if rec.sample_id in dups:
                log.error("%s: sample_id %s is used by %d files, skipped (give each file its own sample_id)",
                          rel(rec.path), rec.sample_id, ids[rec.sample_id])
                skipped.append({"file": rel(rec.path), "sample_id": rec.sample_id, "reason": "duplicate sample_id"})
        candidates = [c for c in candidates if c[0].sample_id not in dups]
    if not candidates:
        log.error("no usable FCS file")
        return 1

    # 2. panel
    panel = P.build([c[1] for c in candidates])
    types = Counter(c[2] for c in candidates)
    main_type = types.most_common(1)[0][0]
    if len(types) > 1:
        log.warning("mixed cytometer types %s; cofactor default for %s", dict(types), main_type)
    default_cof = cofactor_param if cofactor_param is not None else C.DEFAULT_COFACTOR[main_type]
    try:
        for w in P.configure(panel, default_cof, cofactors, overrides):
            log.warning(w)
    except ValueError as exc:
        log.error("bad parameter: %s", exc)
        return 2
    if not panel.markers():
        log.error("no marker channel is common to all samples; nothing to transform")
        return 1
    dropped_channels = sorted({d["channel"] for ds in panel.dropped.values() for d in ds})
    log.info("panel: %d channel(s) kept (%d markers), %d dropped%s; cytometer %s; cofactor %g", len(panel.channels),
             len(panel.markers()), len(dropped_channels),
             (": " + ", ".join(dropped_channels)) if dropped_channels else "", main_type, default_cof)
    write_csv(meta_dir / "panel.csv", ["channel", "marker", "type", "cofactor", "clustering", "state"],
              [c.row() for c in panel.channels])

    # pre-gating plan (run level)
    gate_channels = G.find_channels(panel.channels)
    pregate_on, pregate_reason = G.plan(params["PREGATE"], main_type, gate_channels)
    gating = None
    if pregate_on:
        if main_type != "mass":
            log.warning("pre-gating %s data: the gates are tuned for mass cytometry; check the gate plots", main_type)
        if not gate_channels.dna:
            log.warning("pre-gating without a DNA channel: only %s",
                        " and ".join(k for k in ("beads", "viability") if getattr(gate_channels, k)))
        gating = {"channels": gate_channels, "cofactors": gating_cofactors(panel, gate_channels, main_type),
                  "labels": {c.channel: c.marker for c in panel.channels}, "plot_dir": meta_dir / "gates",
                  "mass": main_type == "mass"}
    log.info("pregate: %s", pregate_reason)

    # 3. one sample at a time
    sample_rows: list[dict] = []
    prep_samples: dict[str, Any] = {}
    prov_samples: list[dict] = []
    for rec, header, cyt_type, version, tot in candidates:
        text = header.text
        details = {
            "file": rec.file, "input": rel(rec.path), "in_sample_sheet": rec.in_sheet if sheet.path else None,
            "fcs_version": version, "cytometer": text.get("cyt"), "cytometer_serial": text.get("cytsn"),
            "cytometer_type": cyt_type, "creator": text.get("creator"), "date": text.get("date"),
            "batch": rec.batch, "condition": rec.condition,
            "channels_dropped": panel.dropped.get(rec.sample_id, []),
            "markers_renamed": panel.renamed.get(rec.sample_id, []),
            "cofactors": {c.marker: c.cofactor for c in panel.markers()},
            "notes": rec.notes,
        }
        try:
            result = process_sample(rec, header, cyt_type, panel, params, fcs_dir, app_version, gating)
        except Exception as exc:  # noqa: BLE001 - one sample's failure must not sink the run
            log.error("%s: FAILED: %s\n%s", rel(rec.path), exc, traceback.format_exc())
            skipped.append({"file": rel(rec.path), "sample_id": rec.sample_id, "reason": f"{type(exc).__name__}: {exc}"})
            prov_samples.append({"sample_id": rec.sample_id, "input": rel(rec.path), "status": "failed",
                                 "error": str(exc)})
            continue
        details.update(result)
        prep_samples[rec.sample_id] = details
        sample_rows.append({"sample_id": rec.sample_id, "file": rec.file, "batch": rec.batch,
                            "condition": rec.condition, "n_events": result["events_out"], **rec.extra})
        prov_samples.append({"sample_id": rec.sample_id, "input": rel(rec.path), "status": "ok",
                             "output": f"fcs/{rec.sample_id}.fcs", "events_in": result["events_in"],
                             "events_out": result["events_out"], "seconds": result["seconds"]})
        log.info("%s: %s -> %d events%s, %d channels (%d markers), compensation %s, %d dropped, %.1fs",
                 rec.sample_id, rec.file, result["events_out"],
                 (f" (of {result['pregate']['events_out']} pre-gated from {result['events_in']})"
                  if "pregate" in result else f" (of {result['events_in']})") if result["subsampled"]
                 else (f" (pre-gated from {result['events_in']})" if "pregate" in result else ""), len(panel.channels),
                 len(panel.markers()), "applied" if result["compensation"]["applied"] else "none",
                 len(details["channels_dropped"]), result["seconds"])

    write_csv(meta_dir / "samples.csv", list(I.REQUIRED_COLUMNS) + sheet.columns, sample_rows)
    if gating is not None:
        write_csv(meta_dir / "gating.csv", ["sample_id", "events_in", "beads_removed", "dna_removed", "dead_removed",
                                            "events_out", "thresholds"],
                  [{"sample_id": sid, **{k: d["pregate"][k] for k in ("events_in", "beads_removed", "dna_removed",
                                                                      "dead_removed", "events_out")},
                    "thresholds": threshold_text(d["pregate"]["thresholds"])}
                   for sid, d in prep_samples.items() if "pregate" in d])
    finished = datetime.now(timezone.utc)
    prep_doc = {
        "app": APP_NAME, "appVersion": app_version, "generatedAt": finished.isoformat(timespec="seconds"),
        "cytometerType": main_type, "defaultCofactor": default_cof,
        "panel": [c.row() | {"spellings": c.spellings} for c in panel.channels],
        "channelsDropped": dropped_channels,
        "pregate": {"mode": params["PREGATE"], "applied": pregate_on, "reason": pregate_reason,
                    "channels": gate_channels.as_dict(),
                    **({"cofactors": gating["cofactors"], "method": G.METHOD, "settings": G.settings()}
                       if gating else {})},
        "samples": prep_samples, "skipped": skipped,
    }
    with open(meta_dir / "prep.json", "w") as fh:
        json.dump(prep_doc, fh, indent=2)

    provenance = {
        "app": APP_NAME, "appVersion": app_version, "tools": tool_versions(),
        "startedAt": started.isoformat(timespec="seconds"), "finishedAt": finished.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2),
        "parameters": {k: v for k, v in params.items() if k.isupper()},
        "inputsFound": {"fcs": [rel(p) for p in fcs_files], "sampleSheet": rel(sheet_path) if sheet_path else None},
        "samples": prov_samples, "skipped": skipped,
    }
    for d in (fcs_dir, meta_dir):
        with open(d / "provenance.json", "w") as fh:
            json.dump(provenance, fh, indent=2)

    log.info("%d of %d file(s) prepared, %d skipped, %.1fs; wrote %s and %s", len(sample_rows), len(fcs_files),
             len(skipped), time.monotonic() - t_start, fcs_dir, meta_dir)
    if not sample_rows:
        log.error("every sample failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
