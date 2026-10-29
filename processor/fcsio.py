"""
Reading (FlowKit/FlowIO) and writing (FlowIO) FCS files, and compensation.
"""
from __future__ import annotations

from array import array
from pathlib import Path

import numpy as np


def read_header(path: Path) -> tuple[dict, list[str], list[str], int, str]:
    """(text keywords, $PnN list, $PnS list, event count, FCS version) without reading events."""
    import flowio

    fd = flowio.FlowData(str(path), only_text=True)
    text = fd.text
    n = int(text["par"])
    pnn = [text.get(f"p{i}n", f"P{i}").strip() for i in range(1, n + 1)]
    # '/' is the delimiter of the TEXT segment we write: FlowIO escapes it as '//' (valid FCS),
    # which flowCore's read.FCS takes for an empty value and refuses. Marker names keep their
    # meaning with '|' ('CD154/CD40L' -> 'CD154|CD40L'; marker_aliases splits on both).
    pns = [fcs_safe((text.get(f"p{i}s") or "").strip()) for i in range(1, n + 1)]
    if len(set(pnn)) != len(pnn):
        raise ValueError("duplicate $PnN channel names")
    return text, pnn, pns, int(text.get("tot", 0) or 0), fd.version


def parse_spill(text: dict) -> tuple[list[str], np.ndarray] | None:
    """The $SPILLOVER / $SPILL / SPILL matrix as (detectors, matrix), or None."""
    raw = text.get("spillover") or text.get("spill") or text.get("comp")
    if not raw or not raw.strip():
        return None
    parts = [p.strip() for p in raw.split(",")]
    try:
        n = int(float(parts[0]))
    except ValueError:
        return None
    if n <= 0 or len(parts) < 1 + n + n * n:
        return None
    names = parts[1:1 + n]
    values = np.array([float(v) for v in parts[1 + n:1 + n + n * n]]).reshape(n, n)
    return names, values


def is_identity(matrix: np.ndarray, atol: float = 1e-4) -> bool:
    """True when the spillover matrix does nothing (spectral exports write identity, give or take 1e-6)."""
    return bool(np.allclose(matrix, np.eye(matrix.shape[0]), atol=atol))


def is_offset_off_by_one(exc: BaseException) -> bool:
    """FlowIO's error for files whose $ENDDATA points one byte past the DATA segment."""
    from flowio.exceptions import FCSParsingError

    return isinstance(exc, FCSParsingError) and "off by 1" in str(exc)


def open_sample(path: Path):
    """
    FlowKit Sample of all events, and whether ignore_offset_error was needed.

    Some writers (the CyTOF software and the Nolan lab bead normaliser among them)
    record $ENDDATA one byte past the DATA segment; FlowIO refuses those files unless
    told to ignore that one error. Only that error is retried; any other offset
    problem still fails the file.
    """
    import warnings

    import flowkit as fk

    try:
        return fk.Sample(str(path), subsample=0), False
    except Exception as exc:  # noqa: BLE001
        if not is_offset_off_by_one(exc):
            raise
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # FlowIO's own warning; the caller logs one line per file
        return fk.Sample(str(path), subsample=0, ignore_offset_error=True), True


def read_events(path: Path, spill: tuple[list[str], np.ndarray] | None):
    """
    All events as float64 (events x channels), via FlowKit's Sample, which applies
    $PnG/$PnE scaling for integer data. With `spill`, the events of the detectors in
    the matrix are compensated (FlowKit Matrix: events @ inverse(spillover)).
    Returns (events, $PnN list, offset_off_by_one: bool).
    """
    import flowkit as fk

    sample, off_by_one = open_sample(path)
    if spill is not None:
        names, values = spill
        matrix = fk.Matrix(values, detectors=names)
        sample.apply_compensation(matrix)
        events = sample.get_events(source="comp")
    else:
        events = sample.get_events(source="raw")
    # FlowIO's preprocessing (gain and log decoding, needed for integer data) also
    # multiplies the Time channel by $TIMESTEP. Downstream tools (flowCore) expect
    # Time as acquired with $TIMESTEP alongside, so undo that one step.
    ti = getattr(sample, "time_index", None)
    timestep = (sample.metadata.get("timestep") or "").strip()
    if ti is not None and timestep:
        try:
            step = float(timestep)
        except ValueError:
            step = 0.0
        if step > 0:
            events = events.copy() if not events.flags.writeable else events
            events[:, ti] = events[:, ti] / step
    return events, list(sample.pnn_labels), off_by_one


def fcs_safe(value: str) -> str:
    """A keyword value without the TEXT delimiter '/', so every reader parses it (see read_header)."""
    return value.replace("/", "|")


def write_fcs(path: Path, events: np.ndarray, pnn: list[str], pns: list[str], keywords: dict[str, str],
              ranges: list[str]) -> None:
    """Write FCS 3.1, float32, little-endian, with $PnN/$PnS/$PnR and extra keywords (FlowIO create_fcs)."""
    import flowio

    data = array("f")
    data.frombytes(np.ascontiguousarray(events, dtype="<f4").tobytes())
    meta = {k: fcs_safe(str(v)) for k, v in keywords.items() if v not in (None, "")}
    pns = [fcs_safe(p) for p in pns]
    for i, r in enumerate(ranges, start=1):
        if r:
            meta[f"p{i}r"] = r
    tmp = path.with_suffix(".fcs.part")
    with open(tmp, "wb") as fh:
        flowio.create_fcs(fh, data, pnn, opt_channel_names=pns, metadata_dict=meta)
    tmp.replace(path)
