"""
Finding the FCS files and the sample sheet under INPUT_DIR, and building one
sample record per file (sample_id, batch, condition, extra sheet columns).
"""
from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

log = logging.getLogger("processor-fcs-prep")

SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")
REQUIRED_COLUMNS = ("sample_id", "file", "batch", "condition", "n_events")
NA = "NA"


def walk_files(root: Path) -> list[Path]:
    """Every file under root, following symlinked folders and files (INPUT_DIR entries are symlinks)."""
    found: list[Path] = []
    seen_dirs: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen_dirs:  # symlink loop
            dirnames[:] = []
            continue
        seen_dirs.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        for name in sorted(filenames):
            if name.startswith("."):
                continue
            path = Path(dirpath) / name
            if path.is_file():
                found.append(path)
    return found


def find_inputs(input_dir: Path, sheet_name: str) -> tuple[list[Path], Path | None]:
    """(*.fcs files, the sample sheet or None). Several sheets with the name: the shallowest wins."""
    files = walk_files(input_dir)
    fcs = [p for p in files if p.suffix.lower() == ".fcs"]
    sheets = [p for p in files if p.name == sheet_name] if sheet_name else []
    sheet = None
    if sheets:
        sheets.sort(key=lambda p: (len(p.relative_to(input_dir).parts), str(p)))
        sheet = sheets[0]
        if len(sheets) > 1:
            log.warning("%d files named %s; using %s", len(sheets), sheet_name, sheet.relative_to(input_dir))
    return fcs, sheet


def sanitise(stem: str) -> str:
    """A sample_id from a file stem: anything outside [A-Za-z0-9_-] becomes '_'."""
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", stem).strip("_")
    s = re.sub(r"_+", "_", s)
    return s or "sample"


def batch_from_date(date: str | None) -> str:
    """$DATE (FCS 3.x 'dd-mmm-yyyy', and a few variants seen in the wild) as an ISO date."""
    if not date or not date.strip():
        return NA
    d = date.strip()
    for fmt in ("%d-%b-%Y", "%d-%B-%Y", "%Y-%m-%d", "%d-%b-%y", "%m/%d/%Y", "%d.%m.%Y", "%Y/%m/%d", "%d %b %Y"):
        try:
            return datetime.strptime(d, fmt).date().isoformat()
        except ValueError:
            pass
    return sanitise(d)


@dataclass
class SampleSheet:
    path: Path | None = None
    columns: list[str] = field(default_factory=list)      # sheet columns other than the required ones
    rows: dict[str, dict[str, str]] = field(default_factory=dict)  # file name -> row

    def lookup(self, file_name: str) -> dict[str, str] | None:
        if file_name in self.rows:
            return self.rows[file_name]
        low = {k.lower(): v for k, v in self.rows.items()}
        return low.get(file_name.lower())


def read_sheet(path: Path | None) -> SampleSheet:
    """The sample sheet; needs a `file` column. Raises ValueError when it cannot be used."""
    if path is None:
        return SampleSheet()
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(fh, dialect=dialect)
        header = [h.strip() for h in (reader.fieldnames or [])]
        if "file" not in header:
            raise ValueError(f"{path.name} has no 'file' column (columns: {', '.join(header)})")
        rows: dict[str, dict[str, str]] = {}
        for raw in reader:
            row = {(k or "").strip(): (v or "").strip() for k, v in raw.items() if k is not None}
            name = os.path.basename(row.get("file", ""))
            if not name:
                continue
            if name in rows:
                log.warning("sample sheet lists %s twice; using the first row", name)
                continue
            rows[name] = row
    extra = [h for h in header if h not in REQUIRED_COLUMNS]
    return SampleSheet(path=path, columns=extra, rows=rows)


@dataclass
class SampleRecord:
    path: Path
    file: str
    sample_id: str
    batch: str
    condition: str
    extra: dict[str, str]
    in_sheet: bool
    notes: list[str] = field(default_factory=list)


def make_record(path: Path, text: dict, sheet: SampleSheet) -> SampleRecord:
    """Sample identity for one file: the sheet row when there is one, defaults otherwise."""
    row = sheet.lookup(path.name) if sheet.path else None
    notes: list[str] = []
    default_id = sanitise(path.stem)
    default_batch = batch_from_date(text.get("date"))
    if row is None:
        if sheet.path:
            notes.append(f"not in {sheet.path.name}: defaults used")
        return SampleRecord(path, path.name, default_id, default_batch, NA, {c: "" for c in sheet.columns}, False, notes)
    sample_id = row.get("sample_id") or default_id
    if not SAMPLE_ID_RE.match(sample_id):
        fixed = sanitise(sample_id)
        notes.append(f"sample_id {sample_id!r} sanitised to {fixed!r}")
        sample_id = fixed
    batch = row.get("batch") or default_batch
    condition = row.get("condition") or NA
    return SampleRecord(path, path.name, sample_id, batch, condition, {c: row.get(c, "") for c in sheet.columns}, True,
                        notes)
