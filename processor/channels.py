"""
Channel classification, cytometer detection, default cofactors and the
clustering/state flags of panel.csv.

Everything here works on FCS TEXT keywords only (lower-cased keys without
the leading '$', as FlowIO returns them), so it can run on headers before
any events are read.
"""
from __future__ import annotations

import re
from typing import Iterable

# --- channel types -----------------------------------------------------------

SCATTER_RE = re.compile(r"^(FSC|SSC)", re.IGNORECASE)
TIME_RE = re.compile(r"^time$|^time[\s_-]", re.IGNORECASE)
# Aurora / ID7000 autofluorescence "channels" written by unmixing (AF-A, AF 1-A, ...).
AUTOFLUORESCENCE_RE = re.compile(r"^AF(\s*\d+)?(-[AHW])?$|autofluor", re.IGNORECASE)
# Mass cytometry technical parameters (CyTOF/Helios/XT): not markers.
CYTOF_TECHNICAL = {"event_length", "center", "offset", "width", "residual", "beaddist", "sample_id"}
# A mass channel name: Yb171Di, 171Yb, (Yb171)Dd, Yb171 ...
METAL_RE = re.compile(r"^\(?([A-Z][a-z]?)(\d{2,3})\)?(Di|Dd)?$|^(\d{2,3})([A-Z][a-z]?)(Di|Dd)?$")
# Fluorochrome names as they appear in $PnN of conventional and spectral cytometers.
FLUOR_RE = re.compile(
    r"^("
    r"FITC|PE|APC|PerCP|PacBlue|Pacific|BV\d|BUV\d|BB\d|BYG\d|BV|BUV|"
    r"Alexa|AF\d{3}|A\d{3}|eFluor|eF\d|Spark|cFluor|RB\d|NovaFluor|StarBright|Super\s*Bright|SB\d|"
    r"Real\s*Blue|RY\d|Brilliant|Vio|Qdot|Cy\d|DyLight|KIRAVIA|iFluor|NFluor|"
    r"DAPI|Hoechst|7-?AAD|PI\b|Zombie|LIVE|L/D|Viability|GFP|EGFP|YFP|RFP|BFP|CFP|mCherry|tdTomato|"
    r"(UV|V|B|YG|R|Y|G)\d{1,2}(-[AHW])?$"
    r")",
    re.IGNORECASE,
)


def norm(name: str | None) -> str:
    """Case- and punctuation-insensitive key: 'HLA-DR' -> 'HLADR', 'PD-1' -> 'PD1'."""
    return re.sub(r"[^A-Z0-9]", "", (name or "").upper())


def is_metal(channel: str) -> bool:
    return bool(METAL_RE.match(channel.strip()))


def classify(channel: str, marker: str | None) -> str:
    """'scatter' | 'time' | 'marker' | 'other' for one channel ($PnN, $PnS)."""
    ch = channel.strip()
    if TIME_RE.match(ch):
        return "time"
    if SCATTER_RE.match(ch):
        return "scatter"
    if AUTOFLUORESCENCE_RE.search(ch) or (marker and AUTOFLUORESCENCE_RE.search(marker.strip())):
        return "other"
    if norm(ch).lower() in {norm(t).lower() for t in CYTOF_TECHNICAL}:
        return "other"
    if marker and marker.strip():
        return "marker"
    if is_metal(ch) or FLUOR_RE.match(ch):
        return "marker"
    return "other"


# --- cytometer type ----------------------------------------------------------

SPECTRAL_CYT = ("aurora", "northern lights", "id7000", "sp6800", "sa3800", "facsdiscover", "a5 se", "symphony a5se",
                "bigfoot", "mosaic", "spectral")
CYTOF_CYT = ("cytof", "helios", "fluidigm", "standard biotools", "cytof xt", "dvs")


def cytometer_type(text: dict, channels: Iterable[str]) -> str:
    """'mass' (CyTOF), 'spectral' or 'flow' (conventional), from $CYT/$CREATOR and channel names."""
    cyt = (text.get("cyt") or "").lower()
    creator = (text.get("creator") or "").lower()
    channels = list(channels)
    if any(k in cyt for k in CYTOF_CYT):
        return "mass"
    metals = sum(1 for c in channels if is_metal(c))
    if channels and metals >= max(3, len(channels) // 3):
        return "mass"
    if any(k in cyt for k in SPECTRAL_CYT) or "spectroflo" in creator:
        return "spectral"
    return "flow"


# arcsinh(x / cofactor) defaults per cytometer type; see README "Transformation".
DEFAULT_COFACTOR = {"mass": 5.0, "spectral": 6000.0, "flow": 150.0}


# --- clustering / state flags --------------------------------------------------

# Functional / activation ("state") markers: state=true, clustering=false. They
# change with stimulation, so clustering on them would split one cell type into
# condition-specific clusters and confound differential abundance; diffcyt tests
# them within clusters instead (differential state). See README "Panel flags".
STATE_MARKERS = {
    # activation
    "CD25", "IL2RA", "CD69", "HLADR", "CD38", "CD71", "CD137", "41BB", "CD134", "OX40", "CD154", "CD40L",
    "ICOS", "CD278", "GITR", "CD357",
    # inhibitory receptors / exhaustion
    "PD1", "CD279", "TIGIT", "CTLA4", "CD152", "LAG3", "CD223", "TIM3", "CD366", "CD39", "2B4", "CD244",
    # cytotoxic granules and degranulation
    "GZMA", "GZMB", "GZMK", "GZMH", "GRANZYMEA", "GRANZYMEB", "GRANZYMEK", "PFN", "PRF1", "PERFORIN",
    "GNLY", "GRANULYSIN", "CD107A", "LAMP1",
    # death receptor, proliferation
    "CD95", "FAS", "KI67", "MKI67",
    # cytokines (intracellular staining)
    "IFNG", "IFNGAMMA", "TNF", "TNFA", "TNFALPHA", "IL2", "IL4", "IL5", "IL6", "IL10", "IL13", "IL17",
    "IL17A", "IL21", "IL22", "GMCSF", "MIP1B", "CCL4",
}
# Signalling (phospho-epitopes, e.g. pSTAT5, pERK1_2, pp38) are state markers too.
PHOSPHO_RE = re.compile(r"^P(STAT|ERK|P38|S6|AKT|CREB|NFKB|ZAP|SLP|PLC|BTK|SYK|LCK|MAPK|JNK|HH3|RB)", re.IGNORECASE)

# Gating, viability, DNA, beads and barcodes: kept and transformed, but neither
# clustering nor state markers (they carry no cell-type information after gating).
GATING_MARKERS = {"CD45", "DNA", "DNA1", "DNA2", "IR191", "IR193", "CISPLATIN", "VIABILITY",
                  "LIVEDEAD", "LD", "ZOMBIE", "BEAD", "BEADS", "BC", "BARCODE"}
GATING_RE = re.compile(r"(DNA|VIABILITY|LIVE.?DEAD|ZOMBIE|CISPLATIN|BEAD|BARCODE|^BC\d|^CD45$)", re.IGNORECASE)


def marker_key(marker: str) -> str:
    """Key for the built-in lists: strips mass tags ('171Yb_CD3', 'CD3_Yb171') and punctuation."""
    m = (marker or "").strip()
    m = re.sub(r"^\(?\d{2,3}[A-Z][a-z]?\)?[_\s-]+", "", m)          # 171Yb_CD3
    m = re.sub(r"[_\s-]+\(?[A-Z][a-z]?\d{2,3}(Di|Dd)?\)?$", "", m)    # CD3_Yb171Di
    m = re.sub(r"[_\s-]+\(?\d{2,3}[A-Z][a-z]?\)?$", "", m)            # CD3_171Yb
    return norm(m)


def default_flags(marker: str) -> tuple[bool, bool]:
    """(clustering, state) for a marker channel."""
    key = marker_key(marker)
    if key in GATING_MARKERS or GATING_RE.search(marker or ""):
        return False, False
    if key in STATE_MARKERS or PHOSPHO_RE.match(key):
        return False, True
    return True, False
