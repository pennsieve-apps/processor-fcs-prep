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
    # immune checkpoint ligand, induced by IFN-gamma and activation
    "PDL1", "CD274", "PDL2", "CD273",
    # DNA damage (gamma-H2AX: H2A.X phosphorylated at S139)
    "H2AX", "PH2AX", "GH2AX", "GAMMAH2AX", "YH2AX",
}
# Signalling (phospho-epitopes, e.g. pSTAT5, pERK1_2, pp38, pH2A.X) are state markers too.
PHOSPHO_RE = re.compile(
    r"^P(STAT|ERK|P38|S6|AKT|CREB|NFKB|ZAP|SLP|PLC|BTK|SYK|LCK|MAPK|JNK|HH3|RB|H2AX|4EBP|IKB|SRC|CBL|LAT|BLNK|"
    r"CD3Z|MEK|GSK|MTOR|AMPK|TBK|IRF|SMAD|CDK|CHK|ATM|ATR|P53|RPS6)",
    re.IGNORECASE)

# Gating, viability, DNA, beads and barcodes: kept and transformed, but neither
# clustering nor state markers (they carry no cell-type information after gating).
GATING_MARKERS = {"CD45", "BC", "BARCODE"}
DNA_RE = re.compile(r"DNA|IRIDIUM|^IR19[13]|INTERCALATORIR|IRINTERCALATOR", re.IGNORECASE)
BEAD_RE = re.compile(r"BEAD", re.IGNORECASE)
BARCODE_RE = re.compile(r"BARCODE|^BC\d*$", re.IGNORECASE)
# Viability: whole words (Dead, LD, L/D, Live, Viab..., Cisplatin, Zombie, 7-AAD, DAPI), so
# that LDLR, LDHA, OLIG2 or "Deadbox" style names are not caught.
VIABILITY_WORDS = {"DEAD", "LD", "LIVE", "LIVEDEAD", "DEADLIVE", "VIABILITY", "VIAB", "VIA", "CISPLATIN", "CISPT",
                   "CIS", "ZOMBIE", "7AAD", "DAPI", "PI", "PROPIDIUM", "FVS", "FVD", "VIOBILITY"}
VIABILITY_RE = re.compile(r"^(VIAB|CISPLAT|ZOMBIE|LIVEDEAD|LIVE?DEAD)", re.IGNORECASE)
# Unlabelled channels of the usual viability reagents: cisplatin (195Pt, 198Pt) and Rh103 intercalator.
VIABILITY_METALS = {("Pt", 195), ("Pt", 198), ("Rh", 103)}
METAL_PARTS = re.compile(r"^\(?([A-Z][a-z]?)(\d{2,3})\)?(?:Di|Dd)?$|^(\d{2,3})([A-Z][a-z]?)(?:Di|Dd)?$")


def metal_of(channel: str | None) -> tuple[str, int] | None:
    m = METAL_PARTS.match((channel or "").strip())
    if not m:
        return None
    return (m.group(1), int(m.group(2))) if m.group(1) else (m.group(4), int(m.group(3)))


# Element symbols of the metal tags used in mass cytometry (lanthanides, In, Pt, Ir, Bi, Y, Rh, Cd, Pd ...).
ELEMENTS = ("La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|In|Pt|Ir|Bi|Rh|Cd|Pd|Ba|Cs|Au|Ag|Os|Sn|Te|Xe|"
            "Y|I")
_TAG = rf"(?:\d{{2,3}}(?:{ELEMENTS})|(?:{ELEMENTS})\d{{2,3}}(?:Di|Dd)?)"


def tag_metal(marker: str | None) -> tuple[str, int] | None:
    """The metal of a marker's mass tag: '103Rh_Intercalator' -> ('Rh', 103)."""
    m = re.search(rf"(?:^|[_\s(-])\(?(\d{{2,3}})({ELEMENTS})(?![a-z])|(?:^|[_\s(-])\(?({ELEMENTS})(\d{{2,3}})",
                  marker or "")
    if not m:
        return None
    return (m.group(2), int(m.group(1))) if m.group(1) else (m.group(3), int(m.group(4)))


def strip_mass_tag(marker: str) -> str:
    """'171Yb_CD3' / 'Yb171Di_CD3' / 'CD3_Yb171Di' / 'CD3 (171Yb)' -> 'CD3' (punctuation otherwise kept)."""
    m = (marker or "").strip()
    m = re.sub(rf"^\(?{_TAG}\)?[_\s-]+", "", m)        # 171Yb_CD3, Yb171Di_CD3
    m = re.sub(rf"[_\s-]+\(?{_TAG}\)?$", "", m)        # CD3_Yb171Di, CD3_171Yb
    m = re.sub(rf"\s*\({_TAG}\)$", "", m)              # CD3(171Yb)
    return m.strip() or (marker or "").strip()


def marker_key(marker: str) -> str:
    """Key for the built-in lists: strips mass tags ('171Yb_CD3', 'CD3_Yb171') and punctuation."""
    return norm(strip_mass_tag(marker))


def marker_aliases(marker: str) -> list[str]:
    """
    Keys of a marker and its aliases: 'CD154 (CD40L)', 'CD154/CD40L' and Fluidigm's
    'CD154__CD40L_' (spaces and brackets written as '_') give ['CD154CD40L', 'CD154', 'CD40L'].
    The first entry is the whole name, the second the primary name.
    """
    m = strip_mass_tag(marker)
    keys = [norm(m)]
    alt = re.sub(r"_{2,}|[()\[\]/]", "|", m)
    alt = re.sub(r"_([^_|]+)_$", r"|\1", alt)  # 'CD154_CD40L_' -> 'CD154|CD40L'
    parts = [norm(p) for p in alt.split("|")]
    if len([p for p in parts if p]) > 1:
        keys += [p for p in parts if p]
    return keys


def _words(marker: str) -> list[str]:
    return [w for w in re.split(r"[^A-Z0-9]+", strip_mass_tag(marker).upper()) if w]


def is_viability(marker: str, channel: str | None = None) -> bool:
    """Dead/live/viability stains, cisplatin and the Rh103 intercalator (never a real marker name)."""
    name = strip_mass_tag(marker)
    words = _words(marker)
    key = norm(name)
    if any(w in VIABILITY_WORDS for w in words) or key in {"LD", "LIVEDEAD", "DEADLIVE"} or VIABILITY_RE.match(key):
        return True
    metal = metal_of(channel) or metal_of(marker) or tag_metal(marker)
    if "INTERCALATOR" in key and (metal and metal[0] == "Rh" or "RH" in words or "RH103" in key or "103RH" in key):
        return True
    if metal in VIABILITY_METALS and (not marker or key == norm(channel or "") or key in {"PT", "RH", "CISPLATIN"}):
        return True
    return False


def gating_role(marker: str, channel: str | None = None) -> str | None:
    """'dna' | 'beads' | 'viability' | 'barcode' | 'gating' (CD45) or None for a real marker."""
    key = marker_key(marker)
    if is_viability(marker, channel):
        return "viability"
    metal = metal_of(channel)
    if DNA_RE.search(key) or (metal and metal[0] == "Ir" and (not marker or key == norm(channel or ""))):
        return "dna"
    if "INTERCALATOR" in key and ((metal and metal[0] == "Ir") or (tag_metal(marker) or ("",))[0] == "Ir"):
        return "dna"
    if BEAD_RE.search(key):
        return "beads"
    if BARCODE_RE.search(key):
        return "barcode"
    if key in GATING_MARKERS:
        return "gating"
    return None


def is_state(marker: str) -> bool:
    """A state marker when the marker or any of its aliases is one (CD154 (CD40L), pH2A.X ...)."""
    return any(k in STATE_MARKERS or PHOSPHO_RE.match(k) for k in marker_aliases(marker))


def default_flags(marker: str, channel: str | None = None) -> tuple[bool, bool]:
    """(clustering, state) for a marker channel."""
    if gating_role(marker, channel):
        return False, False
    if is_state(marker):
        return False, True
    return True, False
