"""
Pre-gating of raw mass cytometry (CyTOF) data: beads out, DNA+ singlet cells,
live cells (PREGATE).

The standard CyTOF clean-up sequence (Bagwell et al. 2020, Cytometry A,
doi:10.1002/cyto.a.23926; the manual gating of Kowli et al. on Ir DNA, 140Ce
beads and 115In live/dead), automated with deterministic, data-driven
thresholds on arcsinh-transformed values (cofactor 5 for mass cytometry):

  a) beads:  events high in the bead channel(s) are removed. Low-DNA ones are
             free beads; DNA+ ones are bead-cell doublets, removed as well
             because EQ bead metals (151Eu, 153Eu, 165Ho, 175Lu) would show as
             false marker signal. The threshold is the density valley between
             the cells (near 0) and the bead peak.
  b) DNA:    events in the main Ir191/Ir193 population. In rotated coordinates
             (s = mean DNA along the diagonal, r = DNA2 - DNA1 across it), the
             main mode of s is found on a smoothed histogram; the low cut
             (debris, DNA-) and the high cut (aggregates) are set from the
             half-widths of that mode, with bounds; the high cut is generous,
             since doublets are left to cleanet. Events far off the diagonal
             (the two Ir isotopes disagree) are removed too (MAD bound).
  c) live:   events bright on the viability channel (cisplatin, 115In, 103Rh,
             live/dead dyes) are removed. Dead cells are the brightest
             population: the threshold is the brightest valley (else shoulder)
             of the density that holds a plausible dead fraction, else Tukey's
             upper fence (outliers only); it is never below the median + 0.75
             (nor, for mass cytometry, arcsinh 3: the live cells' dim
             background) and never removes more than MAX_DEAD_FRACTION.

Nothing here is random; the gate plots draw a seeded subsample.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from processor import channels as C

log = logging.getLogger("processor-fcs-prep")

# Gate bounds (arcsinh units; with cofactor 5, 0.69 = a factor 2 in signal).
DNA_LOW_SIGMAS, DNA_LOW_MIN, DNA_LOW_MAX = 4.0, 0.6, 1.5      # low cut: mode - clip(4 sigma_left, 0.6, 1.5)
DNA_HIGH_SIGMAS, DNA_HIGH_MIN, DNA_HIGH_MAX = 5.0, 1.0, 2.0   # high cut: mode + clip(5 sigma_right, 1.0, 2.0)
DNA_DIAG_MADS, DNA_DIAG_MIN = 6.0, 0.5                        # |r - median r| <= max(6 MAD, 0.5)
BEAD_MIN_SIGNAL = 3.0          # a bead population is searched above this (arcsinh)
BEAD_THRESHOLD_FLOOR = 2.0     # the bead threshold is never below this
DEAD_MIN_GAP = 0.75            # the viability threshold is at least the median + 0.75
MIN_DEAD_FRACTION = 0.005      # a dead-cell population holds at least this much of the gate's input
DEAD_MIN_SIGNAL = 3.0          # mass cytometry: no viability threshold below arcsinh(x/5) = 3 (~50 counts)
DEAD_SLOPE_RISE = 0.5          # rise of the log-density slope (within DEAD_SLOPE_WINDOW) that marks
DEAD_SLOPE_WINDOW = 1.5        # a brighter (dead-cell) mode or shoulder
MAX_DEAD_FRACTION = 0.30       # ... and the viability gate never removes more than this
WARN_FRACTION = 0.60           # a gate removing more than this of its input is flagged
HWHM_TO_SIGMA = 1.0 / np.sqrt(2 * np.log(2))
BIN = 0.02

METHOD = ("beads (bead-channel valley; free beads and bead-cell doublets), then DNA+ (main mode of mean Ir "
          "with half-width bounds, off-diagonal MAD bound), then live (brightest valley/shoulder of the viability "
          "density, else Tukey upper fence; floor and cap); arcsinh values; after Bagwell et al. 2020")


def settings() -> dict:
    """The gate constants, for prep.json."""
    names = ("DNA_LOW_SIGMAS", "DNA_LOW_MIN", "DNA_LOW_MAX", "DNA_HIGH_SIGMAS", "DNA_HIGH_MIN", "DNA_HIGH_MAX",
             "DNA_DIAG_MADS", "DNA_DIAG_MIN", "BEAD_MIN_SIGNAL", "BEAD_THRESHOLD_FLOOR", "DEAD_MIN_GAP", "DEAD_MIN_SIGNAL",
             "MIN_DEAD_FRACTION", "MAX_DEAD_FRACTION", "DEAD_SLOPE_RISE", "DEAD_SLOPE_WINDOW", "WARN_FRACTION")
    return {n.lower(): globals()[n] for n in names}


@dataclass
class GateChannels:
    dna: list[str] = field(default_factory=list)
    beads: list[str] = field(default_factory=list)
    viability: list[str] = field(default_factory=list)

    def any(self) -> bool:
        return bool(self.dna or self.beads or self.viability)

    def as_dict(self) -> dict:
        return {"dna": self.dna, "beads": self.beads, "viability": self.viability}


# EQ / normalisation bead metals: 140Ce, 151Eu, 153Eu, 165Ho, 175Lu. Only labelled bead
# channels are used, plus an unlabelled 140Ce (the others usually carry antibodies).
def find_channels(chans) -> GateChannels:
    """DNA, bead and viability channels among the panel's marker channels (panel.Channel objects)."""
    found = GateChannels()
    for ch in chans:
        if ch.type != "marker":
            continue
        role = C.gating_role(ch.marker, ch.channel)
        unlabelled = not ch.marker or ch.marker == ch.channel
        if role == "dna":
            found.dna.append(ch.channel)
        elif role == "beads" or (unlabelled and C.metal_of(ch.channel) == ("Ce", 140)):
            found.beads.append(ch.channel)
        elif role == "viability":
            found.viability.append(ch.channel)
    # Ir191 before Ir193, at most two
    found.dna = sorted(found.dna, key=lambda c: (C.metal_of(c) or ("", 0))[1])[:2]
    return found


def plan(mode: str, cyt_type: str, gc: GateChannels) -> tuple[bool, str]:
    """(apply?, why) for PREGATE = auto | true | false at run level."""
    if mode == "false":
        return False, "PREGATE=false"
    if not gc.any():
        return False, "no DNA, bead or viability channel in the panel"
    names = ", ".join(f"{k}: {'/'.join(v)}" for k, v in gc.as_dict().items() if v)
    if mode == "true":
        return True, f"PREGATE=true ({names})"
    if cyt_type != "mass":
        return False, f"PREGATE=auto: {cyt_type} cytometer, not mass cytometry; not gated"
    return True, f"PREGATE=auto: mass cytometry with gating channels ({names})"


# --- density helpers -----------------------------------------------------------

def smoothed_density(x: np.ndarray, lo: float, hi: float, sigma: float = 0.06) -> tuple[np.ndarray, np.ndarray]:
    """(bin centres, Gaussian-smoothed counts) of x on [lo, hi]."""
    from scipy.ndimage import gaussian_filter1d

    hi = max(hi, lo + 10 * BIN)
    edges = np.arange(lo, hi + BIN, BIN)
    counts, _ = np.histogram(x, edges)
    return (edges[:-1] + edges[1:]) / 2, gaussian_filter1d(counts.astype(float), sigma / BIN, mode="constant")


def mode_and_halfwidths(x: np.ndarray, lo: float, hi: float) -> tuple[float, float, float]:
    """The highest mode of x within [lo, hi] and its left / right half-widths at half maximum."""
    centres, dens = smoothed_density(x, lo, hi)
    i = int(np.argmax(dens))
    half = dens[i] / 2
    left = i
    while left > 0 and dens[left] > half:
        left -= 1
    right = i
    while right < len(dens) - 1 and dens[right] > half:
        right += 1
    return float(centres[i]), float(centres[i] - centres[left]), float(centres[right] - centres[i])


# --- the gates -----------------------------------------------------------------

def bead_threshold(b: np.ndarray) -> tuple[float | None, str]:
    """Threshold on the bead signal, or None when there is no bead population."""
    high = b[b > BEAD_MIN_SIGNAL]
    if len(high) < max(10, 0.0002 * len(b)):
        return None, f"no bead population ({len(high)} events above {BEAD_MIN_SIGNAL:g})"
    centres, dens = smoothed_density(b, 0.0, float(np.percentile(high, 99.9)) + 0.5)
    above = centres > BEAD_MIN_SIGNAL
    peak = float(centres[above][np.argmax(dens[above])])
    window = (centres >= 1.0) & (centres <= peak - 0.3)
    if not window.any():
        return max(BEAD_THRESHOLD_FLOOR, peak - 1.0), "bead peak close to the cells; peak - 1"
    # the valley between the cells and the bead peak: the middle of the longest
    # (nearly) empty stretch, else the density minimum
    w_c, w_d = centres[window], dens[window]
    empty = w_d <= max(0.5, w_d.min())
    runs, start = [], None
    for k, e in enumerate(np.append(empty, False)):
        if e and start is None:
            start = k
        elif not e and start is not None:
            runs.append((k - start, start, k - 1))
            start = None
    _, a, b = max(runs)
    valley = float((w_c[a] + w_c[b]) / 2)
    thr = float(np.clip(valley, BEAD_THRESHOLD_FLOOR, peak - 0.5))
    return thr, f"valley between cells and the bead peak at {peak:.2f}"


def dna_gate(d: list[np.ndarray]) -> tuple[np.ndarray, dict]:
    """Keep mask and thresholds of the DNA+ gate on one or two (arcsinh) DNA channels."""
    s = d[0] if len(d) == 1 else (d[0] + d[1]) / 2
    top = float(np.percentile(s, 99.9)) + 0.5
    m, hl, hr = mode_and_halfwidths(s, 0.5, top)
    sl, sr = hl * HWHM_TO_SIGMA, hr * HWHM_TO_SIGMA
    low = m - float(np.clip(DNA_LOW_SIGMAS * sl, DNA_LOW_MIN, DNA_LOW_MAX))
    high = m + float(np.clip(DNA_HIGH_SIGMAS * sr, DNA_HIGH_MIN, DNA_HIGH_MAX))
    keep = (s >= low) & (s <= high)
    info = {"dna_mode": round(m, 3), "dna_sigma_left": round(sl, 3), "dna_sigma_right": round(sr, 3),
            "dna_low": round(low, 3), "dna_high": round(high, 3)}
    if len(d) == 2:
        r = d[1] - d[0]
        core = r[keep] if keep.any() else r
        r0 = float(np.median(core))
        mad = float(np.median(np.abs(core - r0))) * 1.4826
        width = max(DNA_DIAG_MADS * mad, DNA_DIAG_MIN)
        keep &= np.abs(r - r0) <= width
        info.update({"dna_diagonal_center": round(r0, 3), "dna_diagonal_halfwidth": round(width, 3)})
    return keep, info


def viability_threshold(v: np.ndarray, min_signal: float = 0.0) -> tuple[float | None, dict]:
    """
    Threshold on the viability signal (dead = high), from the data, with a floor and a cap.

    Dead cells are the brightest population on a viability channel. Mass cytometry
    channels are zero-inflated, so the density is estimated on the events with
    signal (> 0.1), and live cells often show as a zero spike plus a dim mode. Within
    one population the slope of the log density keeps falling towards the bright
    side; a brighter population makes it rise again, either through a valley (two
    modes) or a shoulder. Every rise of at least DEAD_SLOPE_RISE (from the lowest
    slope in the DEAD_SLOPE_WINDOW below a local maximum of the slope) gives a candidate
    threshold: the density minimum of a valley, or for a shoulder the point where
    the slope is half-way back up. Of the candidates that remove between
    MIN_DEAD_FRACTION and MAX_DEAD_FRACTION of the gate's input, the brightest valley
    is used, else the brightest shoulder. Without one
    (a unimodal channel) only outliers go: Tukey's upper fence, Q3 + 1.5 IQR of the
    events with signal. The threshold is never below the median + DEAD_MIN_GAP or
    `min_signal`, and never removes more than MAX_DEAD_FRACTION.

    `min_signal` (DEAD_MIN_SIGNAL for mass cytometry) keeps the threshold out of the
    live cells' background: on CyTOF viability channels live cells are a zero spike
    plus a dim mode (arcsinh ~2-3.5, i.e. 20-80 counts, e.g. 115In in Kowli et al.),
    and the valley between the two is not a live/dead boundary.
    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks

    n = len(v)
    pos = v[v > 0.1]
    if len(pos) < 50:
        return None, {"dead_method": "too few events with signal"}
    sorted_v = np.sort(v)

    def removed(t: float) -> float:
        return (n - np.searchsorted(sorted_v, t, side="right")) / n

    top = float(np.percentile(pos, 99.9)) + 0.5
    centres, dens = smoothed_density(pos, 0.1, top, sigma=0.1)
    slope = gaussian_filter1d(np.gradient(np.log(dens + 1.0), BIN), 0.1 / BIN)
    end = int(np.flatnonzero(dens >= 0.002 * dens.max()).max())
    candidates: list[tuple[float, str]] = []
    window = int(round(DEAD_SLOPE_WINDOW / BIN))
    for p in find_peaks(slope[:end + 1])[0]:
        start = max(0, p - window)
        lo = start + int(np.argmin(slope[start:p + 1]))
        if slope[p] - slope[lo] < DEAD_SLOPE_RISE:
            continue
        if slope[p] > 0:  # a valley: the density minimum between the two modes
            k = lo + int(np.argmin(dens[lo:p + 1]))
            kind = "valley"
        else:             # a shoulder: half-way up the slope's rise
            mid = (slope[lo] + slope[p]) / 2
            k = lo + int(np.argmax(slope[lo:p + 1] >= mid))
            kind = "shoulder"
        t = float(centres[k])
        if t >= min_signal and MIN_DEAD_FRACTION <= removed(t) <= MAX_DEAD_FRACTION:
            candidates.append((t, kind))
    if candidates:
        # a valley (two modes) is stronger evidence than a shoulder; the brightest of the strongest kind
        valleys = [c for c in candidates if c[1] == "valley"]
        thr, method = max(valleys or candidates)
    else:
        q1, q3 = np.percentile(pos, [25, 75])
        thr, method = float(q3 + 1.5 * (q3 - q1)), "upper fence"
    floor = max(float(np.median(v)) + DEAD_MIN_GAP, min_signal)
    cap = float(np.quantile(v, 1 - MAX_DEAD_FRACTION))
    final = max(thr, floor, cap)
    if final > thr:
        method += " (floor)" if floor >= cap else " (cap)"
    return float(final), {"dead_method": method, "dead_candidates": [round(t, 3) for t, _ in candidates]}


@dataclass
class GateResult:
    keep: np.ndarray
    events_in: int
    beads_removed: int = 0
    beads_free: int = 0
    bead_cell_doublets: int = 0
    dna_removed: int = 0
    dead_removed: int = 0
    thresholds: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    # arcsinh values for the plot
    plot: dict = field(default_factory=dict)

    @property
    def events_out(self) -> int:
        return int(self.keep.sum())

    def summary(self) -> dict:
        return {"events_in": self.events_in, "beads_removed": self.beads_removed,
                "beads_free": self.beads_free, "bead_cell_doublets": self.bead_cell_doublets,
                "dna_removed": self.dna_removed, "dead_removed": self.dead_removed,
                "events_out": self.events_out, "thresholds": self.thresholds, "notes": self.notes,
                "warnings": self.warnings}


def _check(res: GateResult, name: str, removed: int, entering: int) -> None:
    if entering and removed / entering > WARN_FRACTION:
        res.warnings.append(f"{name} gate removed {removed} of {entering} events ({removed / entering:.0%}), "
                            f"more than {WARN_FRACTION:.0%}; applied anyway, check the gate plot")


def gate(columns: dict[str, np.ndarray], gc: GateChannels, cofactors: dict[str, float],
         mass: bool = True) -> GateResult:
    """
    Apply the gates to one sample. `columns` maps channel -> raw values (all events),
    `cofactors` channel -> arcsinh cofactor; `mass` for mass cytometry (DEAD_MIN_SIGNAL
    applies). Returns the keep mask and the counts.
    """
    t = {ch: np.arcsinh(columns[ch] / cofactors.get(ch, 5.0)) for ch in gc.dna + gc.beads + gc.viability}
    n = len(next(iter(columns.values())))
    res = GateResult(keep=np.ones(n, dtype=bool), events_in=n)
    dna = [t[c] for c in gc.dna]
    s = None if not dna else (dna[0] if len(dna) == 1 else (dna[0] + dna[1]) / 2)

    # a) beads
    bead_out = np.zeros(n, dtype=bool)
    if gc.beads:
        b = np.min(np.column_stack([t[c] for c in gc.beads]), axis=1)
        thr, why = bead_threshold(b)
        res.thresholds["bead"] = None if thr is None else round(thr, 3)
        res.notes.append(f"beads: {why}")
        if thr is not None:
            bead_out = b >= thr
        res.beads_removed = int(bead_out.sum())
        _check(res, "bead", res.beads_removed, n)
        res.keep &= ~bead_out
        res.plot["bead"] = b

    # b) DNA+ (thresholds from the bead-free events)
    if dna:
        entering = int(res.keep.sum())
        sub_keep, info = dna_gate([x[res.keep] for x in dna])
        res.thresholds.update(info)
        keep_dna = np.zeros(n, dtype=bool)
        keep_dna[np.flatnonzero(res.keep)[sub_keep]] = True
        res.dna_removed = entering - int(sub_keep.sum())
        _check(res, "DNA", res.dna_removed, entering)
        res.keep = keep_dna
        res.beads_free = int((bead_out & (s < info["dna_low"])).sum())
        res.bead_cell_doublets = res.beads_removed - res.beads_free
        res.plot["dna"] = dna
        res.plot["dna_in"] = ~bead_out
    else:
        res.beads_free = res.beads_removed

    # c) live
    if gc.viability:
        v = t[gc.viability[0]]
        entering = int(res.keep.sum())
        res.plot["viability_in"] = res.keep.copy()
        thr, info = viability_threshold(v[res.keep], DEAD_MIN_SIGNAL if mass else 0.0)
        res.thresholds.update(info)
        res.thresholds["dead"] = None if thr is None else round(thr, 3)
        if thr is not None:
            dead = res.keep & (v > thr)
            res.dead_removed = int(dead.sum())
            _check(res, "viability", res.dead_removed, entering)
            res.keep &= ~dead
        res.plot["viability"] = v
        if len(gc.viability) > 1:
            res.notes.append(f"viability gated on {gc.viability[0]} only (also found: {', '.join(gc.viability[1:])})")
    return res


# --- plot ----------------------------------------------------------------------

def plot(res: GateResult, gc: GateChannels, labels: dict[str, str], title: str, path: Path, seed: int = 42,
         max_points: int = 100_000) -> None:
    """A small PNG: bead vs DNA, DNA1 vs DNA2 with the gate, viability histogram with its threshold."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    n = res.events_in
    idx = np.arange(n)
    if n > max_points:
        idx = np.sort(np.random.default_rng(seed).choice(n, max_points, replace=False))
    panels = [k for k in ("bead", "dna", "viability") if k in res.plot]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 3.0), dpi=80, squeeze=False)
    axes = axes[0]
    th = res.thresholds
    dna = res.plot.get("dna")
    s = None if dna is None else (dna[0] if len(dna) == 1 else (dna[0] + dna[1]) / 2)
    cmap = "viridis"
    for ax, key in zip(axes, panels):
        ax.tick_params(labelsize=7)
        if key == "bead":
            b = res.plot["bead"]
            if s is not None:
                ax.hist2d(b[idx], s[idx], bins=120, range=[[0, max(8, b.max())], [0, max(8, s.max())]], cmap=cmap,
                          norm=LogNorm(), rasterized=True)
                ax.set_ylabel("DNA (mean Ir)", fontsize=8)
            else:
                ax.hist(b[idx], bins=150, log=True, color="#3b6ea5")
            if th.get("bead") is not None:
                ax.axvline(th["bead"], color="red", lw=1)
            ax.set_xlabel(" / ".join(labels[c] for c in gc.beads), fontsize=8)
            ax.set_title(f"beads out: -{res.beads_removed}", fontsize=9)
        elif key == "dna":
            k = idx[res.plot["dna_in"][idx]]
            if len(dna) == 2:
                ax.hist2d(dna[0][k], dna[1][k], bins=120, range=[[0, 8], [0, 8]], cmap=cmap, norm=LogNorm(),
                          rasterized=True)
                lo, hi = th["dna_low"], th["dna_high"]
                r0, w = th["dna_diagonal_center"], th["dna_diagonal_halfwidth"]
                # the gate in (DNA1, DNA2): s in [lo, hi], r in [r0 - w, r0 + w]; d1 = s - r/2, d2 = s + r/2
                corners = [(lo, r0 - w), (hi, r0 - w), (hi, r0 + w), (lo, r0 + w), (lo, r0 - w)]
                ax.plot([c[0] - c[1] / 2 for c in corners], [c[0] + c[1] / 2 for c in corners], color="red", lw=1)
                ax.set_xlabel(labels[gc.dna[0]], fontsize=8)
                ax.set_ylabel(labels[gc.dna[1]], fontsize=8)
            else:
                ax.hist(dna[0][k], bins=150, log=True, color="#3b6ea5")
                for x in (th["dna_low"], th["dna_high"]):
                    ax.axvline(x, color="red", lw=1)
                ax.set_xlabel(labels[gc.dna[0]], fontsize=8)
            ax.set_title(f"DNA+ cells: -{res.dna_removed}", fontsize=9)
        else:
            v = res.plot["viability"]
            vv = v[idx][res.plot["viability_in"][idx]]  # the gate's input: after beads and DNA
            ax.hist(vv, bins=150, range=(0, max(8, float(vv.max()) if len(vv) else 8)), log=True, color="#3b6ea5")
            if th.get("dead") is not None:
                ax.axvline(th["dead"], color="red", lw=1)
            ax.set_xlabel(labels[gc.viability[0]], fontsize=8)
            ax.set_title(f"live: -{res.dead_removed}", fontsize=9)
    fig.suptitle(f"{title}: {res.events_in} -> {res.events_out} events", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=80)
    plt.close(fig)
