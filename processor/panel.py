"""
The common panel across samples, its cofactors and flags (panel.csv).

A channel ($PnN) is kept when every sample has it AND every sample says the
same marker on it ($PnS compared case- and punctuation-insensitively, so
'HLADR' and 'HLA-DR' agree). A channel that carries different markers in
different batches (a panel change) is dropped from every file: comparing two
markers as one would be wrong. The marker name written to panel.csv and to
every output file's $PnS is the most common spelling.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field

from processor import channels as C


@dataclass
class Header:
    """What the panel needs from one file's TEXT segment."""
    sample_id: str
    pnn: list[str]
    pns: list[str]
    text: dict

    @property
    def markers(self) -> dict[str, str]:
        return dict(zip(self.pnn, self.pns))


@dataclass
class Channel:
    channel: str
    marker: str
    type: str
    cofactor: float | None = None
    clustering: bool = False
    state: bool = False
    spellings: list[str] = field(default_factory=list)

    def row(self) -> dict[str, str]:
        cof = "" if self.cofactor is None else (f"{self.cofactor:g}")
        return {"channel": self.channel, "marker": self.marker, "type": self.type, "cofactor": cof,
                "clustering": str(self.clustering).lower(), "state": str(self.state).lower()}


@dataclass
class Panel:
    channels: list[Channel]
    dropped: dict[str, list[dict]]  # sample_id -> [{channel, marker, reason}]
    renamed: dict[str, list[dict]]  # sample_id -> [{channel, from, to}]

    def by_channel(self) -> dict[str, Channel]:
        return {c.channel: c for c in self.channels}

    def markers(self) -> list[Channel]:
        return [c for c in self.channels if c.type == "marker"]


def build(headers: list[Header]) -> Panel:
    """Intersect channels (by $PnN, with agreeing $PnS) across the given headers."""
    if not headers:
        return Panel([], {}, {})
    n = len(headers)
    present: Counter[str] = Counter()
    spellings: dict[str, Counter[str]] = {}
    keys: dict[str, set[str]] = {}
    for h in headers:
        for ch, mk in h.markers.items():
            present[ch] += 1
            spellings.setdefault(ch, Counter())[mk] += 1
            keys.setdefault(ch, set()).add(C.norm(mk))
    order = list(headers[0].pnn) + [c for h in headers[1:] for c in h.pnn if c not in headers[0].pnn]
    seen: set[str] = set()
    kept: list[Channel] = []
    reasons: dict[str, str] = {}
    for ch in order:
        if ch in seen:
            continue
        seen.add(ch)
        if present[ch] < n:
            reasons[ch] = f"missing in {n - present[ch]} of {n} samples"
            continue
        if len(keys[ch]) > 1:
            reasons[ch] = "different markers in different samples: " + ", ".join(
                f"{m or '(none)'} ({k})" for m, k in spellings[ch].most_common())
            continue
        marker = spellings[ch].most_common(1)[0][0]
        kept.append(Channel(ch, marker or ch, C.classify(ch, marker), spellings=sorted(spellings[ch])))
    dropped: dict[str, list[dict]] = {}
    renamed: dict[str, list[dict]] = {}
    kept_map = {c.channel: c for c in kept}
    for h in headers:
        dropped[h.sample_id] = [{"channel": ch, "marker": mk, "reason": reasons[ch]}
                                for ch, mk in h.markers.items() if ch in reasons]
        renamed[h.sample_id] = [{"channel": ch, "from": mk, "to": kept_map[ch].marker}
                                for ch, mk in h.markers.items()
                                if ch in kept_map and mk and mk != kept_map[ch].marker]
    return Panel(kept, dropped, renamed)


def parse_json_object(value: str | None, name: str) -> dict:
    if value is None or not value.strip():
        return {}
    try:
        obj = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} is not valid JSON: {exc}") from exc
    if not isinstance(obj, dict):
        raise ValueError(f"{name} must be a JSON object")
    return obj


def _lookup(overrides: dict, ch: Channel) -> tuple[str, object] | None:
    """An override keyed by channel or marker (case- and punctuation-insensitive)."""
    for key, value in overrides.items():
        k = C.norm(str(key))
        if k and k in (C.norm(ch.channel), C.norm(ch.marker), C.marker_key(ch.marker)):
            return key, value
    return None


def configure(panel: Panel, cofactor_default: float, cofactors: dict, overrides: dict) -> list[str]:
    """Set cofactors and flags on the panel's channels. Returns warnings (unused override keys)."""
    used: set[str] = set()
    for ch in panel.channels:
        if ch.type != "marker":
            ch.cofactor, ch.clustering, ch.state = None, False, False
            continue
        ch.cofactor = cofactor_default
        hit = _lookup(cofactors, ch)
        if hit:
            used.add("COFACTORS:" + hit[0])
            try:
                value = float(hit[1])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"COFACTORS[{hit[0]!r}] is not a number: {hit[1]!r}") from exc
            if value <= 0:
                raise ValueError(f"COFACTORS[{hit[0]!r}] must be > 0")
            ch.cofactor = value
        ch.clustering, ch.state = C.default_flags(ch.marker)
        hit = _lookup(overrides, ch)
        if hit:
            used.add("PANEL_OVERRIDES:" + hit[0])
            if not isinstance(hit[1], dict):
                raise ValueError(f"PANEL_OVERRIDES[{hit[0]!r}] must be an object like {{\"clustering\": true}}")
            for flag in ("clustering", "state"):
                if flag in hit[1]:
                    v = hit[1][flag]
                    if isinstance(v, str):
                        v = v.strip().lower() in ("true", "1", "yes")
                    setattr(ch, flag, bool(v))
    warnings = []
    for name, obj in (("COFACTORS", cofactors), ("PANEL_OVERRIDES", overrides)):
        for key in obj:
            if f"{name}:{key}" not in used:
                warnings.append(f"{name}: {key!r} matches no marker channel in the panel")
    return warnings
