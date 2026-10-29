# processor-fcs-prep

The first step of the two immune cytometry workflows (Immune Profiling & Differential
Abundance; Multi-batch Harmonization & QC). It takes the FCS files of a study as they come off
the cytometer (or out of FlowJo after pre-gating), plus an optional sample sheet, and turns
them into one consistent set every later step can use without glue:

- **one panel**: only the channels that every sample measured, with the same marker, are kept,
  so a FlowSOM model, a CytoNorm reference or a diffcyt test never compares CD163 in one batch
  with CD28 in another;
- **compensated when needed**: conventional flow files are compensated with their own
  `$SPILLOVER`; spectral (already unmixed) and mass cytometry files are not;
- **cleaned up when raw**: raw mass cytometry (CyTOF) files are pre-gated to DNA⁺ live cells,
  beads out (see *Pre-gating*); pre-gated files are left alone;
- **on one scale**: marker channels are transformed with arcsinh(x / cofactor), the variance-
  stabilising transform used by FlowSOM, CytoNorm, cyCombine and diffcyt; scatter, Time and
  autofluorescence stay as acquired;
- **described**: `samples.csv` (who is who: batch, condition, donor…), `panel.csv` (which
  markers define cell types and which describe their state) and `prep.json` (what was done to
  each file and why).

FCS reading and compensation use FlowKit / FlowIO (White et al. 2021); the writer is FlowIO's
FCS 3.1 writer. Pre-gating is this app's own code (NumPy/SciPy), following the published
CyTOF clean-up sequence (Bagwell et al. 2020); the gate plots use Matplotlib.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `fcs` | input | root | `*.fcs` (FCS 2.0/3.0/3.1) anywhere under `INPUT_DIR` (folders and symlinks followed), and the optional sample sheet `SAMPLE_SHEET` (default `samples.csv`). |
| `fcs` | output | `OUTPUT_DIR/fcs/` | `<sample_id>.fcs` per sample (FCS 3.1, float32, panel channels only, markers transformed) and `provenance.json`. |
| `metadata` | output | `OUTPUT_DIR/metadata/` | `samples.csv`, `panel.csv`, `prep.json`, `provenance.json`; when pre-gated, `gating.csv` and `gates/<sample_id>.png`. |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | Name of the sample sheet, looked for anywhere in the input (the shallowest one wins). A CSV (`,`, `;` or tab) with a `file` column holding FCS file names; `sample_id`, `batch` and `condition` are used when present, and every other column (e.g. `donor`, `timepoint`, `control`) is carried into `samples.csv` unchanged. |
| `COMPENSATE` | `auto` | `auto`: apply the file's `$SPILLOVER`/`$SPILL` when it is not the identity and the cytometer is neither spectral nor mass. `true`: apply any non-identity matrix (also on spectral data). `false`: never. |
| `PREGATE` | `auto` | Pre-gating of raw mass cytometry: beads out, DNA⁺ cells, live cells (see *Pre-gating*). `auto`: only for mass cytometry whose panel has a DNA (Ir191/Ir193), bead (140Ce…) or viability (115In_Dead, cisplatin, Rh103…) channel; spectral and conventional flow data, and panels without such channels, are not touched. `true`: gate whenever any of those channels is present (also on flow data, with a warning: the gates are tuned for CyTOF). `false`: never (for files already gated, e.g. in Cytobank or FlowJo). |
| `COFACTOR` | `auto` | Cofactor *c* in arcsinh(x / *c*) for every marker channel. `auto`: 5 mass cytometry, 6000 spectral flow, 150 conventional flow (see *Transformation*). A number sets it for all markers. |
| `COFACTORS` | `{}` | JSON object of per-marker overrides, keyed by marker or channel, e.g. `{"CD3": 3000, "BV421-A": 2000}`. |
| `PANEL_OVERRIDES` | `{}` | JSON object changing `panel.csv` flags, e.g. `{"CD25": {"clustering": true}, "HLADR": {"clustering": true, "state": false}}`. Keys match marker or channel names, ignoring case and punctuation (`HLA-DR` = `HLADR`). Unknown keys are logged. |
| `MAX_EVENTS_PER_SAMPLE` | `0` | Keep at most this many events per sample: a seeded random subset, kept in acquisition order (Time stays monotonic for PeacoQC). `0` keeps every event. For quick test runs. |
| `SEED` | `42` | Seed of that subsample (combined with the sample_id, so each sample gets its own reproducible subset). |

## What it does, step by step

1. **Finds** every `*.fcs` and the sample sheet under `INPUT_DIR`, reads only the headers.
   A file whose `$ENDDATA` points one byte past its DATA segment (written so by the CyTOF
   software and the Nolan lab bead normaliser, among others) fails FlowIO's check; that one
   error is retried with FlowIO's `ignore_offset_error=True`, logged as a WARNING per file and
   recorded as `prep.json → samples.<id>.ignore_offset_error`. Any other offset problem still
   fails the file.
2. **Names the samples.** With a sheet, a file's row (matched by file name, then case-
   insensitively) gives `sample_id`, `batch`, `condition` and the extra columns. A file that is
   not in the sheet is **processed with defaults and a warning** (so one forgotten row does not
   lose a sample; check the log and `prep.json → notes`); a sheet row without a file is logged.
   Defaults: `sample_id` = the file stem with anything outside `[A-Za-z0-9_-]` replaced by `_`,
   `batch` = the file's `$DATE` as an ISO date, `condition` = `NA`.
3. **Validates.** Unreadable files and files with no events are skipped and listed in
   `prep.json → skipped`. Two files that end up with the same `sample_id` are both skipped
   (the step cannot tell which one is meant); give them distinct ids in the sheet.
4. **Builds the panel.** A channel (`$PnN`) is kept when every sample has it **and** every
   sample says the same marker on it (`$PnS`, compared ignoring case and punctuation). Channels
   missing from some samples, or carrying different markers in different samples (a panel
   change between batches), are dropped from every output file and listed per sample in
   `prep.json → channels_dropped` with the reason. The marker name written to `panel.csv` and
   to every output `$PnS` is the most common spelling; renames are listed in
   `prep.json → markers_renamed`.
5. **Classifies channels**: `scatter` (FSC\*/SSC\*), `time` (Time), `other` (autofluorescence
   `AF-A` on the Aurora, CyTOF technical channels such as `Event_length`, and anything without a
   marker that is not a fluorochrome or metal), `marker` (has a `$PnS`, or is named after a
   fluorochrome or a metal isotope, e.g. `FITC-A`, `Yb171Di`).
6. **Compensates** (see `COMPENSATE`) with FlowKit (events × inverse(spillover)), per file.
7. **Pre-gates** raw mass cytometry (see `PREGATE` and *Pre-gating*): beads, then DNA⁻ debris and
   aggregates, then dead cells are removed, on all events of the file, before
   `MAX_EVENTS_PER_SAMPLE` subsamples the cells that are left.
8. **Transforms** marker channels with arcsinh(x / cofactor); the rest are written unchanged.
9. **Writes** `<sample_id>.fcs`: FCS 3.1, float32, the panel's channels in the first file's
   order, `$PnN`/`$PnS` kept (harmonised marker names), `$DATE`, `$BTIM`, `$ETIM`, `$CYT`,
   `$CYTSN`, `$TIMESTEP` and other descriptive keywords copied, `$FIL` set to the new name, plus
   `PENNSIEVE_PREP` (e.g. `transform=arcsinh;cofactor=6000;compensation=none;events=559685;
   app=processor-fcs-prep@0.1.0`; pre-gated files add `pregate=beads+dna+live;pregated_from=N`),
   `PENNSIEVE_SAMPLE_ID`, `PENNSIEVE_SOURCE_FILE` and
   `PENNSIEVE_COFACTORS` (channel=cofactor list). `$SPILLOVER` is not copied: the output is
   either compensated already or was never meant to be. `$PnR` is the source's range, raised
   when the data exceed it (SpectroFlo writes `$PnR=129` for a Time channel that runs to 10⁶,
   and flowCore's `read.FCS` truncates at `$PnR` by default).

Files are handled one at a time; only the headers of all files are in memory together.

## Outputs

**`metadata/samples.csv`**: one row per prepared sample.

| Column | |
|---|---|
| `sample_id` | `[A-Za-z0-9_-]+`, the output file is `fcs/<sample_id>.fcs` |
| `file` | original FCS file name |
| `batch` | from the sheet, else the file's `$DATE` (ISO) |
| `condition` | from the sheet, else `NA` |
| `n_events` | events written (after `MAX_EVENTS_PER_SAMPLE`) |
| … | every other sample sheet column, unchanged (empty for files not in the sheet) |

**`metadata/panel.csv`**: one row per channel kept: `channel` (`$PnN`), `marker` (`$PnS` without its
mass tag, e.g. `154Sm_CD3` → `CD3`, `CD154__CD40L_` → `CD154 (CD40L)`, unless two channels
would then share a name; or the channel name), `type` (`marker` / `scatter` / `time` / `other`), `cofactor` (empty unless
transformed), `clustering` and `state` (`true`/`false`, markers only; see *Panel flags*).

**`metadata/prep.json`**: run-level (`cytometerType`, `defaultCofactor`, the panel with the
spellings seen, `channelsDropped`, `pregate` (`mode`, `applied`, `reason`, the DNA / bead /
viability `channels`, and when applied the gating `cofactors`, `method` and `settings`),
`skipped`) and per sample: input path, FCS version, cytometer, serial, cytometer type,
`$CREATOR`, `$DATE`, batch, condition, `channels_dropped` (channel, marker, reason),
`markers_renamed`, `cofactors`, `compensation` (`applied`, `reason`, detectors),
`events_in` (events in the file), `events_out` (written), `subsampled`,
`ignore_offset_error`, `nonfinite_values`, seconds, notes, and when pre-gated `pregate`
(`events_in`, `beads_removed` = `beads_free` + `bead_cell_doublets`, `dna_removed`,
`dead_removed`, `events_out`, `thresholds`, `notes`, `warnings`, `plot`).

**`metadata/gating.csv`** (only when pre-gating applied): one row per sample, `sample_id`,
`events_in`, `beads_removed`, `dna_removed`, `dead_removed`, `events_out` (before
`MAX_EVENTS_PER_SAMPLE`) and `thresholds` (`bead=…;dna_low=…;dna_high=…;
dna_diagonal_center=…;dna_diagonal_halfwidth=…;dead=…`, arcsinh units).

**`metadata/gates/<sample_id>.png`** (only when pre-gating applied): a small plot per sample
(about 770 × 240 px): DNA vs bead channel with the bead threshold, Ir191 vs Ir193 with the DNA
gate, and the viability histogram (log counts) with the dead threshold.

**`provenance.json`** (in both folders): app and version, tool versions (FlowKit, FlowIO,
FlowUtils, NumPy, pandas, Python), parameters, inputs found, timings and per-sample status.

## Transformation

arcsinh(x / *c*) is linear around zero (so the negative spread of unmixed or compensated
data stays symmetric) and logarithmic for bright events; the cofactor *c* sets where one turns
into the other and has to match the instrument's scale:

| Cytometer (detected from `$CYT`, `$CREATOR` and channel names) | `COFACTOR=auto` |
|---|---|
| mass (CyTOF, Helios, XT; or mostly metal channels such as `Yb171Di`) | 5 |
| spectral (Cytek Aurora / Northern Lights, Sony ID7000 / SP6800, BD FACSDiscover / A5 SE, …; or `$CREATOR` SpectroFlo) | 6000 |
| conventional flow (anything else) | 150 |

These are the values cyCombine's authors recommend and use ("recommended cofactors are 5 for
CyTOF, 150 for flow cytometry, and 6000 for spectral flow cytometry"; Pedersen et al. 2022,
Methods), and 5 is the CyTOF convention since Bendall et al. 2011. On the demo Aurora data,
after arcsinh(x/6000), negative populations are centred near 0 and the 99th percentile of the
markers lies between 1.5 and 5.3. A cofactor that is right for one fluorochrome can be too small or too large
for a dim or very bright one; check the biaxial plots and adjust with `COFACTORS`.

## Panel flags

`clustering` markers define cell types (FlowSOM, UMAP); `state` markers describe what those
cells are doing and are tested *within* clusters for differential state (diffcyt), the
"type" vs "state" split of Nowicka et al. 2017 and Weber et al. 2019. A marker that changes
with the experimental condition should not drive clustering: in a stimulated-vs-unstimulated
design it would split one cell type into a "stimulated" and an "unstimulated" cluster and turn
a change of state into a spurious change of abundance.

Defaults: every marker `clustering=true, state=false`, except

- **state markers** (`clustering=false, state=true`): activation — CD25, CD69, HLA-DR, CD38,
  CD71, CD137 (4-1BB), CD134 (OX40), CD154 (CD40L), ICOS, GITR; inhibitory receptors — PD-1,
  TIGIT, CTLA-4, LAG-3, TIM-3, CD39, 2B4; cytotoxic granules — granzymes A/B/K/H, perforin,
  granulysin, CD107a; CD95 (Fas); Ki-67; intracellular cytokines (IFN-γ, TNF, IL-2, IL-4, IL-17…);
  phospho-epitopes (pSTAT…, pERK…, pS6…); DNA damage (γH2AX: `pH2A.X`, `pH2AX`, `gH2AX`);
  PD-L1 (CD274) and PD-L2 (CD273), induced by IFN-γ and activation.
- **neither** (`clustering=false, state=false`): CD45 (a gating marker, constant after the
  CD45⁺ gate), viability stains (`Dead`, `LD`, `L/D`, `Live`, `LIVE/DEAD …`, `Viab…`,
  `Cisplatin`, unlabelled 195Pt/198Pt, the Rh103 intercalator, Zombie, 7-AAD, DAPI, FVS…;
  matched as whole words, so `LDLR` or `LDHA` are markers), DNA intercalators (Ir191/Ir193,
  `DNA1`/`DNA2`), beads and barcodes.

Names are matched after the mass tag is stripped (`171Yb_CD3`, `CD3_Yb171Di` → `CD3`) and
ignoring case and punctuation. A marker written with an alias — `CD154 (CD40L)`,
`CD154/CD40L`, or Fluidigm's `168Er_CD154__CD40L_` (spaces and brackets become `_`) — is a
state marker when any of its names is one.

Choices an immunologist may want to revisit, and why the defaults are what they are:

- **FOXP3 clusters, CD25 does not.** FOXP3 (with HELIOS) is the lineage transcription factor
  of regulatory T cells, so it defines the Treg cluster. CD25 also marks Tregs, but every
  conventional T cell up-regulates it on activation (the demo's anti-CD3/CD28 + IL-2
  stimulation does exactly that), so using it for clustering would merge activated
  conventional T cells with Tregs. Tregs are still found through FOXP3, CD45RA, CD4 and
  HELIOS, and CD25 is tested as a state marker. In a resting-only study,
  `PANEL_OVERRIDES={"CD25": {"clustering": true}}` restores the classic CD25⁺CD127^lo gate.
- **HLA-DR is a state marker.** On T and NK cells it is an activation marker (with CD38, the
  "activated" T cell phenotype of Mathew et al. 2020). On B cells, monocytes and dendritic
  cells it is lineage; for myeloid or DC-focused panels set it to `clustering: true`.
- **Granzyme B, perforin, granulysin are state markers.** They change with activation and
  differentiation, but in NK- or CTL-focused studies they also define subsets; override if so.
- **CD95** is listed as state (it rises on activation); it is also how stem-cell memory T cells
  (CD45RA⁺CCR7⁺CD95⁺) are told from naïve ones, so turn on clustering when that matters.
- Differentiation markers (CD45RA, CCR7, CD27, CD28, CD57, CD127, chemokine receptors) and
  lineage markers cluster by default.

## Pre-gating (mass cytometry)

Raw CyTOF files hold every event the instrument pushed: normalisation beads, debris, cell
aggregates and dead cells next to the cells. The standard clean-up before any analysis is
beads out → DNA⁺ (intercalator) cells → live cells, then singlets (Bagwell et al. 2020 describe
and automate this sequence; Kowli et al. 2025 gated the test data below by hand in Cytobank on
Ir DNA, 140Ce beads, event length and 115In live/dead). With `PREGATE=auto` this step does the
first three automatically for mass cytometry; doublets are left to cleanet, which handles them
properly, and acquisition anomalies to peacoqc.

The gates work on arcsinh(x / 5) of all events of a file, in this order, each on what the
previous one kept. Every threshold comes from the file's own data, by deterministic rules (no
random step; the plots draw a seeded subsample):

1. **Beads out.** Bead channels: those labelled as beads (`140Ce_Beads`, or 151Eu, 153Eu, 165Ho,
   175Lu when they say "bead"), and an unlabelled 140Ce. On the lowest of them (beads are high
   in all), the threshold is the middle of the empty stretch between the cells (near 0) and the
   bead peak (searched above 3). Events above it are removed: free beads (low DNA) and
   bead–cell doublets (DNA⁺), the latter because EQ bead metals would show as false signal on
   the 151Eu/153Eu/165Ho/175Lu marker channels. `prep.json` counts both (`beads_free`,
   `bead_cell_doublets`).
2. **DNA⁺ cells.** On Ir191 and Ir193, rotated to s = (DNA1 + DNA2) / 2 along the diagonal and
   r = DNA2 − DNA1 across it: the main mode of s on a smoothed histogram (above 0.5), its
   left and right half-widths at half maximum turned into σ (HWHM / 1.177), then
   low = mode − clip(4 σ_left, 0.6, 1.5) (debris, DNA⁻ events) and
   high = mode + clip(5 σ_right, 1.0, 2.0) (aggregates; at least a factor e ≈ 2.7 in DNA, so
   doublets at ×2 stay for cleanet). Events off the diagonal (the two isotopes disagree) by more
   than max(6 MAD, 0.5) of r are removed too. One Ir channel: the same on it alone.
3. **Live cells.** On the viability channel (dead = bright: cisplatin, 115In, 103Rh
   intercalator, live/dead dyes), the density of the events with signal (> 0.1; CyTOF channels
   are zero-inflated) is smoothed and the slope of its logarithm followed: within one
   population it keeps falling towards the bright side, and a brighter population makes it
   rise again, through a valley (two modes) or a shoulder. Each rise of at least 0.5 (within
   1.5 arcsinh units) is a candidate: the density minimum of a valley, or the point half-way up
   the slope for a shoulder. Candidates must remove between 0.5% and 30% of the gate's input
   and, for mass cytometry, lie above arcsinh 3 (~50 counts): live cells show as a zero spike
   plus a dim mode around 2–3.5, and the valley between those two is not a live/dead boundary.
   The brightest valley is used, else the brightest shoulder; with none (a unimodal channel)
   only outliers go, above Tukey's upper fence (Q3 + 1.5 IQR). The threshold is never below
   the median + 0.75 and never removes more than 30%.

A gate that removes more than 60% of its input is applied anyway, but logged as a WARNING and
recorded in `prep.json → samples.<id>.pregate.warnings`: check that sample's plot. Each sample
gets one log line with the counts per gate and the thresholds, a row in `metadata/gating.csv`
and a plot in `metadata/gates/`. The constants are listed in `prep.json → pregate.settings`.

`PREGATE=auto` gates only mass cytometry (by `$CYT`, or mostly metal channel names) with at
least one DNA, bead or viability channel; a gate whose channels are missing is skipped (e.g.
no viability stain: beads and DNA only). Spectral and conventional flow files are unchanged
(on the demo data the output is identical to the previous release apart from the version in
`PENNSIEVE_PREP` and the `pregate` entry in `prep.json`).

### Mass cytometry test data

Kowli et al. 2025 (renal transplant recipients and healthy donors, pre-transplant and 6
months; Zenodo doi:10.5281/zenodo.17238592, CC BY 4.0): 48 CyTOF (Helios, DVS CyTOF 6.7) files,
bead-normalised with the Nolan lab normaliser, not gated, 64k–355k events. Every file has
`$ENDDATA` off by one (read with `ignore_offset_error`). Run at `MAX_EVENTS_PER_SAMPLE=20000`
(48 files in 24 s in Docker):

| Gate | Removed (median, range over 48 files) |
|---|---|
| beads (of all events) | 2.1% (0.8–4.1%); 89% free beads, 11% bead–cell doublets |
| DNA⁻ / aggregates (of what the bead gate kept) | 5.2% (2.1–17.1%) |
| dead (of what the DNA gate kept) | 2.0% (0.6–4.7%); 38 files by a valley, 10 by a shoulder |
| all three | 9.2% (4.2–21.6%) |

The files were acquired clean (debris is low), and on 115In most live cells sit in the zero
spike and a dim mode (arcsinh ~3) made of the same cell types as the zero spike (checked on
three files); dead cells form a separate bright mode or shoulder at arcsinh 4.5–5.5. Event length was
not gated (the Gaussian parameters `Center`, `Offset`, `Width`, `Residual` and `Event_length`
are kept as `other` channels for cleanet). Panel (42 markers): clustering (18) CD57, CD19,
CD4, CD8, CD20, CD40, CD45RA, CD3, CD28, CD33, CD14, TCRgd, CD127, CD27, CCR7, CD21, CD56,
CD16; state (20) HLA-DR, IL-10, IL-4, CTLA-4, MIP-1b, CD107a, TNFa, CD38, CD25, GM-CSF, IFNg,
CD69, IL-17, IL-2, CD154 (CD40L), PD-1, Granzyme B, pH2A.X, perforin, PD-L1; neither: Dead,
Beads, DNA1, DNA2.

## Runtime

`computeTypes: [standard]`, 2 vCPU / 4 GiB, 1 h timeout. The work is single-threaded and
streams one file at a time; memory scales with the largest file (its events, as float64,
are held a few times over: raw, selected, float32 output).

Measured on the demo data (24 files, 10,212,112 events, 1.46 GB), full event counts, in Docker
(linux/amd64 emulated on an M-series Mac, `--cpus 2 --memory 2g`): **26 s** for the whole run
(1.4 s for the largest file, 602,586 events × 35 channels), cgroup `memory.peak` **1.77 GB**
(this includes the page cache of the 1.3 GB written; the process's own peak RSS is ~1.25 GB
measured natively). Pre-gating holds a few extra float64 columns per file and draws one
small plot (~0.2 s per file); the 48 CyTOF files above take 24 s. 4 GiB leaves room for files of ~1.5 M events; for larger files (e.g. a
5 M-event CyTOF file) give it 8–16 GiB. The image is 190 MB compressed (~590 MB unpacked).

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python        # unit tests on a synthetic fixture
make run PYTHON=.venv/bin/python         # Docker over the fixture -> data/output/
make run-demo DEMO=../workflow-immune-cytometry/data/zenodo-15388817   # demo, 20k events/file
```

Or by hand:

```sh
docker build --platform linux/amd64 -t processor-fcs-prep:dev .
docker run --rm --platform linux/amd64 --memory 4g --env-file dev.env \
  -v /path/to/fcs:/data/input:ro -v "$PWD/data/output":/data/output processor-fcs-prep:dev
```

## Demo data result

The CytoBatchFlagR dataset (Koenig et al. 2025; Cytek Aurora, CD45⁺CD14⁻ PBMC, 6 batches × 4
samples) gives 24 samples, 33 channels kept (6 scatter, AF-A, 25 markers, Time).

- `$SPILLOVER`: batches 01–04 carry SpectroFlo's identity matrix (one 1e-6 entry), batches
  05–06 none; no compensation is applied.
- Dropped: `BV785-A` (CD28, absent in batches 05–06), `PerCP-A` (CD163 in 01–04, CD28 in
  05–06) and `APC-A` (GITR in 01–04, CCR4 in 05–06): the panel changed between batches 04 and
  05. `BV750-A` is `HLADR` in 01–04 and `HLA-DR` in 05–06 and is harmonised to `HLADR`.
- Clustering (14): HELIOS, CCR7, CD4, CD56, CD3, CD19, CD8, CD103, CD16, CD14, CD27, CCR6,
  CD45RA, FOXP3; state (10): GZMB, GNLY, CD95, CD69, PFN, HLADR, CD25, PD1, CD39, TIGIT;
  neither: CD45.

## Limitations

- A marker measured on **different channels** in different batches (CD28 on BV785 then on
  PerCP in the demo) is dropped, not merged: merging would put two fluorochromes with
  different brightness and spillover on one axis. Re-acquire, or prepare the batches
  separately.
- Cytometer type is detected from keywords; files with an unusual `$CYT` fall back to
  "conventional flow" (cofactor 150). Set `COFACTOR` when that is wrong. Mixed instruments in
  one run use the most common type's default (logged).
- `COMPENSATE=auto` trusts the file's `$SPILLOVER`. FlowJo-exported files whose matrix was
  edited in FlowJo carry it in the workspace, not the FCS file; export compensated data or
  supply files with the right `$SPILLOVER`.
- Multi-dataset FCS files (several data sets in one file) are read as their first data set.
- Pre-gating is for raw mass cytometry only (`PREGATE=auto`). Conventional and spectral flow
  data must arrive gated to live single cells (as the demo does); doublets are left to cleanet,
  acquisition anomalies to peacoqc. The gates are automatic: check `metadata/gates/` (and the
  60% warnings), and gate by hand (then `PREGATE=false`) when a sample looks wrong.

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles FlowKit, FlowIO, FlowUtils
(BSD-3-Clause) and their dependencies; see `THIRD_PARTY.md`.

## Citations

FCS reading, compensation (FlowKit, FlowIO, FlowUtils):

> White S, Quinn J, Enzor J, Staats J, Mosier SM, Almarode J, Denny TN, Weinhold KJ, Ferrari G,
> Chan C. **FlowKit: A Python Toolkit for Integrated Manual and Automated Cytometry Analysis
> Workflows.** Front Immunol 12:768541 (2021). doi:10.3389/fimmu.2021.768541, PMID 34804056.

Pre-gating of mass cytometry (beads, DNA, viability):

> Bagwell CB, Inokuma M, Hunsberger B, Herbert D, Bray C, Hill B, Stelzer G, Li S, Kollipara A,
> Ornatsky O, Baranov V. **Automated Data Cleanup for Mass Cytometry.** Cytometry A
> 97(2):184–198 (2020). doi:10.1002/cyto.a.23926.

> Kowli S, Minocherhomji S, Martinez OM, Busque S, Lebrec H, Maecker HT. **Characterization of
> immune phenotypes in peripheral blood of adult renal transplant recipients using mass
> cytometry (CyTOF).** ImmunoHorizons 9(4):vlae013 (2025). doi:10.1093/immhor/vlae013.

Cofactors 5 / 150 / 6000:

> Pedersen CB, Dam SH, Barnkob MB, Leipold MD, Purroy N, Rassenti LZ, Kipps TJ, Nguyen J,
> Lederer JA, Gohil SH, Wu CJ, Olsen LR. **cyCombine allows for robust integration of
> single-cell cytometry datasets within and across technologies.** Nat Commun 13:1698 (2022).
> doi:10.1038/s41467-022-29383-5, PMID 35361793.

> Bendall SC, Simonds EF, Qiu P, Amir ED, Krutzik PO, Finck R, Bruggner RV, Melamed R, Trejo A,
> Ornatsky OI, Balderas RS, Plevritis SK, Sachs K, Pe'er D, Tanner SD, Nolan GP. **Single-cell
> mass cytometry of differential immune and drug responses across a human hematopoietic
> continuum.** Science 332:687–696 (2011). doi:10.1126/science.1198704, PMID 21551058.

Type vs state markers:

> Nowicka M, Krieg C, Crowell HL, Weber LM, Hartmann FJ, Guglietta S, Becher B, Levesque MP,
> Robinson MD. **CyTOF workflow: differential discovery in high-throughput high-dimensional
> cytometry datasets.** F1000Res 6:748 (2017). doi:10.12688/f1000research.11622.3, PMID 28663787.

> Weber LM, Nowicka M, Soneson C, Robinson MD. **diffcyt: Differential discovery in
> high-dimensional cytometry via high-resolution clustering.** Commun Biol 2:183 (2019).
> doi:10.1038/s42003-019-0415-5, PMID 31098416.

> Mathew D, Giles JR, Baxter AE, Oldridge DA, Greenplate AR, Wu JE, et al. **Deep immune
> profiling of COVID-19 patients reveals distinct immunotypes with therapeutic implications.**
> Science 369:eabc8511 (2020). doi:10.1126/science.abc8511, PMID 32669297.

Demo data:

> Koenig Z, Andorf S, Tilburgs T, Eswar S. **CytoBatchFlagR: Spectral flow cytometry dataset.**
> Zenodo (2025). doi:10.5281/zenodo.15388817 (CC BY 4.0).

> Kowli S. **Characterization of immune phenotypes in peripheral blood of adult renal transplant
> recipients using mass cytometry (CyTOF).** Zenodo (2025). doi:10.5281/zenodo.17238592
> (CC BY 4.0). Mass cytometry test data.

All PMIDs were checked against PubMed E-utilities (title, first author, year), the dataset
DOIs against DataCite, and the Bagwell and Kowli article DOIs against Crossref.
