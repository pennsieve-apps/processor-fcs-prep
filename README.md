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
- **on one scale**: marker channels are transformed with arcsinh(x / cofactor), the variance-
  stabilising transform used by FlowSOM, CytoNorm, cyCombine and diffcyt; scatter, Time and
  autofluorescence stay as acquired;
- **described**: `samples.csv` (who is who: batch, condition, donor…), `panel.csv` (which
  markers define cell types and which describe their state) and `prep.json` (what was done to
  each file and why).

FCS reading and compensation use FlowKit / FlowIO (White et al. 2021); the writer is FlowIO's
FCS 3.1 writer.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `fcs` | input | root | `*.fcs` (FCS 2.0/3.0/3.1) anywhere under `INPUT_DIR` (folders and symlinks followed), and the optional sample sheet `SAMPLE_SHEET` (default `samples.csv`). |
| `fcs` | output | `OUTPUT_DIR/fcs/` | `<sample_id>.fcs` per sample (FCS 3.1, float32, panel channels only, markers transformed) and `provenance.json`. |
| `metadata` | output | `OUTPUT_DIR/metadata/` | `samples.csv`, `panel.csv`, `prep.json`, `provenance.json`. |

## Parameters

| Name | Default | What it does |
|---|---|---|
| `SAMPLE_SHEET` | `samples.csv` | Name of the sample sheet, looked for anywhere in the input (the shallowest one wins). A CSV (`,`, `;` or tab) with a `file` column holding FCS file names; `sample_id`, `batch` and `condition` are used when present, and every other column (e.g. `donor`, `timepoint`, `control`) is carried into `samples.csv` unchanged. |
| `COMPENSATE` | `auto` | `auto`: apply the file's `$SPILLOVER`/`$SPILL` when it is not the identity and the cytometer is neither spectral nor mass. `true`: apply any non-identity matrix (also on spectral data). `false`: never. |
| `COFACTOR` | `auto` | Cofactor *c* in arcsinh(x / *c*) for every marker channel. `auto`: 5 mass cytometry, 6000 spectral flow, 150 conventional flow (see *Transformation*). A number sets it for all markers. |
| `COFACTORS` | `{}` | JSON object of per-marker overrides, keyed by marker or channel, e.g. `{"CD3": 3000, "BV421-A": 2000}`. |
| `PANEL_OVERRIDES` | `{}` | JSON object changing `panel.csv` flags, e.g. `{"CD25": {"clustering": true}, "HLADR": {"clustering": true, "state": false}}`. Keys match marker or channel names, ignoring case and punctuation (`HLA-DR` = `HLADR`). Unknown keys are logged. |
| `MAX_EVENTS_PER_SAMPLE` | `0` | Keep at most this many events per sample: a seeded random subset, kept in acquisition order (Time stays monotonic for PeacoQC). `0` keeps every event. For quick test runs. |
| `SEED` | `42` | Seed of that subsample (combined with the sample_id, so each sample gets its own reproducible subset). |

## What it does, step by step

1. **Finds** every `*.fcs` and the sample sheet under `INPUT_DIR`, reads only the headers.
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
7. **Transforms** marker channels with arcsinh(x / cofactor); the rest are written unchanged.
8. **Writes** `<sample_id>.fcs`: FCS 3.1, float32, the panel's channels in the first file's
   order, `$PnN`/`$PnS` kept (harmonised marker names), `$DATE`, `$BTIM`, `$ETIM`, `$CYT`,
   `$CYTSN`, `$TIMESTEP` and other descriptive keywords copied, `$FIL` set to the new name, plus
   `PENNSIEVE_PREP` (e.g. `transform=arcsinh;cofactor=6000;compensation=none;events=559685;
   app=processor-fcs-prep@0.1.0`), `PENNSIEVE_SAMPLE_ID`, `PENNSIEVE_SOURCE_FILE` and
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

**`metadata/panel.csv`**: one row per channel kept: `channel` (`$PnN`), `marker` (`$PnS`, or
the channel name), `type` (`marker` / `scatter` / `time` / `other`), `cofactor` (empty unless
transformed), `clustering` and `state` (`true`/`false`, markers only; see *Panel flags*).

**`metadata/prep.json`**: run-level (`cytometerType`, `defaultCofactor`, the panel with the
spellings seen, `channelsDropped`, `skipped`) and per sample: input path, FCS version,
cytometer, serial, cytometer type, `$CREATOR`, `$DATE`, batch, condition,
`channels_dropped` (channel, marker, reason), `markers_renamed`, `cofactors`,
`compensation` (`applied`, `reason`, detectors), `events_in`, `events_out`, `subsampled`,
`nonfinite_values`, seconds, notes.

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
  phospho-epitopes (pSTAT…, pERK…, pS6…).
- **neither** (`clustering=false, state=false`): CD45 (a gating marker, constant after the
  CD45⁺ gate), viability dyes, DNA intercalators, beads and barcodes.

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

## Runtime

`computeTypes: [standard]`, 2 vCPU / 4 GiB, 1 h timeout. The work is single-threaded and
streams one file at a time; memory scales with the largest file (its events, as float64,
are held a few times over: raw, selected, float32 output).

Measured on the demo data (24 files, 10,212,112 events, 1.46 GB), full event counts, in Docker
(linux/amd64 emulated on an M-series Mac, `--cpus 2 --memory 2g`): **26 s** for the whole run
(1.4 s for the largest file, 602,586 events × 35 channels), cgroup `memory.peak` **1.77 GB**
(this includes the page cache of the 1.3 GB written; the process's own peak RSS is ~1.25 GB
measured natively). 4 GiB leaves room for files of ~1.5 M events; for larger files (e.g. a
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
- No gating, no QC: debris, doublets and dead cells must be gated out before (as in the demo)
  or removed by the QC steps that follow (peacoqc, cleanet).

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles FlowKit, FlowIO, FlowUtils
(BSD-3-Clause) and their dependencies; see `THIRD_PARTY.md`.

## Citations

FCS reading, compensation (FlowKit, FlowIO, FlowUtils):

> White S, Quinn J, Enzor J, Staats J, Mosier SM, Almarode J, Denny TN, Weinhold KJ, Ferrari G,
> Chan C. **FlowKit: A Python Toolkit for Integrated Manual and Automated Cytometry Analysis
> Workflows.** Front Immunol 12:768541 (2021). doi:10.3389/fimmu.2021.768541, PMID 34804056.

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

All PMIDs were checked against PubMed E-utilities (title, first author, year) and the dataset
DOI against DataCite.
