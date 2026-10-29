# sample-input

Real input for this app is a folder of FCS files, optionally with a `samples.csv` sample
sheet. The workflows' demo dataset is public:

> Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry dataset.
> Zenodo 2025. doi:10.5281/zenodo.15388817 (CC BY 4.0)

24 Cytek Aurora files (6 batches × 4 samples, CD45⁺CD14⁻ pre-gated PBMC, 25k–600k events
each, 1.46 GB). Download the FCS files from https://zenodo.org/records/15388817 into one
folder and add a sample sheet, e.g.

```csv
file,sample_id,batch,condition,control
B01_tcu1_CD45+CD14-.fcs,batch01_tcu1,batch01,unstimulated,tcu1
B01_tcs3_CD45+CD14-.fcs,batch01_tcs3,batch01,stimulated,tcs3
...
```

(`tcu1`/`tcu2`: unstimulated technical controls; `tcs3`/`tcs4`: anti-CD3/CD28 + IL-2 for 7
days.) Then:

```sh
make run-demo DEMO=/path/to/that/folder        # 20k events per file, into data/output-demo/
```

For tests and `make run`, `tests/make_fixture.py` writes a small synthetic Aurora-like
fixture instead (`make fixture` -> `data/fixture/linked/fcs/`): two panels across batches
(a channel missing, a channel with a different marker, `HLADR` vs `HLA-DR`), an identity
`$SPILLOVER`, a file not in the sheet, an empty file and a broken one.
