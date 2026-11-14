# Truncation test of the DESI evolving dark energy signal

Code and data products for the paper:

> **"A Truncation Test of the DESI Evolving Dark Energy Signal: Three Supernova Calibrations Compared"**
> Tao Zhang (Independent Researcher, Ningbo, China)
> *Monthly Notices of the Royal Astronomical Society*, 552(4), stag1868 (2026)
> DOI: [10.1093/mnras/stag1868](https://doi.org/10.1093/mnras/stag1868)

A redshift truncation scan of the CPL (`w0`–`wa`) dark energy parametrization:
SNe with `z <= z_cut` are progressively excised and CPL cosmology is re-fitted
at each threshold, for three SN datasets with independent absolute
calibrations (Pantheon+, DESY5, DES-Dovekie), combined with DESI DR2 BAO and
the full Planck 2018 TTTEEE-lite likelihood.

## Contents

| Path | Description |
|------|-------------|
| `truncation_sn.py` | Cobaya likelihood module: truncated Pantheon+/DESY5/DES-Dovekie likelihoods (the exact version used to produce the paper's chains) |
| `configs/` | 39 run configurations (`13 z_cut x 3 datasets`), full Planck 2018 |
| `hotstart/` | Proposal covariance matrices used to (re)start the MCMC chains |
| `data/DES-Dovekie/` | DES-Dovekie distance moduli + full inverse covariance matrix (from the Dovekie data release) |
| `scripts/server_deploy.sh`, `server_run*.sh` | Server deployment and chain-running scripts |
| `scripts/download_data.py` | Fetches the public Pantheon+ and DES-SN5YR data products |
| `scripts/run_scan.py` | Runs the truncation scan (Cobaya MCMC) |
| `scripts/analyze_scan.py` | Turns chains into the summary statistics |
| `scripts/profile_scan.py` | Profile-likelihood cross-check of the 2D significances |
| `scripts/plot_final_fig.py` | Reproduces the paper's Figure 1 |
| `scripts/export_summary_csv.py` | Regenerates `results/cpl_scan_summary.csv` |
| `results/cpl_scan_summary.csv` | The 13 x 3 grid of `w0`, `wa`, and 2D significance used in the paper (data behind Figure 1 and Table 2) |

## Reproducing the analysis

1. Install the sampling environment: `pip install cobaya camb numpy scipy pyyaml matplotlib`
   (tested with Cobaya 3.6.x and CAMB; the MCMC chains were run on Linux servers
   with MPI).
2. Download the SN data: `python scripts/download_data.py` (Pantheon+ and
   DES-SN5YR public releases; DES-Dovekie data is included in `data/`).
   The script installs the data under `packages/`; set `COBAYA_PACKAGES_PATH`
   to this repository's `packages` directory so that Cobaya's built-in SN
   likelihoods find it.
3. Run a single chain, e.g. DESY5 at `z_cut = 0`:
   `cobaya-run configs/fp_dy5_z0p00.yaml`
   (or `bash scripts/server_run.sh` on a server with MPI).
4. Summarize chains and reproduce the figure:
   `python scripts/analyze_scan.py && python scripts/plot_final_fig.py`.

The MCMC chains themselves (~1 GB) are not included in this repository;
`results/cpl_scan_summary.csv` contains the values reported in the paper.

## Data sources

- Pantheon+ (SH0ES Cepheid calibration): Brout et al. 2022, ApJ 938, 110
- DES-SN5YR (DESY5): Abbott et al. 2024, ApJ 973, L14
- DES-Dovekie: Popovic et al. 2025; DES 2026
- DESI DR2 BAO: Abdul-Karim et al. 2025, PRD 112, 083515
- Planck 2018: Planck Collaboration 2020, A&A 641, A6

## License

Code: [MIT](LICENSE). Data products retain their original licenses
(Pantheon+, DES, Dovekie public data releases).
