# Sequencing depth is tissue density in the LUAD invasion sequence

Code, analysis and derived results for the manuscript:

> **Sequencing depth is tissue density in the lung adenocarcinoma invasion sequence:
> niche domains and *in-silico* reversal of the disease state**
> Zhiyang Li, Jiaxuan Yang
> Guangdong Medical University · Monash University Malaysia

This repository is the complete analytical record behind that paper. It is a self-contained
study; the directory it was originally developed in is retained only as a read-only data
source and is not part of this project.

---

## What the study is about

Lung adenocarcinoma (LUAD) progresses through a defined morphological sequence, from
atypical adenomatous hyperplasia (AAH) through adenocarcinoma *in situ* (AIS) and minimally
invasive adenocarcinoma (MIA) to invasive adenocarcinoma (IAC). We profiled paired
single-nucleus RNA sequencing and Visium spatial transcriptomics from the same patients
across all five stages, and asked how the transcriptional state of the tissue changes as
invasion proceeds, and whether that state can be reversed *in silico* by known
perturbations.

The organising finding is that **sequencing depth and tissue density are the same
measurement in this cohort**. Domain identity, marker stability, the apparent
stage-dependence of the spatial programmes, the cohort-level copy-number difference and one
domain's prognostic weight all scale with depth. Several apparently positive results
disappeared once depth was matched, and the paper reports those negatives alongside the
positives. Four results survive; they are listed in the manuscript's Discussion.

## Data

All data are public. Nothing in this repository redistributes controlled-access material.

| Dataset | Accession | Modality | Role |
| :--- | :--- | :--- | :--- |
| snRNA-seq | **GSE308103** | fixed RNA profiling (probe panel) | single-cell reference, 413,697 nuclei analysed |
| Visium | **GSE307534** | CytAssist 11 mm (probe panel) | spatial atlas, 56 sections |
| TCGA-LUAD | archived distribution | bulk RNA + clinical | independent prognostic filter |
| LINCS L1000 | **GSE70138** | compound perturbation | queried through a locally built database |
| SCMG | `xingjiepan/SCMG_data` | gene perturbation (MIT) | second, parallel perturbation library |

23 patients are shared between the two primary modalities. Three further external resources
were evaluated and rejected; they are documented in the manuscript's Methods rather than
listed here.

## Layout

```
├── 00_ingest/          authoritative cohort registry
├── 01_qc/              quality control and doublet detection
├── 02_expression/      expression object assembly
├── 03_cnv/             copy-number inference and the malignancy probes
├── 04_integration/     integration and the traditional (Seurat) pipeline
├── 05_annotation/      marker panels and cell-type annotation
├── 06_scmg/            SCMG subsystem
├── 07_he_pathology/    H&E imaging line, including the slide-registration algorithm
├── 08_spatial_deconv/  RCTD deconvolution and the spatial copy-number arm
├── 09_trajectory/      developmental-trajectory arm
├── 10_niche/           spatial domains, prognosis, and the reversal arm
├── docs/               method provenance, parameter sources, audit records
├── results/            derived outputs (large binaries are git-ignored)
└── figures/            figure source data
```

**Start here if you are new to the project:**
1. `docs/PARAMETERS_AND_SOURCES.md` — every analysis parameter, labelled by whether it comes
   from a source paper, a software default, or from us.
2. `docs/PROJECT_SUMMARY.md` — the verified-facts overview.
3. `10_niche/` — the reversal arm, which is where most of the recent work lives.

## Reproducibility notes

- Scripts are numbered in execution order within each directory.
- Large intermediate objects (`.h5`, `.npz`, `.rds`, `.bin`, model weights) are git-ignored
  and regenerable from the scripts; `.gitignore` records, for each exclusion, why it is
  excluded and what regenerates it.
- Alongside every analysis that produced a negative result, the repository keeps the
  diagnostic that established the failure was real rather than a broken method.

## What is not here

- Structural docking, molecular dynamics and co-folding were **not** run. They appear in the
  planning documents as future work and produced no result.
- No malignant cell label is assigned anywhere in this project. Copy-number output is used
  as a continuous quantity only.

## License

MIT. See [`LICENSE`](LICENSE).
