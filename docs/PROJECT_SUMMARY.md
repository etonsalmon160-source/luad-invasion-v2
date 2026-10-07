# Project Facts Ledger (PROJECT_SUMMARY)

> **Positioning**: this project's **sole up-to-date facts ledger** (authoritative facts ledger).
> **Scope**: only **two paired datasets**—`GSE308103`(snRNA) + `GSE307534`(Visium spatial).
> **Last updated**: 2026-09-17 (v3: paper QC definition 413,697 nuclei + M3-A full dimensionality reduction entering the checkpoint)
> Companion: [`WHITEPAPER.md`](WHITEPAPER.md) (technical roadmap) · [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) (checkpoints)

---

## 1. Data foundation (sample-by-sample verified on GEO)

| Dataset | True identity | Role | Measured scale |
| :--- | :--- | :--- | :--- |
| **GSE308103** | snRNA-seq (FFPE fixed RNA) | single-cell **reference**; the **only single-cell resource containing AAH** | raw **75 samples / 798,100 nuclei** (measured); **analysis definition 413,697 nuclei × 18,069 genes** (paper QC, 2026-09-16) |
| **GSE307534** | Visium CytAssist FFPE spatial | **spatial atlas** (in-situ coordinates); deconvolution target | GEO 56 samples / 25 patients; **local 56 sections** |

**Pairing (23 cases, P3–P25)**: sections exist in both modalities; the local spatial sections **fully cover** these 23 cases (P1/P2 are spatial only, with no snRNA, and are not paired).
> **2026-09-15 correction**: the old value "9 cases" was an artifact of the intersection when **only 19/56 sections had been downloaded**; after the sections were completed, the intersection of the two authoritative GEO tables is in fact 23 cases.
Stage (identical in both modalities): **Normal / AAH / AIS / MIA / IAC**.

**LNM**: for spatial lymph-node metastasis there is **no legitimate LUAD data**. `GSE190811`, verified on GEO, has the series title
*"…paired metastatic lymph node tumors in **breast cancer** patients"* (**breast cancer**),
and the GSM5732148 locked in whitepaper v1 **does not exist** in that repository (the real ones are GSM5732357–5732360) → **that section is void** and is forbidden for any LUAD product.

**Moved out of scope**: GSE131907 / GSE189357 / GSE148071 (three scRNA cohorts) → archived at `/home/eto/luad_invasion/luad_v2_out_of_scope/`.

---

## 2. Measured data characteristics (2026-09-12)

**GSE308103 (snRNA, nuclei)** — cohort-wide quantiles (1/5/25/50/75/95/99%):

| Metric | 1% | 5% | 25% | 50% | 75% | 95% | 99% |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nCount | 500 | 561 | 902 | **1,516** | 2,970 | 11,296 | 31,807 |
| nFeature | 367 | 452 | 697 | **1,082** | 1,786 | 4,234 | 7,085 |
| pct_mt | 0 | 0 | 0.26 | **0.6** | 1.24 | 3.13 | 6.44 |

> ⚠️ **Nucleus data** (low nCount, low mt%) → **whole-cell scRNA thresholds must not be copied over** (blindly using `nCount≥1000` would cut roughly 30% of nuclei;
> `mt<10%` barely filters). M1 itself uses **per-sample MAD outliers** + `pct_mt < 5`.
> ⚠️ **But from 2026-09-16 the final analysis definition has been ruled by the user to be the source paper's definition**: mask = M1 ∩ **paper absolute gate**
> (`nFeature≥500 & nCount≥1000 & pct_mt≤20`) = **413,697 nuclei × 18,069 genes**;
> the old 648,945-nucleus object is retained as a **sensitivity arm**. **The absolute gate is applied nonetheless; "must not be copied over" above refers to M1's own call—do not misread it.**

**GSE307534 (spatial)** — local **56/56** sections (completed 2026-09-15), stage coverage (authoritative GEO table, `IAC` = raw token `LUAD`):

```
Normal ( 1): P4_Normal
AAH    (11): P1_AAH, P2_AAH, P4_AAH, P4_AAH-1, P6_AAH, P9_AAH, P11_AAH, P20_AAH, P22_AAH, P24_AAH, P25_AAH
AIS    (14): P3_AIS, P5_AIS, P8_AIS, P9_AIS, P12_AIS, P14_AIS, P16_AIS, P17_AIS, P19_AIS, P21_AIS, P21_AIS-1, P22_AIS, P23_AIS, P23_AIS-1
MIA    ( 4): P10_MIA, P13_MIA, P15_MIA, P18_MIA
IAC    (26): P1_LUAD, P2_LUAD, P3_LUAD, P4_LUAD, P5_LUAD, P6_LUAD, P7_LUAD, P7_LUAD-1, P8_LUAD, P9_LUAD, P10_LUAD, P11_LUAD, P12_LUAD, P13_LUAD, P14_LUAD, P15_LUAD, P16_LUAD, P17_LUAD, P18_LUAD, P19_LUAD, P20_LUAD, P21_LUAD, P22_LUAD, P23_LUAD, P24_LUAD, P25_LUAD
```

> ⚠️ Sole defect: the tar for `GSM9226176` is **truncated** (56,272,384 B / should be 90,677,930 B; `gzip -t` reports unexpected EOF),
> missing `spatial/scalefactors_json.json` and `spatial/tissue_positions.csv`. Re-downloading was measured to yield the complete 87 MB tar (with only one section root, `P4_AAH2`).

---

## 3. Methodological facts (verified, driving tooling choices)

| Conclusion | Basis |
| :--- | :--- |
| **Visium spot is not single-cell** (55 µm mixes multiple cells) → **deconvolution** is required | platform definition |
| the reference **modality must match the section** (FFPE↔FFPE) → use **GSE308103** as reference | RCTD / benchmark literature |
| **sc↔sn is a system effect**, so `modality as batch` cannot be used | scvi-tools (SysVI was designed for this); Hrovatin 2025 |
| `scvi-tools 0.15.5` scVI parameter defaults are effective; **1.5.x requires Py≥3.10 (cannot install on this machine)** | measured |
| **RCTD `doublet_mode='full'` has no reject category, and `constrain=F`** (weights are not probabilities) | source-code measurement |
| `squidpy.nhood_enrichment` **permutes labels, returns only z, no p-value** | source-code measurement |
| **SCMG has no "state reversal" capability**; `generate_transition_cells` interpolates only between types **already present** in the reference manifold; the reference manifold has no tumor state | source code + measurement |
| **observational single-cell/spatial data cannot establish causality**; legitimate lever = **cis-MR + coloc** or perturbation experiments | methodological consensus |
| single-cell foundation models (scGPT/Geneformer…) **do not beat** scVI/Harmony/PCA baselines | Kedzierska 2025 *Genome Biol*; Ahlmann-Eltze 2025 *Nat Methods* |
| PLIP zero-shot determination of WHO growth pattern is **unvalidated**, and the spot scale is **insufficient** | literature + resolution argument |
| CMap's **Tau is a 0–1 reproducibility metric and cannot be negative**; `Tau≤-90` is a misuse | Subramanian 2017 |

---

## 4. Environment facts (hard constraints)

- **No GPU** (no `/dev/nvidia*`, `torch.cuda.is_available()=False`) → deep-learning/co-folding methods are CPU-only or require an external node;
- **Python 3.8.10, no conda** → `scvi-tools 1.5`, `cellcharter` cannot be installed;
- **shared library owned by root** → packages go into the personal library (`~/.local/lib/python3.8/site-packages`, `~/R/.../4.2`);
  installing `cellcharter` once downgraded torch to 1.12 and broke scvi (rolled back to torch 2.4.1+cu118, pl 1.5.10.post0, torchmetrics 0.7.3);
- **Network**: outbound access is restricted, dependency packages and model weights are obtained via mirrors;
  CRAN ✓, conda-forge ✓, PyPI via the Tsinghua mirror ✓.

**Installed and usable tools**:
- R: `CopyKAT 1.2.5`, `coloc 5.2.3`, `ieugwasr 1.1.0`, `TwoSampleMR 0.7.9`, `scDblFinder 1.12.0`, `spacexr 2.2.1`, `DESeq2`, `scater/scuttle`
- Python: `harmonypy 0.0.10`, `SpaGCN 1.2.7`, `plip` (weights already cached, 1.2 GB), `vina 1.2.7`, `fpocket` (CLI)
- To install: BANKSY, PRECAST, P2Rank, gnina, PoseBusters, OpenMM (GROMACS missing)

---

## 5. Past engineering errors (do not repeat)

See `/home/eto/luad_invasion/luad_v2_out_of_scope/` and the historical record of the v1 whitepaper: TD9 mislabeled AAH; `nLN` recorded as LNM; brain metastases/pleural effusion merged into IAC;
hard-coded values and fake curves; fabricated SCMG reversal factors; misuse of the CMap metric; GSE190811 treated as LUAD LNM.

**Principle**: every label/stage/number must have an authoritative source; **nothing enters the ledger unverified**; anything without a real source → **refused output**.

---

## 6. Related files

- Technical roadmap: [`WHITEPAPER.md`](WHITEPAPER.md)
- Checkpoints: `PLAN_AND_CHECKPOINTS.md`
- Parameter sources: [`PARAMETERS_AND_SOURCES.md`](PARAMETERS_AND_SOURCES.md)
- Target MR/coloc manual: [`M7B_MR_COLOC_TARGET_ANCHORING.md`](M7B_MR_COLOC_TARGET_ANCHORING.md)
- Spatial P0 prohibitions: [`spatial_cohort_and_figure_prohibitions.md`](spatial_cohort_and_figure_prohibitions.md)
- Out-of-scope archive: `/home/eto/luad_invasion/luad_v2_out_of_scope/` (with an explanatory README)
