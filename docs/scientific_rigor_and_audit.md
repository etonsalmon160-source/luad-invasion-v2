# Scientific Rigor & Audit Protocol

> Applies to this project (paired spatial–snRNA atlas). **All scripts, reports, and figures must comply.**
> For scope and milestones see [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) and [`WHITEPAPER.md`](WHITEPAPER.md).

---

## Rule 0 · Data identity and stage ground truth
1. Dataset identity follows **GEO/GSA**; this project has **only two paired datasets**: `GSE308103`(snRNA) + `GSE307534`(spatial);
2. **No silent default for stage** (fallbacks such as `.get(x,'IAC')` are strictly forbidden); any unknown token must raise;
3. **Malignant labels must be confirmed by CNV** (CopyKAT), and **must not** be replaced by the argmax of pan-epithelial markers (EPCAM/KRT…);
4. Doublets must use a **standard algorithm** (scDblFinder), and **must not** be replaced by a nCount+nFeature heuristic.
5. **No inference of lesion ordinal**: a patient's second lesion **must not be guessed from a token** — the two datasets have **different** naming conventions
   (`GSE307534` uses `AAH-1`/`AIS-1` (**with a hyphen**); `GSE308103` uses `AAH1`/`AIS1`/`Normal1`/`LUAD1` (**without a hyphen**)).
   the ordinal is **looked up only in the authoritative GEO table** (`cohort_registry.resolve_lesion_ordinal`); if it is not found, raise.

---

## Rule 1 · Zero tolerance for fabricated and simulated data
1. **Strictly forbidden** to fabricate cell proportions, expression, or spatial abundance with `np.random` / `runif` / `rnorm` / `sample()` / `make_blobs()` / hard-coded dummy values;
2. the input to every figure/statistic must come from a **real sequencing matrix** with a clear source;
3. a fabrication script found → **isolated immediately**, never left on the production trunk.

---

## Rule 2 · Biological specificity (marker cross-validation)

> 🔴 **2026-09-17 rebuilt (user instruction: "markers must come from authoritative papers, don't dig them out of some obscure corner")**.
> The old version of this table (six rows: EPCAM/KRT7/19…) **had not a single source**, flagged as missing in audit. It has now been rebuilt from **primary literature**,
> with per-gene sources in [`05_annotation/marker_panel.py`](../05_annotation/marker_panel.py).
> **The operative definition = that file, not this table**; this table is a human-readable summary.

| Lineage | Positive markers | Expected negative | Source (primary literature) |
| :--- | :--- | :--- | :--- |
| Epithelial | EPCAM, KRT8/18/19, CDH1, NKX2-1, SFTPC, SFTPA1/A2, SFTPB, NAPSA, AGER, CAV1, PDPN, SCGB3A2, SCGB1A1 | PTPRC, CD3D, PECAM1, COL1A1 | Travaglini 2020 *Nature* 587:619 (AT1/AT2/Club/pan-epithelial); Vieira Braga 2019 *Nat Med* 25:1153 (airway); **Peng 2026 *Cancer Cell* (source paper, LUAD lineage NKX2-1)** |
| T/NK | CD3D/E/G, TRAC, CD4, IL7R, CD8A/B, NKG7, GNLY, KLRD1, PRF1, GZMB | EPCAM, COL1A1, CD68 | Travaglini 2020; Vieira Braga 2019; Guo 2018 *Nat Med* 24:978 (NSCLC T cells) |
| B/plasma cells | MS4A1, CD19, CD79A/B, MZB1, JCHAIN, SDC1, IGHG1, IGKC, XBP1, DERL3 | CD3D, EPCAM, ACTA2 | Travaglini 2020; Vieira Braga 2019 |
| Myeloid | LYZ, AIF1, ITGAX, CD68, CD163, MSR1, C1QA/B/C, MARCO, APOE, FCN1, CD14, S100A8/9, SPP1 | CD3D, EPCAM, PECAM1 | Travaglini 2020 (alveolar macrophages); Habermann 2020 *Sci Adv* 6:eaba1972 (SPP1⁺ macrophages); Zilionis 2019 *Immunity* 50:1317 (lung tumor myeloid) |
| Fibroblasts | COL1A1/A2, COL3A1, DCN, LUM, FN1, PDGFRA/B, ACTA2, TAGLN, FAP, CXCL12 | PTPRC, EPCAM, PECAM1 | Travaglini 2020; Habermann 2020; Reyfman 2019 *AJRCCM* 199:1517; Lambrechts 2018 *Nat Med* 24:1277 (tumor stroma) |
| Endothelial | PECAM1, CDH5, KDR, CD34, VWF, EGFL7, RAMP2, EMCN, PLVAP, AQP1, CLDN5, FLT1 | EPCAM, PTPRC, COL1A1 | Travaglini 2020 (EC subtypes); Gillich 2020 *Nature* 586:785 (alveolar capillary aCap/gCap specialization); Lambrechts 2018 |

**Malignant identity is determined by CNV**; markers serve only as **corroboration**, not as a criterion.

**Independence from GP6 Standard B (important)**: this table is the gene source for GP6 **Standard A**, **deliberately not taken from the HLCA integrated atlas**
(Sikkema 2023 *Nat Med* 29:1563)—because CellTypist's `Human_Lung_Atlas.pkl` was trained from exactly that atlas,
so if the two shared a source, κ would become self-validation. **Residual non-independence, truthfully declared**: HLCA integrates several of the primary studies cited in the table above,
so the two are **not statistically independent**; the only guarantee is that the "marker-definition sources differ".

---

## Rule 3 · Parameter transparency and statistical rigor
1. parameters are **explicitly made into constants** and registered in [`PARAMETERS_AND_SOURCES.md`](PARAMETERS_AND_SOURCES.md) (with source and verification status);
2. **once a threshold is fixed it must not be adjusted post hoc** (to prevent p-hacking);
3. single-cell differential expression uses **patient-level pseudobulk** (DESeq2), not single-cell-level Wilcoxon (pseudoreplication);
4. compositional proportion data report **effect size**, not bare FDR;
5. spatial statistics use **category labels** and report **z + empirical p** (`squidpy.nhood_enrichment` **does not return a p-value**; `co_occurrence` **does not permute**);
6. survival analyses must output HR, 95% CI, and Log-rank P, and **tampering with significance is strictly forbidden**.

---

## Rule 4 · Causality and wording boundaries
1. **observational single-cell/spatial data cannot establish causality**; legitimate lever = **cis-MR + coloc** or perturbation experiments;
2. may say: **candidate / genetically supported candidate / prognostic association / computational hypothesis**;
3. **may not say**: `causal`, `driver` (unqualified), `validated target`, "state-reversal factor";
4. docking (M8) conclusions are only **computational hypotheses**; "validated" requires wet-lab experiments.

---

## Rule 5 · Cross-modality (core to this project)
1. spatial (Visium spot, multi-cell mixture) **must be deconvolved**; the reference must be **modality-matched** (FFPE↔FFPE → use GSE308103);
2. cross-modality consistency must pass the **M4 five criteria**; **not passing the gate → AAH may only be labeled a hypothesis**;
3. `modality as batch` is forbidden (sc↔sn is a system effect and is collinear with dataset).

---

## Rule 6 · Freezing and reproducibility
1. artifacts are **reproducible + hashed**; deterministic (explicit seed, no randomness);
2. `patient_id` (true patient) and `sample_id` (section/sample) are **stratified**;
3. inputs and outputs are traceable; the freeze manifest and verification report are submitted with each milestone.

---

## Audit workflow (before artifacts enter the ledger)
1. **Data-source verification**: path valid, source clear;
2. **Count/label consistency verification**: aligned with the freeze table;
3. **Reproducibility verification**: deterministic script + explicit seed + artifact hash;
4. **Wording and methodology comparison**: against this protocol and top-journal standards;
5. **Guards must be falsified**: every new or modified consistency check must be fed a **deliberately corrupted input**, confirming that it really does go red
   (non-zero exit code, and naming the specific item). **Merely "running through" is not validation — a guard that always passes is fake.**
   Example: `00_ingest/03_verify_cohort_consistency.py` self-validates using a bad table with a deleted rule, naming 5 samples and exiting 1.
6. **No fact may have two implementations**: wherever a second derivation of the same fact exists (**especially "the one that is never called"**),
   either delete it or byte-for-byte collide it against the authoritative source. **A defect with no execution path produces no symptom and never exposes itself**,
   and unit tests/pipelines/artifact verification all fail to cover it (example: `lesion_ordinal()`,
   see [`results/03_cnv/GP1_report.md`](../results/03_cnv/GP1_report.md) §6.1).
