# LUAD Early Progression · Paired Spatial–snRNA Atlas
## A paired spatial–snRNA atlas of the pre-invasive→invasive LUAD axis

> **Version**: v2 (rewritten 2026-09-12) · supersedes v1 (contained retracted results and invalid methods)
> **Companion documents**: [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) (milestones and gate conditions) · [`PROJECT_SUMMARY.md`](PROJECT_SUMMARY.md) (verified facts) · [`PARAMETERS_AND_SOURCES.md`](PARAMETERS_AND_SOURCES.md) (parameter sources)

---

## 0. One-sentence positioning

Using paired **Visium spatial** and **snRNA** data from the **same patients, same lesion, adjacent sections**,
we characterize the cell states and spatial niches along the lung adenocarcinoma **pre-invasive→invasive axis (Normal → AAH → AIS → MIA → IAC)**,
and on that basis derive **candidate targets**. Structural docking over those candidates was planned as
milestone M8 and has **not been run**; nothing in this repository reports a docking result.

**Technical positioning**: the contribution of this work is **data uniqueness + methodological rigor (benchmarking and validation)**, **not** the invention of a new causal-inference algorithm.

---

## 1. Scientific questions

1. **AAH is the only actionable pre-cancerous starting point**—how does it differ from the same patient's Normal/AIS in **cell state** and **spatial microenvironment**?
2. How does the **watershed of microinvasion (AIS→MIA)** manifest **spatially** as niche remodeling (malignant–stromal interface, immune exclusion)?
3. Which **candidate targets** simultaneously satisfy: specific to the malignant program/niche, **with genetic causal support**, and **structurally druggable**?

**Why pairing is essential**: cross-patient comparisons cannot distinguish "disease-stage differences" from "individual differences"; pairing (multiple stages from the same patient + dual modality) uses individual variation as its own control.

---

## 2. Data foundation (only two datasets, 23 paired cases)

| | **GSE308103** | **GSE307534** |
| :--- | :--- | :--- |
| Modality | **snRNA** (nuclei; FFPE) | **Visium spot** (FFPE CytAssist; 55 µm, **not single-cell**) |
| Role | single-cell **reference**; the **only single-cell resource containing AAH** | **spatial atlas** (in-situ coordinates); **deconvolution target** |
| Scale | 75 samples; **798,100 nuclei** (measured) | GEO 56 samples / 25 patients; local **56 sections** |
| Stage | Normal / AAH / AIS / MIA / IAC | Normal / AAH / AIS / MIA / IAC |

**Paired patients (23 cases, P3–P25)** — each case has sections in **both modalities**, and the local spatial sections **fully cover** these 23 cases.
> **2026-09-15 correction**: the old value "9 cases" was an artifact of the intersection when **only 19/56 spatial sections had been downloaded**; after the sections were completed, the intersection of the two authoritative GEO tables is in fact 23 cases.

**Key concept**: each Visium spot contains **multiple cells** (mixed signal) → a single-cell reference must be used to **deconvolve** it in order to obtain a spot's cell composition. The role of pairing is to **match the reference modality and match the patient**, thereby making deconvolution trustworthy; it does **not eliminate** the deconvolution step.

**LNM**: for spatial lymph-node metastasis there is **currently no legitimate LUAD data** (`GSE190811` was verified via GEO to be **breast cancer**, retracted). LNM is not treated as a spatial stage; it will be added if real data become available later.

---

## 3. Technical roadmap

```
M0  Input freeze      Freeze the two datasets' samples/patients/stages + hash
M1  QC / doublets     per-sample adaptive MAD + scDblFinder; **the analysis mask additionally overlays the source paper's absolute gate**
                    (`nFeature≥500 & nCount≥1000 & pct_mt≤20`) → 413,697 nuclei (from 2026-09-16, see PLAN §2.0 line 2)
M2  Malignancy confirmation  CopyKAT per-sample CNV (+ infercnv subset cross-validation)
M3  Integration       scVI(batch=sample_id) to build the snRNA atlas (the paper's batch=dataset is retracted, see PLAN §C)
                    + scArches/scANVI cross-modality label transfer; full scIB panel
                    Control arm: pure SCMG zero-shot integration/manifold (orthogonal, produces no conclusions)
M4  Cross-modality AAH  same-patient paired consistency (five criteria)
M5  Spatial deconvolution  RCTD, reference = GSE308103 (modality-matched); gate uses RCTD's native output
M6  Spatial niche     BANKSY (stability-based clustering) + Squidpy (category labels + empirical p)
M7  Targets           a niche prognostic signature (patient-level → TCGA survival)
                    b candidate target pool (cis-MR + coloc genetic anchoring + druggability)
                    c CMap (only when real LINCS data are available)
M8  Structural docking  fpocket+P2Rank → Vina → gnina rescoring → PoseBusters
                    (MD deferred)
```

### M1 · QC: **snRNA thresholds must be set empirically** —— ✅ **Completed** (2026-09-12)
In this data, median nCount 1,516 and median pct_mt 0.6% (**nucleus** characteristics). Copying whole-cell scRNA thresholds
(nCount≥1000) would **cut roughly 30% of nuclei**, whereas mt<10% barely filters. Hence:
**nCount/nFeature per-sample MAD outliers** (`scuttle::isOutlier`, nmads=3, log1p, two-sided) + **pct_mt < 5** (nucleus convention).
Doublets: **scDblFinder per sample** (pooling forbidden).

**Measured results** (`GSE308103`, 75 samples / 798,100 nuclei):

| Metric | Value |
| :--- | :--- |
| QC pass | **767,839 (96.21%)**; removed 30,261 (nCount outliers 8,919 · nFeature outliers 1,077 · mt≥5% 12,918) |
| Doublets | **118,894 (15.48%)**; per-sample median 10.60%, range 4.43–29.03% |

**Robustness (two sensitivities)**:
- **nmads 3 vs 5**: +1.97 pp (per-sample median 2.20 pp) → threshold choice is **insensitive**;
- **doublet-removal fraction**: overlap with a fixed top-10% is about **62.3%** → the call is **fairly sensitive to the rate assumption**, so this result **must be paired with a downstream "remove/keep" sensitivity**.

> ⚠️ **2026-09-16 definition change (the text above in this section still holds, but is no longer the final analysis definition)**: the user ruled that the pipeline be **rebuilt according to the source paper's fixed QC**.
> Final analysis mask = **M1 (`qc_pass & singlet`) ∩ paper absolute gate** (`nFeature≥500 & nCount≥1000 & pct_mt≤20`)
> = **413,697 nuclei × 18,069 genes**. Chain: 798,100 → 767,839 → 648,945 (∩ single cells) → paper gate 555,480 → **intersection 413,697**.
> That is, **the intuition that "absolute thresholds cannot be copied over" is overridden by overlaying the paper's definition**: the absolute gate is applied nonetheless. The old 648,945 object is retained as a **sensitivity arm**.
> The doublet caller differs from the paper by a **known deviation**: the paper uses Scrublet, this project uses **scDblFinder** (in compliance with R2), registered truthfully as such.

**Two independent validations of the doublet call** (this environment cannot perform second-method cross-validation—the detectable fraction of scrublet on sparse nuclei is only ~0.5%;
DoubletFinder requires Seurat 2/3 or 5, and this machine has 4.3.0):
> **A · Count features**: in **75/75 samples** the doublets have significantly higher counts than singlets (nCount ratio **median 2.39**, nFeature median 2.05—consistent with the expected ≈2× for true doublets);
> doublet rate vs sample cell number **r=0.921** (consistent with the 10x loading relationship).
> **B · Cross-lineage co-expression** (independent of scDblFinder): the proportion of doublets co-expressing mutually exclusive lineage markers (EPCAM+ & PTPRC+)
> is **7.2 times** that of singlets (median 3.15% vs 0.49%; **consistent in 72/75 samples**) → the call indeed enriches the true-doublet fingerprint,
> **ruling out "high-count single cells misclassified as doublets".**

> **On the doublet-rate figures**: the pooled 15.48% is **weighted by cell number** (pulled up by large samples); the "typical sample" rate is the **median 10.6%**.
> `dbr=NULL` derives the rate **from cell number** → report the **per-sample presentation** and note it as a model estimate.

**One resolved anomaly**: `P7_LUAD` was once called as 0 doublets; investigation showed this was **the scDblFinder default `xgb` classifier collapsing during training on that sample**
(not a biological cause; the patient's other two samples were normal). Switching to `score="weighted"` yielded 11.69%, consistent with the other samples;
**the main script now has an automatic fallback** (xgb calls 0 → switch to weighted), eliminating this class of silent failure.
**Detailed report**: [`results/01_qc/M1_validation_report.md`](../results/01_qc/M1_validation_report.md).

### M2 · Malignancy confirmation: **CNV is authoritative**
`CopyKAT` per sample; `infercnv` demoted to **optional cross-validation** on 5–10k cells/sample (officially unmaintained, requires JAGS).
**Must not** be replaced by the argmax of pan-epithelial markers (EPCAM/KRT…).

### M3 · Integration: **do not treat "modality as batch"**
sc↔sn is a **system** effect as defined by scvi-tools, and modality is collinear with dataset → a conditional VAE would **under-correct**.
Correct approach: **within snRNA** use `scVI(batch=sample_id)` to build the atlas (**`patient_id` must never be used as batch**—stage is nested within patient);
**across modalities** use **scArches/scANVI label transfer** (not "zero-shot").
Evaluation uses the **full scIB panel**: batch removal (kBET + iLISI + graph-connectivity + PCR) **and** biological conservation (cLISI/ARI/NMI/ASW) **reported together**;
do not judge success on iLISI↑ alone (it can be inflated by over-integration).
**SCMG control arm**: performs only zero-shot cross-dataset integration + manifold + state characterization; **outputs no reversal/causality** (that capability does not exist).

### M4 · Cross-modality AAH (**all five criteria are required**)
① each stage is distinguishable; ② overlapping stages are consistent sn↔spatial; ③ the platform offset δ(stage) is stable; ④ **pairing within sn** (AAH vs the same patient's Normal/AIS) is concordant; ⑤ AAH identity is confirmed by **CNV/markers**.
**Not passing the gate → AAH may only be labeled a hypothesis.**

### M5 · Spatial deconvolution: **the gate must use RCTD's native output**
`RCTD(doublet_mode='full')`, reference = **GSE308103** (modality-matched).
**Forbidden** to treat Σ=1 after `normalize_weights()` as a gate—that normalization is imposed by the project itself; the `full` mode has **no** reject category.
The gate moves to: non-negative weights, a reasonable per-spot distribution, modality-matched reference, and reporting the number of removed spots by an **explicit threshold** on `UMI_min`/`counts_MIN`.

### M6 · Niche: **stability defines clusters, not a single BIC**
**BANKSY** (peer-reviewed, *Nat Genet* 2024) performs spatial domain/niche detection; the number of clusters is decided by **stability/consensus** (not a BIC minimum).
Spatial statistics use **category labels** and report **z + empirical p=(b+1)/(n+1)**, recording `n_perms`.
**Forbidden**: writing `nhood_enrichment` as `P<0.001` (that function **does not return a p-value**; it permutes the labels); claiming permutation was done for `co_occurrence` (it **does not** do it); a distance step smaller than the Visium ~100 µm spot pitch. Enrichment tests report **effect size** (Fisher FDR is inflated by compositionality + large N).

### M7 · Targets: **genetic statistical anchoring is the legitimate lever for "causality"**
- **M7a Prognosis**: niche signature **aggregated at patient/slide level** → TCGA-LUAD bulk scoring → **multivariable Cox (adjusted for stage/age/sex) + KM**; guard against overfitting (penalization/cross-validation). Wording = **prognostic association**. Spatial pairing has only **23 cases** → **no claim of patient subtyping**.
- **M7b Candidate target pool**: sources = malignant program/regulon **+** niche signature → **cis-MR + coloc** (LUAD GWAS × lung eQTL) → **Open Targets druggability + DepMap selective dependency + clinical-stage drug matching**. See [`M7B_MR_COLOC_TARGET_ANCHORING.md`](M7B_MR_COLOC_TARGET_ANCHORING.md).
  Output `genetic_support ∈ {supported, not_supported, not_testable}`; **no post-hoc parameter tuning**, and `not_testable` is **truthfully marked as missing**.
- **M7c CMap**: executed only when **real LINCS data** are obtained and the correct metric (**NCS**, not negative Tau) is used; otherwise **no output**.

### M8 · Structural docking
`fpocket + P2Rank` consensus pocket detection → `AutoDock Vina 1.2.7` → `gnina` rescoring → `PoseBusters` physical-plausibility filtering;
when no crystal structure exists, use co-folding (Boltz-2, **requires GPU**). **MD (100 ns) deferred** (this machine has no GPU/GROMACS).
Docking conclusions are only **computational hypotheses**; "validated" requires wet-lab experiments.

---

## 4. Honest capability boundaries (language conventions)

| May claim | May not claim |
| :--- | :--- |
| cell states / niche architecture at each stage (`we characterize`) | causality (`causal`) |
| spatial localization to the invasion front | `driver` (unqualified) |
| **candidate regulators / candidate targets** (`candidate`) | `validated target` |
| **genetically supported candidate** (`genetically supported candidate`) | "state-reversal factor" (method does not exist) |
| prognostic **association** | prognostic prediction / clinically actionable |
| docking = **computational hypothesis** | docking = validation |

**Legitimate levers for causality**: only **human genetics (cis-MR + coloc)** or **perturbation experiments**; observational single-cell/spatial data **cannot** establish causality.

---

## 5. Validation strategy

1. **Cross-modality paired consistency** (M4): same-patient AAH vs Normal/AIS, cross-validated in both the sn and spatial modalities;
2. **Cross-cohort**: if an independent LUAD cohort is added later, it serves only as **external validation** and is not merged into the main atlas;
3. **Internal consistency**: malignant-cell CNV positivity ↔ canonical markers ↔ spatial localization all agree;
4. **Negative controls / sensitivity analyses**: label permutation, threshold sensitivity (nmads, mt, clustering stability);
5. **Primary/exploratory analyses pre-specified**: to avoid p-hacking.

---

## 6. Explicitly removed / demoted (relative to v1)

| Item | Disposition | Reason |
| :--- | :--- | :--- |
| SCMG "state reversal / reversal factor / causal gene" | **Removed** (only the integration/manifold control is retained) | the source code has no such capability; the reference manifold has no tumor state |
| three scRNA cohorts (GSE131907/189357/148071) | **moved out of scope** (archived at `/home/eto/luad_invasion/luad_v2_out_of_scope/`) | unrelated to the paired main line |
| spatial LNM (GSE190811) | **Removed** | verified to be **breast cancer** |
| CMap (Tau≤-90) | **demoted to optional** | wrong metric (Tau cannot be negative) + no LINCS data |
| WES / TMB / TCGA molecular subtyping | **Removed** | no data / retracted |
| PLIP zero-shot determination of WHO growth pattern | **Removed** | the spot scale (45–90 px ≈ 45 µm) is **insufficient** to determine mm-scale architecture |
| old figures (14 clones, 107,796 spots, `Z<-5.8`, 316,689 cells…) | **void** | old engineering artifacts |
| CellCharter | demoted (cannot install on py3.8) | switched to BANKSY |

---

## 7. Environment constraints (affecting tooling choices)

See [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) §5b. Key points:
**No GPU** (Boltz-2/cell2location/SysVI infeasible or require an external node); **Python 3.8** (scvi-tools 1.5/CellCharter cannot be installed);
shared library is owned by root (packages always go into the personal library); outbound network is restricted, dependencies are obtained via mirrors.

---

## 8. Risk register

| Risk | Mitigation |
| :--- | :--- |
| sn↔spatial cross-modality not transferable | same-patient pairing provides **internal consistency** validation (M4); low-confidence labeled `unassigned` |
| limited number of spatial patients (23 local paired cases) | patient-level aggregation; large TCGA cohort as a fallback; **no claim of subtyping** |
| RCTD errors for rare types <25 cells | merge rare types first or lower `CELL_MIN_INSTANCE` and file a record |
| false positives in compositional proportion statistics | effect size + permutation tests, not bare FDR |
| GWAS needed for MR may be access-controlled | switch to an OpenGWAS lung-cancer GWAS and annotate the source and sample size |
| no GPU limits structure prediction | docking uses Vina/gnina (runs on CPU); co-folding awaits a GPU node |

---

## 9. Deliverables (figure planning)

1. **F1** Paired design + data foundation (**23** paired patients × dual modality × stage)
2. **F2** snRNA atlas and malignant program (CNV confirmation → state/trajectory → regulon)
3. **F3** Spatial deconvolution and niche (RCTD → BANKSY niche → invasion front)
4. **F4** Cross-modality AAH consistency (five criteria)
5. **F5** Candidate target pool (MR/coloc + druggability) and structural docking

> **Principle**: **nothing fabricated that was not computed**; every figure carries its parameters, n, and hash.

---

## 10. Mapping to PLAN

| This whitepaper | PLAN milestone |
| :--- | :--- |
| §3 M0–M1 | M0 / M-1 / M1 |
| §3 M2 | M2 |
| §3 M3 | M3 |
| §3 M4 | M4 |
| §3 M5 | M5 |
| §3 M6 | M6 |
| §3 M7a/b/c | M7 |
| §3 M8 | M8 |
