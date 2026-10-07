# LUAD v2 · Planning Guide and Strict Checkpoints (PLAN & STRICT CHECKPOINTS)

> **Nature**: This file is the **redone planning guide**. Each milestone has **hard pass-through conditions (gate)**,
> and **the next step must not be started without passing the gate**.
> **Companion documents**: [`docs/WHITEPAPER.md`](docs/WHITEPAPER.md) (technical route) · [`docs/PROJECT_SUMMARY.md`](docs/PROJECT_SUMMARY.md) (master statement of facts) · [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) (sources of parameters).

---

## 0. Progress Reset Declaration (PROGRESS RESET)
- **All result-producing artifacts are void** for the old project `luad_invasion` (figures/values/tables/TMB/CMap/docking); see `PROJECT_SUMMARY` for the reasons.
- **This project's progress is reset to zero**, restarting from `M0`. The old directory serves only as a **read-only data source and historical reference**.
- The progress board is in §6 of this file (M0 already complete; from M-1 onward at 0%).

---

## 1. Iron Laws (hard constraints for every step)
1. Data identity is governed by **GEO/GSA**; **no silent default for stage** (code must not contain fallbacks like `.get(x,'IAC')`).
2. **Malignancy labels must be corroborated by CNV**; **doublets use scDblFinder**.
3. **No real source = do not compute** (prefer to refuse/report missing; never fabricate, hard-code, or np.random).
4. Outputs **reproducible + hashed**; `patient_id` (true patient) and `sample_id` (tissue/slice) are **layered**.
5. **The SCMG branch does not mix in traditional algorithms**; the standard branch and the SCMG branch run **in parallel** and are compared at the end.
6. **Capability boundary of the SCMG branch** (verified 2026-09-12): it does only **zero-shot cross-dataset scRNA integration + manifold + state characterization**;
   it **must not** output "state reversal / reversal factors / causal genes" (the source code has no such capability; the reference manifold has no tumor state). Targets are **candidates** and must not be called causal.

---

## 2. Milestones and Hard Checkpoints (Gates)

### 2.0 Execution Order (restructured 2026-09-16 around the **traditional single-cell workflow**)

> **Background**: This file originally placed **M2 (CNV) before M3 (dimensionality reduction and clustering)** — treating CNV as a "leading, main corroboration step".
> The actual order in the traditional single-cell workflow is the opposite: **normalize → dimensionality reduction → clustering → marker annotation, then CNA refinement on the annotated epithelial clusters**.
> After this restructuring, **the M numbers are downgraded to "module identifiers"** (`M2` denotes only the "CNV module", no longer implying order of execution);
> **the execution order is governed by the table below**, and GP numbers are the stable checkpoint numbers.

> 🔴 **2026-09-16 second restructuring (user decisions "rebuild with the paper's QC" + "align with the paper's full set")**: rows 2–5 fully rewritten.
> **GP0 redone with the paper's QC caliber** (413,697 nuclei, not the original 648,945); **GP4a–4c switched to the source paper's Seurat recipe**
> (SCTransform v2 → HVG 3000 → PCA 50 → **Harmony as the main caliber** → FindClusters/Louvain).
> The old scran / HVG-2000 / PCA-30 / kNN-15 / Leiden and the `{20,30,50}` ARI guardrail are **all void** (registration in PARAMETERS §M3-A.1).
> See [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) §M3-A.0–A.3 for details.

> 🔴 **2026-09-16 third restructuring (user decision: CNV is downgraded to "fine malignant calling of epithelial subclusters")**: rows 7–10 rewritten, GP8 split into **GP8a/GP8b**.
> **CNV's role shrinks from an "independent corroboration gate" to a "fine discriminator that labels subclusters as malignant"**, and therefore **must be placed after epithelial subclustering**
> (previously placed after GP6, it could only produce per-cell readouts and could not be assigned to subclusters). **Subclustering of the immune/stromal/endothelial lineages (GP8b) is no longer blocked by CNV scheduling**.
> **GP1 (CopyKAT smoke test) consequently sinks into an execution detail of GP2** and is no longer listed separately.

| No. | Checkpoint | Content | Sub-gate |
| :---: | :--- | :--- | :--- |
| 1 | ~~M1~~ | QC / doublets (per-sample adaptive MAD) | ✅ Passed (🔶 downgraded to **sensitivity arm**) |
| 2 | **GP0** | Expression object rebuild — **🔴 redone 2026-09-16 with the paper QC**: `nFeature≥500 & nCount≥1000 & pct_mt≤20` ∩ scDblFinder singlets → **413,697 nuclei × 18,069 genes** | 🔶 Redone; **gate validation completed by the 2026-09-17 audit (13 PASS / 0 FAIL)**, awaiting your sign-off |
| 3 | **GP4a** | **Preprocessing (paper recipe)**: `SCTransform(vst.flavor="v2", variable.features.n=3000)` | ✅ Full run done (2026-09-17) |
| 4 | **GP4b** | `RunPCA(npcs=50)` — **no ARI guardrail** (the old guardrail is void along with the old grid) | ✅ Full run done (2026-09-17) |
| 5 | **GP4c** | Batch: **Harmony on `sample_id`, `dims.use=1:50` (main caliber, paper recipe)** / no correction (**sensitivity arm, not the paper recipe**) | ✅ **2026-09-17 user ruling: adopt the paper's Harmony main caliber** (with a non-attributability caveat) |
| 6 | **GP5** | Resolution selection: **`{0.5,0.6,0.7,0.8}` × 5 seeds** + metrics **1/2/3/4** (including the **AAH absorption guardrail**) + pre-registered `r*` rule (manually signed **before annotation**) | ✅ **2026-09-17 sign-off: L1 `r*=0.6`**; metric 4 acknowledged as ineffective at full scale. **2026-09-21 checkpoint closed**. 🔴 **2026-09-22 L1 caliber switched to `A_adjudicated` ⇒ six lineages re-clustered and `r*` re-selected**: current values **epithelial 0.7 (relaxed) / T·NK 0.8 / myeloid 0.5 / endothelial 0.7 / B·plasma 0.6 / fibroblast 0.8 (relaxed)**, six JSONs **all signed by the user personally**; metric 3 **after recomputation only myeloid fails the line** (user ruled to keep it, registered as a defect); **after seed recomputation T·NK→3, fibroblast→4, the rest remain 0**. See `results/05_annotation/GP5_report.md` §3.1 / §6.8 / §12.1 for details |
| 7 | **GP6** | **Dual-standard annotation** (Rule 2 marker ⊕ CellTypist) + Cohen's κ — **responsible only for identifying "which cells are epithelial"** | 🔶 Run, **not signed off** (`GP6_report.md`, 413,697 nuclei) |
| 8 | **GP8a** | **Epithelial-lineage subclustering** (restart the full pipeline within the subset, same caliber as GP4a–c) → **epithelial subclusters** | ✅ **Completed 2026-09-22** (report named `GP8c`; 27 subclusters / 133,384 nuclei / 0 removed) — ⚠️ **not signed off** |
| 9 | **GP2** | ~~**CNV fine calling**: per-sample CopyKAT, input = the **epithelial cells** determined by GP6; conclusions aggregated by GP8a's **epithelial subclusters** → per-subcluster malignant calling~~ | ⛔ **Voided 2026-09-24** (single-cell CNV failed) → switch to **spatial SC0–SC4**, placed after M5 |
| 10 | **GP8b** | Subclustering of the remaining **5 lineages** (non-epithelial, unrelated to CNV) + hierarchical self-consistency | ✅ **Completed 2026-09-22** (same batch as GP8a, see `GP8c_report.md`; six lineages mutually exclusive and complete, 413,697 nuclei) — ⚠️ **not signed off** |
| 11 | **GP3 / GP7** | Metric backend freeze (GP3) → scIB dual-panel evaluation (GP7) | ⬜ |
| 12 | **GP9** | **SCMG control arm** (M3-B): zero-shot integration + global manifold + state characterization, **with no mixing of traditional algorithms** (iron law 5) | ⬜ |
| 13 | M4… | Cross-modal AAH → M5 → M6 → M7 → M8 | ⬜ |

**Key Ordering Constraints**:
- **GP2 (CNV) must come after GP8a (epithelial subclustering)** — third adjustment, 2026-09-16. CNV's **sole purpose** in this project
  is **to finely call epithelial subclusters as malignant**, and the object "subcluster" **does not yet exist** before GP8a. Running it right after GP6
  yields only per-epithelial-**cell** readouts that cannot be attributed to subclusters, effectively leaving the most needed step for downstream guessing.
  ⇒ **Order: GP6 defines epithelium → GP8a defines subclusters → GP2 defines subcluster malignancy**.
- **GP2 must not give "malignant clone"-level conclusions** — nuclei have low UMI and contain intronic/ambient RNA, so CNV can support only
  **subcluster-level (aneuploid vs not)** calls. See §M2 for the wording ceiling.
- **GP5 must be signed before GP6** — resolution must not be fitted after the fact.
- **Why GP8a and GP8b are separated** — only the epithelium **needs** to wait for CNV (R2 requires malignancy to be corroborated by CNV);
  subclustering of the immune/stromal/endothelial lineages is unrelated to CNV and **must not** be blocked by CNV scheduling.
- **GP4c main caliber = Harmony (paper recipe)**; the uncorrected arm serves only as a **robustness control**; **GP9's SCMG is a separate branch**.
  The three are not controls at the same level, and **must not be tabulated together** in reports.
- **GP3 must come before GP7** — any **gate-specification change** for the single-dataset-without-ground-truth case must first receive written approval at GP3.

### M0 · Input Freeze Gate
- **Do**: per [`00_ingest/cohort_registry.py`](00_ingest/cohort_registry.py), include the **two paired datasets**
  (`GSE308103` snRNA + `GSE307534` spatial), producing a sample/patient/stage freeze table (including the **23 cases** of pairing relationships).
- **Gate conditions (all must be satisfied)**:
  1. The sample table has complete fields: `dataset, sample_id, patient_id, stage, modality`;
  2. **Stage agrees with the GEO ground truth** (strict token mapping, raise on unknown); **spatial LNM must not be fabricated** (currently missing, honestly annotated);
  3. The **23 paired patients** (P3–P25) all have slices in **both modalities**;
  4. Code review shows **no silent default**; produces a **freeze manifest + SHA-256 + validation report**.
- **Gate not passed → stop** (must not proceed to M1).

### M-1 · Plan Remediation Gate — **before M1, mandatory**

> **Background**: On 2026-09-12 a **stage-by-stage feasibility/rigor audit** of the whitepaper was performed (data identity + 8 methodological checks),
> finding **3 🔴 (method/data not valid) + 7 🟡**. **Without passing this gate, M1 and any later stage must not be entered.**

**A. Data Layer**
- **A1** The spatial atlas is limited to **Normal→AAH→AIS→MIA→IAC** (all GSE307534, same repository, same platform); **GSE190811 is removed** (verified via GEO as **breast cancer**, not LUAD; the original GSM5732148 does not exist in that repository). The rules document has been changed ([`docs/spatial_cohort_and_figure_prohibitions.md`](docs/spatial_cohort_and_figure_prohibitions.md)); the whitepaper is still to be changed.
- **A2** **LNM is out of scope** (no legal LUAD spatial data; the single-cell layer is moved out together with the three scRNA cohorts). To be added later if genuine LUAD LNM spatial data is obtained.
- **A3** Optional: `GSE305258` (ALK+ NSCLC lymph-node/brain metastasis spatial, 10 LNT, **GeoMx ROI, not Visium**) serves only as **orthogonal LNM validation** and is **not merged into the main Visium matrix**.
- **A4** **WES/TMB permanently disabled** (the TCGA mutation files are 0 bytes).

**B. Tools/Environment**
- Install **CopyKAT** (`MCMCpack`+`transport` from CRAN; pure R, no root) → M2 **main workhorse**.
- **infercnv downgraded to an optional cross-validation** (the official README declares it *no longer supported*; requires a system JAGS; Bioconductor is unreachable from this machine).
- **Upgrade scvi-tools to 1.5.x** (currently 0.15.5, a 2022 version) → only then can **SysVI / scArches(scANVI surgery)** be used.
  ⚠️ **This cannot be done in the shared environment** (Py3.8 + root-owned libraries; see §5b) → an **isolated env (micromamba, Py≥3.10)** must be created; otherwise fall back to the 0.15.5 capabilities and **honestly annotate** (no SysVI/scArches surgery).
- Install ~~**spaGCN**~~ (**2026-10-01 user ruling: M6 switches to BANKSY, spaGCN not installed and not used as a control arm**), **harmonypy**, **plip(py)** (weights already cached, 1.2 GB); ~~**CellCharter** (if adopted)~~ → **CellCharter moved to the isolated-env candidate (pending the S7 ruling at pre-registration, see §M6)**.
- **GROMACS missing** → install it, or **drop the 100 ns MD of Stage-8** (keeping only the docking).

**C. Methodological Corrections** (written into the M3/M5/M6/M7 gates below)
- **M3**: Drop "modality as batch" (sc↔sn is a *system* effect, and modality is collinear with dataset) → **scVI(`batch=sample_id`)**.
  **⚠️ 2026-09-15 correction**: In the original text, under the **single dataset** (GSE308103), the `batch=dataset` key is **constant and meaningless**; the 75 sample libraries are the only genuine technical batch axis.
  **Moreover, `patient_id` must never be used as batch** — stage is nested within patient, and correcting for it would erase the within-patient paired comparison needed by M4 gate criterion ④. scIB uses the **full panel**.
- **M5**: Drop self-normalization Σ=1 as a gate (full mode has no reject, so that gate is empty) → the gate moves to **RCTD native output + explicit QC**; the reference modality must match the slice.
- **M6**: GMM/BIC is replaced by **stability/consensus selection**; Fisher is replaced by reporting the **effect size**; Squidpy statistics are replaced by **category labels + empirical p**.
- **M7**: The CMap metric is corrected to **NCS**, or **cut entirely** (no output without LINCS data).
- **SCMG branch contraction**: does only **zero-shot cross-dataset scRNA integration + global manifold + state characterization**, compared against the standard branch via **scIB**.
  **Delete all statements of "state reversal / reversal factors / CausalGenePredictor causality"** (verified: SCMG has no such capability; the reference manifold has no tumor state).
  > 🔴 **2026-10-03 correction (user ruling: "**this was insufficient understanding before**")**: the **factual premise of the above item is wrong**, now checked —
  > ① `CausalGenePredictor` **really exists** (`tools/SCMG/scmg/model/causal_prediction.py`, MIT License,
  >    Copyright 2025 Xingjie Pan), and it **does output negative genes** (`perturbation_sign`); ② the pseudobulk perturbation library
  >    **is already on this machine** (`~/scmg_workspace/hf_data/pseudo_bulk_perturbation_database.h5ad`, 4.43 GB,
  >    20,345 perturbations × 18,108 genes, `obs` contains `perturbed_gene` / `perturbation_sign`).
  > ⇒ **the original text is not deleted** (Rule 3.2), but **that restriction no longer applies to the "reverse the signature against the perturbation library" usage**;
  > the **wording constraint "must not be called causal / must not be called a drug" still holds**, and **a null model is mandatory**.
  > For the implementation plan see `10_niche/TARGET_REVERSAL_PREREG.md` (**awaiting sign-off**).
  Target candidates are instead produced by the **standard branch** (SCENIC+/regulon + CellRank fate association + TCGA survival), with wording limited to **"candidate regulators/targets"**; they must not be called causal.

**D. Document Corrections**
- Whitepaper: delete voided results (Figure 1–5 "completed", 14 clones, 107,796 spots, `Z<-5.8`, 316,689, etc.); change LNM to pending; delete the SCMG reversal statements; fix the CellRank contradiction.
- PARAMETERS: fix the `infercnv denoise` default (**FALSE**, not TRUE), the Tau definition contradiction (**0–1** vs ±100); add references for the QC thresholds or mark them "project convention".

**Gate conditions (all must be satisfied)**:
1. Required tools **installed and smoke-test passed** (CopyKAT runs through on 1 sample; `scvi.external.SysVI` / scArches importable);
2. Whitepaper / PLAN / PARAMETERS revisions **committed to disk**, and a grep for the voided values (`14 clones`, `107,796`, `Z<-5.8`, `316,689`, `PT_3_LNM`/`GSM5732148`) shows no residue (except historical erratum records);
3. GSE190811 **fully purged** from code / documents / rules (except erratum notes);
4. The M3 / M5 / M6 / M7 gates **rewritten and committed**;
5. **No 🔴 remaining**.

### M1 · QC / Doublet Gate

> 🔴 **2026-09-16 caliber reversal (user decision "rebuild with the paper's QC")**. This section originally wrote "**strictly forbidden** to copy whole-cell scRNA thresholds
> (measured: `nCount≥1000` would cut ~30% of nuclei)". That conclusion **still holds numerically** (measured, it does cut ~30%),
> but **rejecting the paper's thresholds on that basis was wrong**: the paper's thresholds come from **this same batch of samples, the same platform (10x Fixed RNA/FFPE)**,
> and are the **native caliber**; this project's per-sample MAD was the sourced-from-nowhere self-chosen value. Hence **the main caliber is changed to the paper's fixed thresholds**,
> and the original MAD caliber is **downgraded to a sensitivity arm** (artifact retained: `results/02_expression/gse308103_counts.h5ad`, 648,945 nuclei).

- **Do**: explicit QC; **scDblFinder per sample** (iron law R2). Data is **GSE308103 (snRNA / nuclei)**.
- **Main caliber (paper, verbatim)**: drop `nFeature<500` or `nCount<1000` or `pct_mt>20%`; keep genes "detected in ≥3 cells".
- **Gate**: report cell counts before and after QC, doublet rate, parameters (with sources, see PARAMETERS); **no heuristic substitution**.
- **★ Measured (2026-09-16)**: new analysis set = **413,697 nuclei × 18,069 genes** (old M1 mask 648,945; the paper reports 401,635, **+12,062 = +3.00%**).
  The paper gate removes **235,248** nuclei within the M1 single-cell set: `nFeature<500` removes **68,694**, a further **166,554** fail only on `nCount<1000`,
  and **`pct_mt>20%` removes 0** (M1's adaptive MAD had already removed them) ⇒ the paper's mitochondrial gate is **effectively vacuous** on this data.
  Boundary: **0** nuclei have `pct_mt` exactly 20.000 ⇒ the two readings `<=20` and `<20` do not differ.
- **⚠️ Registered deviation**: the doublet caliber is **scDblFinder** (R2) vs the paper's **Scrublet**. R2 takes priority over "aligning with the paper"; a Scrublet sensitivity arm already exists.
- **⚠️ Must not** treat `pct_mt` as a fraction in order to "align with the paper" — this project's `pct_mt` is a **percentage** (0–100),
  see `01_qc/00_metrics_gse308103.R`: `pct <- 100 * colSums(m[mt,]) / pmax(cs,1)`.
- **AAH fragility**: the empirical test was **redone on the new cell set** (`results/01_qc/stage_fragility_report_paperqc.json`):
  median `nFeature` is **1.081×** that of Normal, AAH retention 65.48% vs Normal 62.05% (**+3.43 pp**) ⇒ **not fragile at the count level**.
  Therefore AAH **must not** get individually relaxed QC (that would be an unsourced silent default); fragility is instead handled by a falsifiable guardrail at the **GP5 clustering-absorption level** (metric 4).

### M2 · Fine Malignant Calling of Epithelial Subclusters (CNV) — **execution position: GP2, after GP8a epithelial subclustering**

> **⚠️ 2026-09-16 role narrowed twice**: ① from a "**leading corroboration gate before dimensionality reduction and clustering**" to "epithelial refinement after annotation";
> ② further narrowed to a **"fine discriminator that labels epithelial subclusters as malignant"**, with its execution position **moved back from after GP6 to after GP8a**.
> For the reason see §2.0: CNV's **sole purpose** in this project is to call **subclusters** malignant, and subclusters do not exist before GP8a.

- **Do**: **`CopyKAT` as the workhorse** (per sample, **pooling forbidden**), input = cells annotated by GP6 as the **epithelial lineage** and partitioned into subclusters by **GP8a**;
  **`infercnv` serves only as an optional cross-validation on 5–10k cells/sample** (officially unmaintained).
- **Output granularity = subcluster**: after obtaining CNV± per sample, aggregate by **GP8a's epithelial subclusters** (per-subcluster `frac_cnv_pos`, cross-sample consistency),
  and on that basis label **each epithelial subcluster** as `malignant` / `non-malignant` / `not_testable`.
- **Gate**: the malignant labels of epithelial subclusters must be **complete across all three states and traceable at sample level** (no case where "0 aneuploid across all samples" is taken as evidence of normality).
  **Without passing the gate, "malignant clone" must not be produced.** Wording ceiling = **subcluster level (aneuploid vs not)**; clone/subclone/phylogenetic-tree statements are **forbidden**.
- **Downstream use**: M4 criterion ⑤, and every "malignant epithelium" filter, **must all defer to GP2's subcluster labels**; a generic epithelium argmax must not be used as a substitute.

- **★★ 2026-09-16 measured finding: `CopyKAT` has two silent failure paths that must be blocked by guards**
  (source `tools/copykat/R/copykat.R`, copykat 1.2.5):
  1. **`norm.cell.names=""` ⇒ falls back to a self-guessed baseline** (`copykat.R:166-170` `baseline.norm.cl()`). At low confidence
     `WNS="unclassified.prediction"` and it **silently falls back** to `baseline.GMM(max.normal=5, mu.cut=0.05, Nfraq.cut=0.99)`,
     which **rewrites `WNS`** ⇒ the log retains only `"low confidence in classification"`.
     **Guard**: must pass `norm.cell.names=` (immune/stromal reference, `LUAD_NORM_REF=nonepi`),
     and assert that the log contains **`baseline is from known input`** (the falsifiable marker at `copykat.R:136-145`).
  2. **The `:456` all-diploid escape hatch**: `if (cor(conses.diploid, conses.aneuploid) >= 0.6) com.preN[] <- "diploid"`
     — **there is no "uncertain" category**. Pure-epithelial input most easily hits this hatch.
     **Measured**: `P11_LUAD` pure-epithelial input returns **0/160 = 0.0%** aneuploid ⇒ **that arm is blind**,
     and any "0 aneuploid" readout from it **carries no information** and must not be used as evidence of a normal sample.
     **Guard**: every tumor sample must **also** run an "all-cells" control arm; if the control arm is also 0, that sample's CNV conclusion is marked `not_testable`.
  - **Side-effect note**: `copykat.R:18-20` filters only when the **number of failing cells > 1** (with exactly 1 it does not filter) —
    any pipeline replicating copykat's cell filtering must copy this semantics verbatim, otherwise the `n_judged` guard will misfire.

- **Note**: this data is **snRNA (GSE308103)** → CNV requires relaxed parameters and careful interpretation (nuclei, low UMI);
  the reference uses immune/stromal cells from the same sample (i.e. `norm.cell.names` from guard 1 above).
- **Pre-registered handling of degenerate samples** (project convention, no literature threshold): `n_retained<200` → `not_testable`, **never** filled in as "normal";
  `n_cnv_pos==0` → reported as a **genuine result**, and `KS.cut` **must not** be relaxed and rerun (Rule 3.2 forbids post-hoc parameter tuning);
  `frac_cnv_pos>0.95` in Normal/AAH or `<0.01` in IAC → `implausible`.
  **Abort criterion**: `not_testable` + `implausible` totaling > 30% → **stop and report** (R2 is unsatisfiable for that dataset),
  and a generic-epithelium argmax substitute is **disabled**.
- **Subcluster-level propagation rule (added 2026-09-16, fixed before running)**: a sample marked `not_testable` ⇒ **all epithelial subcluster cells** contributed by that sample
  are also marked `not_testable` (**must not** be given a label just because "other samples are normal"); if a subcluster is `not_testable` across **all** samples
  ⇒ that subcluster's malignant status = **`unknown`**, and it **must not** enter any "malignant vs non-malignant" comparison. The report must provide
  a `subcluster × sample` `not_testable` coverage table. **"unknown subcluster" proportion > 30% ⇒ stop and report at GP2.**
- **Must report**: the list of genes silently discarded because `annotateGenes.hg20()`'s internal position table does not cover them (Rule 0's silent-default trap).

### M3-A · Traditional/Mainstream Branch (execution position: GP4a–GP8a → GP2 → GP8b)

> **⚠️ 2026-09-16 split note**: the original text crammed two things at different levels into a single M3 gate; they are now split:
> **the two arms within M3-A** (GP4c) = **Harmony on `sample_id` (main caliber, paper recipe)** vs **no correction (sensitivity arm, not the paper recipe)**;
> **M3-B's SCMG** (GP9) = **a separate branch**. **The three are not controls at the same level and must not be tabulated together in reports.**
>
> 🔴 **2026-09-16 main/control reversal**: originally "no correction = **main**, Harmony = control". The paper writes Harmony into its **main recipe**,
> so the main caliber is changed to Harmony; the uncorrected arm is **retained** but downgraded to a robustness control (the `sample_id`↔`stage` collinearity risk remains).

- **Do** (in the order of §2.0): GP4a preprocessing → GP4b dimensionality reduction → GP4c two batch arms → GP5 resolution → GP6 annotation
  → **GP8a epithelial subclustering** → **GP2 CNV defines subcluster malignancy** (M2 gate, see above) → **GP8b subclustering of the remaining 5 lineages**.
- **The "authoritative parameters and order" for the traditional branch are governed by [`docs/PARAMETERS_AND_SOURCES.md`](docs/PARAMETERS_AND_SOURCES.md) M3-A.0–A.5**,
  which is a pre-execution registration (Rule 3.1); during execution it **must not be deviated from without registration**.
- **Gate** (GP4a/b/c): `SCTransform(vst.flavor="v2")` / HVG `3000` / `npcs=50` / `Harmony dims=1:50` /
  `k.param=20` all **registered with hashes**; **parameter sources traceable item by item** (verbatim paper quote or Seurat source default).
  ⚠️ **The original `{20,30,50}` ARI guardrail is void along with the old PCA-30 caliber** and does not apply to the paper's fixed 50 PCs.
- **Two-arm reporting rule**: per-cell `ARI(Harmony, no correction)` < 0.7 ⇒ **must not report only one**, both arms must be presented together and the matter escalated for a ruling.
- **Execution script**: [`04_integration/10_seurat_traditional.R`](04_integration/10_seurat_traditional.R)
  (`smoke --n_cells N` / `full`, two modes; already passed 3,000- and 60,000-cell smoke tests).
- **Disabled**: `modality as batch` (sc↔sn is a system effect and collinear with dataset); **`patient_id` as batch**
  (stage is nested within patient, and correcting for it would erase the within-patient pairing needed by M4 criterion ④).

### M3-B · SCMG Control Branch (execution position: GP9, after M3-A)

- **Do**: **zero-shot cross-dataset scRNA integration + global manifold + state characterization**
  (**no mixing of traditional algorithms** — iron law 5; **no reversal/causal output** — iron law 6).
- **Gate**: compare against M3-A using the **full scIB caliber** — batch removal (kBET + iLISI + graph-connectivity + PCR)
  **and** biological conservation (cLISI/ARI/NMI/ASW) **reported together**; the batch panel and the bio panel are **reported separately**, with raw per-metric values attached.
- **Criterion rewritten (2026-09-15)**: no longer "maximize iLISI" but rather **batch metrics acceptable _and_ stage separability preserved**
  (demonstrating **resolving power**, not everything smeared together).
- **Gate-specification change**: the single dataset has **no ground truth** ⇒ the biological-conservation ARI/NMI can only be computed **between uncorrected vs corrected cluster labels**,
  and may be computed **only after explicit approval at GP3**.
- **Ruling**: if Arm B's **ARI < 0.7** relative to Arm A ⇒ Arm A is primary, with the methodological limitation reported honestly.
- **Disabled**: judging integration success on iLISI↑ alone (it can be inflated by over-integration).

### M4 · Cross-Modal AAH Gate (**all five criteria required**)
- ① Each stage **distinguishable**; ② overlapping stages sn↔sc **consistent**; ③ platform offset **δ(stage) stable**; ④ **within-sn pairing** (AAH vs the same patient's Normal/AIS) concordant; ⑤ AAH identity **CNV/marker**-corroborated.
- **Gate not passed → do not draw an "AAH conclusion"** (it may only be labelled a hypothesis).
- Data: `GSE308103` (snRNA, main reference); if `HRA001130` (whole-cell scRNA, controlled access) is obtained later it may serve as cross-modal AAH validation.

### M5 · Spatial Deconvolution Gate
- **Do**: RCTD (`doublet_mode='full'`, recommended for Visium), with the reference using the signatures frozen in M3.
- **Gate**: governed by the **RCTD native output** — weights non-negative, **per-spot weight distribution reasonable**, **reference modality matching the slice** (FFPE↔sn, cross-modal DEGs removed first), reporting the **number of low-quality/removed spots** (defined by the explicit `UMI_min`/`counts_MIN` thresholds, **not** using a self-imposed Σ=1 as the gate); signatures and reference **hash-aligned**.
- **Disabled**: treating Σ=1 after `normalize_weights()` as a property of RCTD (that is imposed by the project itself); `full` mode has **no** reject category.

### M6 · Niche Gate
- 🔴 **2026-10-01 tool correction**: this line originally wrote **SpaGCN**, conflicting with the **BANKSY** registered in `docs/PARAMETERS_AND_SOURCES.md` §M6;
  **the user ruled on 2026-10-01 to switch to BANKSY** (SpaGCN is not used as a control arm). **"CellCharter" cannot be installed in the shared Python (3.8.10)** (requires ≥3.10),
  but the **micromamba isolated environment in the plan board §5b is a viable solution** (**not "cannot be installed on this machine"**) ⇒ the cross-slice step **retains CellCharter as a candidate sensitivity arm**,
  pending the **S7** ruling in pre-registration §8 (plan A per-slice→consensus / plan B joint clustering). The number of clusters is always determined by **Leiden multi-resolution + cross-seed stability**.
  For pre-registration see [`10_niche/NICHE_PREREG.md`](10_niche/NICHE_PREREG.md) (🟢 **signed 2026-10-01**: effective values for all of S1–S13 + §12 limitations; 🔴 **2026-10-02 appended the §4.1 Gate 1 FAIL ruling + §12-⑪ ⇒ limitation count nine→ten→eleven**).
- **Do**: BANKSY per-slice domain partitioning → neighborhood aggregation → **data-driven clustering (stability/consensus selection, not a single BIC minimum)** → Squidpy spatial statistics.
- **Gate**: **no hard-coded cluster count/labels**; the cluster count is determined by **stability/consensus** and archived; spatial statistics use **category labels** (not deconvolution proportions), reporting **z-score + empirical p ((b+1)/(n+1))**, with the permutation `n_perms` recorded; enrichment tests report the **effect size** (Fisher FDR is inflated by compositionality and large N).
- **Disabled**: writing `nhood_enrichment`'s output as `P<0.001` (that function **does not return p values**); claiming a permutation was done for `co_occurrence` (it **does not**); a distance step smaller than the Visium ~100 µm spot spacing.

### M7 · Target Gate (traditional track)
**M7a · Niche prognostic signature**: BANKSY niches (**corrected from CellCharter on 2026-10-01**, same reason as §M6) → **patient/slice-level aggregation** into a signature → TCGA-LUAD bulk scoring (ssGSEA/deconvolution) → **multivariable Cox (adjusted for stage/age/sex) + KM**.
- **Gate**: the signature is produced by **patient-level** aggregation (not spot-level DE); overfitting prevention (penalty/cross-validation); wording = **prognostic association** (`niche poor-prognosis signature`). Spatial has only 11–25 cases → **no claim of patient subtyping**.

**M7b · Candidate target pool (genetic-statistical anchoring)**: candidate sources = malignant programs/regulons **＋** niche signatures → **cis-MR + coloc** (ILCCO/TRICL LUAD GWAS × lung eQTL) → **Open Targets tractability + DepMap selective dependency + clinical-stage drug matching**.
- **Operations manual**: [`docs/M7B_MR_COLOC_TARGET_ANCHORING.md`](docs/M7B_MR_COLOC_TARGET_ANCHORING.md)
- **Gate**: each candidate is given `genetic_support ∈ {supported, not_supported, not_testable}` (thresholds in the manual, **no post-hoc parameter tuning**); `not_testable` **honestly marked as missing, not substituted**; wording = **genetically supported candidate target**, and it **must not be called causal**.

**M7c · CMap (optional)**: executed only when **real LINCS data** is obtained and the correct metric is used (**NCS**, not negative Tau); otherwise **refuse to produce a table** (currently no LINCS data → not produced by default).
  > 🟡 **2026-10-03 progress**: **real LINCS is downloading** (GSE70138 Phase II Level5 5.00 GB ＋ GSE92742 Phase I Level5 19.86 GB
  > ＋ metadata, ~25 GB total → `~/lincs_data/`). Once it lands, this section's **gate opens**.
  > Also: the SCMG perturbation library is **unblocked as a second perturbation library** (see above), used for **cross-validation**, and its metric caliber must be fixed separately in writing.

- **Overall gate**: DESeq2 design/truncation has a source; target wording constrained; MR inputs/outputs have hashes and provenance records.

### M8 · Docking Gate
- **Do**: fpocket + AutoDock Vina (targets from M7's real results); **MD is done only after GROMACS is installed**.
- **Gate**: target/structure files **genuinely exist**; parameters have sources; **nothing not computed is fabricated**.

---

## 3. General Checkpoints (checked at every step)
- [ ] No `np.random` / hard-coded values / fake curves;
- [ ] No silent default (falling back to some stage/label);
- [ ] Scripts reproducible (deterministic + explicit seeds);
- [ ] Outputs hashed; inputs and outputs traceable;
- [ ] Parameters have a source in `docs/PARAMETERS_AND_SOURCES.md`.

---

## 4. Dual-Branch Architecture (spanning M3–M8)
- **Standard/mainstream branch**: scVI (`batch=sample_id`)/scANVI + scArches, Harmony, Seurat, RCTD, SpaGCN, Squidpy, DESeq2, CellRank, SCENIC+…
- **Pure SCMG branch**: **zero-shot cross-dataset scRNA integration → global manifold → cell-state characterization** (**no mixing of traditional algorithms**).
  ~~conditional diffusion trajectory → CausalGenePredictor causality~~ — **deleted** (capability does not exist, see iron law 6).
  > 🔴 **2026-10-03 correction**: the "capability does not exist" this deletion was based on **was verified to be erroneous**
  > (`CausalGenePredictor` really exists, see the 2026-10-03 correction block above).
  > **The conditional diffusion trajectory is still not done** (that is a limitation of the manifold route, which genuinely lacks a tumor state);
  > but **"reverse-ranking the signature against the perturbation library" is a different usage and has been unblocked**, see `10_niche/TARGET_REVERSAL_PREREG.md`.
- The two branches are **compared** (full scIB caliber); tool source code is in [`tools/`](tools/).
- Target/regulator candidates are **produced only by the standard branch**, and worded as "candidates".

> **⚠️ 2026-09-16 level clarification**: in this document "arm" and "branch" are **two different levels** and must not be used interchangeably:
> - **Branch** = the top-level dichotomy (**standard/mainstream branch** vs **pure SCMG branch**); the two branches **each run the full pipeline independently** and are then **compared against each other**.
> - **Arm** = the batch-handling comparison **within the standard branch** (**Arm A no correction** vs **Arm B Harmony on `sample_id`**),
>   sharing the same downstream outside preprocessing, diverging only at GP4c.
> Therefore **the Arm A/B comparison conclusion belongs only to M3-A** and **must not** be merged with M3-B's SCMG comparison into a single table.

---

## 5. Data Sources and Acquisition (**scope already narrowed: only two paired datasets**)

> **2026-09-15 correction**: paired-patient count **9 → 23** (P3–P25). **The old value 9 was an artifact of "only 19/56 spatial slices downloaded"**;
> once the slices were fully downloaded, intersecting the two GEO-authoritative tables actually gives **23 cases**. See `PAIRED_PATIENTS_MIN` in `00_ingest/cohort_registry.py`.
> **2026-09-12 narrowing**: use only **`GSE308103`(snRNA) + `GSE307534`(Visium spatial)** — the same study, **modality-matched (FFPE↔FFPE)**.
> The three scRNA cohorts (GSE131907/189357/148071) are **moved out of scope**, archived at `/home/eto/luad_invasion/luad_v2_out_of_scope/`.

- **Single-cell/reference**: `GSE308103` (snRNA, 75 samples / 798,100 verified measured) — the **only AAH-containing** single-cell resource;
- **Spatial**: `GSE307534` (Visium CytAssist FFPE; GEO 56 samples / 25 patients, **all 56 slices present locally**, covering all **23 cases** of paired patients; only P1/P2 have no snRNA and are not paired);
- **LNM**: spatial currently missing (~~GSE190811~~ verified as **breast cancer**, voided);
- (out of scope, archived) `HRA001130` (sc, controlled access) as a reserved interface → `/home/eto/luad_invasion/luad_v2_out_of_scope/`;
- TCGA-LUAD: expression + clinical only (**mutation files 0 bytes, WES/TMB permanently disabled**).

---

## 5b. Environment Constraints and Minimal Decisions (Environment Constraints, verified 2026-09-12)

> These are **hard constraints** of this machine's environment; they determine the tool selection and "which methods simply cannot be installed". **Any plan conflicting with them is void.**

1. **No GPU** (no `/dev/nvidia*`, `torch.cuda.is_available()=False`) → GPU-dependent methods (Boltz-2, cell2location, SysVI…) are **infeasible CPU-only or at ~420,000 cells**.
   ⚠️ **2026-10-01 correction**: this originally listed **CellCharter** among GPU-dependent tools, which is **inaccurate** — per its repository docs, **GPU is not required**; CellCharter's real blocker is **Python ≥3.10** (see item 2), not GPU.
2. **Python 3.8.10 only, no conda** → the modern DL stack cannot be installed:
   - `cellcharter` requires ≥3.9; `scvi-tools` 1.5.x (SysVI / scArches surgery) requires ≥3.10 → **not installable in the shared environment**.
   - **Solution**: when modern DL is needed, use an **isolated environment** (micromamba from conda-forge); PyPI only via mirror (`https://pypi.tuna.tsinghua.edu.cn/simple`), conda-forge via direct connection.
3. **Shared libraries are root-owned** (`/usr/local/lib/R/site-library`, `/usr/local/lib/python3.8/dist-packages`) → **piping community DL packages into the shared environment is forbidden**:
   - Installing `cellcharter` once **silently downgraded torch to 1.12.1**, breaking `scvi`/`pytorch-lightning` (manually rolled back torch→2.4.1+cu118, pytorch-lightning→1.5.10.post0, torchmetrics→0.7.3).
   - `scvi-tools 0.15.5` is **kept** (pin `pytorch-lightning>=1.5,<1.6`); packages always go into **personal libraries** (`~/.local/lib/python3.8/site-packages`, `~/R/.../4.2`).
4. **R installable**: **CopyKAT** (first install `RcppEigen`(Eigen 4.0) into the personal library, then `transport`), `coloc`, `ieugwasr`, `TwoSampleMR`.

**Direct impact on this plan**:
- "Upgrade scvi-tools to 1.5.x" in M-1 §B **cannot be done in the shared environment** → an **isolated env (micromamba)** must be created, or fall back to **scvi 0.15.5's available capabilities** (no SysVI/scArches surgery) with honest annotation.
- M2's workhorse **CopyKAT** and M6's **CellCharter**: CopyKAT is installable (see above); **CellCharter needs an isolated env**.
- Any "GPU-accelerated" statement must be deleted or downgraded to CPU.

---

## 6. Milestone Progress Board

> **2026-09-16 restructuring**: this table was originally ordered by M number, inconsistent with the actual execution order in §2.0 (M2 was placed before M3).
> It is now **reordered by execution order**, and the previously entirely missing GP1–GP9 rows have been added.

| No. | Checkpoint | Module | Status | Gate |
| :---: | :--- | :--- | :---: | :---: |
| — | M0 Input freeze | Input | ⚠️ Gate not passed | ☐ |
| — | M-1 Plan remediation | Plan | 🔶 In progress | ☐ |
| 1 | M1 QC/doublets (🔶 downgraded to sensitivity arm) | QC | ✅ 100% | ✅ |
| 2 | **GP0** Expression object rebuild | Expression | 🔶 **redone 2026-09-16 with the paper QC**; 09-17 audit supplemented validation **13/0** | ☐ |
| 3 | **GP4a** Preprocessing (SCTransform v2 → HVG 3000) | M3-A | ✅ **full run completed 2026-09-17** | ✅ |
| 4 | **GP4b** Dimensionality reduction (PCA 50, no ARI guardrail) | M3-A | ✅ **full run completed 2026-09-17** | ✅ |
| 5 | **GP4c** Two batch arms (**Harmony main** / no-correction control) | M3-A | ✅ **2026-09-17 ruling: adopt the paper's Harmony main caliber** | ✅ |
| 6 | **GP5** Resolution selection (range 0.5–0.8 × 5 seeds + AAH guardrail) | M3-A | ✅ **2026-09-17 sign-off L1 `r*=0.6`**; metric 4 acknowledged ineffective; **2026-09-21 checkpoint closed**. 🔴 **2026-09-22 rerun**: six lineages' `r*` currently **0.7/0.8/0.5/0.7/0.6/0.8** (epithelial, fibroblast carry `relaxed`), six JSONs signed by the user personally; metric 3 after recomputation **only myeloid fails the line** (kept, registered as a defect); seeds **T·NK→3, fibroblast→4** | ✅ |
| 7 | **GP6** Dual-standard annotation + κ (**identifies only "which cells are epithelial"**) | M3-A | 🔶 **run, not signed off**: see `results/05_annotation/GP6_report.md`, artifact `gp6_cell_labels.csv.gz` (413,697 nuclei) | ☐ |
| 8 | **GP8a** **Epithelial** subclustering → epithelial subclusters | M3-A | ✅ **Completed (2026-09-22 night run), report named `GP8c` (see the numbering-fork note below)**. Epithelial `r*=0.7`/seed 0 → **27 subclusters / 133,384 nuclei / 0 clusters removed**; subtypes AT2 88,666 · AT1 30,320 · Ciliated 6,697 · Goblet/Mucous 3,587 · Basal 2,386 · Serous 1,728 (Club/Ionocyte/Neuroendocrine/Tuft did not win out, honestly reported rather than padded to fill). ⚠️ **The gate column is still empty = not formally signed off** | ☐ |
| 9 | ~~**GP2** CNV **subcluster** malignant fine calling~~ | **M2** | ⛔ **voided 2026-09-24** (single-cell CNV failed and exited) → switch to **spatial SC0–SC4** (see §M2 and `08_spatial_deconv/SPATIAL_CNV_PREREG.md`). ⚠️ The original "wait for GP8a" dependency **is void** (GP8a was in fact long completed) | ⛔ |
| 10 | **GP8b** Subclustering of the remaining 5 lineages (not blocked by CNV) | M3-A | ✅ **Completed (2026-09-22 night run, same batch as GP8a, see `GP8c_report.md`)**: T·NK `r*=0.8`/seed 3 (38 subclusters, 63,078 nuclei, 2 clusters removed), B·plasma `r*=0.6`/seed 0 (22 subclusters, 22,915 nuclei, 9 clusters removed), myeloid `r*=0.5`/seed 0 (27 subclusters, 63,869 nuclei, 1 cluster removed, **S4 panel 12 types won out**), endothelial `r*=0.7`/seed 0 (**S4 panel 8 types won out**), fibroblast `r*=0.8`/seed 4 (33 subclusters). Six lineages **mutually exclusive and complete = 413,697 nuclei**. ⚠️ Same as GP8a: **not signed off** | ☐ |
| 11 | **GP3** Metric backend freeze → **GP7** scIB dual panel | M3-A | ⬜ 0% | ☐ |
| 12 | **GP9** SCMG control arm | **M3-B** | ⬜ 0% | ☐ |
| 13 | M4 Cross-modal AAH | — | ⬜ 0% | ☐ |
| 14 | M5 Spatial deconvolution | — | ⬜ 0% | ☐ |
| 15 | M6 Niche | — | 🟢 **Pre-registration signed (2026-10-01)**; tool = R `Banksy` 0.1.6 (already installed, zero install). **2026-10-01 running**: HVG-3000 frozen (§3.1.1 = per-slice vst + Seurat `SelectIntegrationFeatures`, measured **42.4 minutes**) → smoke test **2/2 passed** (~1 minute each) → **full grid 112/112** (56 slices × two AGF versions; concurrency ruled by the user on 2026-10-01 **3 → 8** — a single task measured **~2 GB / 1 core**, concurrency is only a resource knob and **changes no result**) → downstream **03→04→06→05 relay**. First-slice measurement: the smallest P25 (3,776 spots) ran all 40 cells in **10.83 minutes** (≈14.9 s/cell) ⇒ 112 runs ≈ 58 h CPU, ≈ **8 hours** at -P 8. 🔴 **The first run died across the board, two root causes, both fixed** (remounted after the fix at 2026-10-01 19:36): ① `run_banksy_grid.sh` used `export -f` + `xargs bash -c`, but in this machine's PATH `/home/eto/.local/bin/bash` is a `#!/bin/sh` forwarding shell, and **dash drops `BASH_FUNC_*`** ⇒ functions do not propagate (see memory `env_bash_name_resolves_to_dash_wrapper`) ⇒ switched to **self-dispatch** `xargs ... /bin/bash "$SELF" __run_one`; ② the same script's memory guard was written as `if free_gb -lt "$MIN_FREE_GB"`, feeding `-lt 80` to the function as an **argument**, which is **always true** ⇒ all 112 tasks would be misjudged as "insufficient memory" and skipped (measured: 252 GB available yet it still printed "252 GB < 80") ⇒ changed to fetch the value first then compare, `fg="$(free_gb)"; (( fg < MIN_FREE_GB ))`. 🟢 **2026-10-01 20:00 added an unattended supervisor** (user order: checkpoint + memory adaptivity + no loss of intermediate results + session-detached): **① `10_niche/07_supervisor.sh`** (resident with PPID=1, the sole parent process, launching the main chain and downstream itself) — every 30 seconds writes memory/progress into `results/10_niche/banksy/memory_timeline.tsv`; if the chain or downstream **exits unexpectedly it automatically remounts** (up to 10 times each, with backoff); once downstream finishes it writes `results/10_niche/_DOWNSTREAM_DONE` and the supervisor exits on its own. **② Driver given checkpoint-resume** — on startup it reads the completion table and skips (slice, AGF) **unique pairs** already DONE (reruns leave multiple rows, so progress is always counted by unique pair). **③ The memory gate changed from "skip if insufficient" to "wait"** — no more SKIP_MEM quietly dropping work; when the supervisor is tight on memory it sets a pause flag `banksy/_PAUSE`, which the driver also respects; if the supervisor is absent the flag auto-expires and is cleared, so it will not deadlock permanently; only after waiting over 1 hour is `SKIP_MEM_TIMEOUT` recorded. **④ Memory adaptivity** — available < **60 GB** sets the pause flag (no new tasks started, running ones untouched); < **40 GB** kills the **most recently started** task (least invested, a rerun will recompute it). **⑤ `02_banksy_grid.R` writes to disk atomically** (write `.tmp` first, then rename) — being killed midway leaves only a half `.tmp`, and the official file is either the complete old version or the complete new version, never mistaken by downstream for a complete result. ⚠️ The smoke test was switched to an **independent completion table** `banksy/_status_smoke.tsv` (cleared before starting) — if its 2 rows landed in the official table, a resume would misjudge both P25 versions as "already complete" and skip them. P25's 40 official cells (the -P3 round, finished 19:48, 10.83 minutes, gz intact, manifest 40 rows, 5 seeds) are **retained** and back-filled into the completion table, saving a rerun. 🔴 **§5.1 depth-balancing guard (672 reruns): the script is not written and must be attached after 03 produces the configuration**. 🟢 **2026-10-02 09:23:52 full grid finished: 112/112 unique pairs all DONE, zero FAIL** (started 20:09, wall clock **13 hours 14 minutes**; measured **memory was in surplus throughout, 215–252 GB**, the 40 GB hard floor and 60 GB dispatch line **never triggered**). 🔴 **Gate result: Gate 1 (cross-seed ARI ≥ 0.90) measured "0 surviving cells / 40 total cells", identical for both `use_agf` versions** ⇒ **the M6-3 criterion is unmet**; `03_consensus.R` printed "under rule=frac_ge there is no surviving cell in the high-λ range; falling back to all cells and selecting again" and still selected parameters (agfT `k_geom=18/λ=1.0/r=1.0`, agfF `18/1.0/0.5`), and **this fallback behavior was not written in the pre-registration — it is an implementation-precedes-sign-off deviation**. The selected cell's cross-slice median ARI is only **0.721 / 0.711**, pass rate **0.02 / 0.04**. 🟢 **User ruling "A" on 2026-10-02 (already written into `NICHE_PREREG.md` §4.1＋§12-⑪＋the signature slot)**: **keep Gate 1 as is, report FAIL honestly, do not relax the threshold and do not re-sign after the fact (Rule 3.2)**; this set of **16 domains (K\*=16) is downgraded to an "exploratory / sensitivity arm", and the main text must not write it as a "pre-registration-passing main result"**; the report must show "Gate 1 FAIL" on the same page as the result, and **must not** replace or mask it with the §3.5 K\* stability (0.6548/0.6526) (that is cross-slice profile merging, not per-slice cross-seed reproducibility). ⇒ **This arm currently has no main result**; to produce one, a **separate pre-registration** is required. Corroborating evidence (already on disk, must be reported alongside): §3.5 consensus **K\*=16**; §3.5 RCTD composition-alignment sensitivity **ARI = 0.225 / 0.247 ⇒ triggers the preset branch "niche definition is sensitive to method"**. 🔴 **The downstream "silent death" has been identified as a misjudgment**: the downstream relay was not killed externally; rather in `09_postrun_chain.sh` **`06_coverage_guard.R` reported rc=1 ⇒ `exit 3` stopped it in place** (that `!! … 失败` line was in the log all along; my earlier grep missed it) ⇒ each round 03(≈5 min)→04(≈4 min)→**06 crashed in seconds** ⇒ about 9 minutes per round, which was misread as "silently vanishing in 7–8 minutes". **Root cause (my own bug)**: `median()` for **integer** vectors **returns a different type depending on element-count parity** (odd⇒integer, even⇒double) ⇒ when `data.table` groups by `slide`, the slice types differ, directly erroring with `Column 3 of result for group 2 is type 'integer' but expecting type 'double'` ⇒ added `as.numeric()`, which **only unifies the storage type, values unchanged**. Both 03/04 product versions **landed in full** across the repeated reruns; **06 has since run through separately after the fix**; `05_spatial_stats.py` has not yet run. ⚠️ The supervisor side added a `_run_postrun_traced.sh` post-mortem shell (recording rc/signal/lifetime seconds into `logs/_postrun_exit.log`), **retained but without changing the caliber** — it is only insurance, since the real culprit was identified separately. 🔴 **2026-10-02 the entire arm re-signed under the platform-compliant cell (`NICHE_PREREG.md` §13 S4-REV/S6-REV, issued by the user's "commence work")**: it was discovered that **the λ semantics were applied to the wrong platform tier** — our data is **10x Visium v1/v2, 55 µm spots**, and BANKSY's official docs recommend **λ=0.2** for its **domain partitioning** (the general rule "low λ for typing / high λ for domain partitioning" **holds only for high-resolution technologies**); when S4 was signed, 0.2 was deleted from the grid and the main candidate took the high end, and after Gate 1 FAILed across all cells it fell — via an **un-pre-registered fallback rule** — to **λ=1.0** (= zeroing this spot's own expression weight) ⇒ the platform-compliant cell **was never run from start to finish**. New grid **λ {0, 0.2, 0.5, 0.8} (1.0 withdrawn) / res {0.5, 0.6, 0.8} / k_geom {6, 18}**; **λ=0.2 and k_geom=18 are fixed as the main cell, no longer selected by the grid** (§13.5); ranking switched to the **C1 cross-seed median ARI**; **PAC reports the value only, no hard line** (`PAC < 0.1` could not be re-verified from the primary literature ⇒ no line set); **C3 spatial coherence is the sole hard gate**, per-slice `coh_mean > null-distribution q95` (one-sided α=0.05 permutation), cross-slice **user ruling "≥90% of slices pass"** (§13.6; on FAIL report the error honestly, **no silent fallback**). The old λ=1.0 products are archived by **rename, not delete** (`banksy_STALE_lam1.0_20261002` / `consensus_STALE_lam1.0_20261002` / `spatial_stats_STALE_lam1.0_20261002` / `depth_guard_STALE_lam1.0_20261002`). **2026-10-02 11:00 full-pipeline session-detached remount** (`07_supervisor.sh` PPID=1), order **02 → 03 → 04 → 06 → 05 → §5.1 depth guard** (the depth guard has been attached to the end of `09_postrun_chain.sh`, since `08_depth_guard.R` must read the configuration produced by 03). 🔴 **While mounting I poked a hole myself, since fixed**: while clearing old products, `rm -rf $NICHE/banksy` deleted the **in-flight** 8 worker output directories plus the `_status.tsv` header — the driver writes the header only when the file is missing, so a subsequent `>>` creates a **headerless** file ⇒ `awk 'NR>1'` undercounts the first DONE ⇒ **it can never reach 112**; when discovered the table **did not yet exist** (⇒ zero complete, **no data loss**), and I **restored the 56 slice directories + rebuilt the header**, with no need to restart running tasks. Separately a **chain with PPID=1 started at 10:52:30** (not accompanied by a supervisor) was found, and after checking it was **running the new grid** (real-time log λ 0/0.2/0.5/0.8, res 0.5/0.6/0.8), so it was kept running; a `ps` count confirms **there is only one chain and one supervisor on the whole machine** (no risk of a duplicate launch). 🟢 **2026-10-02 12:52 healthy while running**: **37/112 DONE, zero FAIL**; 8 workers; memory **222 GB available** (the dispatch line 60 / hard floor 40 never triggered, supervisor heartbeat `paused=0`); measured rate **0.31 pairs/minute** ⇒ **ETA about 4 hours (finish about 17:00)**. **Parameter check passed**: λ{0,0.2,0.5,0.8}/res{0.5,0.6,0.8}/k_geom{6,18}, npcs=20, Leiden k.neighbors=50, CP10K+log1p, seeds 0:4, N_PERM=200, HVG=3000, the C3 gate, `K_GRID`, and 05's `n_perms=1000/seed=20261001/n_neighs=6/DIST_UM` all match the signed caliber (the only item still awaiting a ruling remains the K-stability threshold of §3.4 supplement-3). 🔴 **Two new registrations (2026-10-02; `NICHE_PREREG.md` §12-⑫/⑬ ＋ the §3.2.2 cross-note ＋ the §13.6 addendum; both "register only, change no number and no gate")**: ⑫ **The C3 hard gate has weak discriminating power on this cohort** — real-time `coh_mean` 0.66–0.90 vs null-distribution q95 0.08–0.19 (z≈140–370), and **even the `λ=0.0` cell (using no spatial information at all) passes** ⇒ C3 can only rule out "clearly anti-spatial" and **cannot prove that the partition used spatial information**; on this cohort it **plays no filtering role** (the actual cell choice rests on C1); **the gate is unchanged**, but the report **must not** write "C3 passed" as "spatial structure established". ⑬ **The scale caliber is self-contradictory within this file** — §3.2.2 records the measured spot spacing **91.1 µm** (395.3 px × 0.2304), §13.0 records **99.03 µm** (which at the same µm/px should correspond to 429.8 px), **the two numbers are mutually exclusive**; the `05_spatial_stats.py` comment carries 91 over; `UM_PER_PX` is **the sole coefficient for the distance-shell conversion** ⇒ if the true value is 99.03, the shell's **actual physical scale differs from the label by about 9%** (the 05 comment's claim "only affects reported numbers, not the figure" **awaits verification**); **no number has been changed**, and the true value must be ruled on and unified across the three places. 🔴 **2026-10-02 13:25 provenance addendum (the user asked "does this gating parameter come from the paper")**: checked item by item, **C3 is not a gate recommended by any source paper — it is project-built** — ① **the source paper (Peng 2026 *Cancer Cell*, the origin of this data) has no spatial domain-partitioning gate**; its use of Visium is to **transfer** snRNA labels back, its own convention is "delete mixed clusters directly + argmax", and it **sets no gate and leaves no "undefined"**; ② **nor does BANKSY** — this machine's `Banksy` 0.1.6 `NAMESPACE` exports only `getARI`/`plotARI`/`ConnectClusters`/`SmoothLabels` and the like, and **no "spatial coherence" function exists**; BANKSY's original evaluation of domain partitioning uses **ARI against ground truth** (only ≈0.35 even on DLPFC) ⇒ it **neither provides nor recommends** any coherence threshold; ③ item-by-item provenance: the **motivation** = SpatialARI 2025 ("ordinary ARI completely ignores spatial location"), the **`6 neighbors`** = a **geometric constant of one ring** of the Visium hexagonal grid, not a recommended value, the **`label-permutation null distribution + q95`** = a generic permutation-test construction, **project-built** (computed by hand in `03_consensus.R`, not a package function), the **`≥90% of slices pass`** threshold's **0.90 borrowed from this project's GP5** (§12-③ already registered as "borrowed, uncalibrated at spot scale"), and **the very act of "using a proportion across slices" was the user's call on 2026-10-02**. ⇒ **The report must not write C3 as a "criterion with literature support"**; it is a **project-built guard** (its "weak discriminating power" is the previous item ⑫; the two are the same source). ⚠️ **Do not conflate**: what genuinely has "source-paper recommendation" are the **parameters** (λ=0.2 / `k_geom`=18 / `npcs`=20 / res≈0.55, official verbatim in §13.0), **not this gate**; and those parameters come from the **BANKSY paper**, **not Peng 2026** (Peng 2026 gave only res 0.5–0.8 for its snRNA clustering). **Already written into `NICHE_PREREG.md` §13.2 table's C3-row "external basis" cell ＋ the new §13.6 provenance-addendum block; the table pipe count was re-checked (6 per row, not broken).** 🟢 **2026-10-02 13:22 healthy while running**: **43/112 DONE, zero FAIL**; 8 workers in flight (P6_AAH / P10_MIA / P14_LUAD / P12_AIS, each two AGF versions); memory **229 GB available** (`free` reading 170 GB free / 36 GB used, dispatch line 60 / hard floor 40 not triggered); measured per-job duration **31–38 minutes** (the recent 6 official jobs average ≈34 minutes; exclude the 8 smoke ones of 1–15 minutes) ⇒ **ETA about 5 hours (finish about 18:00–18:30)**, cumulative rate ≈ **0.29 pairs/minute**. 🔴 **2026-10-02 18:31 during inspection a gap was caught that would silently stall at 111/112 (an after-effect of my earlier `rm -rf`, invisible unless you look)**: in `_status.tsv` **P25_LUAD has a DONE row only for `agfT`, and not a single row for `agfF`** — neither DONE nor FAIL, **the table has no such pair**. Looking at the tail of `logs/nice_grid_GSM9226223_P25_LUAD_agfF.log` shows the truth: that job **11:01:35 rc=1**, reporting that `atomic_write -> gzfile` could not open `domains_agfF.tsv.gz.tmp` (**the output directory had just been deleted by my 10:54 `rm -rf`**), and its `FAIL(rc=1)` row was written by `>>` into the **already-unlinked old inode** ⇒ it could not land in the new table. Consequence chain: the driver **runs each unique pair only once** and never loops back ⇒ **it stalls at 111 forever**; the supervisor requires `d >= 112` to `launch_post` ⇒ **the downstream 03→04→06→05 never starts**, and only after the **72-hour fallback** does it exit (a silent stall, worse than reporting FAIL). **Disposition (completed)**: rerun this pair through the existing mechanism — `run_banksy_grid.sh __run_one GSM9226223_P25_LUAD FALSE` (session-detached, PPID=1, starting 18:32:18, memory 243 GB in surplus) ⇒ on success DONE goes to **112** and downstream relays automatically as usual; the process name also contains `run_banksy_grid.sh`, which incidentally makes the supervisor's `grid_alive` test true, **avoiding mounting an extra chain**. 🔴 **Lesson (folded into memory)**: deleting **the directory containing the completion table** loses not only the directory — **it also loses "failure rows", leaving some pair in neither DONE nor FAIL**; **per-pair progress checking must use `comm` (desired 112 pairs vs pairs already DONE); looking only at the DONE count cannot reveal such a gap**. | ☐ |
| 16 | M7 Target (traditional track) | — | 🟢 M7a caliber signed (merged into `10_niche/NICHE_PREREG.md` §6); **TCGA-LUAD six-step cleaning has produced numbers: 483 cases / 177 deaths (36.6%)**, then after **waiting for the M6 niche signature**, run ssGSEA + multivariable Cox + KM; M7b/M8 not started | ☐ |
| 17 | M8 Docking | — | ⬜ 0% | ☐ |

> 🔴 **2026-10-01 board back-fill correction (GP8 numbering fork)**: §2.0 split GP8 into the two numbers **GP8a** (epithelial) / **GP8b** (the remaining 5 lineages),
> but in actual execution **the six lineages were done together in one batch**, and the completion report was named **`GP8c`** ⇒ the `⬜` in the board's original rows 8/10 is a **back-fill omission of an old state, not "not done"**.
> **Two grounds**: ① `results/05_annotation/GP8c_report.md` (2026-09-22 night run, all exit codes 0, 413,697 nuclei mutually exclusive and complete);
> ② `results/08_spatial_deconv/reference_d.manifest.json` (signed by the user 2026-09-25, `lineage_layer` explicitly says "six-lineage reclustering subset",
> `cell_type_counts` = **39 L2 subtypes**) ⇒ **M5's RCTD reference uses L2, not the first-level six lineages**.
> ⚠️ **"Done" ≠ "signed off"**: the artifacts are complete, but the **gate column of both rows is still empty**, and the L2 result **still lacks a formal sign-off**.
> ⚠️ The same stale text was also carried in `STARTUP_PROMPT.md` ("next step = GP8a epithelial
> subclustering"). That file was a session hand-off note rather than project documentation, and
> has since been **removed**; its rules and environment constraints live in this file (§5b and §3).

> **GP1 (CopyKAT smoke test) is merged into GP2** and no longer listed separately: it is an **execution detail** of GP2 (extrapolate from 3 points before scheduling), not an independent gate.
> Reference: `results/03_cnv/smoke/P19_LUAD/` — that smoke test **failed to converge** (13h CPU, stalled at step 7),
> and is one of the measured grounds that drove "CNV narrowed + input changed to epithelial subclusters".

> Each time a gate is passed → update this table + record one entry of "gate evidence" (artifact path + hash).
>
> ⚠️ **GP4a/4b/4c are merged into a single execution unit**: `04_integration/10_seurat_traditional.R` runs through in one go
> (SCTransform → PCA → Harmony → two-arm kNN → clustering → UMAP), but **gate evidence is recorded separately for the three steps**,
> because they are three independently vetoable criteria and must not be signed off jointly just because "one script ran through".

### Gate Evidence

**GP4a / GP4b / GP4c / GP5 — ✅ all signed (2026-09-17)** (GP4a, GP4b full runs completed at 00:43; GP4c, GP5 signed by user ruling)
- Script: `04_integration/10_seurat_traditional.R` (`full` mode) · log `logs/seurat_full.log`
- Input: `results/02_expression/gse308103_counts_paperqc.h5ad` sha256 `a276cd1a…` (**verified on site to be consistent with the manifest**)
  = **413,697 nuclei × 18,069 genes**, nnz 655,222,534
- Recipe: `SCTransform(vst.flavor="v2", variable.features.n=3000, rv.th=1.3, ncells=5000, seed=1448145)`
  → `RunPCA(npcs=50, seed=42)` → `Harmony(group.by.vars="sample_id", dims.use=1:50)`
  → `FindNeighbors(k.param=20, prune.SNN=1/15)` → Louvain × `{0.5,0.6,0.7,0.8}` × seeds `{0,1,2,3,4}` + no-correction arm
- Elapsed **14,716.6 s = 4 h 05 m**; **peak VmHWM throughout 252.5 GB / 256 GB (98.5%), no swap, available memory once down to only 9 GB**
- Versions: R 4.2.2 · Seurat 4.3.0 · sctransform 0.3.5 · harmony 2.0.5 (paper used 1.2.0, **no claim of numerical equivalence**) · Matrix 1.5.3
- Artifact hashes: `clusters.csv.gz` `4e6e6977…` · `umap.csv.gz` `c33995dd…` · `resolution_metrics.csv` `58c35608…` ·
  `run_manifest.json` `a57c977c…` · `timings.csv` `7d8724f2…` · `embeddings/harmony_f32.bin` `93f2cf96…` · `embeddings/pca_f32.bin` `f529b2f2…`
- Report: [`results/04_integration/seurat_trad/full/GP4_report.md`](results/04_integration/seurat_trad/full/GP4_report.md)

- ✅ **GP4c signed — user ruling (2026-09-17): adopt the source paper's Harmony main caliber.**
  The escalation-clause triggering process was reported honestly, and the user ruled the main caliber = `RunHarmony(group.by.vars="sample_id", dims.use=1:50)` (consistent with the main caliber pre-registered in §M3-A.2);
  **the no-correction arm is retained as a sensitivity/limitation statement and must not be treated as a main result**.
  ⚠️ **This is an explicit ruling on the "two-arm divergence", not the divergence disappearing — the following caveats are binding on all downstream results:**
  Per-cell `ARI(Harmony, no correction)` = r0.6 **0.6575 (5/5 seeds <0.7, robust)**, r0.7 **0.6960 (1/5, marginal)**;
  kNN(k=15) predicting stage from the embedding: **before correction 0.5388 (baseline 0.4546, +8.42 pp) → after correction 0.4866 (+3.19 pp) ⇒ erases about 60% of the recoverable stage signal**;
  `sample_id`↔`stage` is one-to-one nested ⇒ **the erased portion is mathematically non-attributable**.
  ⇒ Any downstream stage-related conclusion **must not claim to have excluded batch confounding**; M4 must write this non-attributability into its interpretation boundary.
  The version discrepancy remains: paper Harmony 1.2.0 vs this machine's 2.0.5 ⇒ **no claim of numerical equivalence**.
- ✅ **GP5 signed, `r* = 0.6`** (the only legal candidates were r=0.6 / r=0.7):
  Metric 1 (cross-seed ARI ≥ 0.90 hard constraint) measured 0.8963–0.9230, **eliminating r=0.5 (0.8963) and r=0.8 (0.8999)** ⇒ this guardrail **genuinely bears weight** at full scale (at 3,000 cells it is always = 1.0, idling; at 60,000 cells 0.9346–0.9834 — **small-scale behavior must not be used to infer this scale**).
  Per the pre-registered tie-break rule (difference <0.01, take the lower) 0.9123−0.9041=0.0082 ⇒ **take r=0.6**.
  ⚠️ **Original line 305 of the script was a bare `which.max`, while the inline comment claimed a tie-break had been applied — the comment does not match the code**, so the full run reported 0.7.
  The script was fixed into an explicit tie-break implementation, and the `is_rstar` column of `resolution_metrics.csv` and the manifest's `rstar_candidate`/`metrics_table[*].is_rstar` were corrected in sync,
  and a **`rstar_correction`** block was added to the manifest for the record; **the original metrics were not changed by a single byte**. **The user accepted this correction through the process.**
- ✅ **Metric 4 (AAH absorption guardrail) — user signed off: acknowledges it is ineffective at full scale.** Across the 4 resolutions **no cluster satisfies "highest Normal proportion and >50%"**,
  so the denominator is empty ⇒ `absorption_rate` is identically 0 by definition and **provides no information for `r*`**.
  ⇒ `r*` **is in practice borne by metric 1 (cross-seed ARI) and metric 2 (cross-resolution stability)**; metric 3 **was computed later on 2026-09-21** (also idling, see the next item); metric 4 does not count.
  The guardrail is **registered but marked "ineffective at full scale"**, neither deleted nor redesigned; if it is ever restarted on other data/granularity, it must be **re-pre-registered**.
- ✅ **Metric 3 (lineage coverage) — computed later on 2026-09-21; 🔴 recomputed under the new L1 caliber on 2026-09-22.**
  🔴 **After recomputation the conclusion changed: only myeloid fails the line at `r*`** (R_mean **0.8929** < 0.90; all five seeds fail; a stricter caliber gives 0.8571).
  Of the 8,321 nuclei (13.0%) that fail the line, **98.9% are cells adjudicated into myeloid** (raw17 mast cells 6,360 + raw26 DC 1,866);
  the root cause is a **registered known limitation** (`marker_panel.py`'s myeloid panel does not include mast-cell/DC markers, and the user has ruled not to add them) ⇒ **structurally unreachable**.
  The three failing clusters are **the same set of cells** across the four resolutions, and r=0.8's "passing the line" (28/31) is **purely a denominator artifact**.
  ⇒ **User ruling 2026-09-22: keep myeloid `r*=0.5` and register the failing of the line as a known defect**, **explicitly overriding** the "must stop at the checkpoint and escalate" pre-registered in §M3-A.3.
  ⇒ The remaining six objects all still pass the line ⇒ **no discriminating power for their `r*`** ⇒ **must not** be written as "metric 3 passing provides support for `r*`".
  **Downstream constraint**: if the myeloid L2 annotates to mast cells/DC, it **requires manual interpretation**, and must not rely on the panel argmax alone.
  The caliber is "mean per-gene detection rate" (R_mean; of the two extreme readings one is always = 1.0 and the other structurally unreachable, both provably degenerate); details in `PARAMETERS §M3-A.3` and `results/05_annotation/GP5_report.md` **§6.8** (the old §6.7 is void).
  ⚠️ **Correcting my earlier statement**: the old entry's "can only be computed after GP6 annotation" is **wrong** — the registered caliber ("≥1 marker set" = any one set suffices) **does not depend on annotation**.
  ⚠️ The old `metric3_absent` text in `results/04_integration/seurat_trad/*/run_manifest.json` is a **2026-09-17 point-in-time record**, superseded by this item.
  🔴 The `caliber` field of the same batch of manifests is hard-coded to the old caliber name `A_frozen` (root cause = the string constant at `04_integration/10_seurat_traditional.R:440`).
  **2026-09-22 user ruling: do not change the script, do not rerun, register honestly** — 7 of the 9 misstate their own input (checked individually by hash),
  and when cited, `caliber` / `rstar_status` / `rstar_candidate` are **all to be treated as invalid**, with the caliber taken from `<tag>_rstar.json`.
  See `results/05_annotation/GP5_report.md` §15.1–§15.2 for details.
- Honest record: `input.n_genes` 18,069 → transformed matrix 18,047, a difference of **22** genes, produced by sctransform v2's internal filtering (checkable in the log, **not silent**).

**Step 0 · GP0 Expression Object Rebuild (paper QC caliber) — 🔶 redone (2026-09-16), validation added by the 2026-09-17 audit**
- Script: `02_expression/04_rebuild_expression_paperqc.py`; source = 75 dense text count matrices → sparse AnnData
- Mask: **M1 `qc_pass & singlet` ∩ paper QC gate** (`nFeature≥500 & nCount≥1000 & pct_mt≤20`), with the inverse's `gene_min_cells=3` applied when rebuilding the matrix
- Result: **(413,697 nuclei × 18,069 genes)**, nnz **655,222,534**; stages IAC 188,087 / Normal 94,506 / AIS 81,555 / AAH 35,283 / MIA 14,266 (**23 patients / 75 samples**)
- Artifacts: `results/02_expression/gse308103_counts_paperqc.h5ad` `a276cd1a…` (3.13 GB) · `results/01_qc/gse308103_analysis_mask_paperqc.csv.gz` `ebc74c1e…` · `rebuild_paperqc_manifest.json`
- ⚠️ **Gap found in this audit (now filled)**: GP0-redo **initially ran no gate validation and produced no report** — because `02_verify_expression_build.py`'s `EXPECT_N_OBS/N_VARS` were **hard-coded to the old object** (648945 / 18082) and never adapted to the redo. **This gap is precisely the root cause of "§6 once recorded the old 648,945 build as ✅ PASS".**
  The 2026-09-17 audit **reran all hard checks** on the paperqc object, with independent recomputation + on-site hash verification, results:
  V1 shape ✅ · V2 per-sample deviation **0** ✅ · V4 barcodes parseable & all `sample_id` in the authoritative table ✅ · V5 `stage == resolve_stage(token)` with no silent default ✅ · V6 all elements non-negative/integer, all-zero cells **0** ✅
  V3 corrected form: `nnz == Σ nFeature − 11 = 655,222,545 − 11 = 655,222,534` ✅ — **the difference of 11 is constructive**, namely the 11 nonzero counts left inside the 13 dropped genes with `<3 cells` (the old script's V3 "exact equality" form **does not hold** for the redo and needs this term)
  **Independent three-way reconciliation of mask ↔ h5ad ↔ manifest**: barcode sets fully equal (`S1 == S2`, bidirectional difference 0), `stage/patient_id/sample_id/stage_token` consistent row by row, no NaN, barcodes unique, `stage_counts` equal
  **Downstream seurat_io export** (what the R side actually reads in): `cell_names`/`gene_names` row counts and **order** both equal to the h5ad, `cell_meta` barcode sets equal, the 6 artifacts' sha256 **matching the manifest one by one**
- **Unverified items (honestly marked missing)**: the redo **has no** old A5 independent re-extraction (3-sample element-by-element recomputation) or A4b all-element integrality dedicated-script trace; this V6 already covers integrality in full, and A5-type re-extraction **has still not been done for the paperqc object**.
- 🔎 **The old 648,945 object is retained as a sensitivity arm, not deleted**: `gse308103_counts.h5ad` `f9dbe382…` · `gse308103_analysis_mask.csv.gz` `2f8bb0f6…` · `GP0_report.md` (that report describes only the old object)
- Old-object history (it did pass the gate then, **now superseded by the paper caliber**): V1–V6 12/12 + adversarial audit 22/22 all green, including A5 independent re-extraction, A2b per-cell row sum == nCount (deviation 0), A6b/c R1 closure 75/75 & 23/23; residual blind spot (row sum, nonzero count) duplicates 252,677/648,945
- Honest record: the old object's **first run crashed on** `KeyError ['cell_barcode']` → located and fixed; **the crash did not contaminate the data**. Log `logs/Step0_build_expression.run1_FAILED.stdout`

**M1 QC / doublets — ✅ PASS (2026-09-12) ｜ 🔶 Current status: per §2.0 row 1, already **downgraded to a sensitivity arm**
- Dataset: `GSE308103` (snRNA) **75 samples / 798,100 nuclei**
- Script: `01_qc/00_metrics_gse308103.R` → `01_qc/01_qc_doublets_gse308103.R` → `02_annotate_doublet_qc.R` → `03_sensitivity_nmads.R` → `06_sensitivity_doublet_rate.R` → `08_validate_doublet_calls.R`
- Result: pre **798,100** → pass **767,839 (96.21%)**; doublets **118,894 (15.48% of pass)**
- Thresholds: `nCount/nFeature` per-sample **MAD outlier** (nmads=3, log1p) + `pct_mt<5` (**set on the merits for nuclear data**, not copied from scRNA)
- ⚠️ **Caliber stacking (do not misread)**: analysis mask = **M1 `qc_pass & singlet` ∩ paper QC gate** (`nFeature≥500 & nCount≥1000 & pct_mt≤20`).
  That is, **the paper's absolute thresholds were applied all the same**; "not copied from scRNA" refers to **M1's own** decision caliber, and **not** to the absolute thresholds being rejected. Chain: 798,100 → 767,839(M1) → 648,945(∩single-cell) → paper gate on the full set 555,480 → **intersection 413,697**
- Sensitivity: nmads 3 vs 5 = +1.97 pp (**insensitive**); doublet rate vs fixed top-10% overlaps by about **62.3%** (**fairly sensitive** → a downstream "remove/keep" sensitivity must be done)
- Positive validation: **A** count features (75/75 samples, median doublet nCount ratio **2.39**; doublet rate vs cell count **r=0.921**); **B** cross-lineage co-expression (EPCAM+PTPRC+ in doublets is **7.2×** that in singlets, consistent in 72/75 samples)
- Cross-validation: **infeasible in this environment** (scrublet does not fit sparse nuclei; DoubletFinder requires Seurat 2/3 or 5) — honestly recorded
- **Resolved anomaly**: `P7_LUAD` was once called 0 doublets → root cause was an **xgb classifier collapse** (not biological) → changed to `score="weighted"` to get 11.69%; the main script now has an **automatic fallback**
- Report: [`results/01_qc/M1_validation_report.md`](results/01_qc/M1_validation_report.md) (v3)
- Artifact hashes: `gse308103_per_cell_qc.csv.gz` `44c890bb…` · `gse308103_qc_per_sample.csv` `2d0bd7df…`
- ⚠️ The v1 per-cell table has the `cell_barcode` failure defect (fread autostart skipping the barcode row) → fixed and fully rerun; the v1/v2 reports are void

**M0 Input Freeze (paired) — ⚠️ Gate not passed (redone 2026-09-15)**
- Script: `00_ingest/01_freeze_paired.py`; registry `00_ingest/cohort_registry.py`
- Result: **131 samples / 25 patients / 23 dual-modality pairs** (P3–P25); unique key = `sample_key` (`dataset:sample_id`)
- Artifacts: `results/00_ingest/paired_{samples,patients,source_files}.csv` + `M0_paired_validation_report.md` + `paired_manifest.json`
- Validation: C1–C7 green; **C8 red on 1 item** → `paired_manifest.json`'s `gate_pass=false`
- **Reason the gate was not passed (the only one)**: the tar of `GSE307534/GSM9226176` on disk is **truncated** (56,272,384 B, should be 90,677,930 B;
  `gzip -t` reports `unexpected end of file`), missing `spatial/scalefactors_json.json` and `spatial/tissue_positions.csv`.
  Measured: re-downloading yields the **complete 87 MB** tar (the manifest contains all required files, and there is **only one slice root**, `P4_AAH2`).
- **Does not block M2/M3-A**: the gap is in the **spatial** dataset (GSE307534), while M2/M3-A run only **snRNA** (GSE308103). The re-download is an M5 prerequisite,
  and the source directory `/home/eto/luad_invasion` is **read-only**, so it requires separate authorization and must be done separately.
- ⚠️ The old M0 (three scRNA cohorts) **is void with the scope narrowing**, and its scripts/artifacts are moved to `/home/eto/luad_invasion/luad_v2_out_of_scope/`.

> **Historical record (old scope, now void)**: 3 scRNA cohorts were once frozen (GSE131907 / GSE189357 / GSE148071, 420,766 cells / 109 samples / 95 patients),
> and GEO cross-validation (stage and origin) was completed. That artifact is archived and **no longer a basis for this project**.


---

## 7. Developmental-Trajectory Arm (`09_trajectory/`) — Closed Routes + Follow-up Task Register

> **This section was added 2026-09-27.** The developmental-trajectory arm previously **did not exist in any tracked document** (§6's board had no row for it);
> its code/pre-registration was committed this time via **PR #14**. This section **attaches it to the plan**: the three closed routes are recorded by status, and the un-done follow-up tasks are registered one by one.
> Two red lines come first: single-cell CNV has exited ⇒ **the whole project has no per-cell malignancy label** (**iron law 2**);
> therefore any "Tumor cell / KAC" wording **may only be reported as an expression state, and must not be treated as a malignancy call**.

**Arm status: all three main-axis conclusions are "negative or constrained", all committed honestly (not "ran successfully")**

| Route | Pre-registration | Conclusion |
| :--- | :--- | :--- |
| CellRank spectral mechanism | `10`/`11`/`12` | **Abandoned** — the assembly yields a block upper-triangular matrix, so the transport values never enter the eigenvalues at all; "endpoint = IAC" is hard-wired into the assembly, not computed from the data |
| WOT direct transport readout | `13`/`14`/`15` | §4.1 passed (small margin), §4.2 partly passed; **§4.3 gene-trend criterion degenerate and void** ("54 candidates" = the genes whose Δ is all 0) |
| CytoTRACE stage axis | `16`/`17`/`18`/`19` | **C1 monotonicity FAIL, direction opposite to the paper** (stage-median highest in IAC); the reversal is robust (22/23 patients positive, IAC highest in every depth bin) |
| Source-paper MP-panel KAC axis | `20`/`21`/`22`/`23`/`24` | D1 hits 8/27 clusters but **all fall within classic AT2**; D2 **the two sources disagree** (MP6 supports / MP9 reverses) |

### Follow-up Task Register (ordered by "can it start now")

| # | Task | Why do it | Can it start now | What is needed before starting |
| :---: | :--- | :--- | :---: | :--- |
| D1 | **Why MP9 reverses** | D2 reported MP6 ρ=+0.722, MP9 ρ=−0.376, and "the two sources disagree" is this arm's hardest falsifying conclusion; clarify what MP9's 50 genes are and whether they are inflated by the Ciliated/Club clusters | ✅ **Completed 2026-09-27** | Conclusion: **MP9 is not "another KAC definition"** — its 50 genes are an immediate-early/stress/inflammatory program (FOSB·FOSL1·EGR1·ATF3·GADD45A/B·CDKN1A·SOCS3·PTGS2·AREG·CCN1·PLAUR·LAMC2), with **zero overlap** with AT2/Ciliated/Club, and **sharing 12 genes with MP7 (`Tumor cell (stress/inflammatory)`)** (the highest among the seven MPs); MP9's highest scores fall in c4/c6 (**AT1**) and c26/c9 (whose argmax is precisely **MP7**), and its lowest fall in AT2/Ciliated. At the per-cell level corr(MP9, KAC signature) = **−0.041** (≈0) ⇒ **that −0.376 is at the cluster-mean level, measuring the "AT1/inflammation vs AT2/ciliated" identity axis, not the KAC axis**. ⚠️ An alternative explanation was not ruled out: `score_genes`' control gene pool may induce part of the artificial negative correlation. Artifacts `09_trajectory/25_kac_followup_cheap.py` / `results/09_trajectory/paper_axis/followup_cheap.md` |
| D2 | **Are clusters 16/21/23 subclones of a single patient** | K3 already flagged `top_patient_frac` = 0.561 / 0.673 / 0.862; the three clusters total 4,695 nuclei and might be contributed by just one patient | ✅ **Completed 2026-09-27** | Conclusion: **not private subclones** — the dominant patients of the three clusters are **all different** (c16→P10 56.1%, c21→P24 67.3%, c23→P5 86.2%), and each dominant patient is himself also spread across **24 other clusters**; the three clusters are **95.4% / 99.1% / 98.8% IAC** respectively, whereas each dominant patient's own body is only 71% / 54% / 37% IAC ⇒ it is "**extreme IAC enrichment + uneven patient composition**", not a private subclone. **K3's alert is downgraded accordingly** (not "single-patient oligarchy" but "stage heavily enriched × IAC distribution uneven across patients") |
| D3 | **AIC arm** | The source chain is AT2→**AIC**→KAC→precursor→invasion, and this arm did only the KAC section, missing AIC | ⛔ **No** | ⚠️ The AIC decision caliber **must** be signed first (the source markers? or MPs? this arm has no AIC panel), and the result **must be explicitly labelled "this is our inference, not a label the paper provided"** |
| D4 | **Re-cluster in MP space** | The current 27 clusters are defined by the **classic panel**, not by MPs; the argmax falling within "classic AT2" may be only a projection of the classic panel | ⛔ **No** | Changing clustering = a new caliber ⇒ must **pre-register first** (the scoring object changes from "the current 27 clusters" to "new clusters in MP space") |
| D5 | **monocle3 pseudotime** (the paper's other half) | The source developmental axis uses monocle3 + CytoTRACE; we replicated only the CytoTRACE half | ⛔ **No** | The T3 pre-registered in `16_` already says "not installed for now"; to run it, **monocle3(R)** must first be installed (or a Python replication, which must be proven equivalent), and it must be pre-registered |
| D6 | **"KAC-like state" vs "invasive-stage tumor cell" indistinguishable** | This arm's most fundamental boundary: IAC accounts for 42% of cells and the KAC-side clusters heavily overlap IAC, so the existing data **cannot separate** the two readings | ⛔ **Blocking** | Requires **a label that can give malignant/epithelial assignment** (**iron law 2**). **The intersection point = spatial CNV (SC0–SC4)**; single-cell CNV has exited, and this is the only possible unlocking path |

**Two loose ends of the WOT arm** (no new tasks; hung after the D table for reference):
- The **remedy for the voided degenerate §4.3 criterion awaits sign-off** (per **Rule 3.2**, do not unilaterally change the criterion and rerun);
- §4.4 **patient-stratified re-reporting and figure generation not done**;
- Two deviations registered but not re-signed: `N_PERM` 1000→**200**, **IAC stage excluded from the decision**.

**Intersection with other modules**: **D6 is the sole intersection of the developmental-trajectory arm and the spatial CNV arm** — the sentence the developmental-trajectory side wants to say ("is the KAC-like state an invasive-stage tumor cell")
is exactly what the spatial CNV side wants to deliver (the malignant assignment of epithelial spots). Neither arm can reach it alone; only together might they.

**Progress on the spatial CNV side (2026-09-27)**: one **caliber contradiction of that arm has been signed** — `08_spatial_deconv/SPATIAL_CNV_PREREG.md`
**§13**: `analysis_mode` set to **`samples`** (not the tool's default `subclusters`), `no_plot=TRUE`, **the cost smoke test limited to a single patient for retesting**.
The measurement that drove the ruling: the default mode on patient P4 (2 slices) **did not finish in 24 hours** and was killed by timeout, with **peak memory of only 7.67 GB**
(not a memory problem) — the root cause is that `subclusters` splits the reference into 181 + observation 65 = **246 units**, running HMM/Bayesian network on each in turn,
plus **9.5 hours of pure plotting**. Taking `samples` reduces the unit count to 2.
⇒ That arm **has still not started**, and §13.3 records two things that must be pinned down before starting (which artifact the per-spot decision comes from; SC2 runs a single patient first).
