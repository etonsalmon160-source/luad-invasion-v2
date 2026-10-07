# M7b · Genetic-Statistical Anchoring of the Target Pool (cis-MR + coloc) Operations Manual

> **Positioning**: This file is the operations manual for **M7b** of [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md).
> Purpose: to upgrade the **candidate genes** obtained from single-cell/spatial data, using **human genetics** evidence, into
> **"genetically supported candidate target"** — the highest honest tier achievable **without wet-lab experiments**,
> and it **feeds directly into M8's molecular docking**.
>
> **Status**: methodology verified (2026-09-12); **not yet run** (it belongs to M7 after M-1).
> This file describes the method, data entry points, decision thresholds, and output schema; **all values must come from real data runs** (iron law R3).

---

## 1. What Problem It Solves (Why It Is Needed)

Single-cell/spatial data gives **association**: a gene is highly expressed in MIA/IAC. Association ≠ causation; three pitfalls:
**Confounding** (driven by smoking/inflammation), **reverse causation** (cancer causes the high expression), **co-regulation** (a whole module shifts together, so the true culprit cannot be distinguished).

**MR (Mendelian randomization)** uses germline variants that are **fixed at fertilization and randomly assorted** as instrumental variables, bypassing these three issues —
equivalent to a "natural randomized controlled trial". **cis-MR + coloc** is the industry standard for drug-target validation
(PCSK9 / IL6R / TYK2 / FXI were all confirmed this way); genetically supported targets roughly double the clinical success rate.

**In one line**: MR gives the **causal direction**, single-cell/spatial gives the **cells and loci**, and the two complement each other.

---

## 2. Output Definition

For **each candidate gene G** (from the M7a niche signature and the M7 malignant program/regulon), produce one row:

| Column | Meaning |
| :--- | :--- |
| `gene` / `ensembl_id` | Candidate gene |
| `source` | Source: `malignant_program` / `niche_signature` |
| `spatial_localization` | Single-cell/spatial localization (e.g. "invasion front", "myCAF niche") |
| `n_cis_snp` | Number of cis-SNPs used as instruments |
| `mr_method` | `Wald` (single SNP) / `IVW` (multiple SNPs) |
| `mr_beta` / `mr_se` / `mr_or` / `mr_or_ci_low` / `mr_or_ci_high` / `mr_p` | MR estimates (OR is the risk ratio **per 1 SD of expression**) |
| `mr_fdr` | FDR after BH correction |
| `f_stat_min` | Minimum instrument F statistic (weak-instrument screening) |
| `egger_intercept_p` | Pleiotropy test (when multiple SNPs) |
| `steiger_ok` | Whether the Steiger direction is correct |
| `coloc_pph4` | Posterior probability of a shared causal variant |
| `coloc_verdict` | `shared`(H4) / `distinct`(H3) / `inconclusive` |
| `ot_tractability` | Open Targets tractability grading |
| `known_drug` | Whether a clinical-stage drug exists (ChEMBL/Open Targets) |
| `genetic_support` | **Final verdict**: `supported` / `not_supported` / `not_testable` |

---

## 3. Principles at a Glance

### 3.1 The Three MR Assumptions (each must be checked)
1. **Relevance**: the instrument indeed affects the exposure (cis-eQTL significant) → screen with **F statistic > 10**.
2. **Independence**: the instrument is unrelated to confounders → cis variants are naturally closer to satisfying this.
3. **Exclusion restriction**: affects disease **only through that gene** → the most fragile one, backstopped by **sensitivity analyses**:
   - **MR-Egger intercept** (Bowden 2015): intercept significantly ≠ 0 → directional pleiotropy present;
   - **Weighted median / weighted mode** (Bowden 2016): robust to some invalid instruments;
   - **Cochran's Q heterogeneity**;
   - **Steiger filtering** (Hemani 2017): ensures the direction is "variant→expression→disease" rather than the reverse.

### 3.2 Why cis-MR and coloc Are Used Together
- **cis-MR** (using only variants within gene ±1 Mb as instruments): biologically most credible, least pleiotropy → gives **direction and effect size**.
- **coloc** (Giambartolomei 2014; Wallace 2021): a Bayesian test of whether the GWAS signal and the eQTL signal are **the same causal variant**
  → excluding the illusion of "two different variants that happen to be in LD". Outputs **PP.H4** (shared) and **PP.H3** (distinct).

### 3.3 Judging with the Two Together
- MR significant **and** `PP.H4 ≥ 0.8` → `supported` (genetically supported)
- MR significant but `PP.H3` high → `distinct` (the locus's effect on disease does not go through that gene) → `not_supported`
- No cis-SNP obtainable → `not_testable` (**honestly mark as missing, do not substitute**)

---

## 4. Data Sources (reachability measured on 2026-09-12)

| Role | Source | Entry | Measured | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Exposure·cis-eQTL (first choice)** | **eQTL Catalogue** | https://www.ebi.ac.uk/eqtl/ (direct connection 200) | ✅ | Cross-tissue cis-eQTL in a **unified format**, including lung |
| Exposure·cis-eQTL (alternative) | GTEx v8 lung | GTEx Portal / dbGaP phs000424 | ⚠️ Portal not reachable from this machine | Can substitute the GTEx processed version in eQTL Catalogue |
| Exposure·cis-eQTL (blood, large sample) | eQTLGen | https://www.eqtlgen.org/ (direct connection 200) | ✅ | n=31,684; different tissue, must be noted |
| **Outcome·GWAS (first choice for practice)** | **IEU OpenGWAS** | https://gwas.mrcieu.ac.uk/ · API `api.opengwas.io` (200) | ✅ | Requires a **free token**; TwoSampleMR connects directly |
| Outcome·GWAS (authoritative) | ILCCO / TRICL lung cancer | McKay 2017 *Nat Genet*; Byun 2022 *Nat Genet* | ⚠️ | Many studies in the GWAS Catalog have `fullPvalueSet=False`; complete summary statistics may be controlled-access (dbGaP) → **availability must be confirmed first** |
| LD reference | 1000 Genomes EUR | 1000G / built into OpenGWAS | ✅ | Ancestry must match the GWAS |
| Tractability | Open Targets Platform | GraphQL API | ✅ | tractability buckets |
| Dependency/repurposing | DepMap, DGIdb, ChEMBL | Public | ✅ | Selective dependency, not broadly essential |

> ⚠️ **Largest uncertainty**: the **complete** LUAD GWAS summary statistics from ILCCO/TRICL may require an application.
> If unobtainable → switch to a lung-cancer GWAS in **IEU OpenGWAS**, and **state the source and sample size explicitly** in the report; do not substitute another phenotype.

---

## 5. Environment

```bash
# R 4.2.2 already installed. Fill in the MR/coloc dependencies:
Rscript -e 'install.packages(c("remotes","data.table","ggplot2"))'
Rscript -e 'remotes::install_github("MRCIEU/TwoSampleMR")'
Rscript -e 'remotes::install_github("MRCIEU/ieugwasr")'
Rscript -e 'install.packages("coloc")'
# Version record (outputs must carry it): sessionInfo() written to the run log
```
- `TwoSampleMR` (Hemani 2018 *eLife*) · `ieugwasr` (OpenGWAS client) · `coloc` (Giambartolomei 2014; v5 Wallace 2021)
- OpenGWAS token: `ieugwasr::get_opengwas_jwt()` (first time requires registering at https://api.opengwas.io)

---

## 6. Procedure

> What follows is a **structural skeleton**; a real run must record input hashes, versions, and seeds (project iron laws R4/R5).

```r
# 6.0 Input: candidate gene table (from M7a/M7)
cand <- read.csv("results/M7/candidate_genes.csv")   # gene, source, spatial_localization

for (g in cand$gene) {
  # 6.1 cis-eQTL (lung): take significant SNPs within ±1Mb (eQTL Catalogue / GTEx processed version)
  #     → data.frame(snp, beta_exposure, se_exposure, ea, nea, p, n)
  expo <- fetch_cis_eqtl(g, tissue = "lung")

  # 6.2 Outcome GWAS: take the effects of the same set of SNPs (OpenGWAS / ILCCO)
  out <- fetch_outcome_gwas(expo$snp, gwas_id = LUAD_GWAS_ID)

  # 6.3 Harmonization (allele alignment, avoiding strand flips)
  dat <- harmonise_data(expo, out)

  # 6.4 MR: single SNP→Wald; multiple SNPs→IVW; sensitivity analyses
  res    <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
  pleio  <- mr_pleiotropy_test(dat)      # Egger intercept
  het    <- mr_heterogeneity(dat)        # Cochran's Q
  steig  <- directionality_test(dat)     # Steiger
  fstat  <- min(dat$beta.exposure^2 / dat$se.exposure^2)   # weak instrument

  # 6.5 coloc: same region (eQTL and GWAS)
  col <- coloc::coloc.abf(dataset1 = eqtl_region, dataset2 = gwas_region)
  pph4 <- col$summary["PP.H4.abf"]; pph3 <- col$summary["PP.H3.abf"]

  # 6.6 Verdict (see §7), write one row
}
```

### Key Quality-Control Points
- **Allele harmonization**: `harmonise_data(action=2)`, removing ambiguous/palindromic SNPs;
- **Weak instruments**: `F > 10`, otherwise remove and note;
- **LD/ancestry**: the LD reference must match the GWAS population (EUR);
- **Sample overlap**: exposure and outcome samples must not overlap (eQTL and GWAS populations differ → usually satisfied);
- **Smoking confounding**: report whether a smoking-adjusted GWAS was used; if not adjusted, it must be declared in the limitations.

---

## 7. Decision Thresholds

| Metric | Threshold | Source/Nature |
| :--- | :--- | :--- |
| Number of cis-SNPs | ≥ 1 (single SNP uses Wald) | Methodology |
| Instrument F statistic | **> 10** | Weak-instrument convention (Burgess & Thompson 2011) |
| MR P value | **BH-FDR < 0.05** | Project convention (multiple testing) |
| MR-Egger intercept | P > 0.05 (no directional pleiotropy) | Bowden 2015 |
| Steiger | Direction correct | Hemani 2017 |
| coloc PP.H4 | **≥ 0.8** | Convention (Giambartolomei 2014) |
| coloc PP.H3 | High → judged `distinct` | Convention |

> ⚠️ The thresholds are **conventions/project agreements** and must be registered in `docs/PARAMETERS_AND_SOURCES.md`; **thresholds must not be tuned after the fact**.

---

## 8. Limitations and Red Lines (Wording Rules)

- ✅ Permitted to say: **"genetically supported candidate target"**, **"prioritized target"**, **"repurposing candidate"**
- ❌ Not permitted to say: `causal` (unconditionally), `driver`, `drug target`, `validated`
- MR is **population-level** causation and is **not equivalent** to "this gene drives MIA in a particular cell" — cell-level conclusions still come from single-cell/spatial data;
- Lung eQTL sample size is limited → some genes will be `not_testable`; **honestly mark as missing, do not substitute**;
- Target **structural docking** (M8) remains a **computational hypothesis**; `validated` requires wet-lab experiments (SPR/ITC/organoids).

---

## 9. Reproducibility and Output Requirements (iron laws)

1. Inputs (eQTL table, GWAS summary statistics, candidate gene table) **all carry a SHA-256**;
2. Output table: `results/M7b/mr_coloc_target_anchoring.csv` + run log (`sessionInfo()` + parameters + seed);
3. Intermediate outputs per gene (harmonised table, coloc region) retained for review;
4. All scripts deterministic, with no `np.random`/`runif` and no hard-coded values.

---

## 10. To-Do (executed after M-1)

- [ ] Confirm LUAD GWAS availability: take a usable ID from OpenGWAS **or** obtain ILCCO/TRICL summary statistics (record source and n)
- [ ] Download eQTL Catalogue lung cis-eQTL (record version)
- [ ] Register an OpenGWAS token
- [ ] Register the §7 thresholds into `docs/PARAMETERS_AND_SOURCES.md`
- [ ] Run 1 gene (smoke test) with M7a's candidate genes → then full scale
