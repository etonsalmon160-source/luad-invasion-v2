# Number-traceability audit of the manuscript

**Date:** 2026-10-07 · **Scope:** every quantitative claim in the manuscript checked against
the artifact that produced it, in the repository.

This record exists because the audit found errors. It lists what was wrong, what it was
changed to, and the artifact the new value comes from, so that a reader can re-check any
line of it. The audit was run against the manuscripts in `paper/tex/` and the analysis
outputs under `results/`.

---

## Findings that were corrected

| # | Claim in the manuscript | Was | Now | Artifact |
| :-- | :--- | :--- | :--- | :--- |
| 1 | Agreement between the expression-based domain partition and an RCTD-aligned partition | ARI 0.225--0.247 | **0.352** | `results/10_niche/consensus/agfF/sensitivity_rctd_ari.tsv` |
| 2 | Stability of the archetype clustering, described as *cross-seed* | "cross-seed stability 0.945" | **subsampling stability 0.945**, with the cross-seed ARI (median 0.782) reported separately | `consensus/*/k_stability.tsv` (80% subsample × 100, see `03_consensus.R`) and `consensus/*/selection_report.txt` |
| 3 | Domain enrichment fold changes (Table 5) | D4 2.5 / 1.9, D5 0.480, D6 5.1, D7 9.1 / 4.7 | **3.56 / 2.43, 0.426, 4.26, 4.82 / 3.13** | `results/10_niche/domain_rctd_subtype_enrich.tsv`, `domain_rctd_subtype_prop.tsv` |
| 4 | Perturbations reaching the top 100 in ≥1 subtype axis | 1,748 | **1,466** (gene and direction) | `results/10_niche/tr_singlecell/PERCELLTYPE_causal.tsv`, recomputed |
| 5 | D5 copy-number positive in precursor sections | 17 of 20 | **21 of 25** | `results/10_niche/domain_cnv_delta_perslide.tsv`, recomputed |
| 6 | Nuclei "retained after quality control" | 399,579 | **413,697** after QC; 399,579 is the six-lineage reference subset | `results/01_qc/regate_paper_qc_report.json`; `results/08_spatial_deconv/reference_d.manifest.json` |
| 7 | Registration residual | within 1 px on 53 of 56 | **within 1 px on 46 of 56** | `results/07_he_pathology/resolution_sweep/registration_final.csv`, `shift_native` |
| 8 | Per-section spot-removal range | 0.88% to 29.6% | **0.02% to 29.6%** (median 2.6%) | `results/08_spatial_deconv/spot_mask.tsv.gz`, recomputed per slide |
| 9 | D5 score-versus-UMI correlation | r = 0.740 | **r = 0.755** (median across sections) | `results/10_niche/kstar_diag/d9_d5_depth_control.tsv` |
| 10 | Split-half top-500 compound overlap | 438 on average | **423 on median** | `results/10_niche/sigsearch/splits/cor_split*_h{1,2}.tsv`, recomputed |
| 11 | Query axis length | "a 16,104-gene axis" | **18,069 genes, of which 16,104 map to SCMG** | `results/06_scmg/scmg_branch_summary.json` (`n_genes_mapped`) |

Two of these matter more than the others. **#1** used a value from an analysis that had been
superseded and archived under `consensus_STALE_lam1.0_20261002/`, so the manuscript was
quoting a retracted configuration. **#2** labelled a subsampling stability as a cross-seed
one, which made the Results contradict the Limitations a few pages later, where the true
cross-seed figure (below the 0.90 gate) is reported. The two quantities are now named
separately wherever they appear.

## Verified as correct, no change needed

- Cohort counts: 798,100 called; 767,839 passing QC; 118,894 doublets (15.48%); 92,337
  normal nuclei; 56 sections / 25 patients / 23 paired; 14,336 spots per section; 18,082
  probe panel and 18,085 spatial probes; sections per stage 1/11/14/4/26.
- 39 L2 subtypes, 67 lineage–subtype combinations, Basophil/Mast 1 at 6,360 epithelial and 8
  myeloid cells.
- Domain structure: 434 domains; K-stability values 0.945 / 0.966 / 0.814; the three
  partitions resolving 7, 6 and 5 domains; domain shares and stage occupancy.
- Depth matching: the 1.00–7.95× range; D3 with 16 of 150 markers surviving; Spearman
  ρ = −0.89; the PLIP density correlation ρ = 0.754.
- ECM programme: 20 of 22 patients; +0.402; P = 0.0033; the 14-gene and 8-gene subsets; all
  22 leave-one-out variants.
- Single-cell copy number: the 0% / 13.6% / 49.8% instability; AUROC 0.5000 on 19,748 cells;
  AUROC 0.652 with the 0.0826 control.
- Cohort-level copy number: the 71% shrinkage; p = 0.092; D5 Δ = +0.044 in 40 of 46; D3 in
  12 of 12; the per-domain significance pattern.
- Prognosis: 483 patients / 177 deaths and 381 / 143; the D7 and D3 hazard ratios; the D3
  adjusted HR moving 2.301 → 1.303 at p = 0.090.
- The compound-level reversal in full: 34 compounds; chance 5.48; 16 PI3K/mTOR; the
  correlation-scheme agreement values; median ρ = 0.962 with one split at 0.752; cross-modal
  ρ = 0.7345 and 52 shared; the positive control ranking TBXT 1st of 8,450 through EVX1 15th;
  the causal-score cross-check; 9,878 gated entries and the tier counts; 31 axes with two
  survivors; SFTPC at 22 and STAT1 at 21; the NicheNet outcome.
- Imaging: the 0.729 ceiling; AUC 0.972 with ρ = 0.574 falling to 0.759 when density is
  projected out; the density-neutral prompt at 0.833 / 0.129 / 0.791; lesion score 0.718 and
  −0.794; AAH-versus-Normal AUCs 0.2836 / 0.3047 / 0.3228; 1,344 prompt combinations; the
  four image-source correlations.
- Registration: 4.6233 µm/px; 56 of 56 passing; reflection applied to 53 of 56; the affine
  grid residual ≤ 0.18 px; the NCC gate values.
- Spatial QC cohort-wide: 4.97% removed against the source study's 1.05%.

## Registered as weak provenance, not corrected

Three numbers appear in the manuscript but could not be traced to an artifact in the
repository. They are not contradicted by anything; they simply have no on-disk source, and
are recorded here rather than silently kept.

- **"19,819 entries"** in the causal-score numerical cross-check. The cross-check itself is
  traceable for the weighted-KS and connectivity statistics (2,000 signatures in
  `results/10_niche/target_reversal/xcheck/`), but no artifact states the 19,819 figure.
- **"237 of 2,439 genes at the minimum attainable $p$"** for the cosine projection. The 2,439
  matches the row count of `results/10_niche/tr_singlecell/scmg_k562_control.tsv`, but the
  237 has no artifact behind it.
- **"20,345 perturbations"** for the SCMG library. This is documented in the project's prose
  records but no in-repository artifact states the count.

These should be re-derived from the scripts before the manuscript is finalised, or removed.

## Method

Every DOI in the reference list was resolved through Crossref and its returned title compared
against the cited title; all 46 resolved and matched. Each in-text citation was read against
the claim it supports. Every software tool named in the Methods was checked for its presence
in the analysis scripts. The numeric checks above were run by reading the named artifact or
by recomputing from it.
