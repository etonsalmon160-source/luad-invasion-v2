# Cover letter — Computational Biology and Chemistry

*Draft. To be adjusted once the journal's current guide for authors is confirmed.*

---

Dear Editors,

We submit for your consideration our manuscript, **"Paired single-nucleus and spatial
transcriptomics of the lung adenocarcinoma invasion sequence reveal niche domains and
candidate reversal compounds"**, for publication as a research article in *Computational
Biology and Chemistry*.

**What the paper does.** We assembled paired single-nucleus RNA sequencing and Visium
spatial transcriptomics from 23 patients spanning the full lung adenocarcinoma invasion
sequence, and used the two modalities together to ask two questions: how the
transcriptional state of the tissue changes as invasion proceeds, and whether that state
can be reversed *in silico* by known perturbations.

**Why we believe it suits this journal.** The study is, at its centre, an evaluation of
what the computational methods now standard in spatial transcriptomics can and cannot
establish on real data — and we report that evaluation quantitatively.

- **A confound that propagates and can be measured.** Sequencing depth and tissue density
  are the same measurement in these sections. We show that this single coupling moves
  domain identity, marker-gene stability, the apparent stage-dependence of spatial
  programmes, the cohort-level copy-number difference and one domain's prognostic weight —
  and we quantify how far each one moves. Several apparently positive results disappeared
  once depth was matched. The finding is not specific to this cohort; it applies to any
  spatial study in which tumour and non-tumour regions differ in cellularity.

- **Method concordance tested where it can be tested.** Domain structure is
  method-sensitive: three defensible partitions of the same data return 7, 6 and 5
  domains, with an adjusted Rand index of 0.352 between the expression-based and
  RCTD-aligned partitions. We also report that a spatial-coherence criterion we built
  ourselves is passed by every configuration we tried, including one that uses no spatial
  information at all — a check that therefore carries no weight, and which we say so.

- **Verification against reference implementations.** Our reimplementation of the CMap
  weighted-KS statistic reproduces the official `signatureSearch` implementation to
  6.6 × 10⁻¹⁴ with 2000/2000 signs agreeing, and our vectorised causal score reproduces
  the authors' Python at Spearman 1.00000000. A positive control with a known answer
  (the official epiblast-to-nascent-mesoderm benchmark, TBXT ranked 1st of 8,450)
  establishes that an empty result from this pipeline reflects the query rather than a
  broken method.

**What we deliberately do not claim.** The manuscript contains an extensive negative
section, and we would rather it be read as rigour than as incompleteness. Six analyses
failed, and for each we identified why: single-cell copy-number inference could not
produce a usable malignant/non-malignant separator across six tool families; a
developmental trajectory ran opposite to the source study's and our measurement cannot
distinguish two readings of that; a spectral trajectory read-out turned out to be a
property of how the matrix was assembled rather than of the data; and a ligand-level
analysis produced no usable pair. We report these because they constrain what the positive
results can mean, and because the diagnoses are reusable.

**Reproducibility.** All raw data are public (GSE308103, GSE307534). All analysis code,
including the scripts that generate every figure, is openly available, together with a
per-parameter provenance record that labels each analysis setting as coming from a source
paper, a software default, or from us. We have also published an audit of every
quantitative claim in the manuscript against the artifact that produced it, including the
errors that audit found and corrected.

**Declarations.** The corresponding author is a reviewer for *Computational Biology and
Chemistry*. We disclose this for transparency and leave the handling of the review process
entirely to the Editor. The authors declare no competing interests. This work received no
specific funding.

This manuscript has been posted as a preprint on Research Square and on bioRxiv; it has
not been published by a journal and is not under consideration elsewhere. Both postings are
disclosed here for transparency. All authors have approved the submission and agree to its
content.

We hope the manuscript is suitable for *Computational Biology and Chemistry*, and we thank
you for considering it.

Yours sincerely,

**Zhiyang Li**
Guangdong Medical University, Zhanjiang, Guangdong, China
eto-1024@gdmu.edu.cn

*on behalf of both authors*
