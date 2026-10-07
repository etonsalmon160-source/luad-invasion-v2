# M0 / M1 Academic-Grade Review Report

> ⚠️ **Historical document (2026-09-12 snapshot); its M1 portion is no longer the final analysis caliber.**
> On 2026-09-16 the user ruled to rebuild the analysis mask using the **source paper's fixed QC**: final caliber = M1 ∩ paper gate = **413,697 nuclei × 18,069 genes**.
> The M1 described here (767,839 passed, 648,945 singlets) is now a **sensitivity arm**. The M0 conclusions (including the `GSM9226176` truncation) **remain valid**.
> For the current caliber see [`PLAN_AND_CHECKPOINTS.md`](../PLAN_AND_CHECKPOINTS.md) §2.0 / §6.

> **Review date**: 2026-09-12 · **Scope**: M0 input freeze (paired datasets) + M1 QC/doublets
> **Authoritative source**: GEO per-sample metadata (GSE308103 / GSE307534, `targ=gsm&form=text`)
> **Review purpose**: to investigate **data/label mapping confusion** and any issue that could seriously affect downstream stages

---

## Conclusion

**No serious academic problems were found.** Label mapping vs GEO shows **zero per-sample discrepancy**; outputs are joinable, reproducible, and hashed.
All findings are **documented limitations** (§3) and do not constitute a blocker.

---

## 1. Verification Items and Results

| # | Check | Method | Result |
| :--- | :--- | :--- | :---: |
| A1 | **Stage label mapping** | Per-GSM comparison of GEO title vs freeze table (both datasets) | ✅ **0 discrepancy** |
| A2 | **Stage counts** | Compared against GEO title tallies | ✅ GSE308103: AAH 9 / AIS 14 / IAC 24 / MIA 4 / Normal 24 (**fully consistent**) |
| A3 | **Sample coverage** | Local vs full GEO set | ✅ GSE308103 75/75; GSE307534 local **56/56** (completed 2026-09-15; only the `GSM9226176` tar is truncated, see §3-①) |
| B1 | **Patient identity** | Taken from GEO titles `… of patient N` | ✅ The `P*` numbering agrees across both datasets; **GEO has no independent patient field**, and pairing is established on this basis (§3-②) |
| B2 | **Pairing relationship** | Intersection of dual-modality patients | ✅ **23 cases** (P3–P25), each with slices in both modalities. **2026-09-15 correction**: the old value 9 was an artifact of "only 19/56 spatial slices downloaded" |
| C1 | **cell_barcode integrity** | Uniqueness + join with the source matrix | ✅ 798,100 rows / 798,100 unique; spot-check of P3_Normal matches the matrix **exactly (2631=2631)** |
| C2 | **Gene set** | Source-file line count | ✅ All **18,082** (including the previously lost SAMD11) |
| D1 | **QC ↔ doublet-annotation consistency** | Logical cross-check | ✅ QC-fail but tested 0; QC-pass but untested 0; tested with missing score 0 |
| D2 | **per-cell vs per-sample tables** | Three-way count alignment | ✅ **0 inconsistencies** |
| E1 | **Validity of doublet calling** | Count features + cross-lineage co-expression | ✅ Median count ratio 2.39 (75/75); cross-lineage enrichment **7.2×** (72/75) |
| F1 | **Anomalous sample** | All-sample distribution | ✅ `P7_LUAD` root cause located and fixed (§2) |

---

## 2. Located and Fixed Issues

| Issue | Root cause | Impact | Disposition |
| :--- | :--- | :--- | :--- |
| `cell_barcode` all equal to `"0\|<sample>"` | `fread` **autostart skipped the barcode row** (first-row field count ≠ data rows) | Identifier voided, cannot join; also lost 1 gene | Read the barcode row separately + `skip=1` + `stopifnot` assertion; **full rerun** |
| After the fix, one extra data row was dropped | `dt[-1,-1]` dropped rows redundantly | Gene count disagreed with the matrix (immediate error) | Changed to `dt[,-1]` |
| `P7_LUAD` doublet rate = 0 | **scDblFinder `xgb` classifier collapse** (non-biological; the patient's other two samples were normal) | No removal for that sample | Changed to `score="weighted"` → 11.69%; **main script gets an automatic fallback** |
| Normalization lost the original GEO label | No traceback after `LUAD` → `IAC` | Insufficient transparency | Freeze table gains a **`stage_token`** column preserving the original label |
| Spatial samples lacked a GSM | The source-file table did not record the GSM in the directory name | Insufficient traceability | Source-file table gains a **`gsm`** column |

---

## 3. Known Limitations (non-blocking, but downstream stages must declare them)

**① Spatial data local 56/56 slices (completed 2026-09-15), but 1 tar is truncated**
GSE307534 has 56 samples in full on GEO; **56** are decompressed locally. The only defect: the `GSM9226176` tar is only 56,272,384 B (should be 90,677,930 B),
`gzip -t` reports unexpected EOF → missing `spatial/scalefactors_json.json` and `spatial/tissue_positions.csv`.
Measured: re-downloading yields the complete 87 MB tar (**with only one slice root**, `P4_AAH2`). This is an M5 prerequisite; the source directory is read-only, so it requires separate authorization to fill in.

**② Patient identity is inferred from GEO titles**
GEO has no independent `patient id` field; the `P*` numbering is taken from the title `… of patient N` (the title itself is GEO-authoritative).
The pairing relationship (**23 cases**) is built on this basis — it has been committed as the **GEO authoritative sample tables** `00_ingest/geo_metadata/*_samples.tsv` (per-GSM cross-verification, 0 discrepancy).

**③ The stage normalization `LUAD → IAC` is a project convention**
That study names invasive-stage samples `LUAD` (in parallel with AIS/MIA), hence the normalization to `IAC`. The original label is preserved in `stage_token`.

**④ Doublet detection uses a single method**
A second method is infeasible in this environment: **scrublet** can detect only ~0.5% on sparse nuclei; **DoubletFinder** requires Seurat 2/3 or 5, while this machine runs 4.3.0.
→ Replaced with two independent validations (count features, cross-lineage co-expression), and **honestly declared**.

**⑤ The doublet-removal proportion is somewhat sensitive**
The main result overlaps a fixed top-10% by about **62.3%** → at key analyses a downstream "remove/keep doublets" sensitivity **must be performed**.

**⑥ The pooled doublet rate of 15.48% is misleading**
It is **weighted by cell count** (large samples inflate it); the "typical-sample" rate has a median of **10.6%**. It should be **presented per sample** and noted as a model estimate.

**⑦ Spatial spot QC not performed**
M1 covers only GSE308103 (snRNA); spot-level QC for GSE307534 belongs to spatial preprocessing and has not yet been executed.

---

## 4. Review Verdict

| Dimension | Verdict |
| :--- | :--- |
| Label/stage mapping | ✅ **No confusion** (GEO per-sample 0 discrepancy) |
| Data integrity (barcode/gene/join) | ✅ Pass |
| Internal consistency (QC↔doublet, cross-table alignment) | ✅ Pass |
| Outlier handling | ✅ Root cause located and fixed, not merely "flagged and dismissed" |
| Transparency and traceability | ✅ `stage_token` / `gsm` added |
| Known limitations | ⚠️ 7 items, all documented, **constituting no blocker** |

**→ Conclusion at the time: may proceed to the next step (M2).**
> ⚠️ Addendum: M1 was subsequently downgraded to a **sensitivity arm**, and the analysis caliber was switched to the paper QC (413,697 nuclei). "M2" here refers to the CNV in the plan at that time;
> in the current plan, CNV has been repositioned as **GP2**, and **must come after GP8a epithelial subclustering** (see PLAN §Key Ordering Constraints).
