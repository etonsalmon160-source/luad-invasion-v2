# 🛑 Spatial Transcriptomics 6-Stage Slice Mapping and Figure Rigor: Supreme Prohibitions (Spatial Cohort & Figure Prohibitions)

> **Level in force**: supreme research-integrity iron law (P0-level permanent constraint); at every stage, any script, report, or agent behavior must comply unconditionally, and any violation is treated as academic misconduct.

---

## 1. The 6-Stage Spatial Slice Mapping Iron Law (Golden Standard Slice Manifest)

In every script across the entire project that generates Figure 1, spatial matrices, microenvironment deconvolution, and multimodal atlases, **the slice-to-stage mapping must be 100% locked as follows; any form of unauthorized substitution or fallback is strictly forbidden**:

| Stage | Stage Standard Name (Label) | Sole Legal Slice ID | Authoritative Database | GSM / Sample ID | Histological Features and Strict Prohibitions |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Stage 1** | `Stage 1: Normal in PT` | **`P4_Normal`** | GSE307534 | GSM9226174 | Normal lung parenchyma control; use of any cancer-containing slice is strictly forbidden. |
| **Stage 2** | `Stage 2: AAH Precursor` | **`P1_AAH`** | GSE307534 | GSM9226168 | A genuine atypical adenomatous hyperplasia (AAH) precursor slice. |
| **Stage 3** | `Stage 3: AIS In Situ` | **`P3_AIS`** | GSE307534 | GSM9226172 | Adenocarcinoma in situ (AIS) slice with pure lepidic growth. |
| **Stage 4** | `Stage 4: MIA Micro-inv` | **`P10_MIA`** | GSE307534 | GSM9226189 | **【Iron-law prohibition】A genuine microinvasive `P10_MIA` must be used; using `P4_AAH1` / `P4_AAH` as MIA is absolutely forbidden!** |
| **Stage 5** | `Stage 5: IAC Invasive` | **`P3_LUAD`** | GSE307534 | GSM9226173 | **【Iron-law prohibition】Primary invasive adenocarcinoma (IAC), solid/acinar pattern with fibrotic stroma!** |
| **Stage 6** | `Stage 6: LNM Metastasis` | **⛔ No legal slice available for now (awaiting real LUAD data)** | — | — | **【2026-09-12 erratum】The originally locked `PT_3_LNM`(GSE190811) was verified against GEO as a 【breast cancer】 lymph-node metastasis (series title "…breast cancer patients"), and its GSM number GSM5732148 does not exist in that repository (the real ones are GSM5732357–2360). That slice is 【void】 and forbidden for any LUAD output; the Stage 6 spatial anchor is currently missing, to be re-locked once the project provides real LUAD lymph-node-metastasis spatial transcriptomics data.** |

> **⛔ Erratum record (2026-09-12)**: This file originally locked `PT_3_LNM`(GSE190811, GSM5732148) as the permanent Stage 6 slice,
> but per-sample verification against GEO identified it as **breast cancer** data. It is voided. **LNM-stage spatial outputs (deconvolution / niche / PLIP / docking) must not be produced in any form until it is replaced**;
> in the interim the spatial atlas is limited to **Normal → AAH → AIS → MIA → IAC** (all from the single GSE307534 repository on a single platform).
> At the single-cell level, LNM may still use GSE131907 `mLN` (genuine metastatic lymph nodes, GEO-verified among 44 patients).

---

## 2. The Three Red-Line Prohibitions for Algorithms and Figure Rendering (Rendering Prohibitions)

1. 🚫 **Synthetic modulo or fabricated striping is strictly forbidden (Zero Synthetic Modulo Logic)**:
   - The WHO pathological classification (Lepidic, Acinar, Papillary, Solid, Fibrotic Stroma) and the TCGA molecular subtypes (TRU, PP, PI) must be 100% computed from genuine RCTD subclones and microenvironment deconvolution weights.
   - **It is absolutely forbidden to introduce any artificial alternating-stripe code such as `i % 2 == 0`, `i % 3 != 0`, or `np.random`!**

2. 🚫 **Painting false positives across non-tumor background is strictly forbidden (Clean Non-tumor Background Rule)**:
   - Non-tumor cells (normal alveolar epithelium, lymphocytes, macrophages, vascular endothelium, etc.) must uniformly keep a clean light-gray background (`#E2E8F0`) in both the WHO and TCGA columns, strictly consistent with `Non-tumor / NA` in the legend; assigning a malignant-subtype color to normal tissue is strictly forbidden.

3. 🚫 **The ROI box must not degenerate into a pure-yellow tumor block (Multilineage Frontier ROI Law)**:
   - The dashed box in column 5 and the magnified hexagon in column 6 must show a **"multilineage invasion-front breach (Multilineage Frontier)" in which malignant tumor (yellow/orange/red), fibroblast stroma (dark-red myCAF), macrophages (purple TAM), T cells (blue), and endothelial cells (pink) are in close contact and actively remodeling**.
   - **It is absolutely forbidden to center the ROI on a 100% pure-malignant-cell yellow block inside the tumor parenchyma!**
