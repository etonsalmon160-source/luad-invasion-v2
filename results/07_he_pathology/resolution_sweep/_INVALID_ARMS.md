# 本目录里哪些输出作废

`INVALID_` 前缀的四个臂（`INVALID_F1P0` / `INVALID_F1P1` / `INVALID_F2P0` / `INVALID_F2P1`）
以及 `INVALID_resolution_sweep_summary.json`：**作废，不得引用。**

原因：它们的源图 `aligned_tissue_image.jpg` 是**灰度的**（`std(R−G) = 0.000`），
而 PLIP 依赖颜色 —— 灰度裁框被读成「空玻璃片」
（实测 `none_A` = 0.906 / 0.930，`neo_A` = 0.037 / 0.035）。
跑出来的不是分辨率结果，是"喂错输入"的结果。

`H448` 臂（`per_slide/H448/`）**有效**：hires 源、旧基线配置，
单张 P1_AAH 上 `neo_A` = 0.2092（旧基线全量 0.2205）⇒ 代码路径正常。

`INVALID_source_highfreq.png`：**作废。** 它的对照项取「hires 放大 **3 倍**」，
建立在「aligned 与 hires 像素数差 3 倍 ⇒ 同视野」这个**未验证假设**上；
实测视野比是 **1.30** ⇒ 对照项选错。
（「aligned 图上可见方块状拼接伪影」这一观察仍在，但**结论不下**。）

`source_audit_ncc.csv` 里 **`aligned_tissue_image.jpg ↔ tissue_hires_image.png` 那一格作废**——
同一原因（该表按"各自归一化到 800×800"算整幅相关，只对同视野图对有效）。

`source_correspondence.png` / `.csv` / `_summary.json`：**有效**（现行测法，
`07_he_pathology/16b_source_correspondence.py`）。结论：
aligned 与 hires 是**同一块组织、同一朝向**（4/4 张最佳摆法 `as-is`，镜像最差），
**视野比 q = 1.30**（4/4 一致）⇒ aligned 视野 ≈ 14.74 mm、≈ **2.457 µm/px**，比 hires 细 **2.31 倍**。

`aligned_detail.png` / `.csv` / `_summary.json`：**有效**（`07_he_pathology/16c_aligned_detail.py`，
按**实测 2.31 倍**重做，无信息对照换成**理想插值** —— Nyquist 以上功率按构造为 0）。
结论：同位置像素 rho **0.83 / 0.87 / 0.89**；亚像素配准后 hires 能表达的**每个**尺度两图都一致
（≥50 µm **0.92**、22–50 µm **0.84**、11–22 µm **0.70**；整数对齐时细带那个 **−0.39 是配准残差**，
残差 ≈ (+0.8,+1.0) px，**不是数据性质**）；aligned 在 hires 够不到的频段**仍有功率**、
JPEG 8×8 块指纹满格比 **1.10**（无峰）
⇒ **aligned 确实同景、更细、细节是真的**。⚠️ 但**灰度** ⇒ PLIP 不可用；图上另有**拼接缝**（灰度背景不同的矩形块）。

**⇒ 结论**：**比 hires 更细的彩色 H&E 在本 deposit 里不存在。**
（分辨率**是**有的 —— aligned 真细 2.31 倍且细节是真的；**缺的是颜色**。）

详见 `PREREG.md` 顶部横幅与 `PREREG_v2.md` §1；证据图 `source_audit.png`、
`source_correspondence.png`。
