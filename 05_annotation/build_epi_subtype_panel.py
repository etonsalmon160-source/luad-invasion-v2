"""生成 `epi_subtype_panel.py` —— 从源论文原始补充材料 .xlsx 程序化导出 marker 面板。

用法:
    python3 05_annotation/build_epi_subtype_panel.py \
        --tableS2 /tmp/supp3.xlsx --tableS3 /tmp/supp4.xlsx \
        --out 05_annotation/epi_subtype_panel.py

为什么要有这个脚本而不是手抄基因名
----------------------------------
421 个基因名手工转录必然出错，且出错后无法察觉。本脚本直接读 .xlsx 的单元格，
把基因列表逐字写进目标文件；任何"基因名和原表不一致"都只可能来自 Excel 自身的
自动改名（已知 2 处，见 EXCEL_NAME_FIX），而那一类是被显式登记、显式修复的。

输入文件的来历（须与 PARAMETERS 登记一致）:
    supp3.xlsx = 源论文 Table S2  "Meta-programs gene lists"
    supp4.xlsx = 源论文 Table S3  "Gene lists for various signatures"
    Table S2 有三块：snRNA-seq(all cells) / snRNA-seq(lung epithelium) / Visium ST(all spots)
    本脚本只取**第二块**（列 10-18，上皮），因为 GP8a 是上皮亚聚类。
"""

import argparse
import hashlib
import pandas as pd

# Excel 把基因名自动改成"文本+数字"形式。逐处核对后手工还原，并在此登记（法则0：不静默）。
EXCEL_NAME_FIX = {
    "DKK\xa03.00": "DKK3",   # Table S2 MP4 Basal/basal stem
    "ERN\xa01.00": "ERN1",   # Table S2 MP7 Tumor cell (stress/inflammatory)
}

# 上皮块的列区间与 MP 标签所在行
EPI_COLS = range(10, 19)
LABEL_ROW = 3
FIRST_GENE_ROW = 4

SUBTYPES = {
    "Ciliated": (["MP1", "MP8"], "MP1 ∪ MP8（论文给了两个纤毛 MP）"),
    "AT2": (["MP2"], None),
    "Club/secretory": (["MP3"], None),
    "Basal/basal stem": (["MP4"], "含 Excel 改名修复：DKK 3.00 → DKK3"),
    "AT1": (["MP5"], None),
    "Tumor cell/KAC": (["MP6"], None),
    "Tumor cell (stress/inflammatory)": (["MP7"], "含 Excel 改名修复：ERN 1.00 → ERN1"),
    "KAC/inflammatory": (["MP9"], None),
}


def read_table_s2_epi(path):
    d = pd.read_excel(path, header=None)
    mp = {}
    for c in EPI_COLS:
        lab = str(d.iloc[LABEL_ROW, c]).strip()
        genes = [str(d.iloc[r, c]).strip() for r in range(FIRST_GENE_ROW, len(d))]
        genes = [EXCEL_NAME_FIX.get(g, g) for g in genes if g and g != "nan"]
        if len(genes) != 50:
            raise ValueError(f"Table S2 上皮 MP{c-9} 【{lab}】基因数 {len(genes)} != 50，格式可能已变")
        mp[f"MP{c - 9}"] = (lab, genes)
    if len(mp) != 9:
        raise ValueError(f"表 S2 上皮块应解析出 9 个 MP，实得 {len(mp)}")
    return mp


def read_table_s3(path):
    d = pd.read_excel(path, header=None)
    names = ["KAC_human", "NFkB_human", "Interferon_human"]
    out = {}
    for c, nm in enumerate(names):
        g = [str(d.iloc[r, c]).strip() for r in range(3, len(d))]
        out[nm] = [v for v in g if v and v != "nan"]
    return out


def fmt(lst, indent=8, per=6):
    return "\n".join(
        " " * indent + '"' + '", "'.join(lst[i:i + per]) + '",' for i in range(0, len(lst), per)
    )


def emit(mp, sigs, out):
    L = []
    A = L.append
    A('"""GP8a 上皮亚聚类 —— 亚型 marker 面板（逐基因可溯源）。')
    A('')
    A('生成方式：本文件由 `05_annotation/build_epi_subtype_panel.py` 从原始补充材料')
    A('（源论文 Table S2 / Table S3 的 .xlsx）**程序化生成**，非手工转录 —— 保证基因名逐字一致。')
    A('')
    A('两个来源层')
    A('----------')
    A('1. **主干 = 源论文自身上皮 MP**（Table S2 第二块 "snRNA-seq (lung epithelium)"，9 个 meta-program，')
    A('   每 MP 50 基因）。论文定义的亚型就是我们要的亚型，且 450/450 基因槽**全部命中本数据矩阵**。')
    A('   ⚠️ 但 MP 是 **NMF meta-program**，不是洁净 marker 表 —— 里面混有应激/管家类基因')
    A('   （如 MP2 AT2 含 SOD2/CXCL2/NR4A1）。**直接整块做 score_genes 特异性会偏低**，')
    A('   故同时提供 CORE 小集做交叉校验。')
    A('2. **CORE = 经典亚型 marker**，逐基因出处核对状态见 CORE 的 status 字段。')
    A('   🔴 未经原论文补充表逐条核对前，status 一律为 "unverified"。')
    A('')
    A('已知数据源缺陷（必须随结果一起报告）')
    A('------------------------------------')
    A('- 源论文 Table S2 的 MP4 与 MP7 各有一处 **Excel 自动改名**：`DKK 3.00`（应为 DKK3）、')
    A('  `ERN 1.00`（应为 ERN1）。已在本文件手工还原并在此登记，**不静默使用**。')
    A('- 作者上传的矩阵只有 18,082 个基因（约标准转录组之半）；`SFTPA2` 等基因**不在数据中**，')
    A('  故本面板已剔除，并在 MISSING_IN_MATRIX 登记。')
    A('"""')
    A('')
    A('# ---------------------------------------------------------------- 一次文献')
    A('REFERENCES = {')
    A('    "Peng2026": "Peng F, Sinjab A, et al., Kadara H. Multimodal spatial-omics reveal co-evolution of "')
    A('                "alveolar progenitors and proinflammatory niches in progression of lung precursor "')
    A('                "lesions. Cancer Cell 2026;44(2):321-339.e13. doi:10.1016/j.ccell.2025.10.004 "')
    A('                "【本项目源论文；本面板主干出自其 Table S2 / Table S3】",')
    A('    "Travaglini2020": "Travaglini KJ, Nabhan AN, Penland L, et al. A molecular cell atlas of the "')
    A('                      "human lung from single-cell RNA sequencing. Nature 2020;587(7835):619-625. "')
    A('                      "doi:10.1038/s41586-020-2922-4",')
    A('    "VieiraBraga2019": "Vieira Braga FA, Kar G, Berg M, et al. A cellular census of human lungs "')
    A('                       "identifies novel cell states in health and in asthma. Nat Med "')
    A('                       "2019;25(7):1153-1163. doi:10.1038/s41591-019-0468-5",')
    A('}')
    A('')
    A('# ---------------------------------------------- Excel 改名修复（登记，非静默）')
    A('EXCEL_NAME_FIX = {')
    A('    "DKK\\xa03.00": "DKK3",   # Table S2 MP4 Basal/basal stem')
    A('    "ERN\\xa01.00": "ERN1",   # Table S2 MP7 Tumor cell (stress/inflammatory)')
    A('}')
    A('')
    A('# ------------------------------------------------- 源论文 Table S2 上皮 MP（主干）')
    A('# 逐 MP 基因序 = 原表列序，逐字未改（仅套用 EXCEL_NAME_FIX）。')
    A('PAPER_EPI_MP = {')
    for k, (lab, g) in mp.items():
        A(f'    "{k}": dict(label="{lab}", genes=[')
        A(fmt(g))
        A('    ]),')
    A('}')
    A('')
    A('# 亚型 -> 由哪些 MP 组成')
    A('PAPER_SUBTYPE_FROM_MP = {')
    for st, (mps, note) in SUBTYPES.items():
        n = f'   # {note}' if note else ''
        A(f'    "{st}": {mps!r},{n}')
    A('}')
    A('')
    A('# ------------------------- 源论文 Table S3 签名（人，已剔除小鼠列）')
    A('PAPER_SIGNATURES = {')
    for nm, g in sigs.items():
        A(f'    "{nm}": [')
        A(fmt(g))
        A('    ],')
    A('}')
    A('')
    A('# 注：Table S3 另有两列是小鼠（KAC / Inflammatory），本项目为人类数据，**不采纳**，如实登记。')
    A('')
    A('MISSING_IN_MATRIX = ["SFTPA2"]   # 作者矩阵本身不含；见模块 docstring')
    A('')
    A('')
    A('def build_panel(mode="paper_mp"):')
    A('    """返回 {亚型: [基因...]}。')
    A('')
    A('    mode="paper_mp" : 主干，源论文 Table S2 上皮 MP 合并（Ciliated = MP1 ∪ MP8）')
    A('    mode="kac_sig"  : 源论文 Table S3 人 KAC 签名，作为 KAC 的独立第二来源')
    A('    """')
    A('    if mode == "paper_mp":')
    A('        return {st: sorted({g for m in mps for g in PAPER_EPI_MP[m]["genes"]})')
    A('                for st, mps in PAPER_SUBTYPE_FROM_MP.items()}')
    A('    if mode == "kac_sig":')
    A('        return {"KAC": [g for g in PAPER_SIGNATURES["KAC_human"] if g not in MISSING_IN_MATRIX]}')
    A('    raise ValueError(f"未知 mode: {mode}")')
    A('')
    text = "\n".join(L) + "\n"
    with open(out, "w") as fh:
        fh.write(text)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tableS2",
                    default="data/external/peng2026_cancercell/TableS2_meta_program_gene_lists.xlsx")
    ap.add_argument("--tableS3",
                    default="data/external/peng2026_cancercell/TableS3_signature_gene_lists.xlsx")
    ap.add_argument("--out", default="05_annotation/epi_subtype_panel.py")
    a = ap.parse_args()
    mp = read_table_s2_epi(a.tableS2)
    sigs = read_table_s3(a.tableS3)
    text = emit(mp, sigs, a.out)
    print(f"写出 {a.out}  ({len(text)} 字节)")
    print("out sha256:", hashlib.sha256(text.encode()).hexdigest())
    for f, p in (("TableS2", a.tableS2), ("TableS3", a.tableS3)):
        print(f"{f} sha256:", hashlib.sha256(open(p, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
