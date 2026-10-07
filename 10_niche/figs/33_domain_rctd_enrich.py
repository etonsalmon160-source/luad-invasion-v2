#!/usr/bin/env python3
# 33_domain_rctd_enrich.py —— 每个生态位域类型(D1..D7)的 RCTD 组成与富集
import numpy as np, pandas as pd, glob, os
ROOT = "/home/eto/luad_v2"
KH = f"{ROOT}/results/10_niche/kstar_diag/d8_parts"
OUT = f"{ROOT}/results/paper_figures"

map6 = pd.read_csv("/tmp/_rctd2lin.csv")
r2l = dict(zip(map6.rctd, map6.L1))
REF = sorted(map6.rctd.unique())
LIN_ORDER = ["上皮", "成纤维", "髓系", "内皮", "T/NK", "B/浆"]
LIN_EN = {"上皮":"Epithelial","成纤维":"Fibroblast","髓系":"Myeloid",
          "内皮":"Endothelial","T/NK":"T/NK","B/浆":"B/Plasma"}

d7 = pd.read_csv(f"{ROOT}/results/10_niche/kstar_diag/d7_domain_assign.tsv", sep="\t")
D7 = {(s, int(d)): int(a) for s, d, a in zip(d7.slide, d7.domain, d7.archetype)}

rows = []
for f in sorted(glob.glob(f"{KH}/*.tsv")):
    sl = os.path.basename(f)[:-4]
    dm = pd.read_csv(f, sep="\t")
    dm["arch"] = [D7.get((sl, int(d)), np.nan) for d in dm.domain]
    W = pd.read_csv(f"{ROOT}/results/08_spatial_deconv/rctd_d/per_slide/{sl}.weights.tsv.gz",
                    sep="\t", index_col=0)
    j = W.join(dm.set_index("barcode")["arch"], how="inner")
    j["arch"] = j["arch"].astype(int)
    rows.append(j.groupby("arch").mean())

M = pd.concat(rows).groupby(level=0).mean()          # archetype × 39 亚型（mean weight）
N = pd.concat([pd.read_csv(f, sep="\t").assign(sl=os.path.basename(f)[:-4]) for f in
               sorted(glob.glob(f"{KH}/*.tsv"))])
N["arch"] = [D7.get((s, int(d)), np.nan) for s, d in zip(N.sl, N.domain)]
NSPOT = N.groupby("arch").size()

base = M.mean(0)
enr = M.div(base, axis=1)

# 6 谱系汇总
L6 = pd.DataFrame({LIN_EN[L]: M[[r for r in REF if r2l.get(r) == L]].sum(1) for L in LIN_ORDER})
L6b = L6.mean(0)
L6enr = L6.div(L6b, axis=1)

print("=== 每个域类型：spot 数 ===")
for a in NSPOT.index:
    print(f"  D{int(a)}: {NSPOT[a]:7d}  ({NSPOT[a]/NSPOT.sum()*100:4.1f}%)")

print("\n=== 6 谱系平均占比（%）===")
print((L6 * 100).round(1).to_string())

print("\n=== 6 谱系富集倍数（÷ 全队列基线）===")
print(L6enr.round(2).to_string())

print("\n=== 每域 top6 富集 RCTD 亚型 ===")
for a in enr.index:
    top = enr.loc[a].sort_values(ascending=False).head(6)
    print(f"  D{int(a)}: " + " | ".join(f"{k} x{v:.2f}" for k, v in top.items()))

L6.to_csv(f"{OUT}/../10_niche/domain_rctd_lineage_prop.tsv", sep="\t")
L6enr.to_csv(f"{OUT}/../10_niche/domain_rctd_lineage_enrich.tsv", sep="\t")
M.to_csv(f"{OUT}/../10_niche/domain_rctd_subtype_prop.tsv", sep="\t")
enr.to_csv(f"{OUT}/../10_niche/domain_rctd_subtype_enrich.tsv", sep="\t")
pd.DataFrame({"n_spot": NSPOT,
              "pct": (NSPOT / NSPOT.sum() * 100).round(2)}).to_csv(
    f"{OUT}/../10_niche/domain_nspot.tsv", sep="\t")
print("\nsaved: results/10_niche/domain_rctd_{lineage,subtype}_{prop,enrich}.tsv, domain_nspot.tsv")
