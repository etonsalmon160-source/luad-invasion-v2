#!/usr/bin/env python3
# 66_merge_decomposed.py —— 修正版合并：逐细胞型**预算名次表**，再统计全局覆盖
import numpy as np, pandas as pd
OUT="/home/eto/luad_v2/results/10_niche/tr_singlecell"; TOP_N=100
A=pd.read_csv(f"{OUT}/PERCELLTYPE_causal.tsv",sep="\t")
HR=pd.read_csv(f"{OUT}/tcga_pergene_HR.tsv",sep="\t"); hrd=dict(zip(HR.gene,HR.HR))
# 逐型：(基因)->名次（在该型全部扰动里）
RANKS={}
for t,gt in A.groupby("cell_type"):
    gt=gt.sort_values("causal_score",ascending=False).drop_duplicates("perturbed_gene_name")
    RANKS[t]=dict(zip(gt.perturbed_gene_name.astype(str), range(len(gt))))
CT=sorted(RANKS); NCT=len(CT)
print(f"细胞型 {NCT} 个；每型扰动数中位 {int(np.median([len(v) for v in RANKS.values()])):,}")
rows=[]
for (g,s),grp in A.groupby(["perturbed_gene_name","perturbation_sign"]):
    g=str(g); hits=[t for t in CT if RANKS[t].get(g,10**9)<TOP_N]
    if not hits: continue
    rows.append({"gene":g,"sign":s,"intent":"敲低" if s<0 else "过表达",
                 "n_celltypes_top":len(hits),"frac":len(hits)/NCT,
                 "median_causal":grp.causal_score.median(),
                 "median_shift_z":grp.gene_shift_z.median(),
                 "median_rank":int(np.median([RANKS[t].get(g,10**9) for t in CT if RANKS[t].get(g,10**9)<10**9])),
                 "celltypes":",".join(hits[:6])+("…" if len(hits)>6 else "")})
G=pd.DataFrame(rows); G["HR"]=G.gene.map(hrd)
G["qc_keep"]=(((G.intent=="敲低")&(G.HR>1))|((G.intent=="过表达")&(G.HR<1))).fillna(False)
G=G.sort_values(["n_celltypes_top","median_causal"],ascending=False)
G.to_csv(f"{OUT}/GLOBAL_DECOMPOSED_targets_v2.tsv",sep="\t",index=False)
print(f"\n全局榜 {len(G):,} 条；质控后 {int(G.qc_keep.sum()):,}")
print("\nn_celltypes_top 分布:", G.n_celltypes_top.value_counts().sort_index(ascending=False).head(8).to_dict())
print(f"\n=== 在最多细胞型里都进 top{TOP_N} 的（质控过）top 25 ===")
print(G[G.qc_keep].head(25)[["gene","intent","n_celltypes_top","frac","median_rank","median_causal","HR"]].round(3).to_string(index=False))
print(f"\n=== 细胞型特异的（只 1-3 型命中，质控过，按 causal 排）top 12 ===")
sp=G[(G.n_celltypes_top<=3)&G.qc_keep].sort_values("median_causal",ascending=False)
print(sp.head(12)[["gene","intent","n_celltypes_top","median_causal","celltypes","HR"]].round(4).to_string(index=False))
