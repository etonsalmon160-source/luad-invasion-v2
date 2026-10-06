#!/usr/bin/env python3
# 53_target_table.py —— 把所有方法的证据合并成**一张候选靶点表**（交付物）
import pandas as pd, numpy as np, glob
ROOT="/home/eto/luad_v2"; SC=f"{ROOT}/results/10_niche/tr_singlecell"; TRD=f"{ROOT}/results/10_niche/target_reversal"
nz=lambda s: str(s).strip().upper()

# ① 主表：SCMG 系统生物学口径（基因×方向，逐基因经验 p + 跨数据集一致性）
S=pd.read_csv(f"{SC}/scmg_systems_ranked_p.tsv",sep="\t")
S=S.rename(columns={"cos_fib":"cos_state","p":"p_axis","frac_pos":"frac_datasets","n":"n_ds"})
S["gene"]=S.gene.map(nz)
# cos_domain（复合域轴口径）从全量表补回来
AL=pd.read_csv(f"{SC}/scmg_systems_all.tsv",sep="\t"); AL["gene"]=AL.gene.map(nz)
CD=AL.groupby(["gene","sign"]).cos_comp.median().rename("cos_domain").reset_index()
S=S.merge(CD,on=["gene","sign"],how="left")
S["dir_human"]=np.where(S.sign<0,"敲低","过表达")

# ② SCMG CMap 口径（τ）
C=pd.read_csv(f"{TRD}/scL2_A_ECMhi_ranked.tsv",sep="\t").rename(columns={"pert_gene":"gene"})
G=C.groupby(["gene","sign"]).agg(tau_min=("tau","min"), tau_med=("tau","median"),
                                 n_tau=("tau","size")).reset_index()
G["gene"]=G.gene.map(nz); S["gene"]=S.gene.map(nz)
M=S.merge(G,on=["gene","sign"],how="outer")

# ③ 复合域轴（路 1）两个版本
for tag,lab in [("ddom","cos_ddom"),("ddlc","cos_ddlc")]:
    f=f"{SC}/composite_{tag}_ranked.tsv"
    try:
        t=pd.read_csv(f,sep="\t").rename(columns={"cos":lab,"frac":"frac_"+tag})
        t["gene"]=t.gene.map(nz); t=t[["gene","sign",lab,"frac_"+tag]]
        M=M.merge(t,on=["gene","sign"],how="left")
    except Exception as e: print("跳过",tag,e)

# ④ NicheNet（若该基因是配体）
try:
    N=pd.read_csv(f"{SC}/nichenet_ligand_activities.tsv",sep="\t")
    N=N.rename(columns={"test_ligand":"gene","aupr_corrected":"nichenet_aupr"})
    N["gene"]=N.gene.map(nz); M=M.merge(N[["gene","nichenet_aupr"]],on="gene",how="left")
except Exception as e: print("NicheNet 未并入:",e)

# ⑤ 分层
def tier(r):
    ok_axis   = (r.p_axis<0.05) if pd.notna(r.p_axis) else False
    all_ds    = (r.frac_datasets==1.0) if pd.notna(r.frac_datasets) else False
    ok_tau    = (r.tau_min<=-95) if pd.notna(r.tau_min) else False
    both_axes = (pd.notna(r.cos_state) and pd.notna(r.cos_domain) and r.cos_state>0 and r.cos_domain>0)
    if all_ds and ok_axis and ok_tau and both_axes: return "T1_四法一致"
    if all_ds and ok_axis and (ok_tau or both_axes): return "T2_三法支持"
    if all_ds and (ok_axis or ok_tau):               return "T3_两法支持"
    if ok_axis or ok_tau:                            return "T4_单法"
    return "T5_其余"
M["tier"]=M.apply(tier,axis=1)
M=M.sort_values(["tier","cos_state"],ascending=[True,False])
cols=["tier","gene","dir_human","n_ds","cos_state","cos_domain","frac_datasets","p_axis",
      "tau_min","tau_med","cos_ddom","cos_ddlc","nichenet_aupr"]
M[cols].round(4).to_csv(f"{SC}/CANDIDATE_TARGETS.tsv",sep="\t",index=False)

print("=== 分层计数 ===")
print(M.tier.value_counts().sort_index().to_string())
print(f"\n总条目（基因×方向）：{len(M)}")
print("\n=== T1/T2 全部 ===")
top=M[M.tier.isin(["T1_四法一致","T2_三法支持"])]
print(top[cols].head(40).round(4).to_string(index=False))
print(f"\n落盘 {SC}/CANDIDATE_TARGETS.tsv")
