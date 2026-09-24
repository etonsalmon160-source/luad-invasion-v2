#!/usr/bin/env python3
"""从权威 GENCODE hg38 注释建 inferCNV 要的基因位置表。

输入：data/external/gencode/gencode.v44.annotation.gtf.gz（GENCODE v44, hg38, EBI 镜像）
      data/external/hgnc/hgnc_complete_set.txt（HGNC 官方命名表，用于解析旧名/别名）
      results/02_expression/gse308103_counts_paperqc.h5ad 的 var_names（我们的 18,069 个基因）
输出：results/03_cnv/infercnv_smoke/gene_order_hg38.tsv      （4 列无表头：gene chr start end）
      results/03_cnv/infercnv_smoke/gene_order_unmatched_genes.txt

我们的基因名来自较老的 cellranger 参考，GENCODE v44 已给一批基因改名
（YARS→YARS1、SARS→SARS1、KIARA1324→ELAPOR1 等）。故走 HGNC 官方 prev_symbol/alias_symbol
解析，而**不是**自己编映射。

法则 0：**最终仍匹配不上的基因必须逐个列出**，绝不静默丢弃
（CopyKAT 在内嵌表上静默丢了 724 个）。
"""
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict

import anndata as ad

ROOT = "/home/eto/luad_v2"
GTF = f"{ROOT}/data/external/gencode/gencode.v44.annotation.gtf.gz"
HGNC = f"{ROOT}/data/external/hgnc/hgnc_complete_set.txt"
H5AD = f"{ROOT}/results/02_expression/gse308103_counts_paperqc.h5ad"
OUTDIR = f"{ROOT}/results/03_cnv/infercnv_smoke"
OUT = f"{OUTDIR}/gene_order_hg38.tsv"
UNMATCHED = f"{OUTDIR}/gene_order_unmatched_genes.txt"

# 主染色体的规范顺序（inferCNV 按此排序画图）
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def parse_gtf_gene_lines(path):
    """只取 feature == gene 的行，返回 gene_name -> [(chr, start, end, gene_id), ...]"""
    sym2recs = defaultdict(list)
    n_lines = n_gene = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            n_lines += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            n_gene += 1
            attrs = {}
            for kv in parts[8].split(";"):
                kv = kv.strip()
                if not kv or " " not in kv:
                    continue
                k, v = kv.split(" ", 1)
                attrs[k] = v.strip('"')
            name = attrs.get("gene_name")
            if name:
                sym2recs[name].append((parts[0], int(parts[3]), int(parts[4]),
                                       attrs.get("gene_id", "")))
    print(f"[gtf] 扫描 {n_lines} 行，其中 gene 行 {n_gene}，唯一 gene_name {len(sym2recs)}")
    return sym2recs


def load_hgnc_maps(path):
    """HGNC 官方表 → prev_symbol/alias_symbol 到当前 symbol 的映射（只保留无歧义的）。"""
    prev_map, alias_map = defaultdict(set), defaultdict(set)
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        i_sym = header.index("symbol")
        i_prev = header.index("prev_symbol")
        i_alias = header.index("alias_symbol")
        n = 0
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(i_sym, i_prev, i_alias):
                continue
            n += 1
            sym = p[i_sym].strip().strip('"')
            for field, d in ((p[i_prev], prev_map), (p[i_alias], alias_map)):
                for v in field.split("|"):
                    v = v.strip().strip('"')
                    if v:
                        d[v].add(sym)
    print(f"[hgnc] 读 {n} 行；prev_symbol 键 {len(prev_map)}，alias_symbol 键 {len(alias_map)}")
    return prev_map, alias_map


def pick_record(recs):
    """同名基因取哪一条：优先非 PAR_Y、优先主染色体、再取最长 span。"""
    def key(r):
        chrom, start, end, gid = r
        is_par = 1 if gid.endswith("_PAR_Y") else 0
        chrom_rank = CHR_ORDER.index(chrom) if chrom in CHR_ORDER else len(CHR_ORDER)
        return (is_par, chrom_rank, -(end - start))
    return sorted(recs, key=key)[0]


def resolve(g, sym2recs, prev_map, alias_map):
    """返回 (当前 symbol 或 None, 解析方式)。先精确名，再官方 prev_symbol，再官方 alias_symbol。"""
    if g in sym2recs:
        return g, "exact"
    for m, tag in ((prev_map, "prev_symbol"), (alias_map, "alias_symbol")):
        cands = m.get(g)
        if cands and len(cands) == 1:
            c = next(iter(cands))
            if c in sym2recs:
                return c, tag
    return None, "unmatched"


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print(f"[in] 读我们的基因列表 {H5AD}")
    a = ad.read_h5ad(H5AD, backed="r")
    genes = list(a.var_names)
    print(f"[in] {len(genes)} 个基因")
    assert len(genes) == 18069, f"基因数 {len(genes)} != 18069，停下"
    assert len(set(genes)) == len(genes), "我们自己的基因列表有重复，停下"

    print(f"[in] 读 GENCODE {GTF}")
    sym2recs = parse_gtf_gene_lines(GTF)
    print(f"[in] 读 HGNC {HGNC}")
    prev_map, alias_map = load_hgnc_maps(HGNC)

    name_rows, unmatched, multi = [], [], []
    n_by_method = defaultdict(int)
    renamed = []  # (我们的旧名, 解析到的当前名, 方式)
    for g in genes:
        cur, how = resolve(g, sym2recs, prev_map, alias_map)
        n_by_method[how] += 1
        if cur is None:
            unmatched.append(g)
            continue
        if cur != g:
            renamed.append((g, cur, how))
        recs = sym2recs[cur]
        if len(recs) > 1:
            multi.append((cur, len(recs)))
        chrom, start, end, _gid = pick_record(recs)
        # 🔴 位置表第一列必须写**我们的**基因名（旧名）：inferCNV 是按基因名 join 的，
        # 写成 GENCODE 当前名会让 301 个改过名的基因**静默掉**。
        name_rows.append((g, chrom, start, end))

    n_total = len(genes)
    print(f"[match] 精确名 {n_by_method['exact']}")
    print(f"[match] 经 HGNC prev_symbol 解析 {n_by_method['prev_symbol']}")
    print(f"[match] 经 HGNC alias_symbol 解析 {n_by_method['alias_symbol']}")
    print(f"[match] 合计匹配 {len(name_rows)} / {n_total} "
          f"({100.0*len(name_rows)/n_total:.2f}%)")
    print(f"[match] 仍**未匹配** {len(unmatched)}")
    print(f"[match] 多候选记录（已按规则择一）{len(multi)}")

    if renamed:
        print("[match] 被 GENCODE 改名、经 HGNC 解析回来的（前 15）：")
        for old, new, how in renamed[:15]:
            print(f"         {old} -> {new}  ({how})")

    if multi:
        print("[match] 多候选示例（前 5）：")
        for g, n in multi[:5]:
            print(f"         {g}: {n} 条")

    with open(UNMATCHED, "w") as f:
        if unmatched:
            f.write("\n".join(unmatched) + "\n")
    if unmatched:
        print(f"[match] 仍未匹配示例（前 20）：{unmatched[:20]}")
    print(f"[out] 未匹配清单 → {UNMATCHED}")

    off = [r for r in name_rows if r[1] not in CHR_ORDER]
    if off:
        print(f"[warn] {len(off)} 个基因落在规范染色体之外：{sorted(set(r[1] for r in off))}")

    name_rows.sort(key=lambda r: (CHR_ORDER.index(r[1]) if r[1] in CHR_ORDER else len(CHR_ORDER),
                                  r[2]))
    with open(OUT, "w") as f:
        for g, chrom, start, end in name_rows:
            f.write(f"{g}\t{chrom}\t{start}\t{end}\n")
    print(f"[out] {OUT}  {len(name_rows)} 行")

    manifest = {
        "script": "03_cnv/14_build_gene_order.py",
        "gencode_gtf": GTF,
        "gencode_gtf_sha256": sha256(GTF),
        "gencode_version": "v44 (hg38)",
        "hgnc_table": HGNC,
        "hgnc_table_sha256": sha256(HGNC),
        "source_of_our_genes": H5AD,
        "n_genes_ours": n_total,
        "n_matched": len(name_rows),
        "n_unmatched": len(unmatched),
        "match_rate_pct": round(100.0 * len(name_rows) / n_total, 2),
        "n_exact": n_by_method["exact"],
        "n_via_prev_symbol": n_by_method["prev_symbol"],
        "n_via_alias_symbol": n_by_method["alias_symbol"],
        "n_symbol_renamed": len(renamed),
        "renamed_sample": renamed[:50],
        "n_multi_candidate": len(multi),
        "out": OUT,
        "out_sha256": sha256(OUT),
        "unmatched_list": UNMATCHED,
        "chr_order": CHR_ORDER,
        "note": "位置表第一列写**我们数据里的基因名**（可能是旧名）；坐标取自 GENCODE v44，"
                "经 HGNC prev_symbol/alias_symbol 解析改名，映射见 renamed_sample",
    }
    with open(f"{OUTDIR}/gene_order_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in manifest.items() if k != "renamed_sample"},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
