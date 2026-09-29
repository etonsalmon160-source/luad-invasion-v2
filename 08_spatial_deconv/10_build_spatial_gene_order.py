#!/usr/bin/env python3
"""空转（Visium）版基因位置表 —— SPATIAL_CNV_PREREG.md §2.3 要求**重跑、重算哈希**。

为什么不沿用 `results/03_cnv/infercnv_smoke/gene_order_hg38.tsv`：
  那份是给 **snRNA**（18,069 基因）建的，空转是 **18,085 行 / 18,082 唯一 symbol**，
  基因集合不同 ⇒ 位置表内容与哈希都不同。§2.3 明写「不得直接沿用旧哈希」。
  位置解析逻辑与 `03_cnv/14_build_gene_order.py` **逐字相同**（同一 GENCODE v44 + 同一 HGNC 表），
  只换「我们的基因列表从哪来」这一步。

空转的一个额外坑（§2.3）：features 里 `HSPA14` / `TBCE` / `TMSB15B` 各是**双探针**，
两行同名。位置表按 **symbol 去重**（取一次），而**表达矩阵那边必须先把两行按 symbol 求和**
（在 `<NN>_build_spatial_cnv_input.py` 里做）。若位置表按行写、矩阵也按行写却不对齐，
inferCNV 会按名字 join，把这两个探针之一静默丢掉。

输入：data/visium_spatial/*/filtered_feature_bc_matrix/features.tsv.gz（全 56 张，核对只有一套基因集）
      data/external/gencode/gencode.v44.annotation.gtf.gz
      data/external/hgnc/hgnc_complete_set.txt
输出：results/08_spatial_deconv/spatial_cnv/gene_order_spatial_hg38.tsv（4 列无表头：gene chr start end）
      results/08_spatial_deconv/spatial_cnv/gene_order_spatial_manifest.json
      results/08_spatial_deconv/spatial_cnv/gene_order_spatial_unmatched.txt

法则 0：匹配不上的基因**逐个列出**，绝不静默丢弃。
"""
import gzip
import glob
import hashlib
import json
import os
import sys
from collections import defaultdict

ROOT = "/home/eto/luad_v2"
VISIUM = f"{ROOT}/data/visium_spatial"
GTF = f"{ROOT}/data/external/gencode/gencode.v44.annotation.gtf.gz"
HGNC = f"{ROOT}/data/external/hgnc/hgnc_complete_set.txt"
OUTDIR = f"{ROOT}/results/08_spatial_deconv/spatial_cnv"
OUT = f"{OUTDIR}/gene_order_spatial_hg38.tsv"
UNMATCHED = f"{OUTDIR}/gene_order_spatial_unmatched.txt"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# SC0 底座核对用的期望值（SPATIAL_CNV_PREREG.md §4）
EXP_GE_ROWS = 18085
EXP_UNIQUE_SYMBOLS = 18082
EXP_AB_ROWS = 35


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def collect_visium_genes():
    """全 56 张的 features.tsv.gz → 唯一 Gene Expression symbol 集合；顺带做 SC0 的行数核对。"""
    files = sorted(glob.glob(f"{VISIUM}/*/filtered_feature_bc_matrix/features.tsv.gz"))
    assert len(files) == 56, f"期望 56 张切片，实得 {len(files)}"
    sets, n_ge_seen, n_ab_seen = [], set(), defaultdict(int)
    for f in files:
        slide = f.split("/")[-3]
        ge_rows = 0
        ab_rows = 0
        syms = []
        with gzip.open(f, "rt") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                if p[2] == "Gene Expression":
                    ge_rows += 1
                    syms.append(p[1])
                elif p[2] == "Antibody Capture":
                    ab_rows += 1
        # 🔴 逐张硬核对：GE 行数、抗体行数（有则必须恰为 35）
        assert ge_rows == EXP_GE_ROWS, f"{slide}: GE 行 {ge_rows} != {EXP_GE_ROWS}"
        assert ab_rows in (0, EXP_AB_ROWS), f"{slide}: 抗体行 {ab_rows} 既不为 0 也不为 {EXP_AB_ROWS}"
        n_ge_seen.add(ge_rows)
        n_ab_seen[ab_rows] += 1
        sets.append((slide, frozenset(syms)))
    uniq_sets = set(s for _, s in sets)
    assert len(uniq_sets) == 1, f"56 张的基因集合并不相同（{len(uniq_sets)} 种）—— 停下"
    genes = set(next(iter(uniq_sets)))
    assert len(genes) == EXP_UNIQUE_SYMBOLS, \
        f"唯一 symbol {len(genes)} != {EXP_UNIQUE_SYMBOLS}（双探针未按 symbol 折叠？）"
    print(f"[visium] 56 张全部核对通过：GE {sorted(n_ge_seen)} 行、"
          f"抗体分布 {dict(n_ab_seen)}、唯一 symbol {len(genes)}")
    return sorted(genes)


def parse_gtf_gene_lines(path):
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

    genes = collect_visium_genes()
    print(f"[in] 读 GENCODE {GTF}")
    sym2recs = parse_gtf_gene_lines(GTF)
    print(f"[in] 读 HGNC {HGNC}")
    prev_map, alias_map = load_hgnc_maps(HGNC)

    name_rows, unmatched, multi, renamed = [], [], [], []
    n_by_method = defaultdict(int)
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
        # 第一列写**我们数据里的名字**：inferCNV 按名字 join，写成 GENCODE 当前名会让改名基因静默掉
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
    if unmatched:
        print(f"[match] 仍未匹配示例（前 20）：{unmatched[:20]}")

    with open(UNMATCHED, "w") as f:
        if unmatched:
            f.write("\n".join(unmatched) + "\n")

    off = [r for r in name_rows if r[1] not in CHR_ORDER]
    if off:
        print(f"[warn] {len(off)} 个基因落在规范染色体之外：{sorted(set(r[1] for r in off))}")

    name_rows.sort(key=lambda r: (CHR_ORDER.index(r[1]) if r[1] in CHR_ORDER else len(CHR_ORDER),
                                  r[2]))
    with open(OUT, "w") as f:
        for g, chrom, start, end in name_rows:
            f.write(f"{g}\t{chrom}\t{start}\t{end}\n")
    print(f"[out] {OUT}  {len(name_rows)} 行")

    # SC0：匹配率 ≥ 99%
    rate = 100.0 * len(name_rows) / n_total
    assert rate >= 99.0, f"SC0 不过：匹配率 {rate:.2f}% < 99%"

    manifest = {
        "script": "08_spatial_deconv/10_build_spatial_gene_order.py",
        "why_not_reuse_snrna_table": "SPATIAL_CNV_PREREG.md §2.3：空转基因集与 snRNA 不同，须重跑重算哈希",
        "gencode_gtf": GTF,
        "gencode_gtf_sha256": sha256(GTF),
        "gencode_version": "v44 (hg38)",
        "hgnc_table": HGNC,
        "hgnc_table_sha256": sha256(HGNC),
        "source_of_our_genes": f"{VISIUM}/*/filtered_feature_bc_matrix/features.tsv.gz（56 张，基因集唯一）",
        "gene_calling_rule": "Gene Expression 行；抗体行 35 已剔；双探针基因按 symbol 折叠取一次",
        "n_ge_rows_per_slide": EXP_GE_ROWS,
        "n_ab_rows": EXP_AB_ROWS,
        "n_slides": 56,
        "n_genes_ours": n_total,
        "n_matched": len(name_rows),
        "n_unmatched": len(unmatched),
        "match_rate_pct": round(rate, 2),
        "n_exact": n_by_method["exact"],
        "n_via_prev_symbol": n_by_method["prev_symbol"],
        "n_via_alias_symbol": n_by_method["alias_symbol"],
        "n_symbol_renamed": len(renamed),
        "renamed_sample": renamed[:50],
        "n_multi_candidate": len(multi),
        "unmatched_list": UNMATCHED,
        "out": OUT,
        "out_sha256": sha256(OUT),
        "chr_order": CHR_ORDER,
        "note": "第一列写我们数据里的基因名（可能是旧名）；坐标取自 GENCODE v44，经 HGNC "
                "prev_symbol/alias_symbol 解析改名。⚠️ 本哈希与 snRNA 那份**不同**是预期的。",
    }
    with open(f"{OUTDIR}/gene_order_spatial_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print(json.dumps({k: v for k, v in manifest.items() if k != "renamed_sample"},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    sys.exit(main())
