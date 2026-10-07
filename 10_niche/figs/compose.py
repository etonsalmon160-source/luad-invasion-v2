#!/usr/bin/env python3
# compose.py —— 把已有面板图拼成正文大图（统一版式 + 小写字母编号）
#
# 学术版式约定：
#   · 面板编号用小写加粗字母 a, b, c…（Nature/Cell 现行标准），置于面板左上角
#   · 全图统一字体（DejaVu Sans）、统一字号
#   · 面板间留白按"紧但不粘连"取值；整图外侧留 0.2 in
#   · 输出 PDF（矢量容器）+ 400 dpi PNG
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread

SRC = "/home/eto/luad_v2/results/paper_figures"
OUT = f"{SRC}/composite"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,          # TrueType 嵌入，投稿要求
    "ps.fonttype": 42,
})

LET = "abcdefghijklmnopqrstuvwxyz"


def compose(rows, fname, width=7.2, title=None, hspace=0.06, wspace=0.05,
            letters=True, row_gap=0.055, panel_pad=0.0):
    """rows: list[list[(png_name, rel_width)]]，每个内层列表是一行面板。"""
    ims = [[(n, imread(f"{SRC}/{n}.png"), w) for n, w in r] for r in rows]
    # 各行高度 = 该行最高面板高度 / 总宽 × 目标宽
    hs, nr = [], len(ims)
    for r in ims:
        tot_w = sum(w for _, _, w in r)
        arr = [im.shape[0] / im.shape[1] * w for _, im, w in r]
        hs.append(max(arr) * width / tot_w)
    H = sum(hs) + row_gap * width * (nr - 1) + 0.42
    fig = plt.figure(figsize=(width + 0.02, H))
    fig.patch.set_facecolor("white")
    y = 1.0
    for ri, (r, rh) in enumerate(zip(ims, hs)):
        frac_h = rh / H
        tot_w = sum(w for _, _, w in r)
        x = 0.0
        for ci, (name, im, w) in enumerate(r):
            fw = (w / tot_w)
            ax = fig.add_axes([x, y - frac_h, fw, frac_h])
            ax.imshow(im); ax.axis("off")
            if letters:
                k = sum(len(rr) for rr in ims[:ri]) + ci
                ax.text(-0.012, 1.012, LET[k], transform=ax.transAxes, fontsize=12,
                        fontweight="bold", va="bottom", ha="left", color="black")
            x += fw
        y -= frac_h + row_gap / H * (1 if ri < nr - 1 else 0)
    if title:
        fig.text(0.0, 1.0, title, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.subplots_adjust(top=0.99)
    fig.savefig(f"{OUT}/{fname}.pdf", bbox_inches="tight", pad_inches=0.14,
                facecolor="white")
    fig.savefig(f"{OUT}/{fname}.png", dpi=400, bbox_inches="tight", pad_inches=0.14,
                facecolor="white")
    plt.close(fig)
    print(f"  ✓ {fname}   ({nr} 行, {sum(len(r) for r in rows)} 面板)")
