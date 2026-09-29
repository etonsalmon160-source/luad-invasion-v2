#!/usr/bin/env python3
"""
15_make_external_he_qc.py —— 外部对照切片（GSE248082 N1/N3）的 H&E 质量对照图。

只做「看图」这一件事：把外部两张 H&E 与我们自己队列的两张（P4_Normal / P4_LUAD）
并排放，并把 spot 栅格叠上去（用 tissue_positions × tissue_hires_scalef）。
不产生任何判定、不写任何口径。

坐标约定（见 reference_visium_coordinate_conventions）：
    hires_px = fullres_px × tissue_hires_scalef     ← **乘**，不是除
"""
import json
import os
import sys

from PIL import Image, ImageDraw, ImageFont

ROOT = "/home/eto/luad_v2"
OUTD = os.path.join(ROOT, "results", "08_spatial_deconv", "external_he_qc")
os.makedirs(OUTD, exist_ok=True)

EXT_IMG = os.path.join(ROOT, "data", "external", "GSE248082", "images")
EXT_DIR = os.path.join(ROOT, "data", "external", "GSE248082")   # 坐标文件在这里，不在 images/
OURS = os.path.join(ROOT, "data", "visium_spatial")

PANELS = [
    dict(tag="N1  (external, GSE248082)", kind="ext", gsm="GSM8087031",
         png=os.path.join(EXT_IMG, "GSM8087031_N1_tissue_hires_image.png"),
         pos=os.path.join(EXT_DIR, "GSM8087031_N1_tissue_positions_list.csv.gz"),
         sca=os.path.join(EXT_IMG, "GSM8087031_N1_scalefactors_json.json")),
    dict(tag="N3  (external, GSE248082)", kind="ext", gsm="GSM8087033",
         png=os.path.join(EXT_IMG, "GSM8087033_N3_tissue_hires_image.png"),
         pos=os.path.join(EXT_DIR, "GSM8087033_N3_tissue_positions_list.csv.gz"),
         sca=os.path.join(EXT_IMG, "GSM8087033_N3_scalefactors_json.json")),
    dict(tag="P4_Normal  (ours)", kind="ours", gsm="GSM9226174",
         dir=os.path.join(OURS, "GSM9226174_P4_Normal", "spatial")),
    dict(tag="P4_LUAD  (ours)", kind="ours", gsm="GSM9226177",
         dir=os.path.join(OURS, "GSM9226177_P4_LUAD", "spatial")),
]
for p in PANELS:
    if p["kind"] == "ours":
        p["png"] = os.path.join(p["dir"], "tissue_hires_image.png")
        p["pos"] = os.path.join(p["dir"], "tissue_positions.csv")
        p["sca"] = os.path.join(p["dir"], "scalefactors_json.json")

BOX = 760          # 每格画布边长
CAP = 34           # 标题条高度


def read_positions(p):
    """返回 {barcode: (fullres_row, fullres_col, in_tissue)}；兼容带/不带表头。"""
    op = open
    if p.endswith(".gz"):
        import gzip
        op = gzip.open
    out = {}
    with op(p, "rt") as fh:
        for k, line in enumerate(fh):
            f = line.rstrip("\n").split(",")
            if len(f) != 6:
                continue
            if k == 0 and f[0].lower() == "barcode":
                continue
            try:
                out[f[0]] = (float(f[4]), float(f[5]), int(f[1]))
            except ValueError:
                continue
    return out


def font(sz):
    for cand in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        if os.path.exists(cand):
            return ImageFont.truetype(cand, sz)
    return ImageFont.load_default()


def main():
    print(f"{'panel':38s} {'hires px':>12s} {'spots':>7s} {'in_tissue':>9s} "
          f"{'hires_scalef':>12s} {'画布 µm/px(近似)':>16s}")
    rows = []
    for p in PANELS:
        sca = json.load(open(p["sca"]))
        im = Image.open(p["png"]).convert("RGB")
        pos = read_positions(p["pos"])
        p["_im"], p["_pos"], p["_sca"] = im, pos, sca
        n_in = sum(1 for v in pos.values() if v[2] == 1)
        rows.append(dict(p, w=im.size[0], h=im.size[1], n=len(pos), n_in=n_in,
                         sc=sca["tissue_hires_scalef"]))
        print(f"{p['tag']:38s} {im.size[0]:>6d}x{im.size[1]:<5d} {len(pos):>7d} "
              f"{n_in:>9d} {sca['tissue_hires_scalef']:>12.5f}")

    # ---- 拼图：2 × 2，每格 H&E + spot 栅格 ----
    W = 2 * BOX
    H = 2 * (BOX + CAP)
    canvas = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(canvas)
    f_cap, f_leg = font(20), font(15)

    for i, r in enumerate(rows):
        cx, cy = (i % 2) * BOX, (i // 2) * (BOX + CAP)
        im, sc = r["_im"], r["sc"]
        # 等比缩放到 BOX 见方
        s = BOX / max(im.size)
        im2 = im.resize((max(1, int(im.size[0] * s)), max(1, int(im.size[1] * s))),
                        Image.LANCZOS)
        ox = cx + (BOX - im2.size[0]) // 2
        oy = cy + CAP + (BOX - im2.size[1]) // 2
        canvas.paste(im2, (ox, oy))

        # spot 栅格：hires = fullres × scalef
        d2 = ImageDraw.Draw(canvas)
        rad = max(1.6, 0.5 * 238.8 * sc * s) if r["kind"] == "ours" else \
              max(1.6, 0.5 * 91.63 * sc * s)
        hit = 0
        for bc, (fr, fc, it) in r["_pos"].items():
            x = ox + fc * sc * s
            y = oy + fr * sc * s
            if not (ox <= x <= ox + im2.size[0] and oy <= y <= oy + im2.size[1]):
                continue
            hit += 1
            col = (0, 170, 0) if it == 1 else (230, 60, 60)
            d2.ellipse([x - rad, y - rad, x + rad, y + rad], outline=col,
                       width=1 if it != 1 else 2)

        dr.rectangle([cx, cy, cx + BOX, cy + CAP], fill=(28, 28, 32))
        dr.text((cx + 10, cy + 7), r["tag"], fill="white", font=f_cap)
        dr.text((cx + BOX - 250, cy + 9),
                f"{r['w']}x{r['h']} px  ·  n={r['n_in']}", fill=(200, 200, 200),
                font=f_leg)
        dr.text((ox + 6, oy + 6), f"drawn {hit}", fill=(255, 255, 0), font=f_leg)
        print(f"  {r['tag']:38s} 重画 spot {hit}/{r['n']}")

    dr.rectangle([0, 0, W - 1, H - 1], outline=(120, 120, 120), width=2)
    out = os.path.join(OUTD, "external_vs_ours_he_grid.png")
    canvas.save(out, optimize=True)
    print(f"\n图：{out}  ({os.path.getsize(out)/1e6:.2f} MB)  {W}x{H}")
    print("绿圈 = in_tissue==1，红圈 = in_tissue==0；半径按各自 spot_diameter_fullres×scalef")

    # ---- 另一张：仅外部两张的近距离对照（不叠栅格，看染色质量）----
    W2, H2 = 2 * 900, 900 + CAP
    c2 = Image.new("RGB", (W2, H2), "white")
    d2 = ImageDraw.Draw(c2)
    for i, r in enumerate(rows[:2]):
        im = r["_im"]
        s = 900 / max(im.size)
        im2 = im.resize((int(im.size[0] * s), int(im.size[1] * s)), Image.LANCZOS)
        cx = i * 900
        d2.rectangle([cx, 0, cx + 900, CAP], fill=(28, 28, 32))
        d2.text((cx + 10, 7), r["tag"] + f"   hires {r['w']}x{r['h']}",
                fill="white", font=f_cap)
        c2.paste(im2, (cx + (900 - im2.size[0]) // 2, CAP + (900 - im2.size[1]) // 2))
    out2 = os.path.join(OUTD, "external_he_only.png")
    c2.save(out2, optimize=True)
    print(f"图：{out2}  ({os.path.getsize(out2)/1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
