#!/usr/bin/env python3
"""
14_fetch_geo_tar_members.py —— 从 GEO 的**未压缩** RAW.tar 里按需取单个成员。

为什么不用直接下整个 tar：GSE248082_RAW.tar 是 1.4 GB，而我们只要 N1/N3 的
H&E 图（tissue_hires / tissue_lowres / aligned_fiducials / detected_tissue_image），
合计约 31 MB。NCBI 的 FTP 支持 HTTP Range（Accept-Ranges: bytes），
且该 tar 是 **stored**（未压缩）的 ⇒ 可以先读 512 字节的 tar 头拿到成员的
名字与长度，再按偏移量单独把该成员抓下来。

不做的事：不写任何判定、不碰任何口径；这只是取文件。
"""
import gzip
import io
import json
import os
import sys
import time
import urllib.request

URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248082/"
       "suppl/GSE248082_RAW.tar")
OUTDIR = "/home/eto/luad_v2/data/external/GSE248082/images"

WANT = [
    "GSM8087031_N1_tissue_hires_image.png.gz",
    "GSM8087031_N1_tissue_lowres_image.png.gz",
    "GSM8087031_N1_aligned_fiducials.jpg.gz",
    "GSM8087031_N1_detected_tissue_image.jpg.gz",
    "GSM8087031_N1_scalefactors_json.json.gz",
    "GSM8087033_N3_tissue_hires_image.png.gz",
    "GSM8087033_N3_tissue_lowres_image.png.gz",
    "GSM8087033_N3_aligned_fiducials.jpg.gz",
    "GSM8087033_N3_detected_tissue_image.jpg.gz",
    "GSM8087033_N3_scalefactors_json.json.gz",
]

BLK = 512


def range_get(start, length, tries=5):
    """取 [start, start+length) 的字节；失败重试。"""
    end = start + length - 1
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                b = r.read()
            if len(b) == length:
                return b
            last = f"短读 {len(b)}/{length}"
        except Exception as e:  # noqa: BLE001
            last = repr(e)
        time.sleep(2.0 * (k + 1))
    raise RuntimeError(f"range {start}+{length} 失败：{last}")


def parse_header(buf):
    """返回 (name, size)；全零块 ⇒ (None, 0) 表示归档结束。"""
    if buf == b"\0" * BLK:
        return None, 0
    name = buf[0:100].split(b"\0", 1)[0].decode("utf-8", "replace")
    try:
        size = int(buf[124:136].split(b"\0", 1)[0].strip() or b"0", 8)
    except ValueError:
        size = 0
    prefix = buf[345:500].split(b"\0", 1)[0].decode("utf-8", "replace")
    if prefix:                       # GNU 长名扩展
        name = prefix + "/" + name
    return name, size


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    todo = set(WANT)
    found = {}
    off = 0
    n_hdr = 0
    t0 = time.time()

    while todo:
        buf = range_get(off, BLK)
        name, size = parse_header(buf)
        if name is None:
            print(f"[end] 走到归档末尾（{n_hdr} 个头），仍未找到：{sorted(todo)}")
            break
        n_hdr += 1
        # 类型标志：'0' / '\0' = 普通文件
        typeflag = buf[156:157]
        if name in todo and typeflag in (b"0", b"\0"):
            print(f"[取] 偏移 {off:>12,}  {name}  ({size:,} B, "
                  f"{size/1e6:.1f} MB) …", flush=True)
            data = range_get(off + BLK, size)
            dest = os.path.join(OUTDIR, name)
            with open(dest, "wb") as fh:
                fh.write(data)
            # 立刻校验能不能解压（坏文件早发现）
            try:
                if name.endswith(".gz"):
                    with gzip.open(dest, "rb") as g:
                        head = g.read(16)
                    print(f"      gz 可解，前 16 字节 {head[:8].hex()}")
                found[name] = size
                todo.discard(name)
            except Exception as e:  # noqa: BLE001
                print(f"      ⚠️ gz 解压失败：{e!r}；删除该文件")
                os.remove(dest)
        off += BLK + ((size + BLK - 1) // BLK) * BLK
        if n_hdr % 25 == 0:
            print(f"  … 扫过 {n_hdr} 个成员，偏移 {off:,} / 1,479,075,840",
                  flush=True)

    meta = {"url": URL, "outdir": OUTDIR,
            "requested": WANT, "got": found,
            "missing": sorted(todo),
            "members_scanned": n_hdr,
            "bytes_scanned_offsets": off,
            "elapsed_sec": round(time.time() - t0, 1)}
    with open(os.path.join(OUTDIR, "fetch_manifest.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"\n用时 {meta['elapsed_sec']} s；拿到 {len(found)}/{len(WANT)} 个")
    print(f"清单：{os.path.join(OUTDIR, 'fetch_manifest.json')}")
    return 0 if not todo else 1


if __name__ == "__main__":
    sys.exit(main())
