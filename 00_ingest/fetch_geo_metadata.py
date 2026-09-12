#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
00_ingest/fetch_geo_metadata.py — 拉取并冻结 GEO 逐样本权威元数据（M0 输入的一部分）

为什么需要
----
本地 per-cell 注释表不含**患者身份**；GEO `!Sample_characteristics_ch1` 才带
`patient id` / `histolgical type` / `tissue origin abbrevation` 等权威字段。
把 GEO 元数据**冻结为本地文件 + 哈希**，下游即可离线、确定性复现，
杜绝"按样本名猜患者号"这类静默假设。

产出：`00_ingest/geo_metadata/<ACC>_gsm.txt`（原始文本）+ `_manifest.json`（URL + SHA-256）。

用法：
    python3 00_ingest/fetch_geo_metadata.py            # 仅在缺失时拉取
    python3 00_ingest/fetch_geo_metadata.py --force    # 强制重拉

注：本机 curl/urllib 的 http(s)_proxy 指向的代理端口未监听（死链），
    故此处**显式禁用代理**直连 NCBI。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import urllib.request
from datetime import datetime, timezone

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(_HERE, "geo_metadata")

ACCESSION = ["GSE131907", "GSE189357", "GSE148071"]

GEO_URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           "?acc={acc}&targ=gsm&form=text&view=brief")


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def fetch(url: str) -> bytes:
    # 显式禁用代理（环境变量中的代理端口未监听）
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    req = urllib.request.Request(url, headers={"User-Agent": "luad_v2-ingest/1.0"})
    with opener.open(req, timeout=60) as resp:
        return resp.read()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    manifest_path = os.path.join(OUT_DIR, "_manifest.json")
    manifest = {}
    if os.path.exists(manifest_path) and not args.force:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)

    for acc in ACCESSION:
        url = GEO_URL.format(acc=acc)
        path = os.path.join(OUT_DIR, f"{acc}_gsm.txt")
        if os.path.exists(path) and not args.force:
            print(f"[skip] {acc} 已存在：{path}")
        else:
            data = fetch(url)
            if not data.startswith(b"^SAMPLE") and b"^SAMPLE" not in data[:5000]:
                raise RuntimeError(f"{acc}: GEO 返回内容异常（非样本元数据）")
            with open(path, "wb") as fh:
                fh.write(data)
            print(f"[fetch] {acc} -> {path} ({len(data)} bytes)")
            manifest[acc] = {
                "url": url,
                "file": os.path.basename(path),
                "sha256": sha256_bytes(data),
                "size_bytes": len(data),
                "fetched_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            }
        # 补齐/刷新哈希
        with open(path, "rb") as fh:
            b = fh.read()
        manifest[acc] = {
            "url": url, "file": os.path.basename(path),
            "sha256": sha256_bytes(b), "size_bytes": len(b),
            "fetched_utc": manifest.get(acc, {}).get(
                "fetched_utc", datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
        }

    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False, sort_keys=True)
    print(f"[manifest] {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
