#!/usr/bin/env python3
# 68_clue_submit.py —— 把我们的签名提交到 **CLUE 官方引擎**（clue.io API）
# 格式取自 cluequery 包的 clue_api.R（逐字）：GMT 文本 + data_type/dataset/tool_id/ts_version
import json, os, sys, time, urllib.request, urllib.error
UP  = "/home/eto/luad_v2/results/10_niche/clue_upload/UP_genes.txt"
DN  = "/home/eto/luad_v2/results/10_niche/clue_upload/DOWN_genes.txt"
OUT = "/home/eto/luad_v2/results/10_niche/clue_upload"
KEY = os.environ["CLUE_KEY"]
def log(*a): print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)

def gmt(path, name):
    genes = [l.strip() for l in open(path) if l.strip()]
    return f"{name}\tsignature\t" + "\t".join(genes)

body = {
    "name": "LUAD_v2_global_state_QC",
    "up-cmapfile": gmt(UP, "UP"),
    "down-cmapfile": gmt(DN, "DOWN"),
    "data_type": "L1000",
    "dataset": "Touchstone",
    "tool_id": "sig_gutc_tool",
    "ts_version": "1.0",
}
log(f"提交：UP {len(open(UP).read().split())} 个 / DOWN {len(open(DN).read().split())} 个")
req = urllib.request.Request(
    "https://api.clue.io/api/jobs",
    data=json.dumps(body).encode(),
    headers={"Content-Type": "application/json", "user_key": KEY, "Accept": "application/json"},
    method="POST")
try:
    r = urllib.request.urlopen(req, timeout=120)
    d = json.load(r)
    log(f"HTTP {r.status}")
    print(json.dumps(d, ensure_ascii=False)[:900])
    open(f"{OUT}/clue_job.json", "w").write(json.dumps(d, indent=2))
    jid = d.get("_id") or d.get("id") or (d.get("jobs") or [{}])[0].get("_id")
    log(f"job id = {jid}")
    open(f"{OUT}/clue_job_id.txt", "w").write(str(jid))
except urllib.error.HTTPError as e:
    log(f"HTTP {e.code}"); print(e.read().decode()[:900])
