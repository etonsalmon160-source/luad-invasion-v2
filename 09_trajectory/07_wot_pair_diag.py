#!/usr/bin/env python
"""
单对诊断：WOT 在满规模（Normal 32,251 × AAH 11,891）上为什么会中途无声死掉。

为什么这么设计：
    上一次全量跑在「开跑 WOT」之后 15 秒无声消失 —— 没有 Python 报错、没有 manifest、
    ulimit -c 是 0 所以也没 core。这只有 C 层死法（SIGKILL/SIGSEGV）才解释得通。
    但外部内存闸（5 秒采一次）显示 RSS 没超过 60 GB。

    所以这里**把采样搬进进程内部**，1 秒一次写盘并 flush：
    即使进程被 SIGKILL，最后一行采样也留在文件里，峰值无处可藏。
    同时把每个重活步骤单独计时、单独标记，崩溃点就能定位。

只诊断，不产出任何生物学结论、不落任何耦合文件。
"""
import json
import os
import threading
import time
import traceback

import numpy as np

IN = "/home/eto/luad_v2/results/09_trajectory/wot_full/epiA_wot_input.h5ad"
OUT = "/home/eto/luad_v2/results/09_trajectory/smoke"
SAMPLE_FILE = os.path.join(OUT, "wot_pair_diag_rss.tsv")
MARK_FILE = os.path.join(OUT, "wot_pair_diag_marks.json")

STEPS = []


def read_rss_gb():
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024 / 1024
    except Exception:                              # noqa: BLE001
        pass
    return None


def say(m):
    line = f"[{time.strftime('%H:%M:%S')}] rss={read_rss_gb():.2f}GB {m}"
    print(line, flush=True)
    STEPS.append(line)
    with open(MARK_FILE, "w") as fh:
        json.dump({"steps": STEPS, "last_rss_gb": read_rss_gb()}, fh,
                  ensure_ascii=False, indent=1)


def sampler(stop):
    """1 秒一次写盘并 flush —— 被 SIGKILL 也留痕。"""
    with open(SAMPLE_FILE, "w") as fh:
        fh.write("epoch\tiso\trss_gb\tvmhwm_gb\n")
        fh.flush()
        while not stop.is_set():
            r = read_rss_gb()
            h = None
            try:
                with open("/proc/self/status") as f2:
                    for ln in f2:
                        if ln.startswith("VmHWM:"):
                            h = int(ln.split()[1]) / 1024 / 1024
            except Exception:                        # noqa: BLE001
                pass
            fh.write(f"{time.time():.1f}\t{time.strftime('%H:%M:%S')}\t{r}\t{h}\n")
            fh.flush()
            os.fsync(fh.fileno())
            stop.wait(1.0)


def step(name, fn):
    """跑一步，计时并单独捕获异常（MemoryError 也是 Exception，抓得到）。"""
    t = time.time()
    say(f">>> 开始 {name}")
    try:
        r = fn()
    except Exception as exc:                          # noqa: BLE001
        say(f"!!! {name} 抛异常 {type(exc).__name__}: {str(exc)[:300]}")
        raise
    say(f"<<< {name} 完成 {time.time() - t:.1f}s")
    return r


def main():
    os.makedirs(OUT, exist_ok=True)
    import anndata as ad
    import scipy.sparse as sp
    import wot
    import wot.ot

    stop = threading.Event()
    th = threading.Thread(target=sampler, args=(stop,), daemon=True)
    th.start()
    say("采样线程已起（1 秒一次，写盘 + fsync）")

    a = step("读入预处理后的输入 h5ad", lambda: ad.read_h5ad(IN))
    say(f"    {a.n_obs} × {a.n_vars}，X 稀疏={sp.issparse(a.X)}")

    d = a.obs["day"].astype(float).values
    p0 = a[d == 0.0]
    p1 = a[d == 1.0]
    say(f"    Normal {p0.n_obs} × {p1.n_obs} AAH；p0.X 稀疏={sp.issparse(p0.X)}")

    # 第一步：只在 compute_pca 这一步 —— 源码 util.py:242 会 m.toarray() 稠密化
    res = step("compute_pca（这一步源码里 .toarray() 稠密化）",
               lambda: wot.ot.compute_pca(p0.X, p1.X, 30))
    p0x, p1x, pca, mean = res
    say(f"    p0_x {p0x.shape} {p0x.dtype}；p1_x {p1x.shape}")
    del res, pca, mean
    say("    已释放 PCA 中间产物")

    # 第二步：代价矩阵
    C = step("compute_default_cost_matrix", lambda: wot.ot.OTModel.compute_default_cost_matrix(p0x, p1x, None))
    say(f"    C {C.shape} {C.dtype} = {C.nbytes/1024**3:.2f} GB")
    say(f"    C 统计 min={C.min():.4f} median={np.median(C):.4f} max={C.max():.4f}")
    del C, p0x, p1x
    say("    已释放代价矩阵与 PC 坐标")

    # 第三步：整对求解（走 WOT 自己的路径）
    om = step("OTModel 构造", lambda: wot.ot.OTModel(a, day_field="day"))
    tmap = step("compute_transport_map(0.0 -> 1.0)",
                lambda: om.compute_transport_map(0.0, 1.0))
    say(f"    tmap {tmap.shape}，X 总和 {float(tmap.X.sum()):.1f}")

    stop.set()
    th.join(timeout=3)
    say("诊断跑完")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except BaseException as e:                        # noqa: BLE001
        with open(MARK_FILE, "w") as fh:
            json.dump({"steps": STEPS, "fatal": f"{type(e).__name__}: {e}",
                       "trace": traceback.format_exc().strip().splitlines()[-6:],
                       "last_rss_gb": read_rss_gb()}, fh, ensure_ascii=False, indent=1)
        raise
