#!/usr/bin/env bash
# 脱离本地会话的独立队列 —— 2026-09-25 夜
#
# 用户指令（2026-09-25）：「RTCD精细映射和全量CNV都要脱离本地会话，独立运行在我的服务器上」
#
# 本队列只装**当前无阻塞**的两件事：
#   [1/2] 空转 CNV 单患者成本冒烟（P4）—— SPATIAL_CNV_PREREG §9 第 2 项，排期前的硬前置
#   [2/2] RCTD 精细版 caliber d（39 型，全 56 张）
#
# ⛔ **不装**全量 CNV 队列，原因登记在 STATUS 末尾：
#    25 个患者里只有 10 个（P1/P2/P4/P6/P9/P11/P20/P22/P24/P25）自带 Normal 或 AAH 切片，
#    另外 15 个必须按 §11.2「借队列良性端」作 C 主锚，而该借锚在脚本里**尚未实现**
#    （11_run_spatial_cnv_smoke.R 第 ~143 行直接 stop）。一跑就硬停，且 §11.2 自身有一处
#    「观测集相同」的不可行矛盾待签字 ⇒ 不在未签字前占用机时。
#
# 三条规矩（沿用 queue_overnight_20260925.sh）：
#   · **串行**：任何时刻只有一个重活。本机 load 长期 30+，并发只会互相拖慢，且成本冒烟
#     量的是**墙钟**，旁边挂个 4 核 RCTD 会把这个数污染成不可用的。
#   · **可续跑**：每步先看哨兵产物，已有就跳过 —— 中断后重跑本脚本即可继续。
#   · **不吞错**：某步失败就记进 STATUS，继续下一步（两步彼此独立）。
#
# 用法（**必须 setsid**，否则随会话一起被带走）：
#   cp <本脚本> /home/eto/luad_v2/logs/detached_20260925/queue.sh
#   setsid nohup bash /home/eto/luad_v2/logs/detached_20260925/queue.sh \
#       < /dev/null > /home/eto/luad_v2/logs/detached_20260925/queue.log 2>&1 &
#
# 本脚本启动时把 R 源码**冻结**到 logs/detached_20260925/code/ 并记 sha256 ——
# 该目录在 .gitignore 里（logs/），且不在工作树内 ⇒ 与工作树是否被回收无关。

set -u

ROOT=/home/eto/luad_v2
WT="$ROOT/.claude/worktrees/vigilant-hawking-799963"
LOGD="$ROOT/logs/detached_20260925"
CODE="$LOGD/code"
STATUS="$LOGD/STATUS.txt"
RES="$ROOT/results/08_spatial_deconv"
mkdir -p "$LOGD" "$CODE"

## rjags 装载必需；infercnv 在自建库；缺一样都跑不起来
export LD_LIBRARY_PATH="/home/eto/local/jags/lib:${LD_LIBRARY_PATH:-}"
export R_LIBS="/home/eto/Rlibs/infercnv:/home/eto/R/x86_64-pc-linux-gnu-library/4.2:/usr/local/lib/R/site-library:/usr/lib/R/site-library"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$STATUS"; }

echo $$ > "$LOGD/queue.pid"
say "===== detached 队列启动（bash PID $$）====="
say "主机 $(hostname)；cwd $(pwd)；父进程 $(ps -o ppid= -p $$ | tr -d ' ')"
say "本进程 SID=$(ps -o sid= -p $$ | tr -d ' ')、TT=$(ps -o tty= -p $$ | tr -d ' ')（应为自身 SID、TT 为 ?）"

## ---------- 0. 前置检查：缺件就地记录，不猜 ----------
PREFLIGHT_FAIL=0
chk() { if [ -e "$1" ]; then say "  ✔ $2"; else say "  ✘ 缺件：$1（$2）"; PREFLIGHT_FAIL=1; fi; }
say "--- 前置检查 ---"
chk /home/eto/local/jags/lib/libjags.so   "JAGS 动态库（rjags 装载用）"
chk "$RES/spatial_cnv/gene_order_spatial_hg38.tsv" "空转基因位置表（10_build_spatial_gene_order.py 产物）"
chk "$RES/spot_mask.tsv.gz"               "spot 掩码"
chk "$RES/reference_d.h5"                 "RCTD 细版参考"
chk "$RES/reference_d.manifest.json"      "RCTD 细版参考 manifest"
chk "$RES/rctd_a/per_slide"               "RCTD 粗版权重（CNV 门控 argmax 的输入）"
chk "$WT/08_spatial_deconv/01_run_rctd.R" "RCTD 主脚本（工作树）"

## ---------- 1. 冻结源码快照 ----------
if [ "$PREFLIGHT_FAIL" = "0" ]; then
  cp -f "$WT/08_spatial_deconv/01_run_rctd.R"               "$CODE/"
  cp -f "$WT/08_spatial_deconv/calibration.json"            "$CODE/"
  cp -f "$WT/08_spatial_deconv/11_run_spatial_cnv_smoke.R"  "$CODE/"
  # 位置表也冻一份：万一后续重跑 10_ 脚本覆盖它，本次运行的输入指纹仍然可查
  cp -f "$RES/spatial_cnv/gene_order_spatial_hg38.tsv"      "$CODE/"
  say "--- 源码冻结到 $CODE ---"
  ( cd "$CODE" && sha256sum 01_run_rctd.R calibration.json 11_run_spatial_cnv_smoke.R \
      gene_order_spatial_hg38.tsv ) | tee -a "$STATUS"
else
  say "❌ 前置检查未通过 ⇒ 不启动任何计算。"
  exit 2
fi

## 跑一步：失败只记不吞，返回码照传
run_step() {   # $1=名字  $2=日志文件  $3..=命令
  local name=$1 log=$2; shift 2
  say "开始：$name"
  local t0=$SECONDS
  if "$@" >> "$log" 2>&1; then
    say "✅ 完成：$name（用时 $(( (SECONDS - t0) / 60 )) 分钟）"
    return 0
  else
    local rc=$?
    say "❌ 失败：$name（exit $rc，用时 $(( (SECONDS - t0) / 60 )) 分钟）—— 见 $(basename "$log")"
    return $rc
  fi
}

## manifest 里**有** elapsed_min_run 才算跑完（--object-only 那份不算）。
## 🔴 §13.6 修复（2026-09-27）：原来只判「文件存在」，超时后会读到上一轮 --object-only
##    留下的 manifest（不含 elapsed_min_run）⇒ 打印分支照样进 ⇒ KeyError 崩掉失败账目。
##    ⇒ 一律判**键**，不判文件。
has_run() {   # $1 = manifest 路径
  [ -f "$1" ] && python3 -c "
import json,sys
sys.exit(0 if 'elapsed_min_run' in json.load(open('$1')) else 1)
" 2>/dev/null
}
print_manifest() {   # $1 = manifest 路径（调用方须先确认 has_run 为真）
  python3 -c "
import json
d=json.load(open('$1'))
print('      建对象 %.1f min ＋ 跑管线 %.1f min ＝ 总 %.1f min；峰值内存 %.2f GB'
      % (d['elapsed_min_object'], d['elapsed_min_run'], d['elapsed_min_total'], d['peak_gb_final']))
" | tee -a "$STATUS"
}

## ---------- 2. [1/2] 空转 CNV 单患者成本冒烟（P4）----------
M1="$RES/spatial_cnv/smoke_P4/smoke_P4_manifest.json"

if has_run "$M1"; then
  say "[1/2] 成本冒烟：已有完整结果，跳过。"
  print_manifest "$M1"
else
  ## timeout 24h：HMM 段没有先例，给足；超时说明该工具在空转分辨率下不可排期
  run_step "[1/2] 空转 CNV 成本冒烟（患者 P4，全量不抽样）" "$LOGD/step1_cnv_smoke_P4.log" \
      timeout 86400 Rscript "$CODE/11_run_spatial_cnv_smoke.R" --patient P4
  if has_run "$M1"; then print_manifest "$M1"; fi
fi

## ---------- 3. [2/2] RCTD 精细版 caliber d（39 型，全 56 张）----------
## 哨兵 = run_manifest.json（脚本只在 56 张全跑完后才写它）
M2="$RES/rctd_d/run_manifest.json"
if [ -f "$M2" ]; then
  say "[2/2] RCTD caliber d：已有 run_manifest.json，跳过。"
else
  ## max_cores=4 与 rctd_a 那次逐字一致（不是新参数，沿用而已）
  run_step "[2/2] RCTD caliber d（39 型 × 56 张）" "$LOGD/step2_rctd_d.log" \
      Rscript "$CODE/01_run_rctd.R" --caliber d --max-cores 4
fi

## ---------- 4. 阻塞登记：全量 CNV 队列为什么不在这里 ----------
say "--- 未排期的项（不是忘了，是过不了门）---"
say "⛔ 全量空转 CNV（25 患者）未排期，两条独立的阻塞："
say "   ① 借锚未实现：仅 10/25 患者自带 Normal|AAH 切片；其余 15 个按 §11.2 须借队列良性端，"
say "      11_run_spatial_cnv_smoke.R 遇到 n_ref==0 直接 stop ⇒ 跑不了。"
say "   ② §11.2 自身矛盾待签：其「两次运行观测集相同、只有参考集不同」与 infercnv 要求"
say "      参考/观测互斥 直接冲突；且照字面定义良性切片只作参考 ⇒ AAH/AIS/MIA 拿不到恶性判定。"
say "   ⇒ 需先签字定借锚口径，再实现，再排期。成本冒烟（[1/2]）的实测数是排期决策的输入。"
say "===== detached 队列结束（bash PID $$，$(date '+%F %T')）====="
