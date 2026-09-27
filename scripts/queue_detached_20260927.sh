#!/usr/bin/env bash
# 脱离本地会话的独立队列 —— 2026-09-27
#
# 用户指令（2026-09-27）：「去做」—— 对预注册 §14 三条裁定 ＋ 单患者成本重测的授权。
#
# 本队列只装**一件事**：
#   [1/1] 空转 CNV 单患者成本冒烟（P4）—— SPATIAL_CNV_PREREG §9 第 2 项
#
# 相对 2026-09-25 那轮（scripts/queue_detached_20260925.sh）的四处改动，
# **全部是已签字的**，没有一处是新口径：
#   · analysis_mode  subclusters → samples   （§13.1；上轮 246 个单位、24 h 未跑完）
#   · no_plot        FALSE      → TRUE       （§13.1；上轮 9.5 h 纯绘图）
#   · chr_exclude    chrX/chrY/chrM → chrM   （§14.3；原文 Step1/Step3 两处逐字）
#   · 失败上报 bug 已修                       （§13.6；上轮超时后 KeyError 崩掉失败账目）
#
# 三条规矩沿用上轮：**串行**（本机 load 长期 30+；成本冒烟量的是墙钟，旁边挂重活会污染）；
# **可续跑**（先看哨兵，已有就跳过）；**不吞错**（失败记进 STATUS）。
#
# 用法（**必须 setsid**）：
#   cp <本脚本> /home/eto/luad_v2/logs/detached_20260927/queue.sh
#   setsid nohup bash /home/eto/luad_v2/logs/detached_20260927/queue.sh \
#       < /dev/null > /home/eto/luad_v2/logs/detached_20260927/queue.log 2>&1 &

set -u

ROOT=/home/eto/luad_v2
WT="$ROOT/.claude/worktrees/vigilant-hawking-799963"
LOGD="$ROOT/logs/detached_20260927"
CODE="$LOGD/code"
STATUS="$LOGD/STATUS.txt"
RES="$ROOT/results/08_spatial_deconv"
M1="$RES/spatial_cnv/smoke_P4/smoke_P4_manifest.json"
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
chk "$RES/rctd_a/per_slide"               "RCTD 粗版权重（CNV 门控 argmax 的输入）"
chk "$WT/08_spatial_deconv/11_run_spatial_cnv_smoke.R" "冒烟脚本（工作树）"
chk "$WT/08_spatial_deconv/SPATIAL_CNV_PREREG.md"      "预注册（§14 签字的落盘）"

## ---------- 1. 冻结源码快照 ----------
if [ "$PREFLIGHT_FAIL" = "0" ]; then
  cp -f "$WT/08_spatial_deconv/11_run_spatial_cnv_smoke.R" "$CODE/"
  cp -f "$WT/08_spatial_deconv/SPATIAL_CNV_PREREG.md"      "$CODE/"
  ## 位置表也冻一份：万一后续重跑 10_ 脚本覆盖它，本次运行的输入指纹仍然可查
  cp -f "$RES/spatial_cnv/gene_order_spatial_hg38.tsv"     "$CODE/"
  say "--- 源码冻结到 $CODE ---"
  ( cd "$CODE" && sha256sum 11_run_spatial_cnv_smoke.R SPATIAL_CNV_PREREG.md \
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

## 🔴 §13.6 修复（2026-09-27）：原来只判「manifest 文件存在」，而超时后会读到上一轮
##    --object-only 留下的那份（不含 elapsed_min_run）⇒ 打印分支照样进 ⇒ KeyError 崩掉
##    失败账目。⇒ 一律判**键**，不判文件。
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
print('      进管线基因 %d；chr_exclude %s；analysis_mode %s'
      % (d['gene_account']['n_into_infercnv'], d['params_signed_later']['chr_exclude'],
         d['params_signed_later']['analysis_mode']))
" | tee -a "$STATUS"
}

## ---------- 2. [1/1] 空转 CNV 单患者成本冒烟（P4）----------
## 哨兵 = 完整跑的 manifest（含 elapsed_min_run）。只跑过 --object-only 的那份不算数。
if has_run "$M1"; then
  say "[1/1] 成本冒烟：已有完整结果，跳过。"
  print_manifest "$M1"
else
  ## timeout 24 h 只作上界：上轮 subclusters 下 24 h 未完成；本轮 samples 后单位数 246→2，
  ## 预期远低于此。真超时即说明该工具在空转分辨率下不可排期，如实上报。
  run_step "[1/1] 空转 CNV 成本冒烟（患者 P4）" "$LOGD/step1_cnv_smoke_P4.log" \
      timeout 86400 Rscript "$CODE/11_run_spatial_cnv_smoke.R" --patient P4
  if has_run "$M1"; then print_manifest "$M1"; fi
fi

## ---------- 3. 阻塞登记：本队列**不**含的项 ----------
say "--- 未排期的项（不是忘了，是过不了门）---"
say "· §13.3 第 1 项（per-spot 判定取自哪个产物）**仍未签字** ⇒ SC4 不得开始。"
say "· §14.6：借锚（§11.2 跨患者借队列良性端）**尚未实现**，脚本遇 n_ref==0 直接 stop"
say "  ⇒ 全队列 25 患者仍不可排期（自带 Normal|AAH 的只有 10 个）；借锚口径须单独签字。"
say "· 本步的实测数是上面两项排期决策的输入。"
say "===== detached 队列结束（bash PID $$，$(date '+%F %T')）====="
