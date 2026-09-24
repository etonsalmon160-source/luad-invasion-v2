#!/usr/bin/env bash
# 2026-09-25 夜间队列：四臂统计 → 分辨率阶梯（PREREG_v2）→ RCTD 参考（粗 6 / 细 39）
#
# 用户 2026-09-25 拍板：「两个都跑：先 PREREG_v2，再 RCTD」「我们准备关电脑睡觉了」
#
# 规矩（沿用 scripts/queue_after_04.sh 的做法）：
#   · **串行**：任何时刻只有一个重活。本机 load 长期 45+，并发只会互相拖慢。
#   · **可续跑**：每步先看哨兵产物，已有就跳过 —— 中断后重跑本脚本即可继续。
#   · **不吞错**：某步失败就记进 STATUS，继续下一步（各步彼此独立）。
#   · 等待循环**只认可执行文件的绝对形状**（^python3 07_he_pathology/17_...），
#     绝不用会匹配到本脚本自身的 pgrep -f 模式 —— 2026-09-25 我在这上面写出过死循环。
#
# 用法（必须脱离终端，关掉电脑上的会话也不受影响）：
#   setsid nohup bash /home/eto/luad_v2/scripts/queue_overnight_20260925.sh \
#       < /dev/null > /home/eto/luad_v2/logs/overnight_20260925/queue.log 2>&1 &

set -u
ROOT=/home/eto/luad_v2
LOGD=$ROOT/logs/overnight_20260925
STATUS=$LOGD/STATUS.txt
mkdir -p "$LOGD"

say() { echo "[$(date '+%F %T')] $*" | tee -a "$STATUS"; }

# ---------- 0. 等四臂分片跑完 ----------
SHARD_PAT='^python3 07_he_pathology/17_plip_source_arms.py'
MAXWAIT=$((10 * 3600))          # 最多等 10 小时，避免异常情况下无限等
t0=$(date +%s)
say "开始：等四臂分片（4 × 14 张）跑完。最多等 $((MAXWAIT / 3600)) 小时。"
while ps -eo args= | grep -q "$SHARD_PAT"; do
    if [ $(( $(date +%s) - t0 )) -gt "$MAXWAIT" ]; then
        say "⚠️ 等待超时（$((MAXWAIT / 3600)) 小时）——不硬等，继续下一步（下面 18 号会因未跑满而硬停）。"
        break
    fi
    sleep 120
done
say "四臂分片已结束（等了 $(( ($(date +%s) - t0) / 60 )) 分钟）。"

# ---------- 1. 四臂统计（自带防呆：分片不满 14 张会硬停）----------
S1_P=$ROOT/results/07_he_pathology/source_arms/arms_stats.json
if [ -f "$S1_P" ]; then
    say "[1/4] 四臂统计：已有 $S1_P，跳过。"
else
    say "[1/4] 四臂统计：开始。"
    if python3 "$ROOT/07_he_pathology/18_arms_stats.py" >> "$LOGD/18_arms_stats.log" 2>&1; then
        say "[1/4] 四臂统计：完成 ✅（$S1_P）"
    else
        say "[1/4] 四臂统计：❌ 失败或硬停，见 18_arms_stats.log（常见原因：分片没跑满）"
    fi
fi

# ---------- 2. 分辨率阶梯 G0 成本实测（1 张切片 × 六臂）----------
S2_P=$ROOT/results/07_he_pathology/resolution_ladder_g0/ladder_summary.json
if [ -f "$S2_P" ]; then
    say "[2/4] 阶梯 G0：已有 $S2_P，跳过。"
else
    say "[2/4] 阶梯 G0（1 张切片 × 六臂）成本实测：开始。"
    if python3 "$ROOT/07_he_pathology/19_resolution_ladder.py" --limit 1 --out-tag g0 \
            >> "$LOGD/19_ladder_g0.log" 2>&1; then
        PER=$(python3 -c "
import json;d=json.load(open('$S2_P'))
print(d['per_slide_log'][0]['sec'] if d['per_slide_log'] else '?')" 2>/dev/null || echo '?')
        say "[2/4] 阶梯 G0：完成 ✅（每张约 ${PER}s ⇒ 56 张约 $(( ${PER:-0} * 56 / 60 )) 分钟）"
    else
        say "[2/4] 阶梯 G0：❌ 失败，见 19_ladder_g0.log —— **全量不跑**，先看日志"
    fi
fi

# ---------- 3. 分辨率阶梯全量（六臂 × 56 张）----------
S3_P=$ROOT/results/07_he_pathology/resolution_ladder/ladder_summary.json
if [ -f "$S3_P" ]; then
    say "[3/4] 阶梯全量：已有 $S3_P，跳过。"
else
    say "[3/4] 阶梯全量（六臂 × 56 张）：开始。"
    if python3 "$ROOT/07_he_pathology/19_resolution_ladder.py" >> "$LOGD/19_ladder.log" 2>&1; then
        say "[3/4] 阶梯全量：完成 ✅（$S3_P）"
    else
        say "[3/4] 阶梯全量：❌ 失败，见 19_ladder.log（脚本逐张落盘，重跑会跳过已完成切片）"
    fi
    # 判读（G1 管线锚：L448P0 的切片级 rho 必须 ∈ [0.60, 0.78]）
    if [ -f "$ROOT/07_he_pathology/20_ladder_stats.py" ]; then
        say "[3b] 阶梯判读（G1 锚）：开始。"
        if python3 "$ROOT/07_he_pathology/20_ladder_stats.py" >> "$LOGD/20_ladder_stats.log" 2>&1; then
            say "[3b] 阶梯判读：完成 ✅"
        else
            say "[3b] 阶梯判读：❌ 失败，见 20_ladder_stats.log"
        fi
    else
        say "[3b] 20_ladder_stats.py 不存在 ⇒ 判读未跑（原始测量值已在库里）"
    fi
fi

# ---------- 4. RCTD 参考：粗 6 型 + 细 39 型（S4 面板）----------
for CAL in a d; do
    S4_P=$ROOT/results/08_spatial_deconv/reference_${CAL}.manifest.json
    if [ -f "$S4_P" ]; then
        say "[4/$CAL] RCTD 参考 caliber=$CAL：已有 manifest，跳过。"
        continue
    fi
    say "[4/$CAL] RCTD 参考 caliber=$CAL：开始。"
    if python3 "$ROOT/08_spatial_deconv/00_build_rctd_reference.py" --caliber "$CAL" \
            >> "$LOGD/00_build_rctd_reference_${CAL}.log" 2>&1; then
        say "[4/$CAL] RCTD 参考 caliber=$CAL：完成 ✅"
    else
        say "[4/$CAL] RCTD 参考 caliber=$CAL：❌ 失败，见 00_build_rctd_reference_${CAL}.log"
    fi
done

# ---------- 收尾：写清「到哪为止」----------
cat >> "$STATUS" <<'EOF'

—— 本次夜间队列的**边界**（RCTD 线）——
RCTD 只做到「两个参考已构建」。后面两步**都还没解开**，本次**没有**也**不该**代签：
  第 0 步 空转 spot 掩码冻结：SP 三个阈值仍未签字（SIGNED 里是 None），
                            且 00_spatial_qc.py 的 --freeze 尚未实现。
  第 2 步 跑 RCTD：01_run_rctd.R 这个脚本还不存在，且它依赖第 0 步的掩码。
⇒ 任何「RCTD 结果」的说法，都必须先说明这两步没跑。详见 RCTD_PREREG.md §10.5。
EOF
say "队列全部结束。"
