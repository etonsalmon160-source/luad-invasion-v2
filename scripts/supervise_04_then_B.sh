#!/usr/bin/env bash
# 守护器：替 04 扛住"会话被关"——被杀就自动续跑；04 全部跑完后接着跑 B 套。
#
# 为什么不能真"脱会话"：进程的 session 在创建那一刻就定死了，事后改不了
# （reptyr 未安装，且它只是接管终端、不换 session）。但 04 具备"已存在就跳过"，
# 续跑一次的代价只有当前那一张切片（约 10 分钟）⇒ 守住 + 死了续 = 效果上等同脱离。
#
# 启动（必须脱离终端，否则关会话会被一起收走）：
#   cd /home/eto/luad_v2 && setsid nohup bash scripts/supervise_04_then_B.sh \
#       < /dev/null > logs/supervise_04_then_B.log 2>&1 &
#
# 唯一保不住的情况：服务器重启。那需要往 crontab 加开机项才能自动续上（需用户同意）。

set -u
ROOT=/home/eto/luad_v2
POLL=${POLL:-60}
MAX_RESTART=${MAX_RESTART:-10}

S04="$ROOT/07_he_pathology/04_spot_annotation.py"
SB="$ROOT/07_he_pathology/11_spot_annotation_B.py"
SENTINEL="$ROOT/results/07_he_pathology/spot_annotation/spot_annotation_summary.json"
PER_SLIDE="$ROOT/results/07_he_pathology/spot_annotation/per_slide"
B_PER="$ROOT/results/07_he_pathology/spot_annotation_B/per_slide"
B_EMB="$ROOT/results/07_he_pathology/spot_annotation_B/embeds"

log(){ echo "[$(date '+%F %T')] $*"; }
n04(){ ls "$PER_SLIDE"/*.csv.gz 2>/dev/null | wc -l; }
nB(){ ls "$B_PER"/*.csv.gz 2>/dev/null | wc -l; }
nE(){ ls "$B_EMB"/*.npy 2>/dev/null | wc -l; }

cd "$ROOT" || exit 1

# 续跑前验已落盘的切片：被杀在写盘中途会留下"半个 gz"（没有校尾），
# 而 04 只看"文件在不在" ⇒ 不检查就会把坏文件当已完成永久跳过。
check_integrity(){
    local f bad=0
    for f in "$PER_SLIDE"/*.csv.gz; do
        [ -e "$f" ] || continue
        if ! gzip -t "$f" 2>/dev/null; then
            log "   🗑 丢弃半截文件：$(basename "$f")（续跑会重做这张）"
            rm -f "$f"; bad=$((bad + 1))
        fi
    done
    [ "$bad" -gt 0 ] && log "   共丢弃 $bad 个，剩 $(n04) 张有效"
    return 0
}

# 有任何一个 04 在跑就返回 0（按 cmdline 认，防止 PID 复用后误判）
running04(){
    local p
    for p in /proc/[0-9]*; do
        [ -r "$p/cmdline" ] || continue
        tr '\0' ' ' < "$p/cmdline" 2>/dev/null | grep -q '04_spot_annotation' && return 0
    done
    return 1
}

log "守护器启动 pid=$$ sid=$(ps -o sid= -p $$ 2>/dev/null | tr -d ' ')"
log "任务一：保 04 跑完（已完成 $(n04)/56 张）"

runs=0
while [ ! -f "$SENTINEL" ]; do
    if running04; then
        sleep "$POLL"
        continue
    fi
    runs=$((runs + 1))
    if [ "$runs" -gt "$MAX_RESTART" ]; then
        log "❌ 已续跑 $MAX_RESTART 次仍未收尾（per_slide=$(n04)/56），放弃守护；B 套仍继续"
        break
    fi
    log "⚠️ 04 已不在且未写哨兵 → 第 $runs 次续跑（已完成 $(n04)/56 张，会自动跳过已有）"
    check_integrity
    python3 -u "$S04" >> "$ROOT/logs/04_rerun${runs}.log" 2>&1
    log "   第 $runs 次续跑退出 rc=$?"
done

if [ -f "$SENTINEL" ]; then
    log "✅ 04 收尾完成，per_slide=$(n04)/56 张"
else
    log "⚠️ 04 未收尾（per_slide=$(n04)/56 张），照常启动 B 套"
fi

# ——— 任务二：B 套全队列 ———
if [ ! -f "$SB" ]; then
    log "❌ 找不到 $SB，退出"
    exit 1
fi
log "▶ 开始 B 套全队列：$SB"
t0=$(date +%s)
python3 -u "$SB"
rc=$?
log "■ B 套结束 rc=$rc，耗时 $(( ($(date +%s) - t0) / 60 )) 分钟"
log "产物：per_slide=$(nB) 张 / embeds=$(nE) 张"
log "守护器退出"
exit $rc
