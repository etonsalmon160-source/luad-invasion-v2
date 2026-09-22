#!/usr/bin/env bash
# 「全量上皮 · 锚定臂」批跑。读 epi_full_joblist.tsv（由 11_gen_epi_full_joblist.py 生成）。
#
# 用法:
#   bash 03_cnv/run_epi_full_anchor.sh              # 按清单全跑，并发 = LUAD_EPI_CONC（默认 8）
#   bash 03_cnv/run_epi_full_anchor.sh --dry-run    # 只打印将要执行什么，不跑
#   bash 03_cnv/run_epi_full_anchor.sh --one <行>   # 内部用：跑清单第 <行> 条（供 xargs 调用）
#
# 三条不可动的口径（改任何一条都必须另开预注册）:
#   · n.cores 恒为 1 —— 它是**结果相关**参数（实测改核数会改结果），不是加速旋钮。
#     并发只能靠**多进程**（本脚本的 -P），不能靠给单个 copykat 加核。
#   · 地板 840 / 锚定 >= 10 原样透传（见 joblist 生成器）。
#   · 输出目录 = $OUT/<亚型><变体>/<样本>/ —— 按亚型分目录，否则同一样本的
#     AT2 与 AT1 两次跑会互相覆盖。
set -uo pipefail

ROOT="/home/eto/luad_v2"
CNV="$ROOT/results/03_cnv"
JOBLIST="${LUAD_EPI_JOBLIST:-$CNV/epi_full_joblist.tsv}"
OUT="${LUAD_EPI_OUTROOT:-$CNV/epi_full}"
LOG="${LUAD_EPI_LOGDIR:-$ROOT/logs/epi_full_anchor}"
CONC="${LUAD_EPI_CONC:-8}"
FLOOR="${LUAD_COVERAGE_FLOOR:-840}"
DISK_SOFT_GB="${LUAD_EPI_DISK_SOFT_GB:-12}"   # 低于此值：暂停等待
DISK_HARD_GB="${LUAD_EPI_DISK_HARD_GB:-6}"    # 低于此值：整批停止

free_gb() { df -BG --output=avail "$ROOT" | tail -1 | tr -dc '0-9'; }

one_job() {
  local rowno="$1"
  local subtype sample aref an lgn ntot pred variant
  # 必须用 awk 取字段：bash 的 `IFS=$'\t' read` 会把**连续的 IFS 空白**（制表符）
  # 并成一个分隔符，从而吞掉空字段（variant 列为空时会把 note 读成 variant）。
  # awk -F'\t' 保留空字段。
  local f; f="$(awk -F'\t' -v n="$rowno" 'NR==n{for(i=1;i<=9;i++) printf "%s\x01", $i}' "$JOBLIST")"
  IFS=$'\x01' read -r subtype sample aref an lgn ntot pred variant _ <<< "${f%$'\x01'}"

  local sub_lc; sub_lc="$(echo "$subtype" | tr '[:upper:]' '[:lower:]')"
  local root_dir="$OUT/${subtype}${variant}"
  local lg="$LOG/${subtype}${variant}_${sample}.log"
  mkdir -p "$root_dir" "$LOG"

  # 断点续跑：R 脚本自己也有 .done 哨兵，这里先挡一道，免得白起进程
  if [[ -f "$root_dir/${sample}.done" && -f "$root_dir/${sample}.json" ]]; then
    echo "[SKIP] ${sample}|${subtype}${variant} 已完成" | tee -a "$lg"; return 0
  fi

  # 磁盘闸门：低于硬线直接拒跑（不删任何东西，只停手）
  local fg; fg="$(free_gb)"
  if (( fg < DISK_HARD_GB )); then
    echo "[ABORT] ${sample}|${subtype}${variant}: 可用 ${fg} GB < 硬线 ${DISK_HARD_GB} GB" | tee -a "$lg"
    return 3
  fi
  while (( fg < DISK_SOFT_GB )); do
    echo "[WAIT] 可用 ${fg} GB < 软线 ${DISK_SOFT_GB} GB，等 60s" | tee -a "$lg"; sleep 60; fg="$(free_gb)"
  done

  echo "[RUN ] ${sample}|${subtype}${variant}  病灶=${lgn} 锚定=${an}(${aref}) n=${ntot} 预测=${pred}h  可用=${fg}GB" | tee -a "$lg"
  LUAD_ANCHOR_MODE=sameSubtypeNormal \
  LUAD_ANCHOR_SUBTYPE="$sub_lc" \
  LUAD_ANCHOR_REF_SAMPLE="$aref" \
  LUAD_GP2_OUTROOT="$root_dir" \
  LUAD_COVERAGE_FLOOR="$FLOOR" \
    Rscript "$ROOT/03_cnv/01_copykat_gse308103.R" "$sample" >> "$lg" 2>&1
  local rc=$?
  # ⚠️ 判成功只能看 .done 哨兵，**不能看退出码**：R 脚本在 copykat 内部报错时
  # 仍以 0 退出（它会打印「**未完成**」且不落哨兵）。实测 n=34 的跑在 copykat
  # 内部报 "either 'k' or 'h' must be specified"，rc 却是 0。
  if [[ -f "$root_dir/${sample}.done" && -f "$root_dir/${sample}.json" ]]; then
    echo "[ OK ] ${sample}|${subtype}${variant}" | tee -a "$lg"; return 0
  fi
  echo "[FAIL] ${sample}|${subtype}${variant} 未落哨兵 rc=${rc}（详见 $lg）" | tee -a "$lg"
  return 1
}

if [[ "${1:-}" == "--one" ]]; then one_job "$2"; exit $?; fi
if [[ "${1:-}" == "--dry-run" ]]; then
  echo "清单: $JOBLIST   输出根: $OUT   日志: $LOG   并发: $CONC   可用: $(free_gb) GB"
  awk -F'\t' 'NR>1{printf "  %-12s %-4s 病灶%-6s 锚定%-6s n=%-6s %sh  %s\n",$2,$1,$5,$4,$6,$7,$3}' "$JOBLIST"
  exit 0
fi

# ---- 预检 ----
[[ -f "$JOBLIST" ]] || { echo "[FAIL] 清单不存在: $JOBLIST"; exit 1; }
mkdir -p "$OUT" "$LOG"
N=$(( $(wc -l < "$JOBLIST") - 1 ))
FG="$(free_gb)"
echo "=============== 全量上皮 · 锚定臂 ==============="
echo "  作业数 $N   并发 $CONC   输出 $OUT"
echo "  可用磁盘 ${FG} GB（软线 ${DISK_SOFT_GB} 硬线 ${DISK_HARD_GB}）"
echo "  内存余量: $(awk '/MemAvailable/{printf "%.0f", $2/1048576}' /proc/meminfo) GB"
echo "  开始 $(date -Is)"
START=$(date +%s)
(( FG < DISK_HARD_GB )) && { echo "[FAIL] 起始可用磁盘 ${FG} GB 已低于硬线"; exit 1; }

seq 2 $(( N + 1 )) | xargs -P "$CONC" -I{} bash "$0" --one {}
RC=$?
echo "  结束 $(date -Is)  用时 $(( ($(date +%s) - START) / 60 )) 分钟  xargs rc=$RC"
echo "  成功 $(grep -rl '^\[ OK \]' "$LOG" 2>/dev/null | wc -l) / $N"
echo "  磁盘余量 $(free_gb) GB"
