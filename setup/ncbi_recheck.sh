#!/usr/bin/env bash
# 6 örneğin YAYINA UYGUN (🟢) contig'lerini NCBI core_nt'ye REMOTE BLAST ile kesin kimlikle
# güncelle (m18 _ncbi_blast). Dossier + INDEX + analiz güncellenir, toplu özet NCBI ile kurulur.
#
# KAPANMAYA DAYANIKLI:
#   (1) PER-CONTIG KAYIT: m18 her contig'in _blast_NCBI_nt.tsv'sini diske yazar ve dolu olanı
#       ATLAR -> kapanırsa tekrar çalıştırınca biten contig'leri atlayıp kaldığı yerden devam.
#   (2) 10 SN HEARTBEAT: arka plan döngüsü 10 sn'de bir ilerlemeyi (biten/toplam, geçen süre,
#       o anki örnek) setup/ncbi_recheck.STATUS'a yazar -> kapanırsa nerede kaldığı bellidir.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$(pwd)"
SBU="/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj"
SAMPLES=(200225319 2858 21857478 21663260 4188mrsa 19576470psa_001)
STATUS="$ROOT/setup/ncbi_recheck.STATUS"
HB=10  # heartbeat aralığı (sn)
CUR_SAMPLE_FILE="$ROOT/setup/.ncbi_cur_sample"
START=$(date +%s)
echo "baslangic" > "$CUR_SAMPLE_FILE"

# yayına uygun toplam contig + biten (dolu _blast_NCBI_nt.tsv) sayısını say
count_progress() {  # echo "done total"
  local tot=0 done=0
  for s in "${SAMPLES[@]}"; do
    local run; run=$(ls -dt runs/*"$s"*/ 2>/dev/null | head -1); [ -z "$run" ] && continue
    local yu="$run/18_Final_Report/yayina_uygun"; [ -d "$yu" ] || continue
    local d c f
    for d in "$yu"/*/; do
      [ -d "$d" ] || continue; c=$(basename "$d"); tot=$((tot+1))
      f="$d${c}_blast_NCBI_nt.tsv"; [ -s "$f" ] && done=$((done+1))
    done
  done
  echo "$done $tot"
}

# 10 sn'de bir heartbeat yazan arka plan döngüsü
heartbeat() {
  while :; do
    read -r dn tt < <(count_progress)
    local now el cur; now=$(date +%s); el=$((now-START))
    cur=$(cat "$CUR_SAMPLE_FILE" 2>/dev/null)
    printf 'status=running started=%s elapsed=%ss progress=%s/%s sample=%s heartbeat=%s\n' \
      "$START" "$el" "$dn" "$tt" "$cur" "$(date '+%H:%M:%S')" > "$STATUS"
    sleep "$HB"
  done
}
heartbeat & HB_PID=$!
trap 'kill "$HB_PID" 2>/dev/null' EXIT

for s in "${SAMPLES[@]}"; do
  echo "$s" > "$CUR_SAMPLE_FILE"
  RUN=$(ls -dt runs/*"$s"*/ 2>/dev/null | head -1)
  [ -z "$RUN" ] && { echo "[atla] $s"; continue; }
  echo "=== $s : NCBI remote BLAST (yayına uygun contig'ler) -> $RUN ==="
  python3 -c "
from pathlib import Path
from ali_wgs.config_loader import load_config
from ali_wgs.resources import detect_resources
from ali_wgs.tool_runner import ToolRunner
from ali_wgs.orchestrator import RunContext
from ali_wgs.modules.m18_report import ReportModule
cfg=load_config(); res=detect_resources(cfg)
ReportModule(RunContext(cfg,res,Path('$RUN').resolve(),ToolRunner(Path('$RUN')/'19_Logs'),'$s')).run()
print('  $s tamam')
"
  DEST="$SBU/$s/analiz"
  if [ -d "$DEST" ]; then
    rm -rf "$DEST/18_Final_Report"; cp -r "$RUN/18_Final_Report" "$DEST/18_Final_Report"
    cp "$RUN/18_Final_Report/ozet.tsv" "$DEST/ozet.tsv" 2>/dev/null
  fi
done

# Toplu özet: her örneğin yayina_uygun/INDEX.tsv'sinden (NCBI kesin kimlik dahil)
python3 - "$SBU" "${SAMPLES[@]}" <<'PY'
import sys, csv, os
SBU=sys.argv[1]; samples=sys.argv[2:]
rows=[["ornek","yayina_uygun","yayina_uygun_fajlar (NCBI core_nt kesin kimlik)"]]
for s in samples:
    idx=os.path.join(SBU,s,"analiz","18_Final_Report","yayina_uygun","INDEX.tsv")
    if not os.path.exists(idx): rows.append([s,"-","INDEX yok"]); continue
    rd=list(csv.reader(open(idx),delimiter='\t'))
    body=rd[1:] if rd and rd[0][0]=="contig" else rd
    hits=[]
    for r in body:
        if len(r)>=5: hits.append(f"{r[0]}: {r[4]}")
    rows.append([s,str(len(body)),"  |  ".join(hits) if hits else "-"])
with open(os.path.join(SBU,"TUM_ORNEKLER_OZET_NCBI.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter='\t').writerows(rows)
for r in rows: print("  "+r[0]+" ("+r[1]+"): "+r[2][:120])
PY

read -r dn tt < <(count_progress)
echo "status=done $(date) progress=$dn/$tt elapsed=$(( $(date +%s)-START ))s" > "$STATUS"
kill "$HB_PID" 2>/dev/null
echo "=== NCBI RECHECK BİTTİ ($dn/$tt) ==="
