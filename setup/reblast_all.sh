#!/usr/bin/env bash
# nt_viruses hazır olunca: 6 örneğin BLAST'ını TAM viral nt ile YENİDEN çalıştır,
# best_hits + ozet + rapor güncelle, analiz/'e kopyala, toplu özeti yeniden kur.
set -uo pipefail
cd "$(dirname "$0")/.."
SBU="/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj"
SAMPLES=(200225319 2858 21857478 21663260 4188mrsa 19576470psa_001)

for s in "${SAMPLES[@]}"; do
  RUN=$(ls -dt runs/*"$s"*/ 2>/dev/null | head -1)
  [ -z "$RUN" ] && { echo "[atla] $s run yok"; continue; }
  echo "=== $s : yeniden BLAST (nt_viruses) -> $RUN ==="
  rm -f "$RUN"/09_BLAST/best_hits.tsv "$RUN"/09_BLAST/blast.tsv
  python3 -c "
from pathlib import Path
from bacforge.config_loader import load_config
from bacforge.resources import detect_resources
from bacforge.tool_runner import ToolRunner
from bacforge.orchestrator import RunContext
from bacforge.modules.m09_blast import BlastModule
from bacforge.modules.m18_report import ReportModule
cfg=load_config(); res=detect_resources(cfg)
ctx=RunContext(cfg,res,Path('$RUN').resolve(),ToolRunner(Path('$RUN')/'19_Logs'),'$s')
BlastModule(ctx).run()
ReportModule(ctx).run()
print('  $s BLAST+rapor güncellendi')
"
  # analiz'e güncel 09_BLAST + 18_Final_Report + ozet kopyala
  DEST="$SBU/$s/analiz"
  if [ -d "$DEST" ]; then
    rm -rf "$DEST/09_BLAST" "$DEST/18_Final_Report"
    cp -r "$RUN"/09_BLAST "$DEST/09_BLAST"
    cp -r "$RUN"/18_Final_Report "$DEST/18_Final_Report"
    cp "$RUN"/18_Final_Report/ozet.tsv "$DEST/ozet.tsv" 2>/dev/null
    cp "$RUN"/18_Final_Report/report.pdf "$DEST/" 2>/dev/null
  fi
done
echo "status=done $(date)" > setup/reblast.DONE
echo "=== YENİDEN BLAST BİTTİ ==="
