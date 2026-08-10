#!/usr/bin/env bash
# contig_698 hipotetik proteinleri YEREL diamond blastp ile yeniden anote et.
# RefSeq viral proteinleri (101MB, fonksiyonel basliklar) indir -> diamond DB -> blastp -> fonksiyon ata.
set -uo pipefail
cd "$(dirname "$0")/.."
ENV="ali-bakta"   # diamond burada
DB="databases/viral_proteins"; mkdir -p "$DB"
RUN=$(ls -dt runs/*19576470psa*/ | head -1)
R="$RUN/18_Final_Report/yayina_uygun/contig_698/04_Annotation/reannotation"
A="$RUN/18_Final_Report/yayina_uygun/contig_698/04_Annotation"
ST="$R/diamond.STATUS"; echo "status=running $(date)" > "$ST"

# 1) indir + diamond DB (yoksa)
if [ ! -f "$DB/refseq_viral.dmnd" ]; then
  echo "[1] RefSeq viral protein indiriliyor (101MB)..."
  curl -sL "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz" -o "$DB/viral.1.protein.faa.gz"
  gunzip -f "$DB/viral.1.protein.faa.gz"
  echo "[2] diamond makedb..."
  conda run -n "$ENV" diamond makedb --in "$DB/viral.1.protein.faa" -d "$DB/refseq_viral" >/dev/null 2>&1
  echo "  DB hazir: $(grep -c '^>' "$DB/viral.1.protein.faa") protein"
else
  echo "[1-2] DB zaten var: $DB/refseq_viral.dmnd"
fi

# 3) diamond blastp (hizli, hassas)
echo "[3] diamond blastp (52 hipotetik protein)..."
conda run -n "$ENV" diamond blastp -q "$R/hypothetical.faa" -d "$DB/refseq_viral" \
  --very-sensitive -k 5 -e 1e-4 -f 6 qseqid sseqid pident length evalue bitscore stitle \
  -o "$R/hypothetical_diamond.tsv" >/dev/null 2>"$R/diamond.log"

# 4) fonksiyon ata + gff guncelle
python3 - "$R" "$A" <<'PY'
import csv, re, sys
from pathlib import Path
R=Path(sys.argv[1]); A=Path(sys.argv[2])
SKIP=re.compile(r"hypothetical|uncharacterized|unknown|unnamed|putative protein|DUF\d|ORF\d|gp\d+\b", re.I)
best={}
res=R/"hypothetical_diamond.tsv"
if res.exists():
    for ln in res.read_text().splitlines():
        c=ln.split("\t")
        if len(c)<7: continue
        q=c[0]; pid=c[2]; stitle=c[6]
        # stitle: "ACC desc [organism]" -> desc
        desc=re.sub(r"\s*\[[^\]]+\]\s*$","",stitle).strip()
        desc=re.sub(r"^[A-Z]{1,3}_?\d[\w.]*\s+","",desc).strip()  # accession prefix at
        if SKIP.search(desc): continue
        if q not in best: best[q]=(desc,pid,stitle[:90])
rows=[["protein_id","atanan_fonksiyon","identity%","kaynak_hit"]]
for q,(d,p,s) in best.items(): rows.append([q,d,p,s])
with open(R/"reannotation_tablo.tsv","w",newline="") as fh:
    csv.writer(fh,delimiter="\t").writerows(rows)
# gff guncelle
gff=A/"contig_698.gff"
if gff.exists():
    out=[]
    for ln in gff.read_text().splitlines():
        if "\tCDS\t" in ln and "product=hypothetical protein" in ln:
            m=re.search(r"ID=([^;]+)",ln)
            if m and m.group(1) in best:
                ln=ln.replace("product=hypothetical protein", f"product={best[m.group(1)][0]} (diamond:RefSeq_viral)")
        out.append(ln)
    (A/"contig_698_REANNOTATED.gff").write_text("\n".join(out)+"\n")
tot=sum(1 for l in open(R/"hypothetical.faa") if l.startswith(">"))
(R/"OZET.txt").write_text(f"contig_698 yeniden anotasyon (yerel diamond, RefSeq viral)\n"
    f"hipotetik: {tot} | fonksiyon ATANAN: {len(best)} | hala bilinmeyen: {tot-len(best)}\n")
print(f"ATANAN: {len(best)}/{tot}")
for q,(d,p,_) in list(best.items())[:20]: print(f"  {q}: {d} ({p}%)")
PY
echo "status=done $(date)" > "$ST"
echo "=== REANNOTATION BITTI ==="
