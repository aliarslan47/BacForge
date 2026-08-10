#!/usr/bin/env python3
"""VIRIDIC-tarzı GERÇEK intergenomik benzerlik — yeni tür adayları için ICTV doğrulaması.
Referans genomu accession'dan efetch ile çek -> contig vs ref İKİ YÖNLÜ blastn ->
intergenomik benzerlik = 100 * 2*ozdes_bp / (len_contig + len_ref).
ICTV: tür <%95, cins <%70. Proxy'yi (pident*qcov) gerçek değerle değiştirir."""
import csv, subprocess, sys
from pathlib import Path

runs = Path("runs")
CACHE = Path("databases/ref_genomes_cache"); CACHE.mkdir(parents=True, exist_ok=True)
ENV = "ali-blast"

# (ornek, contig, ref_accession, ref_isim)  — en anlamlı adaylar
TARGETS = [
    ("4188mrsa", "contig_88", "OP172804.1", "Enterococcus phage EH93P1"),
    ("200225319", "contig_4", "MW929175.1", "Pseudomonas phage K4"),
    ("2858", "contig_4", "PV844377.1", "Pseudomonas phage CP-p-PA-21037"),
    ("4188mrsa", "contig_36", "PX843246.1", "Pseudomonas phage Kovar531"),
    ("21857478", "contig_228", "PZ179922.1", "Pseudomonas phage vB_PaeP_FmM12"),
    ("19576470psa_001", "contig_698", "NC_048634.1", "Staphylococcus phage P282"),
]

def latest(s):
    c = sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name], key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0] if c else None

def seqlen(fa):
    return sum(len(l.strip()) for l in open(fa) if not l.startswith(">"))

def fetch_ref(acc):
    out = CACHE/f"{acc}.fasta"
    if out.exists() and out.stat().st_size > 0:
        return out
    print(f"  efetch {acc} ...", flush=True)
    r = subprocess.run(["conda","run","-n",ENV,"efetch","-db","nuccore","-id",acc,"-format","fasta"],
                       capture_output=True, text=True, timeout=300)
    if r.returncode == 0 and r.stdout.strip().startswith(">"):
        out.write_text(r.stdout); return out
    print(f"  !! efetch başarısız {acc}: {r.stderr[:120]}")
    return None

def identical_bp(query, subject):
    """blastn query->subject, HSP'lerden contig üzerinde örtüşmesiz özdeş bp (overlap düzeltmeli)."""
    r = subprocess.run(["conda","run","-n",ENV,"blastn","-query",str(query),"-subject",str(subject),
                        "-outfmt","6 qstart qend pident length","-evalue","1e-5"],
                       capture_output=True, text=True, timeout=600)
    ivs = []
    for ln in r.stdout.splitlines():
        c = ln.split("\t")
        if len(c) >= 4:
            qs, qe, pid, alen = int(c[0]), int(c[1]), float(c[2]), int(c[3])
            ivs.append((min(qs,qe), max(qs,qe), pid/100))
    # örtüşmesiz: pozisyon başına en yüksek pident'i say
    if not ivs: return 0
    ivs.sort()
    pos_end = 0; tot = 0.0
    for s,e,fr in ivs:
        s = max(s, pos_end+1)
        if e >= s:
            tot += (e-s+1)*fr; pos_end = e
    return tot

rows = [["ornek","contig","len_contig","ref_isim","ref_acc","len_ref",
         "intergenomik_benzerlik%","ICTV_karar(gercek)"]]
for s, c, acc, name in TARGETS:
    run = latest(s)
    cf = run/"18_Final_Report"/"yayina_uygun"/c/f"{c}.fasta"
    if not cf.exists(): print(f"!! contig fasta yok {s}/{c}"); continue
    ref = fetch_ref(acc)
    if not ref:
        rows.append([s,c,seqlen(cf),name,acc,"-","efetch yok","-"]); continue
    lc, lr = seqlen(cf), seqlen(ref)
    ab1 = identical_bp(cf, ref)   # contig -> ref
    ab2 = identical_bp(ref, cf)   # ref -> contig
    sim = round(100*(ab1+ab2)/(lc+lr), 2)
    verdict = ("🔵 bilinen tür" if sim>=95 else "🟢 YENİ TÜR (aynı cins)" if sim>=70
               else "🟠 YENİ CİNS" if sim>=50 else "🔴 ÇOK IRAKSAK (yeni cins/alt-aile)")
    rows.append([s,c,lc,name,acc,lr,sim,verdict])
    print(f"  {s}/{c} vs {name}: intergenomik %{sim}  -> {verdict}", flush=True)

SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
with open(SBU/"YENI_TUR_VIRIDIC_DOGRULAMA.tsv","w",newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(rows)
print("\nYENI_TUR_VIRIDIC_DOGRULAMA.tsv yazıldı.")
