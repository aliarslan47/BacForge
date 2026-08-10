#!/usr/bin/env python3
"""Tüm örnekleri ELİMİZDEKİ NCBI core_nt kimlikleriyle baştan güncelle (eski veri kalmasın).
- READ-ONLY NCBI: yeni uzak BLAST YOK, sadece mevcut _blast_NCBI_nt.tsv okunur.
- Rapor/INDEX/dossier/ozet yeniden üretilir.
- 08_Taxonomy: geNomad soyağacı + NCBI tür-düzeyi birleştirilir (taxonomy_final.tsv) ve
  geNomad taxonomy.tsv'nin lineage'ı yayına uygun contig'ler için tür adıyla güncellenir (.orig yedek).
- sbu-faj/<id>/analiz/ kopyaları + TUM_ORNEKLER_OZET_NCBI.tsv tazelenir.
"""
import os, shutil, csv, time
os.environ["ALI_WGS_NCBI_READONLY"] = "1"   # KRİTİK: yavaş BLAST tetikleme
from pathlib import Path
from ali_wgs.config_loader import load_config
from ali_wgs.resources import detect_resources
from ali_wgs.tool_runner import ToolRunner
from ali_wgs.orchestrator import RunContext
from ali_wgs.modules.m18_report import ReportModule

ROOT = Path(__file__).resolve().parents[1]
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")
SAMPLES = ["200225319", "2858", "21857478", "21663260", "4188mrsa", "19576470psa_001"]
runs = ROOT / "runs"
cfg = load_config(); res = detect_resources(cfg)

def latest_run(s):
    cands = sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name],
                   key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0] if cands else None

def genomad_lineage(run):
    f = run/"08_Taxonomy"/"genomad"/"contigs.filtered_annotate"/"contigs.filtered_taxonomy.tsv"
    m = {}
    if f.exists():
        for r in csv.DictReader(open(f), delimiter="\t"):
            m[r["seq_name"]] = r.get("lineage", "")
    return m

master = [["ornek","contig","sinif","uzunluk","geNomad_soyagaci","NCBI_tur","identity%","qcov%","accession","yayinlanabilirlik"]]
summary = [["ornek","yayina_uygun","NCBI_tur_kimlikleri"]]

for s in SAMPLES:
    run = latest_run(s)
    if not run:
        print(f"!! run yok: {s}"); continue
    print(f"=== {s} -> {run.name} : rapor yeniden üretiliyor (read-only) ===")
    ctx = RunContext(cfg, res, run.resolve(), ToolRunner(run/"19_Logs"), s)
    ReportModule(ctx).run()   # INDEX + dossier + ozet + report.html (BLAST tetiklemez)

    lin = genomad_lineage(run)
    yu = run/"18_Final_Report"/"yayina_uygun"
    tax_rows = [["contig","sinif","uzunluk_bp","geNomad_soyagaci","NCBI_tur(core_nt)","identity%","qcov%","accession","yayinlanabilirlik","not"]]
    # INDEX'ten yayına uygun listesi + sınıf/uzunluk
    idx = yu/"INDEX.tsv"
    pub = []
    if idx.exists():
        rd = list(csv.reader(open(idx), delimiter="\t"))
        for r in rd[1:]:
            if len(r) >= 3: pub.append(r)
    # ozet.tsv -> yayinlanabilirlik + not eşle
    pubinfo = {}
    ozet = run/"18_Final_Report"/"ozet.tsv"
    if ozet.exists():
        orows = list(csv.reader(open(ozet), delimiter="\t"))
        if orows:
            hd = orows[0]
            yi = next((i for i,h in enumerate(hd) if "Yayınlanabilirlik" in h), None)
            ni = next((i for i,h in enumerate(hd) if h.strip()=="not"), None)
            for r in orows[1:]:
                if r and r[0].startswith("contig_"):
                    yl = r[yi] if (yi is not None and yi < len(r)) else ""
                    nt = r[ni] if (ni is not None and ni < len(r)) else ""
                    pubinfo[r[0]] = (yl, nt)
    hits = []
    for d in sorted([x for x in yu.iterdir() if x.is_dir()], key=lambda x: int(x.name.split("_")[1])):
        c = d.name
        kk = d/f"{c}_NCBI_kesin_kimlik.tsv"
        species=ident=qcov=acc="-"
        raw = d/f"{c}_blast_NCBI_nt.tsv"
        if raw.exists() and raw.stat().st_size > 0:
            rr = raw.read_text().splitlines()[0].split("\t")
            species, ident, qcov, acc = rr[8], rr[2], rr[4], rr[1]
        # sınıf/uzunluk INDEX'ten
        klass=length="-"
        for r in pub:
            if r[0]==c: klass, length = r[1], r[2]
        ylp, note = pubinfo.get(c, ("",""))
        tax_rows.append([c, klass, length, lin.get(c,""), species, ident, qcov, acc, ylp, note])
        master.append([s, c, klass, length, lin.get(c,""), species, ident, qcov, acc, ylp])
        if species != "-":
            hits.append(f"{c}: {species} ({ident}%/{qcov}%)")
        else:
            hits.append(f"{c}: (NCBI bekliyor)")
    # taxonomy_final.tsv yaz
    taxf = run/"08_Taxonomy"/"taxonomy_final.tsv"
    with open(taxf, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(tax_rows)
    # geNomad taxonomy.tsv lineage'ını tür adıyla güncelle (.orig yedek)
    gtax = run/"08_Taxonomy"/"genomad"/"contigs.filtered_annotate"/"contigs.filtered_taxonomy.tsv"
    if gtax.exists():
        orig = gtax.with_suffix(".tsv.orig")
        if not orig.exists(): shutil.copy2(gtax, orig)
        sp = {row[0]: row[4] for row in tax_rows[1:] if row[4] != "-"}
        lines = orig.read_text().splitlines()
        out = [lines[0]]
        for ln in lines[1:]:
            col = ln.split("\t")
            if col and col[0] in sp and len(col) >= 5 and "[NCBI:" not in col[4]:
                col[4] = col[4].rstrip(";") + f";{sp[col[0]]} [NCBI:core_nt]"
            out.append("\t".join(col))
        gtax.write_text("\n".join(out) + "\n")

    summary.append([s, str(len(tax_rows)-1), "  |  ".join(hits)])

    # analiz/ kopyala (eski göstermesin diye 18_Final_Report tamamen yenilenir)
    dest = SBU/s/"analiz"
    if dest.exists():
        if (dest/"18_Final_Report").exists(): shutil.rmtree(dest/"18_Final_Report")
        shutil.copytree(run/"18_Final_Report", dest/"18_Final_Report")
        (dest/"08_Taxonomy").mkdir(exist_ok=True)
        shutil.copy2(taxf, dest/"08_Taxonomy"/"taxonomy_final.tsv")
        if gtax.exists(): shutil.copy2(gtax, dest/"08_Taxonomy"/"contigs.filtered_taxonomy.tsv")
        shutil.copy2(run/"18_Final_Report"/"ozet.tsv", dest/"ozet.tsv")
        # copytree kaynak klasörün ESKİ mtime'ını kopyalar -> Windows'ta bayat görünür.
        # Değişen artefaktların tarihini ŞİMDİ'ye çek (klasör + tüm dosyalar) ki tutarlı/güncel olsun.
        now = time.time()
        for base in (dest/"18_Final_Report", dest/"08_Taxonomy"):
            for p in [base, *base.rglob("*")]:
                try: os.utime(p, (now, now))
                except OSError: pass
        os.utime(dest/"ozet.tsv", (now, now))
        print(f"   analiz/ güncellendi (tarihler tazelendi): {dest}")
    else:
        print(f"   (analiz dizini yok, atlandı: {dest})")

# toplu özetler
with open(SBU/"TUM_ORNEKLER_OZET_NCBI.tsv", "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(summary)
with open(SBU/"TUM_CONTIG_NCBI_TAXONOMI.tsv", "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(master)
print("\n=== BİTTİ. Toplu özet + per-contig taksonomi yazıldı. ===")
for r in summary[1:]:
    print(f"  {r[0]} ({r[1]}): {r[2][:140]}")
