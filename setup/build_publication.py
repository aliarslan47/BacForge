#!/usr/bin/env python3
"""YAYIN MİMARİSİ — her yayına uygun faj contig'i için 01-12 numaralı tam karakterizasyon.
Pilot: tek örnek. Kullanım: python3 setup/build_publication.py <ornek>
Bu dosya ileride m20-m25 modüllerine refaktör edilecek; pilot için tek orkestratör.
"""
import csv, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
runs = ROOT/"runs"
SBU = Path("/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj")

def latest(s):
    c = sorted([p for p in runs.iterdir() if p.is_dir() and s in p.name], key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0] if c else None

def w(path, rows):
    with open(path, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)

def read_kv(f):
    return {r[0]: r[1] for r in csv.reader(open(f), delimiter="\t") if len(r) >= 2} if Path(f).exists() else {}

def run_tool(cmd, env, timeout=1200, cwd=None):
    full = ["conda", "run", "-n", env] + cmd
    return subprocess.run(full, capture_output=True, text=True, timeout=timeout, cwd=cwd)

# ---------- numaralı mimari ----------
def folder01_genome(c, d, dest, stat, checkv_dtr):
    o = dest/"01_Genome"; o.mkdir(parents=True, exist_ok=True)
    shutil.copy2(d/f"{c}.fasta", o/f"{c}.fasta")
    rows = [["metrik","deger"],
            ["contig", c], ["uzunluk_bp", stat.get("uzunluk_bp","")],
            ["GC%", stat.get("GC%","")], ["derinlik_x", stat.get("derinlik_x","")]]
    # DTR / genom uçları (CheckV complete_genomes)
    if checkv_dtr:
        rows += [["uc_tipi(prediction)", checkv_dtr.get("prediction_type","")],
                 ["DTR_tekrar_uzunlugu_bp", checkv_dtr.get("repeat_length","")],
                 ["uc_guven", checkv_dtr.get("confidence_level","")]]
    else:
        rows.append(["uc_tipi(prediction)", "DTR tespit edilmedi (lineer/permüte olabilir)"])
    w(o/f"{c}_genome_stats.tsv", rows)

def folder02_completeness(c, dest, qs):
    o = dest/"02_Completeness"; o.mkdir(exist_ok=True)
    w(o/f"{c}_checkv.tsv", [["metrik","deger"],
        ["CheckV_kalite", qs.get("checkv_quality","")], ["MIUViG", qs.get("miuvig_quality","")],
        ["tamlik%", qs.get("completeness","")], ["tamlik_yontemi", qs.get("completeness_method","")],
        ["kontaminasyon%", qs.get("contamination","")], ["gen_sayisi", qs.get("gene_count","")],
        ["viral_gen", qs.get("viral_genes","")], ["host_gen", qs.get("host_genes","")]])

def folder03_taxonomy(c, d, dest):
    o = dest/"03_Taxonomy"; o.mkdir(exist_ok=True)
    for f in [f"{c}_taksonomi.tsv", f"{c}_NCBI_kesin_kimlik.tsv", f"{c}_blast_NCBI_nt.tsv", f"{c}_NCBI_kaynak.txt"]:
        if (d/f).exists(): shutil.copy2(d/f, o/f)

def folder04_annotation(c, d, dest):
    o = dest/"04_Annotation"; o.mkdir(exist_ok=True)
    P = d/"pharokka"
    for f in ["pharokka.gbk","pharokka.tbl","pharokka.gff","pharokka_cds_final_merged_output.tsv","pharokka_cds_functions.tsv"]:
        if (P/f).exists(): shutil.copy2(P/f, o/f.replace("pharokka", c))
    return read_kv(P/"pharokka_cds_functions.tsv")  # not gerçekten kv değil ama kategoriler

def folder05_trna_crispr(c, d, dest):
    o = dest/"05_tRNA_CRISPR"; o.mkdir(exist_ok=True)
    P = d/"pharokka"
    for f in ["trnascan_out.gff","pharokka_minced.gff","pharokka_minced_spacers.txt"]:
        if (P/f).exists(): shutil.copy2(P/f, o/f"{c}_{f}")

def folder06_lifestyle(c, d, dest, cf):
    o = dest/"06_Lifestyle"; o.mkdir(exist_ok=True)
    fa = o/f"{c}.fasta"; shutil.copy2(d/f"{c}.fasta", fa)
    # BACPHLIP
    res = run_tool(["bacphlip","-i",str(fa),"-f"], "ali-bacphlip", timeout=600)
    bac = o/f"{c}.fasta.bacphlip"
    virulent = temperate = ""
    if bac.exists():
        lines = [ln.split("\t") for ln in bac.read_text().splitlines() if ln.strip("\t").strip() != "" or "\t" in ln]
        # format: satır0 = ['', 'Virulent', 'Temperate'], satır1 = ['0', vir, temp]
        if len(lines) >= 2:
            hdr, data = lines[0], lines[1]
            for i, h in enumerate(hdr):
                if h == "Virulent" and i < len(data): virulent = f"{float(data[i]):.3f}"
                if h == "Temperate" and i < len(data): temperate = f"{float(data[i]):.3f}"
    integ = cf.get("integration and excision", "?")
    lytic = ""
    try:
        lytic = "Lytic (virulent)" if float(virulent) >= 0.5 else "Temperate (ılımlı)"
    except: lytic = "belirsiz"
    w(o/f"{c}_lifestyle.tsv", [["metrik","deger"],
        ["BACPHLIP_virulent_olasilik", virulent], ["BACPHLIP_temperate_olasilik", temperate],
        ["integrase_gen(integration&excision)", integ],
        ["yasam_tarzi", lytic + (" + integrase YOK" if integ in ("0","") else " + integrase VAR")]])
    # ara dosyaları temizle (final .bacphlip + _lifestyle.tsv yeter)
    for junk in list(o.glob("*.6frame")) + list(o.glob("*.hmmsearch")) + list(o.glob("*.hmmsearch.tsv")):
        try: junk.unlink()
        except OSError: pass
    try: fa.unlink()  # 01_Genome'da kopyası var
    except OSError: pass
    return {"virulent": virulent, "lytic": lytic, "integrase": integ}

def folder07_safety(c, d, dest, life):
    o = dest/"07_Safety_AMR_Virulence"; o.mkdir(exist_ok=True)
    P = d/"pharokka"
    n_amr = n_vf = n_card = 0
    # AMRFinder
    af = d/f"{c}_amrfinder.tsv"
    if af.exists():
        shutil.copy2(af, o/f"{c}_amrfinder.tsv")
        n_amr = sum(1 for ln in af.read_text().splitlines() if ln and not ln.startswith(("#","Protein id")))
    # VFDB + CARD (pharokka)
    for src, tag in [("top_hits_vfdb.tsv","vfdb"),("top_hits_card.tsv","card")]:
        if (P/src).exists():
            shutil.copy2(P/src, o/f"{c}_{tag}.tsv")
            n = sum(1 for ln in (P/src).read_text().splitlines() if ln.strip()) - 1
            if tag == "vfdb": n_vf = max(n,0)
            else: n_card = max(n,0)
    integ = str(life.get("integrase", "0"))
    integrase_present = integ not in ("0", "", "?")
    lytic_ok = "Lytic" in (life.get("lytic") or "")
    # Terapötik uygunluk: lytic + integrase YOK + AMR yok + virülans yok (integrase=lizojeni riski)
    safe = lytic_ok and not integrase_present and n_amr == 0 and n_vf == 0 and n_card == 0
    verdict = ("✅ TERAPÖTİK ADAY uygun (lytic + integrase yok + AMR/virülans yok)" if safe else
               "⚠️ İNCELE: " + ", ".join(filter(None, [
                   "" if lytic_ok else "BACPHLIP temperate",
                   f"integrase VAR (lizojeni riski)" if integrase_present else "",
                   f"AMR×{n_amr}" if n_amr else "", f"VFDB×{n_vf}" if n_vf else "", f"CARD×{n_card}" if n_card else ""])))
    w(o/f"{c}_safety_karar.tsv", [["kriter","deger"],
        ["yasam_tarzi", life.get("lytic","")], ["integrase_gen", integ],
        ["AMR_gen_sayisi", n_amr], ["virulans_VFDB", n_vf], ["AMR_CARD", n_card],
        ["TERAPOTIK_KARAR", verdict]])
    return {"safe": safe, "verdict": verdict, "n_amr": n_amr, "n_vf": n_vf, "n_card": n_card}

def folder11_map(c, d, dest):
    o = dest/"11_Genome_Map"; o.mkdir(exist_ok=True)
    for f in [f"{c}_genome_map.png",f"{c}_genome_map.svg",f"{c}_genome_map_clean.png",f"{c}_genome_map_clean.svg"]:
        if (d/f).exists(): shutil.copy2(d/f, o/f)

# ---------- ana ----------
def main(sample):
    run = latest(sample)
    yu = run/"18_Final_Report"/"yayina_uygun"
    qs_all = {}
    qf = run/"17_Completeness"/"checkv"/"quality_summary.tsv"
    if qf.exists():
        rd = list(csv.DictReader(open(qf), delimiter="\t"))
        qs_all = {r["contig_id"]: r for r in rd}
    dtr_all = {}
    cg = run/"17_Completeness"/"checkv"/"complete_genomes.tsv"
    if cg.exists():
        for r in csv.DictReader(open(cg), delimiter="\t"): dtr_all[r["contig_id"]] = r
    contigs = sorted([x for x in yu.iterdir() if x.is_dir()], key=lambda x: int(x.name.split("_")[1]))
    summary = [["contig","uzunluk","tamlik%","yasam_tarzi","AMR","VFDB","terapotik_karar"]]
    for d in contigs:
        c = d.name
        dest = d  # numaralı klasörler dossier içine
        stat = read_kv(d/f"{c}_istatistik.tsv")
        print(f"  {c}: 01-07,11 ...", flush=True)
        folder01_genome(c, d, dest, stat, dtr_all.get(c))
        folder02_completeness(c, dest, qs_all.get(c, {}))
        folder03_taxonomy(c, d, dest)
        cf = {}
        ff = d/"pharokka"/"pharokka_cds_functions.tsv"
        if ff.exists():
            for r in csv.reader(open(ff), delimiter="\t"):
                if len(r) >= 2 and r[0] != "Description": cf[r[0]] = r[1]
        folder04_annotation(c, d, dest)
        folder05_trna_crispr(c, d, dest)
        life = folder06_lifestyle(c, d, dest, cf)
        safe = folder07_safety(c, d, dest, life)
        folder11_map(c, d, dest)
        summary.append([c, stat.get("uzunluk_bp",""), stat.get("CheckV_tamlik%",""),
                        life.get("lytic",""), safe["n_amr"], safe["n_vf"], safe["verdict"]])
    w(run/"18_Final_Report"/"YAYIN_OZET.tsv", summary)
    print("\n=== 01-07,11 BİTTİ ===")
    for r in summary[1:]:
        print(f"  {r[0]:<12} {r[1]:>8}bp tamlık%{r[2]:<6} {r[3]:<22} {r[6]}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "4188mrsa")
