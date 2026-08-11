"""M18 -- Final Report, Scientific References & Complete Export
GERCEK dashboard verisinden HTML rapor uretir. KURAL: sabit-kodlu/uydurma sonuc YOK.
Rapor PIPELINE CALISMA SIRASINA gore (M00->M18) yazilir; en basta pipeline akis semasi (Figure 1);
her cikti numarali Table/Figure olarak verilir (genom haritasi + filogeni agaci gomulu PNG, clinker HTML baglantili).
"""
from __future__ import annotations

import base64
import html
import json
import os
import shutil
import zipfile
from pathlib import Path

from .base import Module


class FinalReportExportModule(Module):
    number = "18"
    name = "final_report_export"
    folder = "M18_REPORT_EXPORT"
    enabled_key = "report"

    TOOL_REFERENCES = [
        {"tool": "FastQC", "version": "v0.12.1", "purpose": "Short read quality control", "repo": "https://github.com/s-andrews/FastQC", "doi": "10.1016/j.biotechadv.2012.01.002"},
        {"tool": "fastp", "version": "0.23.4", "purpose": "Ultra-fast FASTQ preprocessing", "repo": "https://github.com/OpenGene/fastp", "doi": "10.1093/bioinformatics/bty560"},
        {"tool": "NanoPlot", "version": "1.42.0", "purpose": "Long-read QC visualization", "repo": "https://github.com/wdecoster/NanoPlot", "doi": "10.1093/bioinformatics/bty149"},
        {"tool": "Filtlong", "version": "0.2.1", "purpose": "Long read quality filtering", "repo": "https://github.com/rrwick/Filtlong", "doi": "10.5281/zenodo.1037300"},
        {"tool": "Kraken2", "version": "2.1.3", "purpose": "Taxonomic sequence classification", "repo": "https://github.com/DerrickWood/kraken2", "doi": "10.1186/s13059-019-1891-0"},
        {"tool": "Bracken", "version": "2.8", "purpose": "Abundance estimation from Kraken2", "repo": "https://github.com/jenniferlu717/Bracken", "doi": "10.7717/peerj-cs.104"},
        {"tool": "SPAdes", "version": "3.15.5", "purpose": "De novo short-read assembler", "repo": "https://github.com/ablab/spades", "doi": "10.1089/cmb.2012.0021"},
        {"tool": "Flye", "version": "2.9.2", "purpose": "De novo long-read assembler", "repo": "https://github.com/mikolmogorov/Flye", "doi": "10.1038/s41587-019-0072-8"},
        {"tool": "Unicycler", "version": "0.5.0", "purpose": "Hybrid assembly pipeline", "repo": "https://github.com/rrwick/Unicycler", "doi": "10.1371/journal.pcbi.1005595"},
        {"tool": "Medaka", "version": "1.11.3", "purpose": "ONT consensus polishing", "repo": "https://github.com/nanoporetech/medaka", "doi": "10.1038/s41587-021-01147-4"},
        {"tool": "Polypolish", "version": "0.6.0", "purpose": "Short-read assembly polisher", "repo": "https://github.com/rrwick/Polypolish", "doi": "10.1371/journal.pcbi.1009802"},
        {"tool": "QUAST", "version": "5.2.0", "purpose": "Assembly quality assessment", "repo": "https://github.com/ablab/quast", "doi": "10.1093/bioinformatics/btt086"},
        {"tool": "CheckM2", "version": "1.0.2", "purpose": "Genome completeness & contamination", "repo": "https://github.com/chklovski/CheckM2", "doi": "10.1038/s41592-023-01940-w"},
        {"tool": "Mash", "version": "2.3", "purpose": "Fast genome distance estimation", "repo": "https://github.com/marbl/Mash", "doi": "10.1186/s13059-016-0997-x"},
        {"tool": "FastANI", "version": "1.33", "purpose": "Whole-genome Average Nucleotide Identity", "repo": "https://github.com/ParBLiSS/FastANI", "doi": "10.1038/s41467-018-07641-9"},
        {"tool": "Bakta", "version": "1.9.3", "purpose": "Rapid bacterial genome annotation", "repo": "https://github.com/oschwengers/bakta", "doi": "10.1093/bioinformatics/btac228"},
        {"tool": "mlst", "version": "2.23.0", "purpose": "Multi-locus sequence typing", "repo": "https://github.com/tseemann/mlst", "doi": "10.1038/nmicrobiol.2017.88"},
        {"tool": "Kleborate", "version": "2.3.2", "purpose": "Klebsiella genomic profiler", "repo": "https://github.com/klebgenomics/Kleborate", "doi": "10.1038/s41467-021-24448-3"},
        {"tool": "Kaptive", "version": "2.0.7", "purpose": "K and O antigen locus typing", "repo": "https://github.com/klebgenomics/Kaptive", "doi": "10.1099/mgen.0.000102"},
        {"tool": "AMRFinderPlus", "version": "3.12.8", "purpose": "NCBI Antimicrobial Resistance gene finder", "repo": "https://github.com/ncbi/amr", "doi": "10.1038/s41598-021-91456-0"},
        {"tool": "CARD/RGI", "version": "6.0.3", "purpose": "Comprehensive Antibiotic Resistance Database", "repo": "https://github.com/arpcard/rgi", "doi": "10.1093/nar/gkz935"},
        {"tool": "VFDB / ABRicate", "version": "1.0.1", "purpose": "Virulence Factor Database screening", "repo": "https://github.com/tseemann/abricate", "doi": "10.1093/nar/gkh087"},
        {"tool": "MOB-suite", "version": "3.1.9", "purpose": "Plasmid typing and reconstruction", "repo": "https://github.com/phac-nml/mob-suite", "doi": "10.1099/mgen.0.000206"},
        {"tool": "IntegronFinder", "version": "2.0.2", "purpose": "Integron detection in bacterial genomes", "repo": "https://github.com/gem-pasteur/Integron_Finder", "doi": "10.1093/nar/gkac375"},
        {"tool": "ISEScan", "version": "1.7.2.3", "purpose": "Insertion Sequence identification", "repo": "https://github.com/xiezhq/ISEScan", "doi": "10.1093/bioinformatics/btx433"},
        {"tool": "geNomad", "version": "1.7.0", "purpose": "Identification of mobile genetic elements & viruses", "repo": "https://github.com/apcamargo/genomad", "doi": "10.1038/s41587-023-01953-y"},
        {"tool": "clinker", "version": "0.0.28", "purpose": "Gene cluster comparative synteny visualization", "repo": "https://github.com/gamcil/clinker", "doi": "10.1093/bioinformatics/btab007"},
        {"tool": "Panaroo", "version": "1.3.4", "purpose": "Pangenome analysis pipeline", "repo": "https://github.com/gtonkinhill/panaroo", "doi": "10.1186/s13059-020-02090-4"},
        {"tool": "IQ-TREE2", "version": "2.2.2", "purpose": "Maximum likelihood phylogenomics", "repo": "https://github.com/iqtree/iqtree2", "doi": "10.1093/molbev/msaa015"},
    ]

    # (kod, ad, arac) — pipeline calisma sirasi
    PIPELINE_STEPS = [
        ("M00", "Input & Detection", "auto-detect"),
        ("M01", "Read QC", "fastp / Filtlong"),
        ("M02", "Taxonomic QC", "Kraken2"),
        ("M03", "Assembly", "SPAdes / Flye / Unicycler"),
        ("M04", "Polishing & QC", "QUAST + CheckM2"),
        ("M05", "Species & Closest-N", "BLAST + FastANI"),
        ("M06", "Annotation", "Bakta"),
        ("M07", "Strain Typing", "mlst"),
        ("M08", "AMR", "AMRFinderPlus"),
        ("M09", "Virulence", "VFDB/ABRicate"),
        ("M10", "Plasmids", "MOB-suite"),
        ("M11", "Mobile Elements", "ISEScan"),
        ("M12", "Phage/CRISPR", "geNomad"),
        ("M13", "Variants", "Snippy"),
        ("M14", "Genomic Context", "clinker"),
        ("M15", "Comparative", "Panaroo"),
        ("M16", "Phylogenomics", "mash + NJ"),
        ("M17", "Statistics", "aggregate"),
        ("M18", "Report", "export"),
    ]

    def inputs(self):
        return [self.ctx.run_dir / "M17_STATISTICS_VISUALIZATION" / "dashboard_data.json"]

    def outputs(self):
        return [self.out_dir / "report.html"]

    def run(self):
        self.check_inputs()
        run_dir = Path(self.ctx.run_dir)
        std_dir = self.out_dir  # M18 klasoru = SADECE rapor

        # Eski kosulardan kalan kalintlari temizle (M18 = sadece rapor + kendi summary'si)
        for stale in ("PROJECT_COMPLETE", "PROJECT_COMPLETE.zip", "scientific_references.json"):
            p = std_dir / stale
            try:
                if p.is_dir():
                    shutil.rmtree(p)
                elif p.exists():
                    p.unlink()
            except Exception:
                pass

        dash_file = run_dir / "M17_STATISTICS_VISUALIZATION" / "dashboard_data.json"
        dash_data = {}
        if dash_file.exists():
            with open(dash_file, "r", encoding="utf-8") as fh:
                dash_data = json.load(fh)

        # M18 klasorunde YALNIZ rapor
        report_path = std_dir / "report.html"
        with open(report_path, "w", encoding="utf-8") as fh:
            fh.write(self._build_html_report(dash_data, run_dir))

        # Referanslar + tam-proje zip'i RUN KOKUNDE (M18'i kirletme)
        with open(run_dir / "scientific_references.json", "w", encoding="utf-8") as fh:
            json.dump({"tools_and_databases": self.TOOL_REFERENCES}, fh, indent=2, ensure_ascii=False)
        zip_path = run_dir / "PROJECT_COMPLETE.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for sub in sorted(run_dir.iterdir()):
                if sub.is_dir() and sub.name.startswith("M"):
                    for rootd, _d, files in os.walk(sub):
                        for f in files:
                            fp = Path(rootd) / f
                            zf.write(fp, fp.relative_to(run_dir))
            for fn in ["project_manifest.json", "README.txt", "scientific_references.json"]:
                if (run_dir / fn).exists():
                    zf.write(run_dir / fn, fn)

        self.write_summary(
            status="PASS",
            statistics={"references_count": len(self.TOOL_REFERENCES)},
            details={"export_zip": str(zip_path), "report": str(report_path)},
        )

    # ---------- yardimcilar ----------
    @staticmethod
    def _esc(v):
        return html.escape(str(v)) if v is not None else "&mdash;"

    @staticmethod
    def _img_b64(path: Path) -> str | None:
        try:
            if path.exists() and path.stat().st_size > 0:
                return base64.b64encode(path.read_bytes()).decode("ascii")
        except Exception:
            pass
        return None

    def _build_html_report(self, dash: dict, run_dir: Path) -> str:
        e = self._esc
        mstat = dash.get("module_status", {}) or {}
        mods = dash.get("modules", {}) or {}
        st = lambda m: mstat.get(m, "—")
        reason = lambda m: (mods.get(m, {}) or {}).get("details", {}).get("reason", "")
        modstat = lambda m: (mods.get(m, {}) or {}).get("statistics", {}) or {}

        counters = {"t": 0, "f": 0}

        def table(caption, headers, rows):
            counters["t"] += 1
            n = counters["t"]
            if not rows:
                return (f'<div class="cap">Table {n}. {e(caption)}</div>'
                        f'<p class="na">Veri yok / analiz uygulanmadi.</p>')
            head = "".join(f"<th>{e(h)}</th>" for h in headers)
            body = "".join("<tr>" + "".join(f"<td>{e(c)}</td>" for c in r) + "</tr>" for r in rows)
            return (f'<div class="cap">Table {n}. {e(caption)}</div>'
                    f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>')

        def figure_img(caption, path, note=""):
            counters["f"] += 1
            n = counters["f"]
            b64 = self._img_b64(path)
            if not b64:
                return f'<div class="cap">Figure {n}. {e(caption)}</div><p class="na">Figur uretilmedi.</p>'
            nt = f'<div class="note">{e(note)}</div>' if note else ""
            return (f'<figure><img src="data:image/png;base64,{b64}" alt="{e(caption)}"/>'
                    f'<figcaption>Figure {n}. {e(caption)}</figcaption></figure>{nt}')

        def figure_links(caption, items):
            # items: [(label, relpath)]
            counters["f"] += 1
            n = counters["f"]
            if not items:
                return f'<div class="cap">Figure {n}. {e(caption)}</div><p class="na">Uretilmedi.</p>'
            li = "".join(f'<li><a href="{e(rel)}" target="_blank">{e(lbl)}</a></li>' for lbl, rel in items)
            return (f'<div class="cap">Figure {n}. {e(caption)}</div>'
                    f'<ul class="figlinks">{li}</ul>')

        def figure_raw(caption, inner_html):
            counters["f"] += 1
            n = counters["f"]
            return f'<figure>{inner_html}<figcaption>Figure {n}. {e(caption)}</figcaption></figure>'

        # ---- veriler ----
        species = dash.get("species") or "Belirlenemedi"
        gs = dash.get("genome_stats", {}) or {}
        cm = dash.get("checkm2", {}) or {}
        tax = dash.get("taxonomy", {}) or {}
        project_id = dash.get("project_id", "run")
        data_type = dash.get("data_type")
        platform = dash.get("platform")

        color = {"PASS": "#2e7d32", "WARNING": "#ed6c02", "NOT_APPLICABLE": "#757575",
                 "SKIPPED": "#757575", "FAIL": "#c62828"}

        # ---- Figure 1: BACFORGE'a OZGU akis semasi (dallan -> birles -> yelpaze) ----
        smap = {c: (nm, tool) for c, nm, tool in self.PIPELINE_STEPS}

        def nd(code, tool_override=None):
            c = color.get(st(code), "#546e7a")
            nm, tool = smap.get(code, (code, ""))
            tool = tool_override if tool_override is not None else tool
            return (f'<div class="nd" style="--nc:{c}"><b>{e(code)}</b>'
                    f'<span>{e(nm)}</span><em>{e(tool)}</em></div>')

        flow_inner = f"""
<div class="bf">
  <div class="bf-row"><span class="bf-io">FASTQ / FASTA girdi</span></div>
  <div class="bf-row">{nd("M00")}</div>
  <div class="bf-split">veri tipi otomatik dallanma</div>
  <div class="bf-branch">
    <div class="bf-lane" data-l="SHORT"><span class="bf-lh sh">SHORT</span>{nd("M01","fastp")}{nd("M03","SPAdes")}</div>
    <div class="bf-lane" data-l="LONG"><span class="bf-lh lo">LONG</span>{nd("M01","Filtlong")}{nd("M03","Flye")}</div>
    <div class="bf-lane" data-l="HYBRID"><span class="bf-lh hy">HYBRID</span>{nd("M01","fastp+Filtlong")}{nd("M03","Unicycler")}</div>
  </div>
  <div class="bf-merge">&#9679; birlesme &rarr; <b>genome.fasta</b> &nbsp;(+ M02 Taksonomi)</div>
  <div class="bf-row">{nd("M04")}{nd("M05")}{nd("M06")}</div>
  <div class="bf-split">genome.fasta uzerinde paralel karakterizasyon</div>
  <div class="bf-fan">{nd("M07")}{nd("M08")}{nd("M09")}{nd("M10")}{nd("M11")}{nd("M12")}{nd("M13")}</div>
  <div class="bf-split">karsilastirmali & filogenomik</div>
  <div class="bf-row">{nd("M14")}{nd("M15")}{nd("M16")}</div>
  <div class="bf-merge">&#9679; toplama &rarr; rapor</div>
  <div class="bf-row">{nd("M17")}{nd("M18")}</div>
</div>
<div class="note">Renk = modul durumu
(<b style="color:#2e7d32">PASS</b> / <b style="color:#ed6c02">WARNING</b> /
<b style="color:#757575">NOT_APPLICABLE</b>). Sema BacForge mimarisine ozgudur:
tek dikey hat degil; <b>dallan &rarr; genome.fasta'da birles &rarr; paralel yelpaze &rarr; topla</b>.</div>"""
        fig_flow = figure_raw("BacForge pipeline mimarisi (dallanma - birlesme - paralel karakterizasyon)", flow_inner)

        # ---- M01 QC (fastp.json) ----
        fastp = {}
        fp = run_dir / "M01_READ_QC_PREPROCESSING" / "fastp.json"
        if fp.exists():
            try:
                fastp = json.load(open(fp))
            except Exception:
                fastp = {}
        bf = (fastp.get("summary", {}) or {}).get("before_filtering", {})
        af = (fastp.get("summary", {}) or {}).get("after_filtering", {})

        # ---- annotation counts (M06 summary) ----
        m06s = modstat("M06")

        # ---- figur dosya yollari (run icine gore relatif; rapor M18/04_standardized'da) ----
        # rapor std_dir icinde; run koku iki ust
        def rel(p: Path):
            try:
                return os.path.relpath(p, run_dir / "M18_REPORT_EXPORT")
            except Exception:
                return str(p)

        genome_map = run_dir / "M06_GENOME_ANNOTATION" / "genome_map.png"
        tree_png = run_dir / "M16_PHYLOGENOMICS" / "phylogeny_tree.png"
        clinker_htmls = sorted((run_dir / "M14_GENOMIC_CONTEXT").glob("clinker_*.html"))
        fastp_html = run_dir / "M01_READ_QC_PREPROCESSING" / "fastp_report.html"

        # ---- tablolar (calisma sirasi) ----
        t_m00 = table("M00 — Girdi tespiti ve okuma istatistikleri",
                      ["Metrik", "Deger"],
                      [["Veri tipi", data_type], ["Platform", platform],
                       ["Tespit edilen tur (kraken2)", tax.get("dominant_organism")]])

        qc_rows = []
        if bf or af:
            qc_rows = [["Okuma (ham)", bf.get("total_reads")], ["Okuma (temiz)", af.get("total_reads")],
                       ["Baz (temiz)", af.get("total_bases")], ["Ort. okuma uzunlugu",
                        af.get("read1_mean_length")], ["Q30 orani", af.get("q30_rate")]]
        t_m01 = table("M01 — Okuma QC (fastp, once/sonra)", ["Metrik", "Deger"], qc_rows)

        t_m02 = table("M02 — Taksonomik QC (Kraken2)", ["Metrik", "Deger"],
                      [["Baskin organizma", tax.get("dominant_organism")],
                       ["Taxonomy ID", tax.get("taxonomy_id")],
                       ["Baskinlik (siniflanan turler icinde)", tax.get("dominance_percent_of_classified_species")],
                       ["Kontaminasyon tahmini %", tax.get("contamination_percent")]])

        t_m04 = table("M03/M04 — Assembly ve genom kalitesi (SPAdes + QUAST + CheckM2)",
                      ["Metrik", "Deger"],
                      [["Genom boyutu (bp)", gs.get("genome_size_bp")], ["Contig", gs.get("contig_count")],
                       ["N50", gs.get("n50")], ["En uzun contig", gs.get("largest_contig")],
                       ["GC %", gs.get("gc_percent")], ["Completeness % (CheckM2)", cm.get("completeness")],
                       ["Contamination % (CheckM2)", cm.get("contamination")],
                       ["Kalite durumu", gs.get("quality_status")]])

        c5 = dash.get("closest_5_strains", []) or []
        t_m05 = table("M05 — En yakin referanslar (FastANI)",
                      ["Rank", "Organizma", "Accession", "ANI %", "Query cov %"],
                      [[c.get("rank"), c.get("organism"), c.get("assembly_accession"),
                        c.get("ani_percent"), c.get("query_coverage")] for c in c5])

        t_m06 = table("M06 — Genom anotasyonu (Bakta)", ["Ozellik", "Sayi"],
                      [["CDS", m06s.get("cds")], ["tRNA", m06s.get("trna")], ["rRNA", m06s.get("rrna")],
                       ["Toplam ozellik", m06s.get("total_features")],
                       ["Benzersiz locus_tag", m06s.get("unique_locus_tags")]])

        mlst = dash.get("mlst")
        import re as _re
        mlst_scheme = mlst_st = None
        mlst_rows = []
        if mlst and len(mlst) >= 3:
            mlst_scheme, mlst_st = mlst[1], mlst[2]
            for tok in mlst[3:]:
                m = _re.match(r"^(.*?)\(([^)]*)\)\s*$", tok or "")
                if m:
                    locus = m.group(1); allele = m.group(2)
                    locus = locus.split("_", 1)[1] if "_" in locus else locus  # 'Pas_cpn60' -> 'cpn60'
                    mlst_rows.append([locus, allele])
                elif tok:
                    mlst_rows.append([tok, "-"])
        mlst_head = (f'<p>Şema: <b>{e(mlst_scheme)}</b> &nbsp;|&nbsp; Sekans Tipi (ST): '
                     f'<b>{e(mlst_st)}</b></p>' if mlst_scheme else '<p class="na">MLST sonucu yok.</p>')
        t_m07 = mlst_head + table("M07 — MLST alel profili (locus / allele numarası)",
                                  ["Lokus", "Alel no"], mlst_rows)

        amr = dash.get("amr_genes", []) or []
        t_m08 = table("M08 — Antimikrobiyal direnç genleri (AMRFinderPlus)",
                      ["Gen", "Ilac sinifi", "Alt sinif", "Identity %", "Contig"],
                      [[a.get("gene_symbol"), a.get("drug_class"), a.get("subclass"),
                        a.get("identity"), a.get("contig")] for a in amr])

        vir = dash.get("virulence_genes", []) or []
        t_m09 = table("M09 — Virülans genleri (VFDB/ABRicate) — ilk 25",
                      ["Gen", "Identity %", "Coverage", "Contig"],
                      [[v.get("gene_symbol"), v.get("identity"), v.get("coverage"), v.get("contig")]
                       for v in vir[:25]])

        plas = dash.get("plasmids", []) or []
        t_m10 = table("M10 — Plazmidler (MOB-suite)", ["Plasmid", "Contig", "Boyut", "Rep tipi"],
                      [[p.get("plasmid_id"), p.get("contig"), p.get("size"), p.get("rep_type")] for p in plas])

        mge = dash.get("mobile_elements", []) or []
        t_m11 = table("M11 — Mobil genetik elemanlar (ISEScan) — ilk 25",
                      ["Eleman", "Tip", "Contig", "Baslangic", "Bitis"],
                      [[m.get("element_id"), m.get("type"), m.get("contig"), m.get("start"), m.get("end")]
                       for m in mge[:25]])

        prop = dash.get("prophages", []) or []
        t_m12 = table("M12 — Prophage/viral bolgeler (geNomad)",
                      ["Phage ID", "Contig", "Uzunluk", "Topoloji", "Virus skoru"],
                      [[p.get("phage_id"), p.get("contig"), p.get("length"), p.get("topology"),
                        p.get("virus_score")] for p in prop])

        var = dash.get("variants", []) or []
        n_coding = sum(1 for v in var if v.get("gene") or v.get("locus_tag"))
        var_head = (f'<p>Toplam varyant: <b>{len(var)}</b> &nbsp;|&nbsp; kodlayan (CDS) bölgede: '
                    f'<b>{n_coding}</b>. Referans: en yakın suş (GenBank anotasyonlu).</p>')
        # Kodlayan (gen/CDS) varyantlari one al -> anlamli
        var_sorted = sorted(var, key=lambda v: 0 if (v.get("gene") or v.get("locus_tag")) else 1)
        t_m13 = var_head + table(
            "M13 — Varyantlar ve CDS etkisi (Snippy, GenBank-anotasyonlu) — ilk 30",
            ["Contig", "Pozisyon", "Tip", "Ref>Alt", "Gen", "Locus_tag", "CDS etkisi (EFFECT)", "Ürün"],
            [[v.get("contig"), v.get("pos"), v.get("type"),
              f"{v.get('ref')}>{v.get('alt')}", v.get("gene") or "—", v.get("locus_tag") or "—",
              v.get("effect") or "—", (v.get("product") or "")[:40]] for v in var_sorted[:30]])

        t_refs = table("Kullanilan araclar ve bilimsel referanslar (DOI)",
                       ["Arac", "Surum", "Amac", "DOI"],
                       [[t["tool"], t["version"], t["purpose"], t["doi"]] for t in self.TOOL_REFERENCES])

        # ---- figurler ----
        fig_qc = figure_links("Okuma QC raporu (fastp, interaktif)",
                              [("fastp_report.html", rel(fastp_html))] if fastp_html.exists() else [])
        fig_map = figure_img("Dairesel genom haritasi ve anotasyon (Bakta)", genome_map)
        clinker_items = [(p.stem.replace("clinker_", "clinker: "), rel(p)) for p in clinker_htmls]
        fig_clink = figure_links("Hedefli AMR/virülans lokuslari — gen-komsulugu sinteni (clinker, interaktif)",
                                 clinker_items)
        fig_tree = figure_img("Akrabalik agaci — sorgu + en yakin referanslar (mash mesafe + Neighbor-Joining)",
                              tree_png, note="Dal uzunluklari mash genom-geneli mesafe ile orantilidir.")

        def na(m):
            return self._na_note(st(m), reason(m))

        sec = lambda code, title, body: (
            f'<section><h2><span class="mc" style="color:{color.get(st(code), "#546e7a")}">{code}</span> '
            f'{e(title)} <span class="badge" style="background:{color.get(st(code),"#546e7a")}">{e(st(code))}</span></h2>{body}</section>')

        body = f"""
<header>
  <h1>Bakteriyel WGS Analiz Raporu</h1>
  <p class="sub">Proje: <b>{e(project_id)}</b> &nbsp;|&nbsp; Tur: <b>{e(species)}</b>
   &nbsp;|&nbsp; Veri tipi: <b>{e(data_type)}</b> &nbsp;|&nbsp; Platform: <b>{e(platform)}</b></p>
</header>

<section><h2>Pipeline Akis Semasi</h2>{fig_flow}</section>

{sec("M00","Girdi & Otomatik Tespit", t_m00)}
{sec("M01","Okuma QC & Preprocessing", t_m01 + fig_qc + na("M01"))}
{sec("M02","Taksonomik QC", t_m02 + na("M02"))}
{sec("M04","Assembly, Polishing & Genom Kalitesi", t_m04 + na("M04"))}
{sec("M05","Tur & En Yakin Referanslar", t_m05 + na("M05"))}
{sec("M06","Genom Anotasyonu", t_m06 + fig_map + na("M06"))}
{sec("M07","Suş Tiplemesi (MLST)", t_m07 + na("M07"))}
{sec("M08","Antimikrobiyal Direnç", t_m08 + na("M08"))}
{sec("M09","Virülans", t_m09 + na("M09"))}
{sec("M10","Plazmidler", t_m10 + na("M10"))}
{sec("M11","Mobil Genetik Elemanlar", t_m11 + na("M11"))}
{sec("M12","Phage / CRISPR / Defense", t_m12 + na("M12"))}
{sec("M13","Varyantlar & Mutasyonlar", t_m13 + na("M13"))}
{sec("M14","Genomik Baglam (Hedefli clinker)", fig_clink + na("M14"))}
{sec("M15","Karsilastirmali Genomik", na("M15"))}
{sec("M16","Filogenomik (Akrabalik Agaci)", fig_tree + na("M16"))}
{sec("M18","Araclar, Veritabanlari & Bilimsel Referanslar", t_refs)}
"""
        return self._wrap(project_id, body)

    def _na_note(self, status, reason):
        if status in ("NOT_APPLICABLE", "SKIPPED"):
            r = f" — {self._esc(reason)}" if reason else ""
            return f'<p class="na">Analiz uygulanmadi ({self._esc(status)}){r}</p>'
        if status == "FAIL":
            return '<p class="fail">Modul basarisiz (FAIL).</p>'
        return ""

    @staticmethod
    def _wrap(project_id, body):
        pid = html.escape(str(project_id))
        return f"""<!DOCTYPE html><html lang="tr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Bacterial WGS Report — {pid}</title>
<style>
:root{{--bg:#f6f7f9;--card:#fff;--bd:#e2e6ea;--pri:#0d8f86;--tx:#14181d;--mut:#6b7682}}
*{{box-sizing:border-box}}
body{{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--tx);margin:0;padding:32px;line-height:1.55}}
header{{border-bottom:2px solid var(--bd);padding-bottom:16px;margin-bottom:24px}}
h1{{color:var(--pri);margin:0 0 6px;font-size:27px}}
.sub{{color:var(--mut);font-size:13px;margin:0}}
section{{background:var(--card);border:1px solid var(--bd);border-radius:10px;padding:18px 22px;margin-bottom:16px}}
h2{{font-size:17px;margin:0 0 12px;display:flex;align-items:center;gap:10px}}
.mc{{font-family:ui-monospace,monospace;font-weight:700}}
.badge{{margin-left:auto;color:#fff;font-size:11px;font-weight:700;padding:2px 8px;border-radius:999px;font-family:ui-monospace,monospace}}
.cap{{font-weight:600;font-size:13px;margin:14px 0 4px;color:var(--tx)}}
.tw{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse;margin-bottom:4px}}
th,td{{text-align:left;padding:6px 10px;border-bottom:1px solid var(--bd);font-size:12.5px;font-variant-numeric:tabular-nums}}
th{{background:#eef4f3;color:var(--pri)}}
a{{color:var(--pri)}}
.na{{color:var(--mut);font-style:italic;font-size:13px}}
.fail{{color:#c62828;font-weight:600}}
.note{{color:var(--mut);font-size:11.5px;margin-top:8px}}
figure{{margin:12px 0;text-align:center}}
figure img{{max-width:100%;height:auto;border:1px solid var(--bd);border-radius:8px;background:#fff}}
figcaption{{font-size:12.5px;color:var(--tx);font-weight:600;margin-top:6px}}
.figlinks{{margin:6px 0 0;padding-left:20px;font-size:13px}}
/* BacForge'a ozgu akis semasi: dallan -> birles -> yelpaze */
.bf{{display:flex;flex-direction:column;align-items:center;gap:6px;padding:6px 0}}
.bf-row{{display:flex;flex-wrap:wrap;gap:8px;justify-content:center}}
.bf-io{{font-family:ui-monospace,monospace;font-size:12px;color:var(--mut);border:1px dashed var(--line-strong,#b9c1ca);padding:4px 12px;border-radius:20px;background:#eef4f3}}
.bf-split{{font-size:11px;color:var(--mut);text-transform:uppercase;letter-spacing:.06em;position:relative;padding-top:6px}}
.bf-split::before{{content:"";display:block;width:2px;height:12px;background:var(--pri);margin:0 auto 4px}}
.bf-branch{{display:flex;gap:14px;flex-wrap:wrap;justify-content:center}}
.bf-lane{{display:flex;flex-direction:column;gap:6px;align-items:center;padding:8px;border:1px dashed var(--line,#e2e6ea);border-radius:12px;background:#fafcfc}}
.bf-lh{{font-family:ui-monospace,monospace;font-size:11px;font-weight:700;color:#fff;padding:2px 12px;border-radius:20px}}
.bf-lh.sh{{background:#2563eb}} .bf-lh.lo{{background:#7c3aed}} .bf-lh.hy{{background:#db2777}}
.bf-merge{{font-size:12px;color:var(--tx);background:var(--pri);color:#fff;padding:4px 14px;border-radius:20px}}
.bf-merge b{{font-family:ui-monospace,monospace}}
.bf-fan{{display:flex;flex-wrap:wrap;gap:8px;justify-content:center;max-width:760px;border:1px dashed var(--pri);border-radius:12px;padding:10px;background:#f2faf9}}
.nd{{display:flex;flex-direction:column;align-items:center;background:#fff;border:2px solid var(--nc,#546e7a);border-radius:9px;padding:6px 10px;min-width:104px;box-shadow:0 1px 2px rgba(0,0,0,.05)}}
.nd b{{font-family:ui-monospace,monospace;font-size:12px;color:var(--nc,#546e7a)}}
.nd span{{font-size:11px;font-weight:600;text-align:center}}
.nd em{{color:var(--mut);font-style:normal;font-size:10px;margin-top:1px}}
</style></head><body>{body}</body></html>"""
