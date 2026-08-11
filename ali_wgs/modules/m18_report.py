"""M18 -- Final Report, Scientific References & Complete Export
GERÇEK dashboard verisinden HTML rapor üretir. KURAL: sabit-kodlu/uydurma sonuç YOK.
Üretilmemiş veriler 'Analiz yapılmadı (NA)' / 'Bulgu yok' olarak dürüstçe gösterilir.
Araç referansları (DOI) platform tool-registry'sidir; per-run sürümler provenance dosyalarındadır.
"""
from __future__ import annotations

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

    MODULE_NAMES = {
        "M00": "Input & Data Detection", "M01": "Read QC", "M02": "Taxonomic QC",
        "M03": "Genome Assembly", "M04": "Polishing & Genome QC", "M05": "Species & Reference ID",
        "M06": "Genome Annotation", "M07": "Strain Typing", "M08": "AMR", "M09": "Virulence",
        "M10": "Plasmids", "M11": "Mobile Genetic Elements", "M12": "Phage/CRISPR/Defense",
        "M13": "Variants", "M14": "Genomic Context", "M15": "Comparative Genomics",
        "M16": "Phylogenomics", "M17": "Statistics", "M18": "Report",
    }

    def inputs(self):
        return [self.ctx.run_dir / "M17_STATISTICS_VISUALIZATION" / "dashboard_data.json"]

    def outputs(self):
        return [self.out_dir / "report.html"]

    def run(self):
        self.check_inputs()
        run_dir = Path(self.ctx.run_dir)
        std_dir = self.sub_dir("04_standardized")

        with open(std_dir / "scientific_references.json", "w", encoding="utf-8") as fh:
            json.dump({"tools_and_databases": self.TOOL_REFERENCES}, fh, indent=2, ensure_ascii=False)

        dash_file = run_dir / "M17_STATISTICS_VISUALIZATION" / "dashboard_data.json"
        dash_data = {}
        if dash_file.exists():
            with open(dash_file, "r", encoding="utf-8") as fh:
                dash_data = json.load(fh)

        report_path = std_dir / "report.html"
        with open(report_path, "w", encoding="utf-8") as fh:
            fh.write(self._build_html_report(dash_data))

        # PROJECT_COMPLETE bundle + zip
        export_dir = self.sub_dir("02_work") / "PROJECT_COMPLETE"
        if export_dir.exists():
            shutil.rmtree(export_dir)
        export_dir.mkdir(parents=True, exist_ok=True)
        for sub in run_dir.iterdir():
            if sub.is_dir() and sub.name.startswith("M") and sub.name != self.folder:
                shutil.copytree(sub, export_dir / sub.name, dirs_exist_ok=True)
        for fn in ["project_manifest.json", "README.txt"]:
            if (run_dir / fn).exists():
                shutil.copy(run_dir / fn, export_dir / fn)
        zip_path = std_dir / "PROJECT_COMPLETE.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(export_dir):
                for f in files:
                    fp = Path(root) / f
                    zf.write(fp, fp.relative_to(export_dir))

        self.write_summary(
            status="PASS",
            statistics={"references_count": len(self.TOOL_REFERENCES)},
            details={"export_zip": str(zip_path), "report": str(report_path)},
        )

    # ---- yardımcılar ----
    @staticmethod
    def _esc(v):
        return html.escape(str(v)) if v is not None else "&mdash;"

    def _na_or(self, status: str, body_html: str, reason: str = "") -> str:
        """Modül NA/SKIPPED ise dürüst not; değilse gerçek içerik."""
        if status in ("NOT_APPLICABLE", "SKIPPED"):
            r = f" — {self._esc(reason)}" if reason else ""
            return f'<p class="na">Analiz uygulanmadı ({self._esc(status)}){r}</p>'
        if status == "FAIL":
            return '<p class="fail">Modül başarısız (FAIL). Ayrıntı için modül loglarına bakın.</p>'
        return body_html

    def _build_html_report(self, dash: dict) -> str:
        e = self._esc
        mstat = dash.get("module_status", {}) or {}
        mods = dash.get("modules", {}) or {}

        def st(m):
            return mstat.get(m, "—")

        def reason(m):
            return (mods.get(m, {}) or {}).get("details", {}).get("reason", "")

        species = dash.get("species")
        gs = dash.get("genome_stats", {}) or {}
        cm = dash.get("checkm2", {}) or {}
        tax = dash.get("taxonomy", {}) or {}

        # Modül durum ızgarası
        color = {"PASS": "#4ade80", "WARNING": "#fbbf24", "NOT_APPLICABLE": "#94a3b8",
                 "SKIPPED": "#94a3b8", "FAIL": "#f87171"}
        status_cells = ""
        for m in [f"M{i:02d}" for i in range(0, 19)]:
            s = st(m)
            c = color.get(s, "#64748b")
            status_cells += (f'<div class="mcell"><span class="mnum">{m}</span>'
                             f'<span class="mname">{e(self.MODULE_NAMES.get(m, ""))}</span>'
                             f'<span class="mstat" style="color:{c}">{e(s)}</span></div>')

        # CheckM2 (gerçek ya da hesaplanmadı)
        comp = cm.get("completeness")
        cont = cm.get("contamination")
        checkm_html = (f"CheckM2 Completeness: <strong>{e(comp)}%</strong> | Contamination: "
                       f"<strong>{e(cont)}%</strong> | Quality: <strong>{e(gs.get('quality_status'))}</strong>"
                       if comp is not None else
                       f'<span class="na">CheckM2 çalışmadı — completeness/contamination hesaplanmadı. '
                       f'{e(cm.get("note"))}</span>')

        # Closest-5 tablo
        c5 = dash.get("closest_5_strains", []) or []
        if c5:
            c5_rows = "".join(
                f"<tr><td>#{e(c.get('rank'))}</td><td>{e(c.get('organism'))}</td>"
                f"<td>{e(c.get('strain'))}</td><td>{e(c.get('ani_percent'))}</td>"
                f"<td>{e(c.get('query_coverage'))}</td></tr>" for c in c5)
            c5_html = (f'<table><thead><tr><th>Rank</th><th>Organism</th><th>Strain/Accession</th>'
                       f'<th>ANI %</th><th>Query Cov %</th></tr></thead><tbody>{c5_rows}</tbody></table>')
        else:
            c5_html = self._na_or(st("M05"), '<p class="na">Closest-5 hesaplanmadı.</p>',
                                  reason("M05") or (mods.get("M05", {}).get("details", {}).get("closest_5_note", "")))

        # MLST
        mlst = dash.get("mlst")
        if mlst and len(mlst) >= 3:
            mlst_html = f"Şema: <strong>{e(mlst[1])}</strong> | ST: <strong>{e(mlst[2])}</strong>"
        else:
            mlst_html = '<p class="na">MLST sonucu yok.</p>'

        # AMR tablo
        amr = dash.get("amr_genes", []) or []
        if amr:
            amr_rows = "".join(
                f"<tr><td>{e(a.get('gene_symbol'))}</td><td>{e(a.get('drug_class'))}</td>"
                f"<td>{e(a.get('identity'))}</td><td>{e(a.get('contig'))}</td></tr>" for a in amr)
            amr_html = (f'<table><thead><tr><th>Gene</th><th>Drug Class</th><th>Identity</th>'
                        f'<th>Contig</th></tr></thead><tbody>{amr_rows}</tbody></table>')
        else:
            amr_html = ('<p>AMRFinderPlus çalıştı; direnç geni bulunamadı (bulgu yok).</p>'
                        if st("M08") == "PASS" else self._na_or(st("M08"), "", reason("M08")))

        def count_section(m, items, label, tool):
            n = len(items)
            if st(m) == "PASS":
                return f"{tool} çalıştı — {e(label)}: <strong>{n}</strong>" + ("" if n else " (bulgu yok)")
            return self._na_or(st(m), f"{tool}: <strong>{n}</strong>", reason(m))

        vir_html = count_section("M09", dash.get("virulence_genes", []), "virülans geni", "VFDB/abricate")
        plas_html = count_section("M10", dash.get("plasmids", []), "plazmid", "MOB-suite")
        mge_html = count_section("M11", dash.get("mobile_elements", []), "MGE", "ISEScan")
        phage_html = count_section("M12", dash.get("prophages", []), "prophage", "geNomad")
        var_html = count_section("M13", dash.get("variants", []), "varyant (SNP)", "Snippy")

        ref_rows = "".join(
            f"<tr><td><strong>{e(t['tool'])}</strong></td><td>{e(t['version'])}</td>"
            f"<td>{e(t['purpose'])}</td><td><a href='https://doi.org/{e(t['doi'])}' target='_blank'>{e(t['doi'])}</a></td></tr>"
            for t in self.TOOL_REFERENCES)

        project_id = dash.get("project_id", "run")
        data_type = dash.get("data_type")
        platform = dash.get("platform")
        species_disp = species if species else "Belirlenemedi"

        return f"""<!DOCTYPE html>
<html lang="tr"><head><meta charset="utf-8">
<title>Bacterial WGS Report — {e(project_id)}</title>
<style>
:root{{--bg:#0f172a;--card:#1e293b;--bd:#334155;--pri:#38bdf8;--tx:#f8fafc;--mut:#94a3b8}}
body{{font-family:system-ui,sans-serif;background:var(--bg);color:var(--tx);margin:0;padding:28px;line-height:1.5}}
h1{{color:var(--pri);margin:0 0 6px;font-size:26px}}
.sub{{color:var(--mut);font-size:13px;margin:0 0 22px}}
.section{{background:var(--card);border:1px solid var(--bd);border-radius:10px;padding:18px 22px;margin-bottom:18px}}
.section h2{{margin:0 0 10px;color:var(--pri);font-size:17px;border-bottom:1px solid var(--bd);padding-bottom:6px}}
table{{width:100%;border-collapse:collapse;margin-top:8px}}
th,td{{text-align:left;padding:7px 10px;border-bottom:1px solid var(--bd);font-size:13px}}
th{{background:#0f172a;color:var(--pri)}}
a{{color:var(--pri)}}
.na{{color:var(--mut);font-style:italic}}
.fail{{color:#f87171;font-weight:600}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px}}
.mcell{{background:#0f172a;border:1px solid var(--bd);border-radius:7px;padding:8px 10px;display:flex;flex-direction:column;gap:2px}}
.mnum{{font-family:monospace;color:var(--pri);font-weight:600;font-size:12px}}
.mname{{font-size:11px;color:var(--mut)}}
.mstat{{font-family:monospace;font-size:12px;font-weight:600}}
.note{{color:var(--mut);font-size:11px;margin-top:14px}}
</style></head><body>
<h1>Bakteriyel WGS Analiz Raporu</h1>
<p class="sub">Proje: <strong>{e(project_id)}</strong> &nbsp;|&nbsp; Tür: <strong>{e(species_disp)}</strong>
 &nbsp;|&nbsp; Veri tipi: <strong>{e(data_type)}</strong> &nbsp;|&nbsp; Platform: <strong>{e(platform)}</strong></p>

<div class="section"><h2>Modül Durumları (gerçek)</h2><div class="grid">{status_cells}</div>
<p class="note">Her durum modülün kendi çalıştırmasından gelir. NOT_APPLICABLE = analiz koşulları oluşmadı (uydurma değil).</p></div>

<div class="section"><h2>Taksonomik QC (M02)</h2>
{self._na_or(st('M02'), f"Baskın organizma: <strong>{e(tax.get('dominant_organism'))}</strong> (taxid {e(tax.get('taxonomy_id'))}) | Kontaminasyon tahmini: <strong>{e(tax.get('contamination_percent'))}%</strong>", reason('M02'))}</div>

<div class="section"><h2>Assembly İstatistikleri (M03/M04)</h2>
Genom boyutu: <strong>{e(gs.get('genome_size_bp'))} bp</strong> | Contig: <strong>{e(gs.get('contig_count'))}</strong> |
N50: <strong>{e(gs.get('n50'))}</strong> | En uzun: <strong>{e(gs.get('largest_contig'))}</strong> | GC: <strong>{e(gs.get('gc_percent'))}%</strong>
<br><span class="note">Polishing: {e(gs.get('polishing_note'))}</span></div>

<div class="section"><h2>Genom Kalitesi — CheckM2 (M04)</h2>{checkm_html}</div>

<div class="section"><h2>Tür & En Yakın Referanslar (M05)</h2>{c5_html}</div>

<div class="section"><h2>Suş Tiplemesi — MLST (M07)</h2>{mlst_html}</div>

<div class="section"><h2>Antimikrobiyal Direnç (M08)</h2>{amr_html}</div>

<div class="section"><h2>Virülans (M09)</h2><p>{vir_html}</p></div>
<div class="section"><h2>Plazmidler (M10)</h2><p>{plas_html}</p></div>
<div class="section"><h2>Mobil Genetik Elemanlar (M11)</h2><p>{mge_html}</p></div>
<div class="section"><h2>Phage / CRISPR / Defense (M12)</h2><p>{phage_html}</p></div>
<div class="section"><h2>Varyantlar (M13)</h2><p>{var_html}</p></div>

<div class="section"><h2>Araç Referansları & DOI</h2>
<p class="note">Bu tablo platformun araç kayıt defteridir (doğrulanmış DOI'ler). Bu çalıştırmada fiilen
kullanılan araçların sürümleri her modülün <code>*.provenance.json</code> dosyalarındadır.</p>
<table><thead><tr><th>Araç</th><th>Sürüm</th><th>Amaç</th><th>DOI</th></tr></thead><tbody>{ref_rows}</tbody></table></div>
</body></html>"""
