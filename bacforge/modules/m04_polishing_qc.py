"""M04 -- Polishing, Assembly QC & Genome Quality
QUAST + CheckM2 ile assembly kalitesi (completeness/contamination, N50, GC...).
KURAL: Sonuç UYDURULMAZ. CheckM2 çalışmazsa completeness=None + WARNING; sabit değer yazılmaz.
NOT: Bu turda gerçek polishing (Medaka/Racon/Polypolish) YOK; draft doğrudan genome.fasta olur.
     Bu durum gizlenmez, summary'de polishing_performed=false olarak raporlanır (Milestone 2'de eklenecek).
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .base import Module
from .. import util


class PolishingGenomeQCModule(Module):
    number = "04"
    name = "polishing_genome_qc"
    folder = "M04_POLISHING_GENOME_QC"
    enabled_key = "assembly_qc"

    def inputs(self):
        return [self.ctx.run_dir / "M03_GENOME_ASSEMBLY" / "draft_genome.fasta"]

    def outputs(self):
        return [self.out_dir / "genome.fasta"]

    def run(self):
        self.check_inputs()
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        draft = self.ctx.run_dir / "M03_GENOME_ASSEMBLY" / "draft_genome.fasta"

        std_dir = self.sub_dir("04_standardized")
        final_genome = std_dir / "genome.fasta"

        # Polishing durumu M03'ten okunur (dürüst raporlama; LONG_READ'de Medaka çalışır).
        m03_summary = self.ctx.run_dir / "M03_GENOME_ASSEMBLY" / "M03_summary.json"
        polishing_performed = False
        polisher = None
        if m03_summary.exists():
            try:
                m03d = json.loads(m03_summary.read_text()).get("details", {})
                polishing_performed = bool(m03d.get("polishing_performed", False))
                polisher = m03d.get("polisher")
            except Exception:
                pass

        # Cilalama M03'te yapıldı; M04 draft'ı (LONG'da polished, diğerinde assembler çıktısı) kopyalar.
        shutil.copy(draft, final_genome)

        # 1. QUAST
        quast_dir = self.sub_dir("02_work") / "quast"
        r.run("quast", ["quast.py", str(final_genome), "-o", str(quast_dir), "-t", str(t)],
              conda_env=E["quast"], version_cmd=["quast.py", "--version"], check=False)
        quast_report = quast_dir / "transposed_report.tsv"
        if quast_report.exists():
            shutil.copy(quast_report, std_dir / "quast_summary.tsv")

        # 2. CheckM2 -- GERÇEK completeness/contamination (uydurma yok)
        checkm2_dir = self.sub_dir("02_work") / "checkm2"
        # DB diamond dosyası (gerçek konum): databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
        checkm2_db_dir = Path(self.ctx.config["paths"]["db"]) / "checkm2" / "CheckM2_database"
        checkm2_db = None
        if checkm2_db_dir.exists():
            dmnds = sorted(checkm2_db_dir.glob("*.dmnd"))
            if dmnds:
                checkm2_db = dmnds[0]

        completeness = None
        contamination = None
        checkm2_note = None
        if checkm2_db is not None:
            prov = r.run("checkm2", [
                "checkm2", "predict", "--threads", str(t),
                "--input", str(final_genome),
                "--output-directory", str(checkm2_dir),
                "--database_path", str(checkm2_db), "--force"
            ], conda_env=E["checkm2"], version_cmd=["checkm2", "--version"],
               db_version=str(checkm2_db), check=False)

            res_tsv = checkm2_dir / "quality_report.tsv"
            if prov.get("exit_code") == 0 and res_tsv.exists():
                lines = res_tsv.read_text().splitlines()
                if len(lines) > 1:
                    parts = lines[1].split("\t")
                    if len(parts) >= 3:
                        try:
                            completeness = float(parts[1])
                            contamination = float(parts[2])
                        except ValueError:
                            checkm2_note = "CheckM2 çıktısı parse edilemedi."
            else:
                checkm2_note = f"CheckM2 başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"
        else:
            checkm2_note = f"CheckM2 DB bulunamadı: {checkm2_db_dir}. Completeness/contamination hesaplanamadı."

        # 3. Genom istatistikleri (her zaman gerçek, FASTA'dan hesaplanır)
        seqs = util.read_fasta(final_genome)
        contig_count = len(seqs)
        total_bp = sum(len(s) for s in seqs.values())
        gc_content = round(sum(s.count("G") + s.count("C") for s in seqs.values()) / max(total_bp, 1) * 100, 2)
        lengths = sorted((len(s) for s in seqs.values()), reverse=True)
        half, acc, n50, l50 = total_bp / 2, 0, 0, 0
        for i, l in enumerate(lengths, 1):
            acc += l
            if acc >= half:
                n50, l50 = l, i
                break
        largest = lengths[0] if lengths else 0

        # Kalite durumu: SADECE gerçek CheckM2 değeri varsa PASS/WARNING; yoksa UNKNOWN + WARNING
        if completeness is not None and contamination is not None:
            quality_status = "PASS" if (completeness >= 90 and contamination <= 5) else "WARNING"
        else:
            quality_status = "UNKNOWN"

        checkm2_summary = {
            "completeness": completeness,
            "contamination": contamination,
            "quality_status": quality_status,
            "note": checkm2_note,
        }
        with open(std_dir / "checkm2_summary.json", "w", encoding="utf-8") as fh:
            json.dump(checkm2_summary, fh, indent=2, ensure_ascii=False)

        genome_stats = {
            "genome_size_bp": total_bp,
            "contig_count": contig_count,
            "n50": n50,
            "l50": l50,
            "largest_contig": largest,
            "gc_percent": gc_content,
            "completeness_percent": completeness,
            "contamination_percent": contamination,
            "quality_status": quality_status,
            "polishing_performed": polishing_performed,
            "polisher": polisher,
            "polishing_note": (
                f"Polisher: {polisher}" if polishing_performed
                else "Polishing uygulanmadi (short/hybrid assembler dahili cilalar; ya da uzun-okuma polishing basarisiz)."
            ),
        }
        with open(std_dir / "genome_stats.json", "w", encoding="utf-8") as fh:
            json.dump(genome_stats, fh, indent=2, ensure_ascii=False)

        warnings = []
        if checkm2_note:
            warnings.append(checkm2_note)
        # Modül durumu: CheckM2 gerçek değer verdiyse onun kararı; vermediyse WARNING (dürüst)
        module_status = quality_status if quality_status in ("PASS", "WARNING") else "WARNING"
        self.write_summary(status=module_status, statistics=genome_stats, warnings=warnings)
