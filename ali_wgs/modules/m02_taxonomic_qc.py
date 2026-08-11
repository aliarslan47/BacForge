"""M02 -- Taxonomic QC
Kraken2 ile taksonomik sınıflandırma; raporu GERÇEKTEN parse eder, baskın türü ve
kontaminasyon tahminini çıkarır ve türü ctx.detection['ncbi_species'] olarak downstream'e taşır.
KURAL: Sonuç uydurulmaz. DB yoksa NOT_APPLICABLE; kraken2 çalışıp bir şey bulamazsa WARNING.
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


def parse_kraken2_report(report_path: Path) -> dict:
    """Kraken2 raporundan baskın tür + kontaminasyon tahmini.
    Sütunlar: pct, clade_reads, direct_reads, rank_code, taxid, name(indented).
    Baskın tür = en yüksek clade_reads'a sahip S (species) satırı.
    """
    species_rows = []
    total_reads = 0
    unclassified_reads = 0
    with open(report_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                pct = float(parts[0].strip())
                clade_reads = int(parts[1].strip())
            except ValueError:
                continue
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name = parts[5].strip()
            if rank == "U":
                unclassified_reads = clade_reads
            if rank in ("R", "U") and taxid in ("1", "0"):
                total_reads += clade_reads  # root + unclassified toplamı = tüm okumalar
            if rank == "S":
                species_rows.append({"pct": pct, "reads": clade_reads, "taxid": taxid, "name": name})

    if not species_rows:
        return {"found": False}

    species_rows.sort(key=lambda r: r["reads"], reverse=True)
    dominant = species_rows[0]
    species_total = sum(r["reads"] for r in species_rows)
    # Kontaminasyon: baskın tür DIŞINDAKİ tür-düzeyi okumaların, tüm tür-düzeyi
    # okumalar içindeki payı (kraken2 tabanlı tahmin).
    contamination = round(100.0 * (species_total - dominant["reads"]) / species_total, 2) if species_total else 0.0
    return {
        "found": True,
        "dominant_organism": dominant["name"],
        "taxonomy_id": dominant["taxid"],
        "dominant_reads": dominant["reads"],
        "dominance_percent_of_total": dominant["pct"],
        "dominance_percent_of_classified_species": round(100.0 * dominant["reads"] / species_total, 2) if species_total else 0.0,
        "contamination_percent": contamination,
        "species_ranked": species_rows[:10],
    }


class TaxonomicQCModule(Module):
    number = "02"
    name = "taxonomic_qc"
    folder = "M02_TAXONOMIC_QC"
    enabled_key = "taxonomy"

    def inputs(self):
        return [self.ctx.run_dir / "M00_INPUT_AUTO_DETECTION" / "data_type.json"]

    def outputs(self):
        return [self.out_dir / "M02_summary.json"]

    def run(self):
        self.check_inputs()
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        data_type = self.ctx.detection.get("data_type", "SHORT_READ")

        kraken_db = Path(self.ctx.config["paths"]["db"]) / "kraken2"
        clean_r1 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R1.fastq.gz"
        clean_r2 = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "clean_R2.fastq.gz"
        filtered_long = self.ctx.run_dir / "M01_READ_QC_PREPROCESSING" / "filtered_long.fastq.gz"

        kraken_report = self.sub_dir("03_native_outputs") / "kraken2_report.txt"
        kraken_out = self.sub_dir("03_native_outputs") / "kraken2.txt"

        # DB yoksa: sessizce PASS yazma -> NOT_APPLICABLE (dürüst)
        if not kraken_db.exists():
            self.write_summary(
                status="NOT_APPLICABLE",
                details={"reason": f"Kraken2 veritabanı yok: {kraken_db}. Taksonomik QC atlandı."},
            )
            self._write_empty_std(std_dir, reason="Kraken2 DB yok")
            return

        # Girdi okumalarını seç (short/long/hybrid)
        if data_type == "LONG_READ" and filtered_long.exists():
            reads_args = [str(filtered_long)]
        elif clean_r1.exists() and clean_r2.exists():
            reads_args = ["--paired", str(clean_r1), str(clean_r2)]
        elif clean_r1.exists():
            reads_args = [str(clean_r1)]
        else:
            self.write_summary(status="WARNING",
                               warnings=["M01 temiz okuma çıktısı bulunamadı; taksonomi çalıştırılamadı."])
            self._write_empty_std(std_dir, reason="Temiz okuma yok")
            return

        cmd = ["kraken2", "--db", str(kraken_db), "--threads", str(util.threads(self.ctx)),
               "--report", str(kraken_report), "--output", str(kraken_out)] + reads_args
        prov = r.run("kraken2", cmd, conda_env=util.ENV.get("qc", "base"),
                     version_cmd=["kraken2", "--version"], db_version=str(kraken_db), check=False)

        if prov.get("exit_code") != 0 or not kraken_report.exists():
            self.write_summary(status="WARNING",
                               warnings=[f"Kraken2 başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
            self._write_empty_std(std_dir, reason="Kraken2 hata")
            return

        parsed = parse_kraken2_report(kraken_report)
        if not parsed.get("found"):
            self.write_summary(status="WARNING",
                               warnings=["Kraken2 çalıştı ancak tür (S) düzeyinde sınıflandırma yok."])
            self._write_empty_std(std_dir, reason="Tür sınıflandırması yok")
            return

        # Bracken: tür-düzeyi abundansı yeniden tahmin et (küçük kraken DB'nin dağıtımını düzeltir)
        bracken = self._run_bracken(kraken_db, kraken_report, r)

        if bracken:
            species = bracken["name"]
            taxid = bracken["taxonomy_id"]
            dominance = bracken["fraction"] * 100.0          # tüm okumalar içinde (Bracken)
            contamination = round(100.0 - dominance, 2)
            source = "bracken"
        else:
            species = parsed["dominant_organism"]
            taxid = parsed["taxonomy_id"]
            dominance = parsed["dominance_percent_of_classified_species"]
            contamination = parsed["contamination_percent"]
            source = "kraken2"

        # Türü downstream modüllere taşı (M05/M07/M13/M17 bunu kullanır)
        self.ctx.detection["ncbi_species"] = species
        self.ctx.detection["ncbi_taxid"] = taxid

        tax_data = {
            "dominant_organism": species,
            "taxonomy_id": taxid,
            "dominance_percent": round(dominance, 2),
            "contamination_percent": contamination,
            "source": source,
            "hybrid_concordance": "N/A",
        }

        with open(std_dir / "taxonomy.json", "w", encoding="utf-8") as fh:
            json.dump(tax_data, fh, indent=2, ensure_ascii=False)
        # M05'in okuduğu tür dosyası
        with open(std_dir / "species_identification.json", "w", encoding="utf-8") as fh:
            json.dump({"species": species, "taxid": taxid,
                       "source": source}, fh, indent=2, ensure_ascii=False)
        with open(std_dir / "taxonomy.tsv", "w", encoding="utf-8") as fh:
            fh.write("Taxonomy_ID\tOrganism\tPercent_of_total\tClade_reads\n")
            for row in parsed["species_ranked"]:
                fh.write(f"{row['taxid']}\t{row['name']}\t{row['pct']}\t{row['reads']}\n")
        with open(std_dir / "contamination.tsv", "w", encoding="utf-8") as fh:
            fh.write("Metric\tValue\n")
            fh.write(f"Dominant_organism\t{species}\n")
            fh.write(f"Dominance_percent\t{round(dominance, 2)}\n")
            fh.write(f"Contamination_percent\t{contamination}\n")

        # Bracken abundansı ile temizlik kararı: baskın tür >= %90 -> PASS
        clean = dominance >= 90.0
        self.write_summary(
            status="PASS" if clean else "WARNING",
            statistics=tax_data,
            warnings=[] if clean else [f"Olası kontaminasyon: baskın tür ({species}) "
                                       f"okumaların yalnızca %{round(dominance, 2)}'i ({source})."],
            details={"note": f"Baskın tür {source} ile belirlendi (Bracken tür-düzeyi abundans re-estimasyonu)."},
        )

    def _run_bracken(self, kraken_db, kraken_report, r) -> dict | None:
        """Bracken ile tür-düzeyi abundans; en baskın türü döndürür (name, taxid, fraction)."""
        # okuma uzunluğuna en yakın kmer DB'sini seç (150 varsayılan)
        rlen = "150"
        out_tsv = self.sub_dir("03_native_outputs") / "bracken_species.tsv"
        brep = self.sub_dir("03_native_outputs") / "bracken_report.txt"
        prov = r.run("bracken", ["bracken", "-d", str(kraken_db), "-i", str(kraken_report),
                                 "-o", str(out_tsv), "-w", str(brep), "-r", rlen, "-l", "S", "-t", "1"],
                     conda_env=util.ENV.get("qc", "base"), version_cmd=["bracken", "--version"], check=False)
        if prov.get("exit_code") != 0 or not out_tsv.exists():
            return None
        best = None
        with open(out_tsv, encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                i_name = header.index("name"); i_tax = header.index("taxonomy_id")
                i_frac = header.index("fraction_total_reads")
            except ValueError:
                return None
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(i_name, i_tax, i_frac):
                    continue
                try:
                    frac = float(p[i_frac])
                except ValueError:
                    continue
                if best is None or frac > best["fraction"]:
                    best = {"name": p[i_name], "taxonomy_id": p[i_tax], "fraction": frac}
        return best

    def _write_empty_std(self, std_dir, reason: str):
        empty = {"dominant_organism": None, "taxonomy_id": None, "dominance_percent": None,
                 "contamination_percent": None, "reason": reason}
        with open(std_dir / "taxonomy.json", "w", encoding="utf-8") as fh:
            json.dump(empty, fh, indent=2, ensure_ascii=False)
