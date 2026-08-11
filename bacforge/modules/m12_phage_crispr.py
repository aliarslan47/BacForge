"""M12 -- Phage, CRISPR & Defense Systems (geNomad; cctyper/DefenseFinder Milestone 2)
KURAL: 'Bulunamadı' satırı + otomatik PASS YOK. geNomad DB yoksa NOT_APPLICABLE;
geNomad çalışıp faj bulamazsa PASS (count=0); hata -> WARNING.
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class PhageCRISPRDefenseModule(Module):
    number = "12"
    name = "phage_crispr_defense"
    folder = "M12_PHAGE_CRISPR_DEFENSE"
    enabled_key = "genomad"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "prophages.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)
        dbp = self.ctx.config["paths"]["db"]

        genomad_out = self.sub_dir("02_work") / "genomad"
        genomad_db = Path(dbp) / "genomad" / "genomad_db"

        if not genomad_db.exists():
            with open(std_dir / "prophages.tsv", "w", encoding="utf-8") as f:
                f.write("Phage_ID\tContig\tLength\tTopology\tVirus_Score\n")
            self.write_summary(status="NOT_APPLICABLE",
                               details={"reason": f"geNomad DB yok: {genomad_db}"})
            return

        prov = r.run("genomad", ["genomad", "end-to-end", "--cleanup", str(genome), str(genomad_out),
                                 str(genomad_db), "--threads", str(t)],
                     conda_env=E["genomad"], version_cmd=["genomad", "--version"],
                     db_version=str(genomad_db), check=False)
        tool_ran = prov.get("exit_code") == 0

        phages = []
        # geNomad çıktı adı girdi dosya adına göre değişir; virus_summary.tsv'yi ara
        for vsum in genomad_out.rglob("*_virus_summary.tsv"):
            with open(vsum, "r", encoding="utf-8") as fh:
                fh.readline()  # header
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        phages.append({
                            "phage_id": parts[0],
                            "contig": parts[0].rsplit("_", 1)[0] if "_" in parts[0] else parts[0],
                            "length": parts[1],
                            "topology": parts[2],
                            "virus_score": parts[5],
                        })

        with open(std_dir / "prophages.tsv", "w", encoding="utf-8") as f:
            f.write("Phage_ID\tContig\tLength\tTopology\tVirus_Score\n")
            for p in phages:
                f.write(f"{p['phage_id']}\t{p['contig']}\t{p['length']}\t{p['topology']}\t{p['virus_score']}\n")
        with open(std_dir / "prophages.json", "w", encoding="utf-8") as f:
            json.dump({"prophages": phages}, f, indent=2, ensure_ascii=False)

        if tool_ran:
            self.write_summary(status="PASS", statistics={"phage_count": len(phages)},
                               details={"tool": "geNomad", "note": "CRISPR (cctyper) & DefenseFinder: Milestone 2"})
        else:
            self.write_summary(status="WARNING", statistics={"phage_count": len(phages)},
                               warnings=[f"geNomad başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
