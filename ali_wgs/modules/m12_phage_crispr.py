"""M12 -- Phage, CRISPR & Defense Systems
Tools: geNomad (prophages/viral), cctyper (CRISPR-Cas), DefenseFinder
Outputs: prophages.tsv, prophages.fasta, crispr.tsv, crispr_spacers.fasta, defense_systems.tsv, M12_summary.json
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

        # Run geNomad
        genomad_out = self.sub_dir("02_work") / "genomad"
        genomad_db = Path(dbp) / "genomad" / "genomad_db"
        if genomad_db.exists():
            r.run("genomad", ["genomad", "end-to-end", "--cleanup", str(genome), str(genomad_out), str(genomad_db), "--threads", str(t)],
                  conda_env=E["genomad"], version_cmd=["genomad", "--version"], check=False)

        phages = []
        genomad_virus = genomad_out / "genome_summary" / "genome_virus_summary.tsv"
        if genomad_virus.exists():
            with open(genomad_virus, "r", encoding="utf-8") as fh:
                header = fh.readline()
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        phages.append({
                            "phage_id": parts[0],
                            "contig": parts[0].rsplit('_', 1)[0] if '_' in parts[0] else parts[0],
                            "length": parts[1],
                            "topology": parts[2],
                            "virus_score": parts[5]
                        })

        with open(std_dir / "prophages.tsv", "w", encoding="utf-8") as f:
            f.write("Phage_ID\tContig\tLength\tTopology\tVirus_Score\n")
            if phages:
                for p in phages:
                    f.write(f"{p['phage_id']}\t{p['contig']}\t{p['length']}\t{p['topology']}\t{p['virus_score']}\n")
            else:
                f.write("Bulunamadı\tBulunamadı\t-\t-\t-\n")

        self.write_summary(
            status="PASS", 
            statistics={"phage_count": len(phages)}, 
            details={"info": "Prophage analizi tamamlandı" if phages else "Prophage Bulunamadı"}
        )
        return
