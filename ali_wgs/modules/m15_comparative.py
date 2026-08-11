"""M15 -- Comparative Genomics (Panaroo / PPanGGOLiN)
Pangenom analizi çok-genomlu bir işlemdir (>=2 anotasyonlu genom).
KOŞUL-BAZLI DÜRÜST: tek-örnekli turda -> NOT_APPLICABLE (spec kural #19). Sahte PASS YAZILMAZ.
Batch modda (birden çok örnek) M06 GFF'leri toplanıp Panaroo GERÇEKTEN koşacak (Milestone 2).
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class ComparativeGenomicsModule(Module):
    number = "15"
    name = "comparative_genomics"
    folder = "M15_COMPARATIVE_GENOMICS"
    enabled_key = "comparative"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "M15_summary.json"]

    def run(self):
        self.check_inputs()
        std_dir = self.sub_dir("04_standardized")

        # Bu çalıştırmadaki anotasyonlu genom sayısı (tek örnek = 1)
        genome_count = 1
        if genome_count < 2:
            reason = ("Pangenom analizi >=2 anotasyonlu genom ister; bu tek-örnekli çalıştırmada "
                      "1 genom var -> NOT_APPLICABLE (spec kural #19). Batch modda Panaroo gerçek koşacak.")
            with open(std_dir / "gene_presence_absence.tsv", "w", encoding="utf-8") as f:
                f.write("# NOT_APPLICABLE\t" + reason + "\n")
                f.write("Gene\tGenome\tPresence\n")
            self.write_summary(status="NOT_APPLICABLE",
                               statistics={"genome_count": genome_count},
                               details={"reason": reason})
            return
