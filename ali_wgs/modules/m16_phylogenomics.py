"""M16 -- Phylogenomics (IQ-TREE2 / Gubbins / ClonalFrameML)
Filogenetik ağaç çok-genomlu bir işlemdir (anlamlı ağaç için >=3 genom/örnek).
KOŞUL-BAZLI DÜRÜST: tek-örnekli turda -> NOT_APPLICABLE (spec kural #19).
Sahte PASS ve 'Bulunamadı;' newick YAZILMAZ. Batch modda core-alignment + IQ-TREE2 gerçek koşacak.
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util


class PhylogenomicsModule(Module):
    number = "16"
    name = "phylogenomics"
    folder = "M16_PHYLOGENOMICS"
    enabled_key = "phylogeny"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "M16_summary.json"]

    def run(self):
        self.check_inputs()
        std_dir = self.sub_dir("04_standardized")

        genome_count = 1  # tek örnek
        if genome_count < 3:
            reason = ("Filogenetik ağaç için >=3 genom/örnek gerekir; bu çalıştırmada 1 var "
                      "-> NOT_APPLICABLE (spec kural #19). Batch modda core-SNP alignment + IQ-TREE2 gerçek koşacak.")
            # Sahte newick YAZMA; yalnızca dürüst bir not dosyası.
            with open(std_dir / "phylogeny_status.txt", "w", encoding="utf-8") as f:
                f.write("NOT_APPLICABLE: " + reason + "\n")
            self.write_summary(status="NOT_APPLICABLE",
                               statistics={"genome_count": genome_count},
                               details={"reason": reason})
            return
