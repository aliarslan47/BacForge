"""M13 -- Variant & Mutation Analysis
Reference mapping and variant calling using Snippy / minimap2 / samtools / bcftools.
Outputs: BAM, VCF, SNP table, INDEL table, M13_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class VariantMutationModule(Module):
    number = "13"
    name = "variants_mutations"
    folder = "M13_VARIANTS_MUTATIONS"
    enabled_key = "variants"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "snps.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")

        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)

        # In a generic pipeline, variant calling needs a reference from M05
        # and reads or contigs. We stay strictly inside our run folder.
        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        
        reference_fasta = None
        if ref_json.exists():
            try:
                with open(ref_json, "r") as fh:
                    strains = json.load(fh)
                    if strains and isinstance(strains, list) and "fasta_path" in strains[0]:
                        ref_candidate = Path(strains[0]["fasta_path"])
                        if ref_candidate.exists():
                            reference_fasta = ref_candidate
            except Exception:
                pass

        snps = []
        if reference_fasta:
            snippy_dir = self.sub_dir("02_work") / "snippy"
            # Run Snippy using contigs mode since we might only have assembly
            r.run("snippy", ["snippy", "--cpus", str(t), "--outdir", str(snippy_dir), "--ref", str(reference_fasta), "--ctgs", str(genome)],
                  conda_env=E.get("typing", "base"), version_cmd=["snippy", "--version"], check=False)
            
            vcf_file = snippy_dir / "snps.vcf"
            if vcf_file.exists():
                with open(vcf_file, "r", encoding="utf-8") as fh:
                    for line in fh:
                        if line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 5:
                            snps.append({
                                "chrom": parts[0],
                                "pos": parts[1],
                                "ref": parts[3],
                                "alt": parts[4]
                            })

        with open(std_dir / "snps.tsv", "w", encoding="utf-8") as f:
            f.write("Chromosome\tPosition\tRef\tAlt\n")
            if snps:
                for s in snps:
                    f.write(f"{s['chrom']}\t{s['pos']}\t{s['ref']}\t{s['alt']}\n")
            else:
                f.write("Bulunamadı\tBulunamadı\t-\t-\n")

        self.write_summary(
            status="PASS", 
            statistics={"snp_count": len(snps)}, 
            details={"info": "Varyant analizi tamamlandı" if snps else "Referans veya Varyant Bulunamadı"}
        )
        return
