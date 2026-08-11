"""M13 -- Variant & Mutation Analysis (Snippy, contigs modu)
Referans = M05 closest-1 (fasta_path). KURAL: 'Bulunamadı' satırı + otomatik PASS YOK.
Referans yoksa NOT_APPLICABLE; snippy çalışıp varyant bulamazsa PASS (count=0); hata -> WARNING.
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

        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        reference_fasta = None
        if ref_json.exists():
            try:
                strains = json.load(open(ref_json))
                if strains and isinstance(strains, list) and strains[0].get("fasta_path"):
                    cand = Path(strains[0]["fasta_path"])
                    if cand.exists():
                        reference_fasta = cand
            except Exception:
                pass

        if reference_fasta is None:
            with open(std_dir / "snps.tsv", "w", encoding="utf-8") as f:
                f.write("Chromosome\tPosition\tRef\tAlt\n")
            self.write_summary(status="NOT_APPLICABLE",
                               details={"reason": "M05 closest-1 referansı yok -> varyant çağrısı yapılamadı."})
            return

        snippy_dir = self.sub_dir("02_work") / "snippy"
        prov = r.run("snippy", ["snippy", "--cpus", str(t), "--outdir", str(snippy_dir), "--force",
                                "--ref", str(reference_fasta), "--ctgs", str(genome)],
                     conda_env=E.get("typing", "base"), version_cmd=["snippy", "--version"], check=False)
        tool_ran = prov.get("exit_code") == 0

        snps = []
        vcf_file = snippy_dir / "snps.vcf"
        if vcf_file.exists():
            with open(vcf_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 5:
                        snps.append({"chrom": parts[0], "pos": parts[1], "ref": parts[3], "alt": parts[4]})

        with open(std_dir / "snps.tsv", "w", encoding="utf-8") as f:
            f.write("Chromosome\tPosition\tRef\tAlt\n")
            for s in snps:
                f.write(f"{s['chrom']}\t{s['pos']}\t{s['ref']}\t{s['alt']}\n")
        with open(std_dir / "variants.json", "w", encoding="utf-8") as f:
            json.dump({"snps": snps, "reference": str(reference_fasta)}, f, indent=2, ensure_ascii=False)

        if tool_ran:
            self.write_summary(status="PASS", statistics={"snp_count": len(snps)},
                               details={"reference": reference_fasta.stem})
        else:
            self.write_summary(status="WARNING", statistics={"snp_count": len(snps)},
                               warnings=[f"Snippy başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
