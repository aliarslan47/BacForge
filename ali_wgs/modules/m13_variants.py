"""M13 -- Variant & Mutation Analysis (Snippy, GenBank referansli -> CDS/gen efekti)
Referans = M05 closest-1. KRITIK: snippy'ye ANOTASYONLU GenBank verilir -> her varyant icin
GENE / LOCUS_TAG / EFFECT (amino asit degisimi) / PRODUCT cikar (CDS baglami). FASTA ref = anlamsiz.
KATI: referans yoksa NOT_APPLICABLE; snippy calisip varyant bulamazsa PASS(0); hata -> WARNING.
"""
from __future__ import annotations

import json
import shutil
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

    def _ensure_ref_gbk(self, fna: Path, r, dbp, t) -> Path | None:
        """Referansin Bakta GenBank'i (cache: fna yaninda .gbff). Yoksa Bakta --skip-sorf ile uret."""
        gbff = fna.with_suffix(".gbff")
        if gbff.exists() and gbff.stat().st_size > 1000:
            return gbff
        outdir = fna.parent / (fna.stem + "_bakta")
        prov = r.run(f"bakta_ref_{fna.stem[:16]}", [
            "bakta", "--db", f"{dbp}/bakta/db-light", "--output", str(outdir),
            "--prefix", "ref", "--skip-sorf", "--skip-plot", "--threads", str(t), "--force", str(fna)
        ], conda_env=util.ENV["bakta"], check=False)
        src = outdir / "ref.gbff"
        if prov.get("exit_code") == 0 and src.exists() and src.stat().st_size > 1000:
            shutil.copy(src, gbff)
            return gbff
        return None

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        t = util.threads(self.ctx)
        dbp = self.ctx.config["paths"]["db"]

        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        ref_fna = None
        if ref_json.exists():
            try:
                strains = json.load(open(ref_json))
                if strains and Path(strains[0].get("fasta_path", "")).exists():
                    ref_fna = Path(strains[0]["fasta_path"])
            except Exception:
                pass

        if ref_fna is None:
            self._write_empty(std_dir)
            self.write_summary(status="NOT_APPLICABLE",
                               details={"reason": "M05 closest-1 referansı yok -> varyant çağrısı yapılamadı."})
            return

        # ANOTASYONLU GenBank referans (CDS/gen efekti icin sart)
        ref_gbk = self._ensure_ref_gbk(ref_fna, r, dbp, t)
        ref_for_snippy = ref_gbk if ref_gbk else ref_fna  # gbk yoksa fasta (efektsiz) - dürüstçe not düşülür

        snippy_dir = self.sub_dir("02_work") / "snippy"
        prov = r.run("snippy", ["snippy", "--cpus", str(t), "--outdir", str(snippy_dir), "--force",
                                "--ref", str(ref_for_snippy), "--ctgs", str(genome)],
                     conda_env=util.ENV.get("typing", "base"), version_cmd=["snippy", "--version"], check=False)
        tool_ran = prov.get("exit_code") == 0

        # snps.tab: CHROM POS TYPE REF ALT EVIDENCE FTYPE STRAND NT_POS AA_POS EFFECT LOCUS_TAG GENE PRODUCT
        variants = []
        tab = snippy_dir / "snps.tab"
        if tab.exists():
            with open(tab, encoding="utf-8") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(header)}
                def g(p, c):
                    i = idx.get(c)
                    return p[i] if i is not None and i < len(p) else ""
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5:
                        continue
                    variants.append({
                        "contig": g(p, "CHROM"), "pos": g(p, "POS"), "type": g(p, "TYPE"),
                        "ref": g(p, "REF"), "alt": g(p, "ALT"),
                        "effect": g(p, "EFFECT"), "locus_tag": g(p, "LOCUS_TAG"),
                        "gene": g(p, "GENE"), "product": g(p, "PRODUCT"),
                    })

        # standardize TSV (CDS/gen dahil)
        with open(std_dir / "snps.tsv", "w", encoding="utf-8") as f:
            f.write("Contig\tPos\tType\tRef\tAlt\tEffect\tLocus_tag\tGene\tProduct\n")
            for v in variants:
                f.write(f"{v['contig']}\t{v['pos']}\t{v['type']}\t{v['ref']}\t{v['alt']}\t"
                        f"{v['effect']}\t{v['locus_tag']}\t{v['gene']}\t{v['product']}\n")
        with open(std_dir / "variants.json", "w", encoding="utf-8") as f:
            json.dump({"snps": variants, "reference": ref_fna.stem,
                       "reference_annotated": bool(ref_gbk)}, f, indent=2, ensure_ascii=False)

        n_coding = sum(1 for v in variants if v.get("gene") or v.get("locus_tag"))
        if tool_ran:
            warns = [] if ref_gbk else ["Referans GenBank üretilemedi; varyantlar gen/CDS efekti OLMADAN raporlandı."]
            self.write_summary(status="PASS" if ref_gbk else "WARNING",
                               statistics={"variant_count": len(variants), "coding_variants": n_coding,
                                           "reference": ref_fna.stem, "annotated": bool(ref_gbk)},
                               warnings=warns)
        else:
            self.write_summary(status="WARNING", statistics={"variant_count": len(variants)},
                               warnings=[f"Snippy başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])

    def _write_empty(self, std_dir):
        with open(std_dir / "snps.tsv", "w", encoding="utf-8") as f:
            f.write("Contig\tPos\tType\tRef\tAlt\tEffect\tLocus_tag\tGene\tProduct\n")
