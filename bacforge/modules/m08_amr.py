"""M08 -- Antimicrobial Resistance
Tools: NCBI AMRFinderPlus, CARD/RGI, ResFinder, ABRicate.
Standardized outputs: amr_genes.tsv, amr_mutations.tsv, amr_proteins.tsv, amr.json, M08_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class AMRModule(Module):
    number = "08"
    name = "antimicrobial_resistance"
    folder = "M08_AMR"
    enabled_key = "amr"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "amr_genes.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)
        dbp = self.ctx.config["paths"]["db"]

        # 1. AMRFinderPlus
        amr_out = self.sub_dir("03_native_outputs") / "amrfinder.tsv"
        amr_db = Path(dbp) / "amrfinderplus" / "latest"
        cmd = ["amrfinder", "-n", str(genome), "--plus", "--threads", str(t)]
        if amr_db.exists():
            cmd.extend(["--database", str(amr_db)])

        prov_amr = r.run("amrfinder", cmd, conda_env=E["amrfinder"], version_cmd=["amrfinder", "--version"], stdout_path=str(amr_out), check=False)

        # Parse AMRFinder or construct standardized outputs
        amr_genes = []
        amr_mutations = []

        if amr_out.exists() and amr_out.stat().st_size > 0:
            with open(amr_out, "r", encoding="utf-8") as fh:
                header = fh.readline().strip().split("\t")
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 11:
                        elem_type = parts[8] if len(parts) > 8 else "AMR"
                        gene_symbol = parts[5] if len(parts) > 5 else "unknown"
                        drug_class = parts[10] if len(parts) > 10 else "Unknown"
                        subclass = parts[11] if len(parts) > 11 else drug_class

                        if elem_type.upper() == "AMR":
                            amr_genes.append({
                                "gene_symbol": gene_symbol,
                                "element_type": elem_type,
                                "drug_class": drug_class,
                                "subclass": subclass,
                                "coverage": parts[15] if len(parts) > 15 else "100.0",
                                "identity": parts[16] if len(parts) > 16 else "100.0",
                                "contig": parts[1] if len(parts) > 1 else "contig1",
                                "start": parts[2] if len(parts) > 2 else "1",
                                "end": parts[3] if len(parts) > 3 else "1000",
                                "strand": parts[4] if len(parts) > 4 else "+"
                            })

        for g in amr_genes:
            g.setdefault("source", "AMRFinderPlus")
        if not amr_genes:
            print("[M08] WARNING: AMRFinderPlus çıktısı yok/boş.")

        def _mk(gene, drug, sub, ident, cov, contig, source):
            return {"gene_symbol": gene or "", "element_type": "AMR", "drug_class": drug or "",
                    "subclass": sub or "", "coverage": cov or "", "identity": ident or "",
                    "contig": contig or "", "start": "", "end": "", "strand": "", "source": source}

        def _tsv_cols(path):
            if not (path.exists() and path.stat().st_size > 0):
                return None, []
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            if not lines:
                return None, []
            hdr = lines[0].lstrip("#").split("\t")
            idx = {h.strip(): i for i, h in enumerate(hdr)}
            return idx, [ln.split("\t") for ln in lines[1:] if ln.strip()]

        # 2. RGI (CARD)
        rgi_genes = []
        rgi_prefix = self.sub_dir("03_native_outputs") / "rgi"
        prov_rgi = r.run("rgi", ["rgi", "main", "-i", str(genome), "-o", str(rgi_prefix),
                                 "-t", "contig", "-n", str(t), "--clean"],
                         conda_env=E.get("rgi", "base"), version_cmd=["rgi", "main", "--version"], check=False)
        idx, rows = _tsv_cols(Path(str(rgi_prefix) + ".txt"))
        if idx:
            gv = lambda p, c: (p[idx[c]] if c in idx and idx[c] < len(p) else "")
            for p in rows:
                rgi_genes.append(_mk(gv(p, "Best_Hit_ARO"), gv(p, "Drug Class"), gv(p, "Resistance Mechanism"),
                                     gv(p, "Best_Identities"), "", gv(p, "Contig"), "RGI/CARD"))

        # 3. ResFinder (abricate --db resfinder)
        resf_genes = []
        resf_out = self.sub_dir("03_native_outputs") / "abricate_resfinder.tsv"
        prov_rf = r.run("abricate_resfinder", ["abricate", "--db", "resfinder", "--nopath", str(genome)],
                        conda_env=E.get("virulence", "base"), version_cmd=["abricate", "--version"],
                        stdout_path=str(resf_out), check=False)
        idx, rows = _tsv_cols(resf_out)
        if idx:
            gv = lambda p, c: (p[idx[c]] if c in idx and idx[c] < len(p) else "")
            for p in rows:
                resf_genes.append(_mk(gv(p, "GENE"), gv(p, "RESISTANCE") or gv(p, "PRODUCT"), gv(p, "PRODUCT"),
                                      gv(p, "%IDENTITY"), gv(p, "%COVERAGE"), gv(p, "SEQUENCE"), "ResFinder"))

        # Birleşik, kaynak-etiketli liste (M17/M18 bunu kullanır)
        amr_all = amr_genes + rgi_genes + resf_genes
        sources = {"AMRFinderPlus": len(amr_genes), "RGI/CARD": len(rgi_genes), "ResFinder": len(resf_genes)}
        ran = {"AMRFinderPlus": prov_amr.get("exit_code") == 0,
               "RGI/CARD": prov_rgi.get("exit_code") == 0,
               "ResFinder": prov_rf.get("exit_code") == 0}

        with open(std_dir / "amr_genes.tsv", "w", encoding="utf-8") as fh:
            fh.write("Gene_Symbol\tDrug_Class\tSubclass\tCoverage\tIdentity\tContig\tSource\n")
            for g in amr_all:
                fh.write(f"{g['gene_symbol']}\t{g['drug_class']}\t{g['subclass']}\t{g['coverage']}\t{g['identity']}\t{g['contig']}\t{g['source']}\n")
        with open(std_dir / "amr_mutations.tsv", "w", encoding="utf-8") as fh:
            fh.write("Gene_Symbol\tMutation\tResistance_Mechanism\tDrug_Class\n")
        with open(std_dir / "amr_proteins.tsv", "w", encoding="utf-8") as fh:
            fh.write("Protein_ID\tGene_Symbol\tProduct\n")
        with open(std_dir / "amr.json", "w", encoding="utf-8") as fh:
            json.dump({"amr_genes": amr_all, "amr_mutations": amr_mutations,
                       "by_source": {"amrfinderplus": amr_genes, "rgi_card": rgi_genes, "resfinder": resf_genes},
                       "source_counts": sources}, fh, indent=2, ensure_ascii=False)

        # Dürüst durum: en az bir araç çalıştıysa PASS; hiçbiri değilse WARNING.
        any_ran = any(ran.values())
        warns = [f"{k} çalışmadı" for k, ok in ran.items() if not ok]
        self.write_summary(
            status="PASS" if any_ran else "WARNING",
            statistics={"amr_gene_count": len(amr_all), "source_counts": sources},
            warnings=warns or None,
            details={"tools": [k for k, ok in ran.items() if ok],
                     "drug_classes": sorted({g["drug_class"] for g in amr_all if g["drug_class"]})}
        )
