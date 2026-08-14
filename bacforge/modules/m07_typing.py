"""M07 -- Strain Typing & Species-Specific Plugins
General MLST (mlst tool) and cgMLST.
Species-specific plugins: Kleborate & Kaptive (Klebsiella), ECTyper (E. coli), SISTR (Salmonella).
Plugins are additive; general pipeline is never skipped.
Outputs: mlst_summary.tsv, cgmlst_summary.tsv, species_plugins.json, M07_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class StrainTypingModule(Module):
    number = "07"
    name = "strain_typing"
    folder = "M07_STRAIN_TYPING"
    enabled_key = "mlst"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "mlst_summary.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV

        # 1. General MLST
        mlst_file = std_dir / "mlst_summary.tsv"
        r.run("mlst", ["mlst", str(genome)], conda_env=E.get("typing", "base"),
              version_cmd=["mlst", "--version"], stdout_path=str(mlst_file), check=False)

        if not mlst_file.exists() or mlst_file.stat().st_size == 0:
            print("MLST failed or produced no output.")
            mlst_status = "FAIL"
        else:
            mlst_status = "PASS"

        # 3. Species-Specific Plugin Evaluation
        species = util.resolve_species(self.ctx) or "Unknown"  # resume'da M02 dosyasından da çözer
        plugin_results = {}

        species_l = (species or "").lower()
        if "klebsiella" in species_l:
            kleb_file = std_dir / "kleborate_results.tsv"
            prov = r.run("kleborate", ["kleborate", "--st258_subtyping", "-a", str(genome), "-o", str(kleb_file)],
                         conda_env=E.get("kleborate", "base"), version_cmd=["kleborate", "--version"], check=False)
            plugin_results["Kleborate"] = "COMPLETED" if prov.get("exit_code") == 0 and kleb_file.exists() else "FAIL"
        elif "escherichia" in species_l or "coli" in species_l:
            # ECTyper bu turda kurulu değil -> dürüstçe atla (Milestone 2). base'de aramaya çalışıp
            # sahte FAIL üretme; env yoksa NOT_AVAILABLE işaretle.
            if E.get("ectyper"):
                ec_dir = std_dir / "ectyper"
                prov = r.run("ectyper", ["ectyper", "-i", str(genome), "-o", str(ec_dir)],
                             conda_env=E["ectyper"], version_cmd=["ectyper", "--version"], check=False)
                plugin_results["ECTyper"] = "COMPLETED" if prov.get("exit_code") == 0 else "FAIL"
            else:
                plugin_results["ECTyper"] = "NOT_AVAILABLE (Milestone 2)"
        elif "salmonella" in species_l:
            plugin_results["SISTR"] = "NOT_AVAILABLE (Milestone 2)"
        else:
            plugin_results["species_specific_module"] = "NOT_APPLICABLE"

        # Kaptive (K/O locus) — Klebsiella (kpsc) ve Acinetobacter baumannii (ab)
        if "klebsiella" in species_l:
            kap_dbs = [("kpsc_k", "K"), ("kpsc_o", "O")]
        elif "acinetobacter" in species_l and "baumannii" in species_l:
            kap_dbs = [("ab_k", "K"), ("ab_o", "OC")]
        else:
            kap_dbs = []
        kaptive_out = {}
        for dbkw, label in kap_dbs:
            kap_tsv = std_dir / f"kaptive_{dbkw}.tsv"
            prov = r.run(f"kaptive_{dbkw}", ["kaptive", "assembly", dbkw, str(genome), "-o", str(kap_tsv)],
                         conda_env=E.get("kaptive", "base"), version_cmd=["kaptive", "--version"], check=False)
            best = None
            if kap_tsv.exists() and kap_tsv.stat().st_size > 0:
                lines = kap_tsv.read_text(encoding="utf-8", errors="replace").splitlines()
                if len(lines) >= 2:
                    hdr = lines[0].split("\t"); ix = {h.strip(): i for i, h in enumerate(hdr)}
                    p = lines[1].split("\t")
                    gv = lambda c: (p[ix[c]] if c in ix and ix[c] < len(p) else "")
                    best = {"locus": gv("Best match locus"), "type": gv("Best match type"),
                            "confidence": gv("Match confidence")}
            if best:
                kaptive_out[label] = best
                plugin_results[f"Kaptive_{label}"] = f"{best['type'] or best['locus']} ({best['confidence']})"
            else:
                plugin_results[f"Kaptive_{label}"] = f"FAIL (exit {prov.get('exit_code')})"
        if kap_dbs:
            with open(std_dir / "kaptive.json", "w", encoding="utf-8") as fh:
                json.dump(kaptive_out, fh, indent=2, ensure_ascii=False)

        with open(std_dir / "species_plugins.json", "w", encoding="utf-8") as fh:
            json.dump(plugin_results, fh, indent=2)

        final_status = "PASS" if mlst_status == "PASS" else "WARNING"

        self.write_summary(
            status=final_status,
            statistics={"plugins_executed": list(plugin_results.keys())},
            details={"mlst_status": mlst_status, **plugin_results}
        )
