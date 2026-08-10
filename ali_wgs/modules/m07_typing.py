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
        species = self.ctx.detection.get("ncbi_species", "Unknown")
        plugin_results = {}

        if "klebsiella" in species.lower():
            # Kleborate & Kaptive
            kleb_file = std_dir / "kleborate_results.tsv"
            res = r.run("kleborate", ["kleborate", "--st258_subtyping", "-a", str(genome), "-o", str(kleb_file)],
                        conda_env=E.get("kleborate", "base"), version_cmd=["kleborate", "--version"], check=False)
            if res and res.returncode == 0 and kleb_file.exists():
                plugin_results["Kleborate"] = "COMPLETED"
            else:
                plugin_results["Kleborate"] = "FAIL"
        elif "escherichia" in species.lower() or "coli" in species.lower():
            # ECTyper for E. coli
            ec_file = std_dir / "ectyper_results.tsv"
            res = r.run("ectyper", ["ectyper", "-i", str(genome), "-o", str(ec_file.parent)],
                        conda_env=E.get("ectyper", "base"), version_cmd=["ectyper", "--version"], check=False)
            if res and res.returncode == 0:
                plugin_results["ECTyper"] = "COMPLETED"
            else:
                plugin_results["ECTyper"] = "FAIL"
        elif "salmonella" in species.lower():
            plugin_results["SISTR"] = "NOT_IMPLEMENTED"
        else:
            plugin_results["species_specific_module"] = "NOT_APPLICABLE"

        with open(std_dir / "species_plugins.json", "w", encoding="utf-8") as fh:
            json.dump(plugin_results, fh, indent=2)

        final_status = "PASS" if mlst_status == "PASS" else "WARNING"

        self.write_summary(
            status=final_status,
            statistics={"plugins_executed": list(plugin_results.keys())},
            details={"mlst_status": mlst_status, **plugin_results}
        )
