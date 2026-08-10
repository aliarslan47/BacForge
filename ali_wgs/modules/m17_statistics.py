"""M17 -- Statistics & Visualization
Central data aggregator gathering statistics from M00--M16 into dashboard_data.json and preparing visual charts.
Outputs: dashboard_data.json, M17_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module


class StatisticsVisualizationModule(Module):
    number = "17"
    name = "statistics_visualization"
    folder = "M17_STATISTICS_VISUALIZATION"
    enabled_key = "stats"

    def inputs(self):
        return [self.ctx.run_dir / "M00_INPUT_AUTO_DETECTION" / "04_standardized" / "data_type.json"]

    def outputs(self):
        return [self.out_dir / "04_standardized" / "dashboard_data.json"]

    def run(self):
        self.check_inputs()
        run_dir = Path(self.ctx.run_dir)
        std_dir = self.sub_dir("04_standardized")

        # Load summaries from all previous modules
        module_summaries = {}
        for sub in run_dir.iterdir():
            if sub.is_dir() and sub.name.startswith("M"):
                mod_num = sub.name.split("_")[0]
                sum_file = sub / f"{mod_num}_summary.json"
                if sum_file.exists():
                    try:
                        with open(sum_file, "r", encoding="utf-8") as fh:
                            module_summaries[mod_num] = json.load(fh)
                    except Exception:
                        pass

        # Load key datasets
        c5_data = []
        c5_file = run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "04_standardized" / "closest_5_strains.json"
        if c5_file.exists():
            with open(c5_file, "r", encoding="utf-8") as fh:
                c5_data = json.load(fh)

        amr_data = []
        amr_file = run_dir / "M08_AMR" / "04_standardized" / "amr.json"
        if amr_file.exists():
            with open(amr_file, "r", encoding="utf-8") as fh:
                amr_data = json.load(fh).get("amr_genes", [])

        vir_data = []
        vir_file = run_dir / "M09_VIRULENCE" / "04_standardized" / "virulence.json"
        if vir_file.exists():
            with open(vir_file, "r", encoding="utf-8") as fh:
                vir_data = json.load(fh).get("virulence_genes", [])

        plasmid_data = []
        plasmid_file = run_dir / "M10_PLASMID" / "04_standardized" / "plasmids.json"
        if plasmid_file.exists():
            with open(plasmid_file, "r", encoding="utf-8") as fh:
                plasmid_data = json.load(fh).get("plasmids", [])

        mge_data = []
        mge_file = run_dir / "M11_MOBILE_GENETIC_ELEMENTS" / "04_standardized" / "mobile_elements.json"
        if mge_file.exists():
            with open(mge_file, "r", encoding="utf-8") as fh:
                mge_data = json.load(fh).get("mobile_elements", [])

        tree_nwk = ""
        tree_file = run_dir / "M16_PHYLOGENOMICS" / "04_standardized" / "tree.nwk"
        if tree_file.exists():
            tree_nwk = tree_file.read_text().strip()

        dashboard = {
            "project_id": run_dir.name,
            "data_type": self.ctx.detection.get("data_type", "SHORT_READ"),
            "platform": self.ctx.detection.get("platform", "Illumina"),
            "species": self.ctx.detection.get("ncbi_species", "Klebsiella pneumoniae"),
            "modules": module_summaries,
            "closest_5_strains": c5_data,
            "amr_genes": amr_data,
            "virulence_genes": vir_data,
            "plasmids": plasmid_data,
            "mobile_elements": mge_data,
            "phylogenetic_tree_newick": tree_nwk,
        }

        with open(std_dir / "dashboard_data.json", "w", encoding="utf-8") as fh:
            json.dump(dashboard, fh, indent=2, ensure_ascii=False)

        self.write_summary(
            status="PASS",
            statistics={"modules_aggregated": len(module_summaries)}
        )
