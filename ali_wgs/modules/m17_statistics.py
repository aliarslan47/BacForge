"""M17 -- Statistics & Visualization
M00--M16 GERÇEK çıktılarını merkezî dashboard_data.json'a toplar. Hiçbir değer uydurulmaz;
üretilmemiş veriler None/boş kalır. Tür M02'den gelir (sabit varsayım yok).
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module


def _load_json(p: Path, default=None):
    if p.exists():
        try:
            return json.load(open(p, encoding="utf-8"))
        except Exception:
            return default
    return default


class StatisticsVisualizationModule(Module):
    number = "17"
    name = "statistics_visualization"
    folder = "M17_STATISTICS_VISUALIZATION"
    enabled_key = "stats"

    def inputs(self):
        return [self.ctx.run_dir / "M00_INPUT_AUTO_DETECTION" / "data_type.json"]

    def outputs(self):
        return [self.out_dir / "dashboard_data.json"]

    def run(self):
        self.check_inputs()
        run_dir = Path(self.ctx.run_dir)
        std_dir = self.sub_dir("04_standardized")

        # Modül özet durumları
        module_summaries = {}
        for sub in run_dir.iterdir():
            if sub.is_dir() and sub.name.startswith("M"):
                mod_num = sub.name.split("_")[0]
                s = _load_json(sub / f"{mod_num}_summary.json")
                if s is not None:
                    module_summaries[mod_num] = s

        taxonomy = _load_json(run_dir / "M02_TAXONOMIC_QC" / "taxonomy.json", {})
        genome_stats = _load_json(run_dir / "M04_POLISHING_GENOME_QC" / "genome_stats.json", {})
        checkm2 = _load_json(run_dir / "M04_POLISHING_GENOME_QC" / "checkm2_summary.json", {})
        c5 = _load_json(run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json", [])
        amr = (_load_json(run_dir / "M08_AMR" / "amr.json", {}) or {}).get("amr_genes", [])
        vir = (_load_json(run_dir / "M09_VIRULENCE" / "virulence.json", {}) or {}).get("virulence_genes", [])
        plasmids = (_load_json(run_dir / "M10_PLASMID" / "plasmids.json", {}) or {}).get("plasmids", [])
        mges = (_load_json(run_dir / "M11_MOBILE_GENETIC_ELEMENTS" / "mobile_elements.json", {}) or {}).get("mobile_elements", [])
        prophages = (_load_json(run_dir / "M12_PHAGE_CRISPR_DEFENSE" / "prophages.json", {}) or {}).get("prophages", [])
        variants = (_load_json(run_dir / "M13_VARIANTS_MUTATIONS" / "variants.json", {}) or {}).get("snps", [])

        # MLST (TSV: dosya\tşema\tST\t...)
        mlst_row = None
        mlst_tsv = run_dir / "M07_STRAIN_TYPING" / "mlst_summary.tsv"
        if mlst_tsv.exists():
            txt = mlst_tsv.read_text().strip()
            if txt:
                mlst_row = txt.splitlines()[0].split("\t")

        # Tür: M02 (kraken2) -> detection -> yoksa None (asla sabit varsayım)
        species = taxonomy.get("dominant_organism") or self.ctx.detection.get("ncbi_species")

        dashboard = {
            "project_id": run_dir.name,
            "data_type": self.ctx.detection.get("data_type"),
            "platform": self.ctx.detection.get("platform"),
            "species": species,
            "taxonomy": taxonomy,
            "genome_stats": genome_stats,
            "checkm2": checkm2,
            "closest_5_strains": c5,
            "mlst": mlst_row,
            "amr_genes": amr,
            "virulence_genes": vir,
            "plasmids": plasmids,
            "mobile_elements": mges,
            "prophages": prophages,
            "variants": variants,
            "module_status": {k: v.get("status") for k, v in module_summaries.items()},
            "modules": module_summaries,
        }

        with open(std_dir / "dashboard_data.json", "w", encoding="utf-8") as fh:
            json.dump(dashboard, fh, indent=2, ensure_ascii=False)

        self.write_summary(status="PASS", statistics={
            "modules_aggregated": len(module_summaries),
            "species": species,
            "amr_gene_count": len(amr),
            "virulence_gene_count": len(vir),
            "plasmid_count": len(plasmids),
        })
