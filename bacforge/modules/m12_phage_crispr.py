"""M12 -- Phage, CRISPR & Defense Systems (geNomad; cctyper/DefenseFinder Milestone 2)
KURAL: 'Bulunamadı' satırı + otomatik PASS YOK. geNomad DB yoksa NOT_APPLICABLE;
geNomad çalışıp faj bulamazsa PASS (count=0); hata -> WARNING.
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class PhageCRISPRDefenseModule(Module):
    number = "12"
    name = "phage_crispr_defense"
    folder = "M12_PHAGE_CRISPR_DEFENSE"
    enabled_key = "genomad"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "prophages.tsv"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)
        dbp = self.ctx.config["paths"]["db"]

        genomad_out = self.sub_dir("02_work") / "genomad"
        genomad_db = Path(dbp) / "genomad" / "genomad_db"

        if not genomad_db.exists():
            with open(std_dir / "prophages.tsv", "w", encoding="utf-8") as f:
                f.write("Phage_ID\tContig\tLength\tTopology\tVirus_Score\n")
            self.write_summary(status="NOT_APPLICABLE",
                               details={"reason": f"geNomad DB yok: {genomad_db}"})
            return

        prov = r.run("genomad", ["genomad", "end-to-end", "--cleanup", str(genome), str(genomad_out),
                                 str(genomad_db), "--threads", str(t)],
                     conda_env=E["genomad"], version_cmd=["genomad", "--version"],
                     db_version=str(genomad_db), check=False)
        tool_ran = prov.get("exit_code") == 0

        phages = []
        # geNomad çıktı adı girdi dosya adına göre değişir; virus_summary.tsv'yi ara
        for vsum in genomad_out.rglob("*_virus_summary.tsv"):
            with open(vsum, "r", encoding="utf-8") as fh:
                fh.readline()  # header
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6:
                        phages.append({
                            "phage_id": parts[0],
                            "contig": parts[0].rsplit("_", 1)[0] if "_" in parts[0] else parts[0],
                            "length": parts[1],
                            "topology": parts[2],
                            "virus_score": parts[5],
                        })

        with open(std_dir / "prophages.tsv", "w", encoding="utf-8") as f:
            f.write("Phage_ID\tContig\tLength\tTopology\tVirus_Score\n")
            for p in phages:
                f.write(f"{p['phage_id']}\t{p['contig']}\t{p['length']}\t{p['topology']}\t{p['virus_score']}\n")
        with open(std_dir / "prophages.json", "w", encoding="utf-8") as f:
            json.dump({"prophages": phages}, f, indent=2, ensure_ascii=False)

        # İçerik-farkında yönlendirme: geNomad her contig'i kromozom/plazmid/virüs sınıflar.
        # Bu sınıflandırma bakteri/faj kolunu yönlendirmek için tablo olarak çıkarılır.
        classification = []
        counts = {"chromosome": 0, "plasmid": 0, "virus": 0}
        for agg in genomad_out.rglob("*_aggregated_classification.tsv"):
            lines = agg.read_text(encoding="utf-8", errors="replace").splitlines()
            if not lines:
                continue
            ix = {h: i for i, h in enumerate(lines[0].split("\t"))}

            def fnum(p, c):
                try:
                    return float(p[ix[c]]) if c in ix and ix[c] < len(p) else 0.0
                except ValueError:
                    return 0.0
            for ln in lines[1:]:
                if not ln.strip():
                    continue
                p = ln.split("\t")
                scores = {"chromosome": fnum(p, "chromosome_score"),
                          "plasmid": fnum(p, "plasmid_score"), "virus": fnum(p, "virus_score")}
                cls = max(scores, key=scores.get)
                counts[cls] += 1
                classification.append({"contig": p[ix.get("seq_name", 0)] if p else "",
                                       "class": cls, **{k: round(v, 4) for k, v in scores.items()}})
        with open(std_dir / "contig_classification.tsv", "w", encoding="utf-8") as f:
            f.write("Contig\tClass\tChromosome_score\tPlasmid_score\tVirus_score\n")
            for c in classification:
                f.write(f"{c['contig']}\t{c['class']}\t{c['chromosome']}\t{c['plasmid']}\t{c['virus']}\n")
        with open(std_dir / "contig_classification.json", "w", encoding="utf-8") as f:
            json.dump({"classification": classification, "counts": counts}, f, indent=2, ensure_ascii=False)

        if tool_ran:
            self.write_summary(status="PASS",
                               statistics={"phage_count": len(phages), "routing": counts},
                               details={"tool": "geNomad", "routing": counts,
                                        "note": "içerik-farkında yönlendirme: kromozom/plazmid/virüs"})
        else:
            self.write_summary(status="WARNING", statistics={"phage_count": len(phages)},
                               warnings=[f"geNomad başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"])
