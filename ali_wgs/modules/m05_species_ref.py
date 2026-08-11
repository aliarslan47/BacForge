"""M05 -- Species & Reference Identification
Tür = M02 (kraken2) sonucundan. Closest-N = FastANI ile yerel referans setine karşı GERÇEK ANI.
KURAL: 'default strain' UYDURULMAZ. Referans yoksa closest_5 boş kalır ve durum dürüstçe belirtilir.
Referans seti: <db>/references altındaki *.fna/*.fasta/*.fa (yoksa closest-5 NOT_APPLICABLE).
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util


class SpeciesReferenceIdentificationModule(Module):
    number = "05"
    name = "species_reference_identification"
    folder = "M05_SPECIES_REFERENCE_IDENTIFICATION"
    enabled_key = "ani"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "closest_5_strains.json"]

    def _resolve_species(self) -> str | None:
        # 1) ctx.detection (M02 aynı çalıştırmada set etti)
        sp = self.ctx.detection.get("ncbi_species")
        if sp:
            return sp
        # 2) M02'nin yazdığı dosya (resume durumu)
        m02_json = self.ctx.run_dir / "M02_TAXONOMIC_QC" / "species_identification.json"
        if m02_json.exists():
            try:
                return json.load(open(m02_json)).get("species")
            except Exception:
                return None
        return None

    def _reference_fastas(self) -> list[Path]:
        ref_dir = Path(self.ctx.config["paths"]["db"]) / "references"
        if not ref_dir.exists():
            return []
        refs = []
        for ext in ("*.fna", "*.fasta", "*.fa"):
            refs.extend(ref_dir.rglob(ext))
        return sorted(set(refs))

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)

        species = self._resolve_species()
        refs = self._reference_fastas()

        strains = []
        closest_note = None

        if refs:
            work = self.sub_dir("02_work")
            reflist = work / "ref_list.txt"
            with open(reflist, "w", encoding="utf-8") as fh:
                for rp in refs:
                    fh.write(str(rp) + "\n")
            ani_out = work / "fastani.txt"
            prov = r.run("fastani",
                         ["fastANI", "-q", str(genome), "--rl", str(reflist),
                          "-o", str(ani_out), "-t", str(t)],
                         conda_env=E.get("species", "base"),
                         version_cmd=["fastANI", "--version"], check=False)
            if prov.get("exit_code") == 0 and ani_out.exists():
                rows = []
                with open(ani_out, "r", encoding="utf-8") as fh:
                    for line in fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 5:
                            try:
                                ani = float(parts[2])
                                matched = int(parts[3]); total = int(parts[4])
                            except ValueError:
                                continue
                            rows.append({
                                "ref_path": parts[1],
                                "ani_percent": round(ani, 4),
                                "query_coverage_percent": round(100.0 * matched / total, 2) if total else 0.0,
                            })
                rows.sort(key=lambda x: x["ani_percent"], reverse=True)
                for i, row in enumerate(rows[:5], 1):
                    accession = Path(row["ref_path"]).stem
                    strains.append({
                        "rank": i,
                        "organism": species or "unknown",
                        "strain": accession,
                        "assembly_accession": accession,
                        "ani_percent": row["ani_percent"],
                        "query_coverage": row["query_coverage_percent"],
                        "fasta_path": row["ref_path"],
                    })
                if not strains:
                    closest_note = "FastANI çalıştı ancak eşik üstü eşleşme yok."
            else:
                closest_note = f"FastANI başarısız (exit {prov.get('exit_code')}). Log: {prov.get('log')}"
        else:
            closest_note = ("Yerel referans seti yok (<db>/references boş) -> Closest-5 NOT_APPLICABLE. "
                            "Milestone 2: NCBI datasets ile otomatik referans çekilecek.")

        with open(std_dir / "closest_5_strains.json", "w", encoding="utf-8") as fh:
            json.dump(strains, fh, indent=2, ensure_ascii=False)
        with open(std_dir / "closest_5_strains.tsv", "w", encoding="utf-8") as fh:
            fh.write("Rank\tOrganism\tStrain\tAccession\tANI_percent\tQuery_coverage\n")
            for s in strains:
                fh.write(f"{s['rank']}\t{s['organism']}\t{s['strain']}\t{s['assembly_accession']}\t"
                         f"{s['ani_percent']}\t{s['query_coverage']}\n")
        with open(std_dir / "species_identification.json", "w", encoding="utf-8") as fh:
            json.dump({"species": species, "source": "kraken2 (M02)",
                       "closest_reference": strains[0] if strains else None}, fh, indent=2, ensure_ascii=False)

        # Durum: türün bilinmesi çekirdek çıktıdır. Closest-5 alt-parçası koşul-bazlı (NA olabilir).
        if not species:
            status = "WARNING"
            warns = ["Tür belirlenemedi (M02 sonucu yok)."]
        else:
            status = "PASS"
            warns = []
        details = {"species": species, "closest_count": len(strains)}
        if closest_note:
            details["closest_5_note"] = closest_note
        self.write_summary(status=status, statistics={"species": species, "closest_5_count": len(strains)},
                           warnings=warns, details=details)
