"""M06 -- Genome Annotation & Functional Annotation (Bakta)
KATI SÖZLEŞME: Bakta gerçekten çalışıp gerçek çıktı (gff3 + GenBank + faa) üretmeden PASS YOK.
Bakta çökerse -> FAIL (kritik modül; mock/sahte gen YAZILMAZ). --skip-sorf: DIAMOND sORF segfault (exit -11)
adımını atlar, tam CDS anotasyonu + GenBank yine üretilir.
Çıktılar: annotation.gff3, annotation.gbk, annotation.faa, annotation.ffn, annotation.tsv, identifiers.tsv
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

from .base import Module
from .. import util


class AnnotationIntegrityValidator:
    """GFF3'ü zorunlu spec kurallarına göre denetler (locus_tag != gene_symbol vb.)."""
    @staticmethod
    def validate_gff3(gff_path: Path) -> dict:
        errors, warnings = [], []
        locus_tags = set()
        feature_count = 0
        locus_copied_as_gene = 0
        with open(gff_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                feature_count += 1
                contig, source, ftype, start, end, score, strand, phase, attr_str = parts[:9]
                try:
                    s_pos, e_pos = int(start), int(end)
                    if s_pos > e_pos or s_pos < 1:
                        errors.append(f"Invalid coordinates: {contig}:{start}-{end}")
                except ValueError:
                    errors.append(f"Non-numeric coordinates: {start}-{end}")
                if strand not in ("+", "-", "."):
                    warnings.append(f"Unknown strand: {strand}")
                attrs = {}
                for item in attr_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k.strip()] = v.strip()
                locus = attrs.get("locus_tag") or attrs.get("ID")
                gene = attrs.get("gene") or attrs.get("gene_symbol")
                if locus:
                    if locus in locus_tags:
                        errors.append(f"Duplicate locus_tag: {locus}")
                    locus_tags.add(locus)
                if locus and gene and locus == gene:
                    locus_copied_as_gene += 1
                    errors.append(f"Integrity Violation: locus_tag '{locus}' copied to gene_symbol")
        return {
            "passed": len(errors) == 0,
            "feature_count": feature_count,
            "unique_locus_tags": len(locus_tags),
            "locus_copied_as_gene": locus_copied_as_gene,
            "errors": errors[:50],
            "warnings": warnings[:50],
        }


class GenomeAnnotationModule(Module):
    number = "06"
    name = "genome_annotation"
    folder = "M06_GENOME_ANNOTATION"
    enabled_key = "annotation"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        # KATI: hem GFF3 hem GERÇEK GenBank üretilmeli
        return [self.out_dir / "annotation.gff3", self.out_dir / "annotation.gbk"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        r = self.ctx.runner
        E = util.ENV
        t = util.threads(self.ctx)
        dbp = self.ctx.config["paths"]["db"]
        bakta_dir = self.sub_dir("02_work") / "bakta"

        # Bakta -- --skip-sorf ile (DIAMOND sORF segfault'unu önler). check=True: gerçekten çökerse yüksek sesle.
        prov = r.run("bakta", [
            "bakta", "--db", f"{dbp}/bakta/db-light", "--output", str(bakta_dir),
            "--prefix", "genome", "--skip-sorf", "--threads", str(t), "--force", str(genome)
        ], conda_env=E["bakta"], version_cmd=["bakta", "--version"], check=False)

        bakta_gff = bakta_dir / "genome.gff3"
        bakta_gbff = bakta_dir / "genome.gbff"

        # KATI DOĞRULAMA: exit 0 + gerçek gff3 + gerçek GenBank (boş değil)
        ok = (prov.get("exit_code") == 0
              and bakta_gff.exists() and bakta_gff.stat().st_size > 1000
              and bakta_gbff.exists() and bakta_gbff.stat().st_size > 1000)
        if not ok:
            # Mock YAZMA. Kritik modül -> yüksek sesle FAIL (orchestrator durdurur).
            raise RuntimeError(
                f"[M06] Bakta gerçek çıktı üretmedi (exit {prov.get('exit_code')}, "
                f"gff3={bakta_gff.exists()}, gbff={bakta_gbff.exists()}). Log: {prov.get('log')}. "
                f"Sahte anotasyon YAZILMADI."
            )

        # Gerçek çıktıları standardize et (clinker/M14 annotation.gbk kullanır)
        shutil.copy(bakta_gff, std_dir / "annotation.gff3")
        shutil.copy(bakta_gbff, std_dir / "annotation.gbk")
        for ext, dst in [("faa", "annotation.faa"), ("ffn", "annotation.ffn"), ("tsv", "annotation.tsv")]:
            src = bakta_dir / f"genome.{ext}"
            if src.exists():
                shutil.copy(src, std_dir / dst)
        for img in ("genome.png", "genome.svg"):
            src = bakta_dir / img
            if src.exists():
                shutil.copy(src, self.sub_dir("06_visualization") / img.replace("genome", "genome_map"))

        # Integrity denetimi (gerçek gff3 üzerinde)
        audit = AnnotationIntegrityValidator.validate_gff3(std_dir / "annotation.gff3")
        with open(std_dir / "annotation_integrity.json", "w", encoding="utf-8") as fh:
            json.dump(audit, fh, indent=2, ensure_ascii=False)

        # Özellik sayıları (gerçek)
        cds = trna = rrna = 0
        with open(std_dir / "annotation.gff3", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or "\t" not in line:
                    continue
                ft = line.split("\t")[2] if len(line.split("\t")) > 2 else ""
                if ft == "CDS": cds += 1
                elif ft == "tRNA": trna += 1
                elif ft == "rRNA": rrna += 1

        id_map_dir = self.sub_dir("08_metadata") / "identifier_mapping"
        id_map_dir.mkdir(parents=True, exist_ok=True)
        (id_map_dir / "identifiers.tsv").write_text(
            "feature_type\tlocus_tag\tgene_symbol\tproduct\tcontig\tstart\tend\tstrand\n", encoding="utf-8")

        # Dairesel genom haritasi: parcali assembly'de per-contig plot "0.0 Mbp" verir.
        # Contig'leri en yakin referansa gore sirala -> tek psodo-kromozom -> tek daire (dogru Mbp).
        try:
            self._reference_ordered_map(genome, std_dir, r, E, t, dbp)
        except Exception as exc:  # harita basarisiz olsa da M06'yi dusurme
            print(f"[M06] referans-sirali harita uretilemedi ({exc}); per-contig harita kullanilacak.")

        self.write_summary(
            status="PASS" if audit["passed"] else "WARNING",
            statistics={"cds": cds, "trna": trna, "rrna": rrna,
                        "total_features": audit["feature_count"],
                        "unique_locus_tags": audit["unique_locus_tags"],
                        "integrity_passed": audit["passed"]},
            warnings=audit["warnings"],
            details={"annotator": "Bakta (--skip-sorf)", "genbank": str(std_dir / "annotation.gbk")},
        )

    def _reference_ordered_map(self, genome, std_dir, r, E, t, dbp):
        """Contig'leri en yakin referansa (M05) gore sirala -> tek psodo-kromozom -> Bakta plot
        (tek daire, dogru Mbp ekseni). Boylece parcali assembly'nin '0.0 Mbp' hatasi cozulur."""
        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        if not ref_json.exists():
            return
        strains = json.load(open(ref_json))
        if not strains:
            return
        ref = Path(strains[0].get("fasta_path", ""))
        if not ref.exists():
            return
        work = self.sub_dir("02_work")
        paf = work / "contigs_vs_ref.paf"
        r.run("minimap2_order", ["minimap2", "-x", "asm20", "-t", str(t), str(ref), str(genome)],
              conda_env=E.get("core", "base"), stdout_path=str(paf), check=False)
        if not paf.exists() or paf.stat().st_size == 0:
            return
        best = {}
        for line in open(paf, encoding="utf-8"):
            p = line.split("\t")
            if len(p) < 12:
                continue
            q, strand, tname = p[0], p[4], p[5]
            try:
                tstart, aln = int(p[7]), int(p[10])
            except ValueError:
                continue
            if q not in best or aln > best[q][0]:
                best[q] = (aln, tname, tstart, strand)
        seqs = util.read_fasta(genome)
        mapped = sorted(best.keys(), key=lambda q: (best[q][1], best[q][2]))
        order = mapped + [c for c in seqs if c not in best]
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
                "a": "t", "t": "a", "g": "c", "c": "g"}
        parts = []
        for c in order:
            s = seqs[c]
            if c in best and best[c][3] == "-":
                s = "".join(comp.get(b, "N") for b in reversed(s))
            parts.append(s)
        pseudo = work / "pseudochromosome.fasta"
        pseudo.write_text(">pseudochromosome_ref_ordered\n" + ("N" * 100).join(parts) + "\n", encoding="utf-8")
        pdir = work / "bakta_pseudo"
        prov = r.run("bakta_map", ["bakta", "--db", f"{dbp}/bakta/db-light", "--output", str(pdir),
                     "--prefix", "pseudo", "--skip-sorf", "--threads", str(t), "--force", str(pseudo)],
                     conda_env=E["bakta"], check=False)
        png = pdir / "pseudo.png"
        if prov.get("exit_code") == 0 and png.exists():
            shutil.copy(png, std_dir / "genome_map.png")   # tek daire, dogru Mbp
