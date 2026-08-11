"""Modül 08 — İçerik Sınıflandırma (geNomad) = içerik-farkında ROUTER.
INPUT : 07_Contig_Filtering/contigs.filtered.fasta
OUTPUT: 08_Taxonomy/classification.tsv + bacterial.fasta / viral.fasta / plasmid.fasta
        (+ prophages.tsv: kromozom içi entegre profajlar)
Sonraki modüller (10 annotation, 17 completeness) bu sınıf dosyalarına göre dallanır.

ÖNEMLİ: Tüm-contig sınıfı *_aggregated_classification.tsv'den argmax(skor) ile belirlenir.
virus_summary.tsv ENTEGRE PROFAJLARI da listeler (örn. 'contig|provirus_x_y'); bunlar
ana kromozomu virüs YAPMAZ — ayrı 'prophages.tsv' olarak raporlanır.
"""
from __future__ import annotations

from pathlib import Path

from .base import Module
from .. import util


class ContentClassificationModule(Module):
    number = "08"
    name = "content_classification"
    folder = "08_Taxonomy"
    enabled_key = "content_classification"

    def inputs(self):
        return [util.filtered_contigs(self.ctx)]

    def outputs(self):
        return [self.out_dir / "classification.tsv"]

    def _genomad_db(self):
        return str(Path(self.ctx.config["paths"]["db"]) / "genomad" / "genomad_db")

    def run(self):
        E = util.ENV
        t = util.threads(self.ctx)
        r = self.ctx.runner
        contigs = util.filtered_contigs(self.ctx)
        gd = self.out_dir / "genomad"

        # geNomad çıktısı zaten varsa tekrar koşturma (hızlı resume)
        done = gd.exists() and any("provirus" not in p.name
                                   for p in gd.rglob("*_aggregated_classification.tsv"))
        if done:
            prov = {"exit_code": 0}
        else:
            prov = r.run("genomad", ["genomad", "end-to-end", "--cleanup",
                                     str(contigs), str(gd), self._genomad_db(),
                                     "--threads", str(t)],
                         conda_env=E["genomad"], version_cmd=["genomad", "--version"],
                         check=False, db_version="genomad_db")

        seqs = util.read_fasta(contigs)
        classes = self._classify(gd, prov["exit_code"], seqs)

        viral = {k: v for k, v in seqs.items() if classes.get(k) == "virus"}
        plasmid = {k: v for k, v in seqs.items() if classes.get(k) == "plasmid"}
        bacterial = {k: v for k, v in seqs.items() if classes.get(k) == "chromosome"}

        util.write_fasta(viral, self.out_dir / "viral.fasta")
        util.write_fasta(plasmid, self.out_dir / "plasmid.fasta")
        util.write_fasta(bacterial, self.out_dir / "bacterial.fasta")
        self._write_prophages(gd)

        with open(self.outputs()[0], "w") as out:
            out.write("contig\tlength\tclass\n")
            for name, seq in seqs.items():
                out.write(f"{name}\t{len(seq)}\t{classes.get(name, 'chromosome')}\n")

    def _classify(self, gd: Path, exit_code: int, seqs: dict) -> dict:
        """Her tam contig'i aggregated_classification skorlarının argmax'ı ile sınıfla."""
        if exit_code != 0:
            (self.out_dir / "WARN.txt").write_text(
                "geNomad başarısız -> hepsi 'chromosome' (generic annotation)\n")
            return {k: "chromosome" for k in seqs}

        # provirus DOSYASINI dışla; sadece tam-contig sınıflaması
        aggs = [p for p in gd.rglob("*_aggregated_classification.tsv")
                if "provirus" not in p.name]
        classes = {}
        if aggs:
            with open(aggs[0]) as fh:
                fh.readline()  # başlık: seq_name chromosome_score plasmid_score virus_score
                for line in fh:
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 4:
                        continue
                    score = {"chromosome": float(c[1]), "plasmid": float(c[2]),
                             "virus": float(c[3])}
                    classes[c[0]] = max(score, key=score.get)
        # aggregated'da olmayan contig -> chromosome (güvenli varsayılan)
        for k in seqs:
            classes.setdefault(k, "chromosome")
        return classes

    def _write_prophages(self, gd: Path):
        """Kromozom içi entegre profajları ayrı raporla (virus_summary 'contig|provirus_x_y')."""
        out = self.out_dir / "prophages.tsv"
        rows = ["region\thost_contig"]
        for vs in gd.rglob("*_virus_summary.tsv"):
            if "provirus" in vs.name:
                continue
            with open(vs) as fh:
                fh.readline()
                for line in fh:
                    sid = line.split("\t")[0]
                    if "|provirus" in sid or "|provirus_" in sid:
                        rows.append(f"{sid}\t{sid.split('|')[0]}")
        out.write_text("\n".join(rows) + "\n")
