"""M16 -- Phylogenomics (akrabalik agaci)
GERCEK agac: sorgu genomu + M05'in buldugu gercek en-yakin referanslar -> mash genom-geneli mesafe
-> Neighbor-Joining (Biopython) -> newick + PNG. (Cok-genom; tek basina referanssiz -> NOT_APPLICABLE.)
KATI: agac gercekten olusmadan (newick + >=3 yaprak) PASS YOK.
Araclar ali-comparative env'inde: mash + Biopython.Phylo + matplotlib.
"""
from __future__ import annotations

import json
from pathlib import Path

from .base import Module
from .. import util

TREE_HELPER = r'''
import sys, json
tri, labelmap_json, out_nwk, out_png = sys.argv[1:5]
labels = json.load(open(labelmap_json))
# mash triangle: 1. satir toplam sayi; sonra her satir: isim d1 d2 ... (alt ucgen)
lines = [l.rstrip("\n") for l in open(tri) if l.strip()]
n = int(lines[0])
names, rows = [], []
for l in lines[1:]:
    parts = l.split("\t")
    names.append(parts[0])
    rows.append([float(x) for x in parts[1:]])
# Biopython alt-ucgen matris (diagonal dahil 0)
matrix = []
for i in range(len(names)):
    matrix.append(rows[i] + [0.0])
disp = [labels.get(nm, nm) for nm in names]
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
dm = DistanceMatrix(names=disp, matrix=matrix)
tree = DistanceTreeConstructor().nj(dm)
# negatif dal uzunluklarini 0'a cek (NJ artefakti)
for c in tree.find_clades():
    if c.branch_length and c.branch_length < 0:
        c.branch_length = 0.0
from Bio import Phylo
Phylo.write(tree, out_nwk, "newick")
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, max(3, len(disp) * 0.5)))
    ax = fig.add_subplot(111)
    tree.ladderize()
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title("Akrabalik agaci (mash NJ) — sorgu + en yakin referanslar")
    fig.tight_layout(); fig.savefig(out_png, dpi=150);
    print("PNG_OK")
except Exception as e:
    print("PNG_FAIL", e)
print("LEAVES", len(disp))
'''


class PhylogenomicsModule(Module):
    number = "16"
    name = "phylogenomics"
    folder = "M16_PHYLOGENOMICS"
    enabled_key = "phylogeny"

    def inputs(self):
        return [self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"]

    def outputs(self):
        return [self.out_dir / "M16_summary.json"]

    def run(self):
        self.check_inputs()
        genome = self.ctx.run_dir / "M04_POLISHING_GENOME_QC" / "genome.fasta"
        std_dir = self.sub_dir("04_standardized")
        work = self.sub_dir("02_work")
        r = self.ctx.runner
        env = util.ENV.get("comparative", "base")

        ref_json = self.ctx.run_dir / "M05_SPECIES_REFERENCE_IDENTIFICATION" / "closest_5_strains.json"
        refs = []
        if ref_json.exists():
            try:
                refs = [x for x in (json.load(open(ref_json)) or []) if Path(x.get("fasta_path", "")).exists()]
            except Exception:
                refs = []

        if len(refs) < 2:
            reason = (f"Akrabalik agaci icin >=2 referans gerekir; M05 {len(refs)} kullanilabilir referans verdi "
                      f"-> NOT_APPLICABLE.")
            self.write_summary(status="NOT_APPLICABLE", statistics={"reference_count": len(refs)},
                               details={"reason": reason})
            return

        # CESITLILIK: sadece en-yakin 5 (hepsi ~ayni) duz/anlamsiz agac verir.
        # M05'in indirdigi havuzdan daha uzak suslari da ekle -> gercek yapili agac (~12 yaprak).
        labels = {str(genome): "QUERY_sample"}
        paths = [str(genome)]
        close_paths = set()
        for x in refs:
            fp = x["fasta_path"]
            paths.append(fp); close_paths.add(fp)
            labels[fp] = f"{x['assembly_accession']}(ANI{x.get('ani_percent','?')})"
        # havuz kok dizini (closest_5'in fasta yolundan yukari) -> tum cekilen genomlar
        pool = []
        try:
            anchor = Path(refs[0]["fasta_path"])
            for up in anchor.parents:
                if up.name == "unz" or (up / "ncbi_dataset").exists():
                    pool = sorted(up.rglob("*.fna")); break
        except Exception:
            pool = []
        extra = [str(p) for p in pool if str(p) not in close_paths]
        # havuzdan esit araliklarla ~7 uzak sus sec (cesitlilik)
        take = 7
        if extra:
            step = max(1, len(extra) // take)
            for p in extra[::step][:take]:
                paths.append(p)
                labels[p] = Path(p).stem.split("_")[0] + "_" + Path(p).stem.split("_")[1] if "_" in Path(p).stem else Path(p).stem
        (work / "labels.json").write_text(json.dumps(labels), encoding="utf-8")

        # mash sketch + triangle
        sketch = work / "combined"
        r.run("mash_sketch", ["mash", "sketch", "-o", str(sketch), *paths],
              conda_env=env, version_cmd=["mash", "--version"], check=False)
        tri = work / "mash_triangle.tsv"
        prov = r.run("mash_triangle", ["mash", "triangle", str(sketch) + ".msh"],
                     conda_env=env, check=False, stdout_path=str(tri))

        if not tri.exists() or tri.stat().st_size == 0:
            self.write_summary(status="WARNING", warnings=["mash triangle mesafe matrisi uretmedi."])
            return

        # NJ agac + PNG (Biopython, ali-comparative)
        helper = work / "build_tree.py"
        helper.write_text(TREE_HELPER, encoding="utf-8")
        out_nwk = std_dir / "tree.nwk"
        out_png = self.sub_dir("06_visualization") / "phylogeny_tree.png"
        prov2 = r.run("nj_tree", ["python", str(helper), str(tri), str(work / "labels.json"),
                                  str(out_nwk), str(out_png)],
                      conda_env=env, check=False)

        # KATI dogrulama: gercek newick olustu mu (>=3 yaprak)?
        ok = out_nwk.exists() and out_nwk.stat().st_size > 0 and out_nwk.read_text().count(",") >= 2
        if not ok:
            self.write_summary(status="WARNING",
                               warnings=[f"NJ agac uretilemedi (exit {prov2.get('exit_code')}). Log: {prov2.get('log')}"])
            return

        leaf_count = len(paths)
        self.write_summary(
            status="PASS",
            statistics={"leaf_count": leaf_count, "reference_count": len(refs),
                        "method": "mash genome-wide distance + Neighbor-Joining (Biopython)"},
            details={"newick": str(out_nwk), "png": str(out_png) if out_png.exists() else None},
        )
