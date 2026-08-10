# Ali WGS Pipeline — Veri Akışı & Input/Output Sözleşmesi (prokaryotik)

> Kural: her modül INPUT'unu bir önceki modülün OUTPUT'undan alır.
> Eksik girdi = modül çalışmaz, açık hata verir (sessiz geçiş YOK).
> Her OUTPUT kendi klasörüne yazılır; üzerine yazma YOK; resume `is_done()` ile.

| # | Modül | INPUT | TOOL | OUTPUT |
|---|-------|-------|------|--------|
| 01 | Input Detection | ham okuma (FASTQ/FAST5/POD5/BAM) | detect.py + seqkit stats | 01_Input/platform.json |
| 02 | Preprocessing | sinyal (FAST5/POD5) veya BAM | dorado / samtools fastq | 02_Preprocessing/reads.fastq.gz |
| 03 | QC | okumalar | NanoPlot (ONT) / FastQC+fastp (Illumina) / seqkit | 03_QC/*.html, qc_stats.tsv |
| 04 | Filtering | okumalar + QC | chopper/Filtlong (ONT) / fastp (Illumina) | 04_Filtering/filtered.fastq.gz |
| 05 | Assembly | filtreli okuma | Flye(+Autocycler)→Medaka / SKESA/SPAdes / hifiasm / Flye→Medaka→Polypolish | 05_Assembly/assembly.fasta |
| 06 | Assembly QC | assembly.fasta + okuma | QUAST, CheckM2, minimap2/bwa | 06_Assembly_QC/quast/, coverage.tsv, circular.tsv |
| 07 | Contig Filtering | assembly.fasta + coverage | seqkit + eşik | 07_Contig_Filtering/contigs.filtered.fasta |
| 08 | Taxonomy | filtreli contig | skani/Mash/GTDB-Tk (ops. Kraken2) | 08_Taxonomy/taxonomy.tsv |
| 09 | BLAST | filtreli contig | blastn (uzak/hedefli) + reciprocal coverage | 09_BLAST/blast.tsv |
| 10 | Annotation | contigs.filtered.fasta | Bakta | 10_Annotation/*.gff3,.gbff,.faa,.ffn,.tsv |
| 11 | AMR | contig + protein | AMRFinderPlus | 11_AMR/amr.tsv |
| 12 | VFDB | protein/contig | AMRFinderPlus --plus / abricate VFDB | 12_VFDB/virulence.tsv |
| 13 | Plasmid | contig | MOB-suite/PlasmidFinder/platon | 13_Plasmid/plasmid_report.tsv |
| 14 | Phage | contig | geNomad/VIBRANT/PhiSpy + CheckV | 14_Phage/prophages.tsv |
| 15 | MLST | contig | mlst / chewBBACA | 15_MLST/mlst.tsv |
| 16 | ANI | contig | FastANI/skani | 16_ANI/ani.tsv |
| 17 | Completeness | contig/protein | CheckM2 (alt BUSCO) | 17_Completeness/quality.tsv |
| 18 | Final Report | tüm çıktılar + provenance | Jinja2 + WeasyPrint | 18_Final_Report/report.html, report.pdf |
| 19 | Logs | her adım | ToolRunner | 19_Logs/<m>.log, <m>.provenance.json |

## Platforma göre 05_Assembly dallanması
- ONT      : Flye → (Autocycler konsensüs) → Medaka polish
- Illumina : SKESA (varsayılan) / SPAdes --isolate
- HiFi     : hifiasm
- Hybrid   : Flye (long-first) → Medaka → Polypolish (short-read)
- Faj/küçük genom: aşırı coverage → rasusa/Filtlong ile ~100x'e subsample (04 öncesi/sonrası)

## 08 İçerik sınıflandırma + otomatik routing (geNomad)
INPUT: contigs.filtered.fasta → geNomad → her contig sınıfı (chromosome/plasmid/virus)
OUTPUT: 08_Taxonomy/genomad/ + sınıf bazlı FASTA (bacterial.fasta, viral.fasta, plasmid.fasta)

Sınıfa göre OTOMATİK dallanan modüller:
| Modül | chromosome/plasmid (bakteri) | virus (faj) |
|-------|------------------------------|-------------|
| 10 Annotation | Bakta | Pharokka |
| 17 Completeness | CheckM2 | CheckV |
| Taxonomy | GTDB-Tk/skani | taxmyPHAGE/geNomad |
| Yaşam tarzı | — | BACPHLIP |
| 15 MLST | mlst | (atlanır) |
| 11 AMR / 12 VFDB | AMRFinderPlus (ortak — her ikisinde de çalışır) | |

`routing.on_unknown=annotate_generic`: sınıflanamayan contig → Bakta(generic) + uyarı, ASLA sessiz geçiş.

## MVP'de açık modüller
01, 03, 04, 05, 06, 07, 08(geNomad router), 10(otomatik Bakta/Pharokka),
11(AMR), 17(otomatik CheckM2/CheckV), 18 (+19 her zaman). Diğerleri config ile sonra.
