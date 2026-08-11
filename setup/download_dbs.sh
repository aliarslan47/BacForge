#!/usr/bin/env bash
# BacForge — veritabanı indirme (BACFORGE_DB altına, sürüm etiketli)
# KURAL: DB'ler koddan ayrı (taşınabilirlik). Çalıştır: bash setup/download_dbs.sh
set -uo pipefail   # not -e: bir DB inmezse diğerleri devam etsin
DB="${BACFORGE_DB:-$(cd "$(dirname "$0")/.." && pwd)/databases}"
mkdir -p "$DB"
echo "== Veritabanları -> $DB =="

# --- Bakta light (~1.5 GB) ---
if [ ! -d "$DB/bakta/db-light" ]; then
  echo "[+] Bakta light DB"
  conda run -n ali-bakta bakta_db download --output "$DB/bakta" --type light
fi

# --- AMRFinderPlus (~0.1 GB) ---
if [ ! -d "$DB/amrfinderplus" ]; then
  echo "[+] AMRFinderPlus DB"
  conda run -n ali-amrfinder amrfinder_update --database "$DB/amrfinderplus"
fi

# --- CheckM2 (~3 GB) ---
if [ ! -e "$DB/checkm2/CheckM2_database/uniref100.KO.1.dmnd" ]; then
  echo "[+] CheckM2 DB"
  conda run -n ali-checkm2 checkm2 database --download --path "$DB/checkm2"
fi

# --- geNomad DB (~1.6 GB) — hedef klasör önceden VAR olmalı ---
if [ ! -d "$DB/genomad/genomad_db" ]; then
  echo "[+] geNomad DB"
  mkdir -p "$DB/genomad"
  conda run -n ali-genomad genomad download-database "$DB/genomad"
fi

# --- Pharokka DB (PHROGs vb. ~1.5 GB) ---
if [ ! -d "$DB/pharokka" ]; then
  echo "[+] Pharokka DB"
  conda run -n ali-pharokka install_databases.py -o "$DB/pharokka"
fi

# --- CheckV DB (~1.5 GB) ---
if [ ! -d "$DB/checkv"/checkv-db-* ]; then
  echo "[+] CheckV DB"
  conda run -n ali-checkv checkv download_database "$DB/checkv"
fi

# --- Kraken2 DB (Standard-8GB ~8GB) ---
if [ ! -d "$DB/kraken2" ] || [ -z "$(ls -A "$DB/kraken2")" ]; then
  echo "[+] Kraken2 DB (Standard-8GB)"
  mkdir -p "$DB/kraken2"
  # Using AWS S3 mirror which is much more reliable than JHU FTP
  wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz | tar zxvf - -C "$DB/kraken2"
fi

# --- Mash RefSeq DB ---
if [ ! -f "$DB/mash/RefSeqSketchesDefaults.msh" ]; then
  echo "[+] Mash RefSeq DB"
  mkdir -p "$DB/mash"
  # Use HTTP instead of HTTPS to avoid SSL handshake issues on gembox
  wget -qO "$DB/mash/RefSeqSketchesDefaults.msh.gz" http://gembox.cbcb.umd.edu/mash/refseq.genomes%2Bplasmid.k21.s1000.msh.gz
  gunzip -f "$DB/mash/RefSeqSketchesDefaults.msh.gz"
fi

# --- Abricate DB ---
if [ ! -d "$DB/abricate" ]; then
  echo "[+] Abricate DB (setupdb)"
  # Abricate updates its own internal DBs, but we'll just run setupdb
  conda run -n ali-virulence abricate --setupdb
  mkdir -p "$DB/abricate"
  touch "$DB/abricate/installed.txt"
fi

echo "== DB indirme bitti. Her DB klasörüne version.txt eklenecek (rapora yazılır). =="
