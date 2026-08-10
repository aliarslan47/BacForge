#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
DB="${ALI_WGS_DB:-$(pwd)/databases}"
echo "status=running" > setup/blast.STATUS
conda env list | grep -q '^ali-blast ' || conda create -n ali-blast -c conda-forge -c bioconda blast --yes
mkdir -p "$DB/blast_viral"
cd "$DB/blast_viral"
# NCBI TAM viral nt (~68 GB) — doğru tür-düzey kimlik için (temsilci set DEĞİL)
conda run -n ali-blast update_blastdb.pl --decompress nt_viruses \
  && echo "status=done" > "$OLDPWD/setup/blast.STATUS" \
  || echo "status=failed" > "$OLDPWD/setup/blast.STATUS"
