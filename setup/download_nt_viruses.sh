#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
DB="${BACFORGE_DB:-$(pwd)/databases}"
echo "status=running started=$(date)" > setup/nt_viruses.STATUS
mkdir -p "$DB/blast_viral"
cd "$DB/blast_viral"
if conda run -n ali-blast update_blastdb.pl --decompress nt_viruses; then
  echo "status=done finished=$(date)" > "$OLDPWD/setup/nt_viruses.STATUS"
else
  echo "status=failed $(date)" > "$OLDPWD/setup/nt_viruses.STATUS"
fi
