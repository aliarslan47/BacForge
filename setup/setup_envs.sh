#!/usr/bin/env bash
# Ali WGS Pipeline — izole conda env kurulumu (MVP: bakteri + faj zinciri)
# Her araç ayrı env'de -> çakışma yok + Docker'a kolay geçiş.
# Çalıştır: bash setup/setup_envs.sh
set -uo pipefail   # not -e: bir araç kurulamazsa diğerleri devam etsin
CHAN="-c conda-forge -c bioconda --yes"
FAILED=()

have_env () { conda env list | awk '{print $1}' | grep -qx "$1"; }
mk () {
  local name="$1"; shift
  if have_env "$name"; then echo "[=] $name zaten var, atlandı"; return; fi
  echo "[+] $name kuruluyor: $*"
  if conda create -n "$name" $CHAN "$@"; then
    echo "[ok] $name"
  else
    echo "[HATA] $name kurulamadı (devam ediliyor)"; FAILED+=("$name")
  fi
}

# --- ortak / okuma işleme ---
mk ali-ont-qc    nanoplot chopper filtlong rasusa seqkit   # 03,04 + subsample
mk ali-illumina-qc fastqc fastp                            # 03,04 (Illumina)
# --- assembly ---
mk ali-flye      flye                                      # 05 ONT/Hybrid/faj
mk ali-medaka    medaka                                    # 05 ONT polishing
mk ali-assembly-sr spades skesa unicycler                  # 05 Illumina
mk ali-hifiasm   hifiasm                                   # 05 HiFi
# --- QC / sınıflandırma ---
mk ali-quast     quast                                     # 06
mk ali-genomad   genomad                                   # 08 içerik router
# --- annotation (sınıfa göre) ---
mk ali-bakta     bakta                                     # 10 bakteri/arkea
mk ali-pharokka  pharokka                                  # 10 faj
# --- completeness (sınıfa göre) ---
mk ali-checkm2   checkm2                                   # 17 bakteri
mk ali-checkv    checkv                                    # 17 faj
mk ali-bacphlip  bacphlip                                  # faj yaşam tarzı
# --- ortak ---
mk ali-amrfinder ncbi-amrfinderplus                        # 11,12 (bakteri+faj)

if [ ${#FAILED[@]} -gt 0 ]; then
  echo "== UYARI: kurulamayan env'ler: ${FAILED[*]} =="
else
  echo "== Tüm izole env'ler kuruldu =="
fi
echo "== Sıradaki: bash setup/download_dbs.sh =="
