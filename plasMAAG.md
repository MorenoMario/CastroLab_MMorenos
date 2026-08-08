
BASE=/home/mmoreno/00_plasmaag/out_comau_90
PLASMIDS="${BASE}/results/candidate_plasmids"

mkdir -p "${BASE}/derep_plasmids/mmseqs_tmp"

find "${PLASMIDS}" \
    -type f \
    \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) \
    -size +0c \
    -exec cat {} + \
    > "${BASE}/derep_plasmids/all_candidate_plasmids.fna"


grep -c '^>' \
    "${BASE}/derep_plasmids/all_candidate_plasmids.fna"

8651

#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=plasmid_derep_mmseqs
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH -o plasmid_derep_mmseqs.%j.out
#SBATCH -e plasmid_derep_mmseqs.%j.err

set -euo pipefail
shopt -s nullglob

THREADS="${SLURM_CPUS_PER_TASK:-16}"

BASE="${HOME}/00_plasmaag/out_comau_90"
PLASMIDS="${BASE}/results/candidate_plasmids"
OUTDIR="${BASE}/derep_plasmids_mmseqs"

mkdir -p "${OUTDIR}"

if ! command -v mmseqs >/dev/null 2>&1; then
    echo "[ERROR] mmseqs no está en PATH" >&2
    exit 1
fi

if [[ ! -d "${PLASMIDS}" ]]; then
    echo "[ERROR] no existe el directorio: ${PLASMIDS}" >&2
    exit 1
fi

echo "=========================================================="
echo "[INFO] inicio      : $(date '+%F %T')"
echo "[INFO] base        : ${BASE}"
echo "[INFO] plasmids    : ${PLASMIDS}"
echo "[INFO] salida      : ${OUTDIR}"
echo "[INFO] nodo        : $(hostname)"
echo "[INFO] hilos       : ${THREADS}"
echo "=========================================================="

echo "[INFO] concatenando FASTA candidatos plasmídicos"
ALL_FASTA="${OUTDIR}/all_candidate_plasmids.fna"

find "${PLASMIDS}" \
    -type f \
    \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) \
    -size +0c \
    -exec cat {} + \
    > "${ALL_FASTA}"

N0=$(grep -c '^>' "${ALL_FASTA}")
echo "[INFO] secuencias de entrada: ${N0}"

if [[ "${N0}" -eq 0 ]]; then
    echo "[ERROR] no se encontraron secuencias en ${ALL_FASTA}" >&2
    exit 1
fi

run_cluster () {
    local id="$1"
    local cov="$2"
    local tmp="${OUTDIR}/tmp_${id}"
    local prefix="${OUTDIR}/plasmids_nr${id}"

    mkdir -p "${tmp}"

    echo
    echo "=========================================================="
    echo "[INFO] MMseqs2 easy-cluster ${id}"
    echo "[INFO] min-seq-id=${cov}"
    echo "=========================================================="

    mmseqs easy-cluster \
        "${ALL_FASTA}" \
        "${prefix}" \
        "${tmp}" \
        --min-seq-id "${cov}" \
        -c 0.95 \
        --cov-mode 0 \
        --threads "${THREADS}"

    local rep="${prefix}_rep_seq.fasta"
    local clu="${prefix}_cluster.tsv"

    if [[ ! -s "${rep}" ]]; then
        echo "[ERROR] no se generó ${rep}" >&2
        exit 1
    fi

    if [[ ! -s "${clu}" ]]; then
        echo "[ERROR] no se generó ${clu}" >&2
        exit 1
    fi

    local nrep nclu
    nrep=$(grep -c '^>' "${rep}")
    nclu=$(cut -f1 "${clu}" | sort -u | wc -l)

    echo "[INFO] representantes: ${nrep}"
    echo "[INFO] clusters      : ${nclu}"

    cp "${rep}" "${OUTDIR}/plasmids_nr${id}.fna"

    python - <<PY
import pandas as pd
from pathlib import Path

clu = Path("${clu}")
out = Path("${OUTDIR}/cluster_sizes_nr${id}.tsv")

# MMseqs cluster.tsv typically has: cluster_rep  member
df = pd.read_csv(clu, sep="\t", header=None, names=["rep", "member"])
sizes = df.groupby("rep").size().reset_index(name="cluster_size")
sizes = sizes.sort_values(["cluster_size", "rep"], ascending=[False, True])
sizes.to_csv(out, sep="\t", index=False)
print(out)
print(sizes["cluster_size"].describe().to_string())
PY

    echo "[INFO] distribución de tamaños de cluster:"
    python - <<PY
import pandas as pd
from pathlib import Path

sizes = pd.read_csv(Path("${OUTDIR}/cluster_sizes_nr${id}.tsv"), sep="\t")
print(sizes["cluster_size"].value_counts().sort_index().to_string())
PY
}

run_cluster 99 0.99
run_cluster 95 0.95

echo
echo "=========================================================="
echo "[INFO] resumen final"
for id in 99 95; do
    rep="${OUTDIR}/plasmids_nr${id}.fna"
    clu="${OUTDIR}/plasmids_nr${id}_cluster.tsv"
    siz="${OUTDIR}/cluster_sizes_nr${id}.tsv"

    echo
    echo "===== ${id} ====="
    echo -n "representantes: "
    grep -c '^>' "${rep}"
    echo -n "clusters      : "
    cut -f1 "${clu}" | sort -u | wc -l
    echo "[INFO] top 20 cluster sizes:"
    python - <<PY
import pandas as pd
from pathlib import Path
sizes = pd.read_csv(Path("${siz}"), sep="\t")
print(sizes["cluster_size"].value_counts().sort_index().tail(20).to_string())
PY
done

echo
echo "[INFO] fin: $(date '+%F %T')"

/home/mmoreno/00_plasmaag/out_comau_90/derep_plasmids_mmseqs


===== 99% identidad =====
Representantes FASTA: 8275
Clusters únicos:      8275
Secuencias agrupadas: 8651
Distribución de tamaños:
   8021 1
    196 2
     29 3
     15 4
      6 5
      2 6
      3 7
      2 9
      1 10

===== 95% identidad =====
Representantes FASTA: 8222
Clusters únicos:      8222
Secuencias agrupadas: 8651
Distribución de tamaños:
   7932 1
    224 2
     37 3
     13 4
      6 5
      2 6
      3 7
      2 8
      2 9
      1 11

