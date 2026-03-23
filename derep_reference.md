```
less deRep.sh


# -----------------------------
# CONFIG (AJUSTADO A TU RUTA)
# -----------------------------
WORK_BASE="/media/mmorenos/disk3/012_eduCastroLab/03_all_MAGs_indv/04_derep"

ALL_MAGS_DIR="${WORK_BASE}/00_all_bins"   # si tus bins están aquí
DREP_OUTPUT="${WORK_BASE}/095"           # output dentro de tu WORK_BASE

COMPLETENESS=50
CONTAMINATION=10
ANI_THRESHOLD=0.95
THREADS=40

mkdir -p "${DREP_OUTPUT}"
LOG="${DREP_OUTPUT}/drep_dereplicate.log"

echo "============================================================" | tee -a "${LOG}"
echo "[INFO] Start       : $(date)" | tee -a "${LOG}"
echo "[INFO] WORK_BASE   : ${WORK_BASE}" | tee -a "${LOG}"
echo "[INFO] ALL_MAGS_DIR: ${ALL_MAGS_DIR}" | tee -a "${LOG}"
echo "[INFO] DREP_OUTPUT : ${DREP_OUTPUT}" | tee -a "${LOG}"
echo "[INFO] Params      : comp=${COMPLETENESS} con=${CONTAMINATION} sa=${ANI_THRESHOLD} threads=${THREADS}" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"

# -----------------------------
# INPUTS
# -----------------------------
genomes=( "${ALL_MAGS_DIR}"/*.fa )
if (( ${#genomes[@]} == 0 )); then
  echo "[ERROR] No se encontraron .fa en: ${ALL_MAGS_DIR}" | tee -a "${LOG}" >&2
  exit 1
fi
echo "[INFO] Genomas encontrados: ${#genomes[@]}" | tee -a "${LOG}"

# -----------------------------
# RUN
# -----------------------------
# Asegúrate de tener dRep en PATH (activa tu env antes si aplica)
echo "[INFO] Running dRep..." | tee -a "${LOG}"

dRep dereplicate \
  "${DREP_OUTPUT}" \
  -g "${ALL_MAGS_DIR}"/*.fa \
  -p "${THREADS}" \
  -comp "${COMPLETENESS}" \
  -con "${CONTAMINATION}" \
  -sa "${ANI_THRESHOLD}" \
  --ignoreGenomeQuality \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[INFO] End          : $(date)" | tee -a "${LOG}"
echo "[INFO] Derep MAGs in: ${DREP_OUTPUT}/dereplicated_genomes/" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
```
