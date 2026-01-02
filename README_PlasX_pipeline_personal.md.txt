# Protocolo personal: Predicción de contigs putativos plasmídicos con PlasX
> Documento de uso personal (pipeline + scripts completos). Diseñado para pegar tal cual en GitHub.
> Última actualización: 2026-01-02

## Resumen del flujo
1. Exportar/armar insumos desde Anvi’o (gene calls / secuencias según tu estrategia).
2. Ejecutar **PlasX de novo** para construir/ajustar el modelo (si aplica en tu set).
3. Ejecutar **PlasX predict** sobre contigs (o sobre el formato de entrada definido por tu script).
4. Consolidar resultados con el script de análisis en Python.

## Scripts (completos) en secciones desplegables
A continuación se incluyen tus scripts **completos** dentro de bloques desplegables (`<details>`), para mantener el README limpio.

<details>
<summary><strong>01_anvio_export_gene_calls_all_n002.slurm</strong></summary>

**Notas personales:**
- Exporta los insumos desde Anvi’o; revisar rutas y listas de entrada antes de correr.
- Si cambias el *gene caller* o el `contigs.db`, este es el punto donde se refleja.

```bash
#!/usr/bin/env bash
set -euo pipefail

#SBATCH --partition=basic
#SBATCH --job-name=anvio_export_gc
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=n002
#SBATCH -o anvio_export_gc.%j.out
#SBATCH -e anvio_export_gc.%j.err

THREADS="${SLURM_CPUS_PER_TASK:-16}"

BASE="/home/databases/comau/00_anvio_comau"
DB_DIR="${BASE}"
OUT_DIR="${BASE}/plasx/gene_calls"
LOG_DIR="${BASE}/plasx/logs_gene_calls"
mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "============================================================"
echo "[INFO] JobID   : ${SLURM_JOB_ID:-NA}"
echo "[INFO] Host    : $(hostname)"
echo "[INFO] Start   : $(date)"
echo "[INFO] DB_DIR  : ${DB_DIR}"
echo "[INFO] OUT_DIR : ${OUT_DIR}"
echo "============================================================"

# --- conda/anvio ---
if [[ -f /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh ]]; then
  source /opt/miniconda/miniconda-py312_24.3.0-0/etc/profile.d/conda.sh
else
  echo "[WARN] No encontré conda.sh en /opt/miniconda/... Ajusta la ruta si aplica." >&2
fi

conda activate anvio-8
command -v anvi-export-gene-calls >/dev/null 2>&1 || { echo "[ERROR] anvi-export-gene-calls no está disponible en env anvio-8" >&2; exit 1; }

report="${OUT_DIR}/anvio_gene_calls_report_${SLURM_JOB_ID}.tsv"
echo -e "prefix\tdb\tgene_calls\tstatus\tmsg" > "$report"

ok=0
fail=0

shopt -s nullglob
dbs=( "${DB_DIR}"/*.db )
shopt -u nullglob

echo "[INFO] DBs detectadas: ${#dbs[@]}"

for db in "${dbs[@]}"; do
  prefix="$(basename "$db" .db)"
  out="${OUT_DIR}/${prefix}-gene-calls.txt"
  log="${LOG_DIR}/${prefix}.log"

  echo "------------------------------------------------------------" | tee "$log"
  echo "[RUN] ${prefix}" | tee -a "$log"
  echo "[DB ] ${db}"     | tee -a "$log"
  echo "[OUT] ${out}"    | tee -a "$log"

  rm -f "$out"

  set +e
  anvi-export-gene-calls --gene-caller prodigal -c "$db" -o "$out" &>> "$log"
  rc=$?
  set -e

  if [[ $rc -ne 0 || ! -s "$out" ]]; then
    echo "[FAIL] rc=${rc}" | tee -a "$log"
    echo -e "${prefix}\t${db}\t${out}\tFAIL\trc=${rc}" >> "$report"
    ((fail++))
    continue
  fi

  # check aa_sequence column exists (needed by plasx)
  if ! head -n 1 "$out" | tr '\t' '\n' | grep -qx "aa_sequence"; then
    echo "[FAIL] falta columna aa_sequence" | tee -a "$log"
    echo -e "${prefix}\t${db}\t${out}\tFAIL\tno_aa_sequence" >> "$report"
    ((fail++))
    continue
  fi

  echo "[OK]" | tee -a "$log"
  echo -e "${prefix}\t${db}\t${out}\tOK\t-" >> "$report"
  ((ok++))
done

echo "============================================================"
echo "[INFO] End    : $(date)"
echo "[INFO] OK     : ${ok}"
echo "[INFO] FAIL   : ${fail}"
echo "[INFO] Report : ${report}"
echo "============================================================"

```

</details>

<details>
<summary><strong>02_plasx_denovo_all_test.slurm</strong></summary>

**Notas personales:**
- Modo *de novo*: úsalo cuando quieras generar/ajustar el modelo a tus datos.
- Este paso puede ser el más pesado en RAM/CPU dependiendo del tamaño del set.

```bash
#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=plasx_denovo_all
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=n002
#SBATCH -o plasx_denovo_all.%j.out
#SBATCH -e plasx_denovo_all.%j.err

set -euo pipefail
shopt -s nullglob

THREADS="${SLURM_CPUS_PER_TASK:-16}"
SPLITS=32

BASE="/home/databases/comau/00_anvio_comau"
GENES_DIR="${BASE}/plasx/gene_calls"
DENOVO_DIR="${BASE}/plasx/plasx_denovo"
LOG_DIR="${BASE}/plasx/logs_denovo"
TMPBASE="${BASE}/plasx/tmp"

mkdir -p "$DENOVO_DIR" "$LOG_DIR" "$TMPBASE"
export TMPDIR="${TMPBASE}/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

REPORT="${DENOVO_DIR}/plasx_denovo_report.tsv"

echo "============================================================"
echo "[INFO] JobID    : ${SLURM_JOB_ID:-NA}"
echo "[INFO] NodeList : ${SLURM_NODELIST:-NA}"
echo "[INFO] Host     : $(hostname)"
echo "[INFO] Start    : $(date)"
echo "[INFO] Threads  : ${THREADS}"
echo "[INFO] Splits   : ${SPLITS}"
echo "[INFO] GENES    : ${GENES_DIR}"
echo "[INFO] DENOVO   : ${DENOVO_DIR}"
echo "[INFO] LOGS     : ${LOG_DIR}"
echo "[INFO] TMPDIR   : ${TMPDIR}"
echo "============================================================"

# Activa tu env de PlasX (ajusta el nombre si corresponde)
module load Miniconda || true
# Ejemplos (deja SOLO el que uses):
# conda activate plasx
# conda activate vamb
# mamba activate plasx
conda activate plasx || true

PLASX_BIN="$(command -v plasx || true)"
if [[ -z "${PLASX_BIN}" ]]; then
  echo "[ERROR] plasx no está en PATH. Activa el conda env correcto (y evita ~/.local/bin)." >&2
  exit 1
fi
echo "[INFO] plasx    : ${PLASX_BIN}"

gene_files=( "${GENES_DIR}"/*-gene-calls.txt )
if (( ${#gene_files[@]} == 0 )); then
  echo "[ERROR] No hay archivos *-gene-calls.txt en ${GENES_DIR}" >&2
  exit 1
fi

echo -e "prefix\tgene_calls\tdenovo_out\tstatus\tmsg" > "${REPORT}"

for gc in "${gene_files[@]}"; do
  prefix="$(basename "${gc}" -gene-calls.txt)"
  out="${DENOVO_DIR}/${prefix}-de-novo-families.txt"
  log="${LOG_DIR}/${prefix}.denovo.log"

  echo "------------------------------------------------------------"
  echo "[INFO] PREFIX : ${prefix}"
  echo "[INFO] GENE   : ${gc}"
  echo "[INFO] OUT    : ${out}"
  echo "[INFO] LOG    : ${log}"
  echo "------------------------------------------------------------"

  # Validación rápida de columna aa_sequence (PlasX la necesita)
  if ! head -n 1 "${gc}" | tr '\t' '\n' | grep -qx "aa_sequence"; then
    echo "[SKIP] Falta columna aa_sequence en ${gc}" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${out}\tSKIP\tmissing_aa_sequence_column" >> "${REPORT}"
    continue
  fi

  set +e
  "${PLASX_BIN}" search_de_novo_families \
      -g "${gc}" \
      -o "${out}" \
      --threads "${THREADS}" \
      --splits "${SPLITS}" \
      --overwrite \
      > "${log}" 2>&1
  rc=$?
  set -e

  if [[ $rc -eq 0 ]]; then
    echo "[OK] ${prefix}" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${out}\tOK\t-" >> "${REPORT}"
  else
    echo "[FAIL] ${prefix} (exit=${rc})" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${out}\tFAIL\texit_${rc}" >> "${REPORT}"
  fi
done

echo "============================================================"
echo "[INFO] End    : $(date)"
echo "[INFO] Report : ${REPORT}"
echo "============================================================"

```

</details>

<details>
<summary><strong>03_plasx_predict_all.slurm</strong></summary>

**Notas personales:**
- Predicción final (scoring) con el modelo definido en tu workflow.
- Verifica que el script apunte al modelo correcto generado/definido en el paso anterior.

```bash
#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=plasx_predict_all
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=n002
#SBATCH -o plasx_predict_all.%j.out
#SBATCH -e plasx_predict_all.%j.err

set -euo pipefail
shopt -s nullglob

THREADS="${SLURM_CPUS_PER_TASK:-16}"

BASE="/home/databases/comau/00_anvio_comau"
GENES_DIR="${BASE}/plasx/gene_calls"
COGS_PFAM_DIR="${BASE}/functions_cogs_pfam"
DENOVO_DIR="${BASE}/plasx/plasx_denovo"
PREDICT_DIR="${BASE}/plasx/plasx_predict"
LOG_DIR="${BASE}/plasx/logs_predict"
TMPBASE="${BASE}/plasx/tmp"

# Crear directorios ANTES de usarlos
mkdir -p "$PREDICT_DIR" "$LOG_DIR" "$TMPBASE"

REPORT="${PREDICT_DIR}/plasx_predict_report.tsv"

export TMPDIR="${TMPBASE}/${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"

echo "============================================================"
echo "[INFO] JobID      : ${SLURM_JOB_ID:-NA}"
echo "[INFO] NodeList   : ${SLURM_NODELIST:-NA}"
echo "[INFO] Host       : $(hostname)"
echo "[INFO] Start      : $(date)"
echo "[INFO] Threads    : ${THREADS}"
echo "[INFO] GENES      : ${GENES_DIR}"
echo "[INFO] COGS_PFAM  : ${COGS_PFAM_DIR}"
echo "[INFO] DENOVO     : ${DENOVO_DIR}"
echo "[INFO] PREDICT    : ${PREDICT_DIR}"
echo "[INFO] LOGS       : ${LOG_DIR}"
echo "[INFO] TMPDIR     : ${TMPDIR}"
echo "============================================================"

# Activa el ambiente de PlasX (sin errores de conda init)
module load Miniconda 2>/dev/null || true

# Método alternativo para activar conda sin "conda activate"
if [[ -f "${HOME}/.conda/envs/plasx/bin/activate" ]]; then
  source "${HOME}/.conda/envs/plasx/bin/activate"
elif [[ -f "${CONDA_PREFIX}/etc/profile.d/conda.sh" ]]; then
  source "${CONDA_PREFIX}/etc/profile.d/conda.sh"
  conda activate plasx 2>/dev/null || true
fi

PLASX_BIN="$(command -v plasx || true)"
if [[ -z "${PLASX_BIN}" ]]; then
  echo "[ERROR] plasx no está en PATH. Activa el conda env correcto." >&2
  exit 1
fi
echo "[INFO] plasx      : ${PLASX_BIN}"

# Obtener lista de archivos gene-calls
gene_files=( "${GENES_DIR}"/*-gene-calls.txt )
if (( ${#gene_files[@]} == 0 )); then
  echo "[ERROR] No hay archivos *-gene-calls.txt en ${GENES_DIR}" >&2
  exit 1
fi

# Crear reporte
echo -e "prefix\tgene_calls\tcogs_pfam\tdenovo\tscores_out\tstatus\tmsg" > "${REPORT}"

total_samples=0
success_count=0
skip_count=0
fail_count=0

for gc in "${gene_files[@]}"; do
  prefix="$(basename "${gc}" -gene-calls.txt)"
  total_samples=$((total_samples + 1))
  
  # Definir rutas de archivos necesarios
  cogs_pfam="${COGS_PFAM_DIR}/${prefix}-cogs-and-pfams.txt"
  denovo="${DENOVO_DIR}/${prefix}-de-novo-families.txt"
  scores_out="${PREDICT_DIR}/${prefix}-scores.txt"
  log="${LOG_DIR}/${prefix}.predict.log"

  echo "------------------------------------------------------------"
  echo "[INFO] SAMPLE ${total_samples}: ${prefix}"
  echo "[INFO] GENE       : ${gc}"
  echo "[INFO] COGS_PFAM  : ${cogs_pfam}"
  echo "[INFO] DENOVO     : ${denovo}"
  echo "[INFO] OUT        : ${scores_out}"
  echo "[INFO] LOG        : ${log}"
  echo "------------------------------------------------------------"

  # Validar que existan todos los archivos necesarios
  missing_files=""
  if [[ ! -f "${cogs_pfam}" ]]; then
    missing_files="${missing_files}cogs-pfam "
  fi
  if [[ ! -f "${denovo}" ]]; then
    missing_files="${missing_files}denovo "
  fi

  if [[ -n "${missing_files}" ]]; then
    echo "[SKIP] Archivos faltantes para ${prefix}: ${missing_files}" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${cogs_pfam}\t${denovo}\t${scores_out}\tSKIP\tmissing_files:${missing_files}" >> "${REPORT}"
    skip_count=$((skip_count + 1))
    continue
  fi

  # Ejecutar plasx predict
  set +e
  "${PLASX_BIN}" predict \
      -a "${cogs_pfam}" "${denovo}" \
      -g "${gc}" \
      -o "${scores_out}" \
      --overwrite \
      > "${log}" 2>&1
  rc=$?
  set -e

  if [[ $rc -eq 0 ]]; then
    echo "[OK] ${prefix}" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${cogs_pfam}\t${denovo}\t${scores_out}\tOK\t-" >> "${REPORT}"
    success_count=$((success_count + 1))
  else
    echo "[FAIL] ${prefix} (exit=${rc})" | tee -a "${log}"
    echo -e "${prefix}\t${gc}\t${cogs_pfam}\t${denovo}\t${scores_out}\tFAIL\texit_${rc}" >> "${REPORT}"
    fail_count=$((fail_count + 1))
  fi
done

echo "============================================================"
echo "[INFO] End        : $(date)"
echo "[INFO] Report     : ${REPORT}"
echo "============================================================"
echo "[SUMMARY] Total samples   : ${total_samples}"
echo "[SUMMARY] Success         : ${success_count}"
echo "[SUMMARY] Skipped         : ${skip_count}"
echo "[SUMMARY] Failed          : ${fail_count}"
echo "============================================================"

# Mostrar primeras líneas del reporte si existe
if [[ -f "${REPORT}" ]]; then
  echo ""
  echo "[INFO] Primeras 10 líneas del reporte:"
  head -n 10 "${REPORT}"
fi

```

</details>

<details>
<summary><strong>06_analyze_plasmids.py</strong></summary>

**Notas personales:**
- Consolidación y parsing de resultados; aquí defines filtros, umbrales y salidas tabulares/figuras.
- Si quieres “trazar” contigs → MAGs, suele ser el lugar más cómodo para cruzar tablas.

```python
#!/usr/bin/env python3
"""
Analiza los resultados de PlasX y genera resúmenes de contigs plasmídicos
El output de 'plasx predict' ya contiene los scores por contig directamente
"""

import pandas as pd
import glob
import os
from pathlib import Path
import sys

# Directorios
BASE = "/home/databases/comau/00_anvio_comau"
SCORES_DIR = f"{BASE}/plasx/plasx_predict"
GENE_CALLS_DIR = f"{BASE}/plasx/gene_calls"
OUTPUT_DIR = f"{BASE}/plasx/plasx_summary"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Umbrales para clasificar como plásmido (puedes ajustarlo)
PLASMID_THRESHOLD_STANDARD = 0.5  # Recomendado por PlasX
PLASMID_THRESHOLD_HIGH_CONF = 0.9  # Alta confianza

def get_gene_count_for_contigs(prefix):
    """Cuenta el número de genes por contig desde gene-calls.txt"""
    gene_calls_file = f"{GENE_CALLS_DIR}/{prefix}-gene-calls.txt"
    
    if not os.path.exists(gene_calls_file):
        return {}
    
    try:
        df = pd.read_csv(gene_calls_file, sep='\t')
        # Contar genes por contig
        gene_counts = df.groupby('contig').size().to_dict()
        return gene_counts
    except Exception as e:
        print(f"[WARNING] No se pudo leer {gene_calls_file}: {e}")
        return {}

def calculate_contig_length(df_genes):
    """Calcula la longitud aproximada de cada contig basándose en coordenadas de genes"""
    if df_genes.empty:
        return {}
    
    contig_lengths = {}
    for contig in df_genes['contig'].unique():
        contig_genes = df_genes[df_genes['contig'] == contig]
        # La longitud es aproximadamente el máximo entre start y stop
        max_coord = max(contig_genes['start'].max(), contig_genes['stop'].max())
        contig_lengths[contig] = max_coord
    
    return contig_lengths

def analyze_plasmid_scores(threshold=PLASMID_THRESHOLD_STANDARD):
    """Analiza todos los archivos de scores de PlasX"""
    
    all_results = []
    all_plasmids_detail = []
    
    score_files = glob.glob(f"{SCORES_DIR}/*-scores.txt")
    
    if len(score_files) == 0:
        print(f"[ERROR] No se encontraron archivos *-scores.txt en {SCORES_DIR}")
        sys.exit(1)
    
    print(f"[INFO] Encontrados {len(score_files)} archivos de scores")
    print(f"[INFO] Usando umbral de score: {threshold}")
    print("="*60)
    
    for score_file in sorted(score_files):
        prefix = os.path.basename(score_file).replace("-scores.txt", "")
        
        try:
            df = pd.read_csv(score_file, sep='\t')
            
            # Obtener información de genes por contig
            gene_counts = get_gene_count_for_contigs(prefix)
            
            # Añadir número de genes al dataframe
            df['num_genes'] = df['contig'].map(gene_counts)
            df['num_genes'] = df['num_genes'].fillna(0).astype(int)
            
            # Clasificar contigs según el umbral
            plasmid_contigs = df[df['score'] >= threshold].copy()
            chromosome_contigs = df[df['score'] < threshold].copy()
            
            # Clasificar por nivel de confianza
            high_conf = df[df['score'] >= PLASMID_THRESHOLD_HIGH_CONF]
            medium_conf = df[(df['score'] >= 0.7) & (df['score'] < 0.9)]
            low_conf = df[(df['score'] >= 0.5) & (df['score'] < 0.7)]
            ambiguous = df[(df['score'] >= 0.4) & (df['score'] < 0.5)]
            
            result = {
                'sample': prefix,
                'total_contigs': len(df),
                'plasmid_contigs': len(plasmid_contigs),
                'chromosome_contigs': len(chromosome_contigs),
                'plasmid_percentage': (len(plasmid_contigs) / len(df) * 100) if len(df) > 0 else 0,
                'high_confidence_plasmids': len(high_conf),
                'medium_confidence': len(medium_conf),
                'low_confidence': len(low_conf),
                'ambiguous_contigs': len(ambiguous),
                'mean_plasmid_score': plasmid_contigs['score'].mean() if len(plasmid_contigs) > 0 else 0,
                'median_plasmid_score': plasmid_contigs['score'].median() if len(plasmid_contigs) > 0 else 0,
                'mean_chromosome_score': chromosome_contigs['score'].mean() if len(chromosome_contigs) > 0 else 0,
                'total_plasmid_genes': plasmid_contigs['num_genes'].sum() if len(plasmid_contigs) > 0 else 0,
                'avg_genes_per_plasmid': plasmid_contigs['num_genes'].mean() if len(plasmid_contigs) > 0 else 0
            }
            
            all_results.append(result)
            
            # Guardar detalles de plásmidos con categoría de confianza
            if len(plasmid_contigs) > 0:
                plasmid_contigs['sample'] = prefix
                plasmid_contigs['confidence'] = pd.cut(
                    plasmid_contigs['score'],
                    bins=[0.5, 0.7, 0.9, 1.0],
                    labels=['low', 'medium', 'high'],
                    include_lowest=True
                )
                all_plasmids_detail.append(plasmid_contigs[['sample', 'contig', 'score', 'num_genes', 'confidence']])
                
                # Guardar lista de contigs plasmídicos por muestra (con scores y genes)
                plasmid_list_file = f"{OUTPUT_DIR}/{prefix}-plasmid-contigs.tsv"
                plasmid_contigs[['contig', 'score', 'num_genes', 'confidence']].to_csv(
                    plasmid_list_file, sep='\t', index=False
                )
                
                # Guardar solo los nombres de contigs para extracción con seqtk
                plasmid_names_file = f"{OUTPUT_DIR}/{prefix}-plasmid-contigs-names.txt"
                plasmid_contigs['contig'].to_csv(plasmid_names_file, index=False, header=False)
                
                print(f"[OK] {prefix}: {len(plasmid_contigs)} plásmidos totales")
                print(f"     └─ Alta confianza (≥0.9): {len(high_conf)}")
                print(f"     └─ Media confianza (0.7-0.9): {len(medium_conf)}")
                print(f"     └─ Baja confianza (0.5-0.7): {len(low_conf)}")
                
                if len(ambiguous) > 0:
                    print(f"     └─ Ambiguos (0.4-0.5): {len(ambiguous)} [revisar manualmente]")
            else:
                print(f"[INFO] {prefix}: 0 plásmidos detectados")
                
        except Exception as e:
            print(f"[ERROR] {prefix}: {e}")
            continue
    
    # Crear resumen general
    summary_df = pd.DataFrame(all_results)
    summary_file = f"{OUTPUT_DIR}/plasmid_classification_summary.tsv"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    
    # Crear archivo consolidado de todos los plásmidos
    if all_plasmids_detail:
        all_plasmids_df = pd.concat(all_plasmids_detail, ignore_index=True)
        all_plasmids_file = f"{OUTPUT_DIR}/all_plasmids_consolidated.tsv"
        all_plasmids_df.to_csv(all_plasmids_file, sep='\t', index=False)
        print(f"\n[INFO] Archivo consolidado de plásmidos: {all_plasmids_file}")
    
    print("\n" + "="*60)
    print("[RESUMEN GENERAL]")
    print("="*60)
    print(f"Total de muestras analizadas: {len(summary_df)}")
    print(f"Muestras con plásmidos: {len(summary_df[summary_df['plasmid_contigs'] > 0])}")
    print(f"Total de contigs plasmídicos: {summary_df['plasmid_contigs'].sum():.0f}")
    print(f"  └─ Alta confianza (≥0.9): {summary_df['high_confidence_plasmids'].sum():.0f}")
    print(f"  └─ Media confianza (0.7-0.9): {summary_df['medium_confidence'].sum():.0f}")
    print(f"  └─ Baja confianza (0.5-0.7): {summary_df['low_confidence'].sum():.0f}")
    print(f"Contigs ambiguos (0.4-0.5): {summary_df['ambiguous_contigs'].sum():.0f}")
    print(f"Promedio de plásmidos por muestra: {summary_df['plasmid_contigs'].mean():.2f}")
    print(f"Promedio de genes por plásmido: {summary_df['avg_genes_per_plasmid'].mean():.2f}")
    print(f"\nResumen guardado en: {summary_file}")
    print("="*60)
    
    # Mostrar top 10 muestras con más plásmidos
    print("\n[TOP 10 muestras con más plásmidos]")
    top10 = summary_df.nlargest(10, 'plasmid_contigs')[
        ['sample', 'plasmid_contigs', 'high_confidence_plasmids', 'total_contigs', 'plasmid_percentage']
    ]
    print(top10.to_string(index=False))
    
    # Estadísticas de confianza
    print("\n[DISTRIBUCIÓN POR NIVEL DE CONFIANZA]")
    conf_stats = pd.DataFrame({
        'Nivel': ['Alta (≥0.9)', 'Media (0.7-0.9)', 'Baja (0.5-0.7)', 'Ambigua (0.4-0.5)'],
        'Total': [
            summary_df['high_confidence_plasmids'].sum(),
            summary_df['medium_confidence'].sum(),
            summary_df['low_confidence'].sum(),
            summary_df['ambiguous_contigs'].sum()
        ]
    })
    print(conf_stats.to_string(index=False))
    
    return summary_df

def extract_plasmid_sequences():
    """Genera comandos para extraer secuencias de plásmidos desde los contigs originales"""
    
    print("\n" + "="*60)
    print("[COMANDOS PARA EXTRAER SECUENCIAS DE PLÁSMIDOS]")
    print("="*60)
    
    plasmid_lists = glob.glob(f"{OUTPUT_DIR}/*-plasmid-contigs-names.txt")
    
    if len(plasmid_lists) == 0:
        print("[INFO] No hay listas de plásmidos para extraer")
        return
    
    extract_script = f"{OUTPUT_DIR}/extract_plasmid_sequences.sh"
    
    with open(extract_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Script para extraer secuencias de contigs plasmídicos identificados por PlasX\n")
        f.write("# Generado automáticamente por 06_analyze_plasmids.py\n\n")
        f.write("set -e\n\n")
        f.write("CONTIGS_BASE=\"/home/databases/comau/00_anvio_comau/02_contigs\"\n")
        f.write(f"PLASMID_SEQS=\"{OUTPUT_DIR}/plasmid_sequences\"\n")
        f.write("mkdir -p \"$PLASMID_SEQS\"\n\n")
        f.write("echo \"[INFO] Extrayendo secuencias de plásmidos...\"\n")
        f.write("echo \"======================================\"\n\n")
        
        extracted = 0
        skipped = 0
        
        for plasmid_list in sorted(plasmid_lists):
            prefix = os.path.basename(plasmid_list).replace("-plasmid-contigs-names.txt", "")
            
            # Verificar cuántos plásmidos hay
            with open(plasmid_list) as pl:
                num_plasmids = len(pl.readlines())
            
            f.write(f"# Extraer plásmidos de {prefix} ({num_plasmids} contigs)\n")
            f.write(f"if [ -f \"$CONTIGS_BASE/{prefix}-contigs.fa\" ]; then\n")
            f.write(f"  seqtk subseq \"$CONTIGS_BASE/{prefix}-contigs.fa\" \\\n")
            f.write(f"    \"{plasmid_list}\" \\\n")
            f.write(f"    > \"$PLASMID_SEQS/{prefix}-plasmids.fa\"\n")
            f.write(f"  echo \"[OK] {prefix}: {num_plasmids} plásmidos extraídos\"\n")
            f.write("else\n")
            f.write(f"  echo \"[SKIP] No encontrado: $CONTIGS_BASE/{prefix}-contigs.fa\"\n")
            f.write("fi\n\n")
            
            extracted += 1
    
    os.chmod(extract_script, 0o755)
    print(f"Script de extracción generado: {extract_script}")
    print(f"  └─ Muestras a procesar: {extracted}")
    print("\nPara extraer las secuencias FASTA, ejecuta:")
    print(f"  bash {extract_script}")

def generate_high_confidence_summary():
    """Genera resumen adicional de plásmidos de alta confianza"""
    
    print("\n" + "="*60)
    print("[ANÁLISIS DE PLÁSMIDOS DE ALTA CONFIANZA (score ≥ 0.9)]")
    print("="*60)
    
    score_files = glob.glob(f"{SCORES_DIR}/*-scores.txt")
    high_conf_dir = f"{OUTPUT_DIR}/high_confidence"
    os.makedirs(high_conf_dir, exist_ok=True)
    
    total_high_conf = 0
    samples_with_high_conf = 0
    
    for score_file in sorted(score_files):
        prefix = os.path.basename(score_file).replace("-scores.txt", "")
        
        try:
            df = pd.read_csv(score_file, sep='\t')
            high_conf_plasmids = df[df['score'] >= PLASMID_THRESHOLD_HIGH_CONF].copy()
            
            if len(high_conf_plasmids) > 0:
                # Añadir número de genes
                gene_counts = get_gene_count_for_contigs(prefix)
                high_conf_plasmids['num_genes'] = high_conf_plasmids['contig'].map(gene_counts)
                high_conf_plasmids['num_genes'] = high_conf_plasmids['num_genes'].fillna(0).astype(int)
                
                out_file = f"{high_conf_dir}/{prefix}-high-conf-plasmids.tsv"
                high_conf_plasmids[['contig', 'score', 'num_genes']].to_csv(out_file, sep='\t', index=False)
                
                total_high_conf += len(high_conf_plasmids)
                samples_with_high_conf += 1
        except:
            continue
    
    print(f"Muestras con plásmidos de alta confianza: {samples_with_high_conf}")
    print(f"Total de plásmidos de alta confianza: {total_high_conf}")
    print(f"Archivos guardados en: {high_conf_dir}/")

def generate_statistics_report(summary_df):
    """Genera reporte estadístico adicional"""
    
    stats_file = f"{OUTPUT_DIR}/statistics_report.txt"
    
    with open(stats_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("REPORTE ESTADÍSTICO - PlasX\n")
        f.write("="*60 + "\n\n")
        
        f.write("[ESTADÍSTICAS GENERALES]\n")
        f.write(f"Muestras analizadas: {len(summary_df)}\n")
        f.write(f"Total de contigs: {summary_df['total_contigs'].sum():.0f}\n")
        f.write(f"Total de plásmidos: {summary_df['plasmid_contigs'].sum():.0f}\n")
        f.write(f"Porcentaje global de plásmidos: {(summary_df['plasmid_contigs'].sum() / summary_df['total_contigs'].sum() * 100):.2f}%\n\n")
        
        f.write("[DISTRIBUCIÓN POR MUESTRA]\n")
        f.write(f"Promedio de plásmidos: {summary_df['plasmid_contigs'].mean():.2f}\n")
        f.write(f"Mediana de plásmidos: {summary_df['plasmid_contigs'].median():.2f}\n")
        f.write(f"Desviación estándar: {summary_df['plasmid_contigs'].std():.2f}\n")
        f.write(f"Mínimo: {summary_df['plasmid_contigs'].min():.0f}\n")
        f.write(f"Máximo: {summary_df['plasmid_contigs'].max():.0f}\n\n")
        
        f.write("[SCORES]\n")
        f.write(f"Score promedio plásmidos: {summary_df['mean_plasmid_score'].mean():.4f}\n")
        f.write(f"Score promedio cromosomas: {summary_df['mean_chromosome_score'].mean():.4f}\n\n")
        
        f.write("[GENES]\n")
        f.write(f"Genes totales en plásmidos: {summary_df['total_plasmid_genes'].sum():.0f}\n")
        f.write(f"Promedio de genes por plásmido: {summary_df['avg_genes_per_plasmid'].mean():.2f}\n")
    
    print(f"\n[INFO] Reporte estadístico guardado en: {stats_file}")

if __name__ == "__main__":
    print("="*60)
    print("ANÁLISIS DE CLASIFICACIÓN DE PLÁSMIDOS - PlasX")
    print("="*60)
    print()
    
    # Análisis principal con umbral estándar (0.5)
    summary = analyze_plasmid_scores(threshold=PLASMID_THRESHOLD_STANDARD)
    
    # Generar script de extracción
    extract_plasmid_sequences()
    
    # Análisis de plásmidos de alta confianza (0.9)
    generate_high_confidence_summary()
    
    # Generar reporte estadístico
    generate_statistics_report(summary)
    
    print("\n" + "="*60)
    print("[COMPLETADO]")
    print("="*60)
    print("\nArchivos generados:")
    print(f"  1. {OUTPUT_DIR}/plasmid_classification_summary.tsv")
    print(f"     └─ Resumen por muestra con estadísticas")
    print(f"  2. {OUTPUT_DIR}/all_plasmids_consolidated.tsv")
    print(f"     └─ Todos los plásmidos de todas las muestras")
    print(f"  3. {OUTPUT_DIR}/[SAMPLE]-plasmid-contigs.tsv")
    print(f"     └─ Plásmidos por muestra con scores y genes")
    print(f"  4. {OUTPUT_DIR}/[SAMPLE]-plasmid-contigs-names.txt")
    print(f"     └─ Nombres de contigs para seqtk")
    print(f"  5. {OUTPUT_DIR}/extract_plasmid_sequences.sh")
    print(f"     └─ Script bash para extraer FASTAs")
    print(f"  6. {OUTPUT_DIR}/high_confidence/")
    print(f"     └─ Plásmidos de alta confianza (score ≥ 0.9)")
    print(f"  7. {OUTPUT_DIR}/statistics_report.txt")
    print(f"     └─ Reporte estadístico completo")
    print("\n" + "="*60)

```

</details>

## Notas operativas rápidas (personales)
- Mantén consistente la nomenclatura de muestras/contigs entre los pasos (especialmente si hay concatenación o demultiplex).
- Si tu entrada a PlasX es **contigs** vs **genes**, deja explícito qué archivo alimenta cada paso (y por qué).
- Cuando integres con MAGs: guarda un mapeo contig→bin/MAG (TSV) para trazabilidad.
- Versiona (en el README o en un `env.yml`) la versión exacta de PlasX, Anvi’o y dependencias.

## Estructura de directorios sugerida (personal)
```text
plasx_pipeline/
├── 01_anvio_export/
├── 02_plasx_denovo/
├── 03_plasx_predict/
├── 04_analysis/
└── scripts/
    ├── 01_anvio_export_gene_calls_all_n002.slurm
    ├── 02_plasx_denovo_all_test.slurm
    ├── 03_plasx_predict_all.slurm
    └── 06_analyze_plasmids.py
```

## Checklist antes de correr
- [ ] Confirmar rutas (contigs, contigs.db, listas, directorios de salida).
- [ ] Confirmar módulos/conda env en el cluster (PlasX, Anvi’o, Python).
- [ ] Confirmar recursos Slurm (CPU/RAM/tiempo) según tamaño real de tus datos.
- [ ] Validar que el modelo (si aplica) exista antes del paso `predict`.
- [ ] Ejecutar un test mínimo con un subconjunto para validar I/O y formato.

