##1 [plasX](https://github.com/michaelkyu/PlasX)

```
anvi-gen-contigs-database -L 0 -T $THREADS --project-name $PREFIX -f unique_contigs.fasta -o $PREFIX.db
(anvio-8) mmorenos@mmorenos:012_eduCastroLab$ nohup anvi-run-ncbi-cogs   -c unique_contigs_CLAB.db  --cog-data-dir /home/mmorenos/github/anvio/anvio/data/misc/COG/COG24/COG14/ --cog-version COG14 &
```
# Manual COG14 functional annotation pipeline with DIAMOND and anvi'o

Este documento describe el pipeline completo que usaste para anotar un `contigs.db` grande con COG14 usando DIAMOND **por fuera** de `anvi-run-ncbi-cogs`, y luego importar las funciones de vuelta a anvi'o.

El flujo es:

1. Exportar secuencias de aminoácidos desde el `contigs.db`.
2. Ejecutar DIAMOND contra la base COG14 (`COG.dmnd`) con control de threads y `TMPDIR`.
3. Usar los archivos RAW de COG14 de NCBI para mapear `gi → COG_ID → descripción`.
4. Construir un archivo `*.functions.txt` en el formato que exige `anvi-import-functions`.
5. Convertir ese archivo a UTF-8 para evitar errores de encoding.
6. Importar las funciones en el `contigs.db`.

Las rutas y nombres de archivo son totalmente editables, pero están inicializados a tu caso de `unique_contigs_CLAB`.

---

## 0. Variables base del proyecto

Crea un archivo `00_config_cog14.sh` con las variables que usarán todos los scripts.

```bash
#!/usr/bin/env bash
set -euo pipefail

########## CONFIGURACIÓN GLOBAL ##########

# Nombre base del proyecto (contigs.db = PREFIX.db)
export PREFIX="unique_contigs_CLAB"

# Contigs database de anvi'o
export CONTIGS_DB="${PREFIX}.db"

# Fichero FASTA de aminoácidos que generaremos
export AA_FASTA="${PREFIX}-aa.fa"

# Directorio donde vive la instalación local de COG14
# (ajusta a tu ruta real)
export COG14_DIR="/media/mmorenos/disk3/012_eduCastroLab/COG_2014/COG14"

# Base de datos DIAMOND de COG14
export COG_DB="${COG14_DIR}/DB_DIAMOND/COG.dmnd"

# Directorio con los RAW_DATA_FROM_NCBI de COG14
export RAW_DIR="${COG14_DIR}/RAW_DATA_FROM_NCBI"

# Directorio temporal grande para DIAMOND
export TMPDIR_BASE="/media/mmorenos/disk3/012_eduCastroLab/tmp_grande_cog14"

# Archivo de salida de DIAMOND (formato m6)
export M6_FILE="${PREFIX}-COG14.m6"

# Archivo intermedio: GI -> COG_ID + descripción
export GI_COG_DESC="${PREFIX}-COG14_gi_cog_desc.tsv"

# Archivo final de funciones (antes de limpiar encoding)
export FUN_FILE="${PREFIX}-COG14.functions.txt"

# Archivo final de funciones en UTF-8
export FUN_FILE_UTF8="${PREFIX}-COG14.functions.utf8.txt"

########## FIN CONFIGURACIÓN ##########

echo "[INFO] Configuración cargada:"
echo "  PREFIX        = ${PREFIX}"
echo "  CONTIGS_DB    = ${CONTIGS_DB}"
echo "  AA_FASTA      = ${AA_FASTA}"
echo "  COG14_DIR     = ${COG14_DIR}"
echo "  COG_DB        = ${COG_DB}"
echo "  RAW_DIR       = ${RAW_DIR}"
echo "  TMPDIR_BASE   = ${TMPDIR_BASE}"
echo "  M6_FILE       = ${M6_FILE}"
echo "  GI_COG_DESC   = ${GI_COG_DESC}"
echo "  FUN_FILE      = ${FUN_FILE}"
echo "  FUN_FILE_UTF8 = ${FUN_FILE_UTF8}"
```

Uso:

```bash
chmod +x 00_config_cog14.sh
source 00_config_cog14.sh
```

Asegúrate de `source`-ar este archivo en cada sesión antes de correr los scripts siguientes.

---

## 1. Exportar secuencias de aminoácidos desde el contigs.db

Script `01_get_aa_sequences.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Cargar configuración
source 00_config_cog14.sh

# Activar el entorno con anvi'o
conda activate anvio-8

if [[ ! -f "${CONTIGS_DB}" ]]; then
    echo "[ERROR] No se encuentra el contigs.db: ${CONTIGS_DB}" >&2
    exit 1
fi

echo "============================================="
echo "[INFO] Exportando secuencias AA desde ${CONTIGS_DB}"
echo "============================================="

anvi-get-sequences-for-gene-calls \
  -c "${CONTIGS_DB}" \
  --get-aa-sequences \
  -o "${AA_FASTA}"

echo "[OK] FASTA de AA generado: ${AA_FASTA}"
grep -c ">" "${AA_FASTA}" || true
```

Uso:

```bash
chmod +x 01_get_aa_sequences.sh
./01_get_aa_sequences.sh
```

---

## 2. Ejecutar DIAMOND contra COG14 (10 threads, TMP controlado)

Script `02_run_cog14_diamond.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Cargar configuración
source 00_config_cog14.sh

# Activar el entorno donde está DIAMOND
conda activate anvio-8

THREADS=10

if [[ ! -f "${AA_FASTA}" ]]; then
    echo "[ERROR] No se encuentra el FASTA de AA: ${AA_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${COG_DB}" ]]; then
    echo "[ERROR] No se encuentra la base DIAMOND COG14: ${COG_DB}" >&2
    echo "Asegúrate de haber corrido anvio-setup-ncbi-cogs o de que la ruta es correcta." >&2
    exit 1
fi

mkdir -p "${TMPDIR_BASE}"

echo "============================================="
echo "[INFO] Corriendo DIAMOND contra COG14"
echo "[INFO] AA_FASTA    : ${AA_FASTA}"
echo "[INFO] COG_DB      : ${COG_DB}"
echo "[INFO] Threads     : ${THREADS}"
echo "[INFO] TMPDIR_BASE : ${TMPDIR_BASE}"
echo "[INFO] Output (m6) : ${M6_FILE}"
echo "============================================="

diamond blastp \
  -q "${AA_FASTA}" \
  -d "${COG_DB}" \
  -o "${M6_FILE}" \
  -p "${THREADS}" \
  --outfmt 6 qseqid sseqid pident length mismatches gapopen qstart qend sstart send evalue bitscore \
  --max-target-seqs 1 \
  --evalue 1e-5 \
  -t "${TMPDIR_BASE}"

echo "[OK] DIAMOND finalizado. Salida m6: ${M6_FILE}"
head -n 5 "${M6_FILE}" || true
```

Uso:

```bash
chmod +x 02_run_cog14_diamond.sh
nohup ./02_run_cog14_diamond.sh > 02_run_cog14_diamond.log 2>&1 &
tail -f 02_run_cog14_diamond.log
```

---

## 3. Construir tabla GI → COG_ID + descripción COG y unir con la salida m6

Aquí explotamos los archivos RAW de COG14 que ya tienes:

- `cognames2003-2014.tab` (COG_ID, categoría, descripción)
- `cog2003-2014.csv` (gi, ..., COG_ID, ...)

Script `03_make_COG14_functions_from_m6.sh` (hace todo: construye el mapa y genera `*.functions.txt`):

```bash
#!/usr/bin/env bash
set -euo pipefail

# Cargar configuración
source 00_config_cog14.sh

COG_NAMES="${RAW_DIR}/cognames2003-2014.tab"
COG_CSV="${RAW_DIR}/cog2003-2014.csv"

for f in "${M6_FILE}" "${COG_NAMES}" "${COG_CSV}"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] No se encuentra: $f" >&2
        exit 1
    fi
done

echo "============================================="
echo "[INFO] Construyendo anotaciones COG14 (GI → COG → descripción)"
echo "[INFO] M6 DIAMOND   : ${M6_FILE}"
echo "[INFO] RAW_DIR COG14: ${RAW_DIR}"
echo "[INFO] GI_COG_DESC  : ${GI_COG_DESC}"
echo "[INFO] FUN_FILE     : ${FUN_FILE}"
echo "============================================="

########## 1) GI -> COG_ID + descripción ##########

awk '
BEGIN{
    OFS = "\t";
}
FNR==1 && FILENAME ~ /cognames/ { FS = "\t"; }
FNR==1 && FILENAME ~ /cog2003-2014.csv/ { FS = ","; }

# Archivo 1: cognames2003-2014.tab
FILENAME ~ /cognames2003-2014.tab/ {
    # COG0001 \t H \t Glutamate-1-semialdehyde aminotransferase
    cog_id = $1;
    desc   = $3;
    cog_desc[cog_id] = desc;
    next;
}

# Archivo 2: cog2003-2014.csv
FILENAME ~ /cog2003-2014.csv/ {
    # 158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
    gi     = $1;
    cog_id = $7;
    if (cog_id != "" && (cog_id in cog_desc)) {
        print gi, cog_id, cog_desc[cog_id];
    }
    next;
}
' "${COG_NAMES}" "${COG_CSV}" > "${GI_COG_DESC}"

echo "[OK] GI->COG_ID+desc generado: ${GI_COG_DESC}"
head -n 3 "${GI_COG_DESC}" || true

########## 2) Unir M6 + GI_COG_DESC y generar FUN_FILE ##########

{
    printf "gene_callers_id\tsource\taccession\tfunction\te_value\n";

    awk -v OFS="\t" '
    BEGIN{
        FS = "\t";
    }

    # Primer archivo: GI_COG_DESC (gi \t COG_ID \t COG_description)
    NR==FNR {
        gi       = $1;
        cog_id   = $2;
        desc     = $3;
        gi_map[gi] = cog_id "\t" desc;
        next;
    }

    # Segundo archivo: M6_FILE (puede venir separado por tab o espacio)
    {
        # reparseamos con split sobre tabs/espacios
        n = split($0, f, /[ \t]+/);
        if (n < 11) {
            next;
        }

        qseqid = f[1];
        sseqid = f[2];
        evalue = f[11];

        gi = "";
        if (match(sseqid, /^gi\|([0-9]+)\|/, m)) {
            gi = m[1];
        }

        if (gi != "" && (gi in gi_map)) {
            split(gi_map[gi], tmp, "\t");
            cog_id = tmp[1];
            desc   = tmp[2];
            print qseqid, "COG14_FUNCTION", cog_id, desc, evalue;
        }
    }
    ' "${GI_COG_DESC}" "${M6_FILE}"
} > "${FUN_FILE}"

echo "[OK] Archivo de funciones generado: ${FUN_FILE}"
head -n 5 "${FUN_FILE}" || true
```

Uso:

```bash
chmod +x 03_make_COG14_functions_from_m6.sh
./03_make_COG14_functions_from_m6.sh
```

El resultado `unique_contigs_CLAB-COG14.functions.txt` tendrá el formato:

```text
gene_callers_id  source           accession  function                                   e_value
0                COG14_FUNCTION   COG0001    Glutamate-1-semialdehyde aminotransferase 4.4e-136
...
```

---

## 4. Convertir el archivo de funciones a UTF-8

Algunas descripciones de COG tienen comillas o caracteres especiales que vienen en Latin-1 / Windows-1252 y provocan:

> UnicodeDecodeError: 'utf-8' codec can't decode byte 0x91 ...

Por eso convertimos explícitamente a UTF-8.

Script `04_convert_functions_to_utf8.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Cargar configuración
source 00_config_cog14.sh

if [[ ! -f "${FUN_FILE}" ]]; then
    echo "[ERROR] No se encuentra el archivo de funciones: ${FUN_FILE}" >&2
    exit 1
fi

echo "============================================="
echo "[INFO] Convirtiendo funciones a UTF-8"
echo "[INFO] Input : ${FUN_FILE}"
echo "[INFO] Output: ${FUN_FILE_UTF8}"
echo "============================================="

if command -v iconv >/dev/null 2>&1; then
    # Intentar conversión desde Windows-1252 (o Latin-1) a UTF-8
    iconv -f WINDOWS-1252 -t UTF-8 -c \
        "${FUN_FILE}" > "${FUN_FILE_UTF8}"
else
    # Fallback en Python si iconv no está disponible
    python - << EOF
import io

inp  = "${FUN_FILE}"
outp = "${FUN_FILE_UTF8}"

with open(inp, "rb") as f_in, open(outp, "w", encoding="utf-8") as f_out:
    for bline in f_in:
        line = bline.decode("latin-1", errors="ignore")
        f_out.write(line)
EOF
fi

echo "[OK] Archivo en UTF-8 generado: ${FUN_FILE_UTF8}"
head -n 5 "${FUN_FILE_UTF8}" || true
```

Uso:

```bash
chmod +x 04_convert_functions_to_utf8.sh
./04_convert_functions_to_utf8.sh
```

---

## 5. Importar las funciones en el contigs.db de anvi'o

Finalmente, usamos `anvi-import-functions` con el archivo ya limpio en UTF-8.

Script `05_import_functions_into_anvio.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Cargar configuración
source 00_config_cog14.sh

conda activate anvio-8

if [[ ! -f "${CONTIGS_DB}" ]]; then
    echo "[ERROR] No se encuentra el contigs.db: ${CONTIGS_DB}" >&2
    exit 1
fi

if [[ ! -f "${FUN_FILE_UTF8}" ]]; then
    echo "[ERROR] No se encuentra el archivo UTF-8 de funciones: ${FUN_FILE_UTF8}" >&2
    exit 1
fi

echo "============================================="
echo "[INFO] Importando funciones COG14 en anvi'o"
echo "[INFO] CONTIGS_DB : ${CONTIGS_DB}"
echo "[INFO] FUN_FILE   : ${FUN_FILE_UTF8}"
echo "============================================="

anvi-import-functions \
  -c "${CONTIGS_DB}" \
  -i "${FUN_FILE_UTF8}"

echo "[OK] Importación completada."

echo "Puedes revisar las annotation sources en tu contigs.db con, por ejemplo:"
echo "  anvi-display-contigs-stats -c ${CONTIGS_DB}"
```

Uso:

```bash
chmod +x 05_import_functions_into_anvio.sh
./05_import_functions_into_anvio.sh
```

---

## 6. Resumen del pipeline completo

1. Configurar variables globales:

   ```bash
   chmod +x 00_config_cog14.sh
   source 00_config_cog14.sh
   ```

2. Exportar aminoácidos desde el `contigs.db`:

   ```bash
   ./01_get_aa_sequences.sh
   ```

3. Anotar contra COG14 con DIAMOND:

   ```bash
   ./02_run_cog14_diamond.sh
   ```

4. Construir `COG14.functions.txt` en formato anvi'o:

   ```bash
   ./03_make_COG14_functions_from_m6.sh
   ```

5. Convertir a UTF-8:

   ```bash
   ./04_convert_functions_to_utf8.sh
   ```

6. Importar las funciones en el `contigs.db`:

   ```bash
   ./05_import_functions_into_anvio.sh
   ```

Con esto tienes la versión manual pero completamente controlada de la anotación COG14_FUNCTION, compatible con tu instalación de COG14 y evitando los problemas de encoding y de `/tmp` que te estaban frenando `anvi-run-ncbi-cogs`.
