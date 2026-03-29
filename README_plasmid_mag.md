# Análisis de Asociación Plásmido-MAG

## Descripción

Este pipeline determina qué contigs plasmídicos (detectados por geNomad) están insertados en los genomas reconstruidos (MAGs - Metagenome-Assembled Genomes) para cada muestra.

## Lógica de matching

Los contigs tienen IDs con diferente formato en cada archivo:

| Archivo | Ejemplo ID de contig |
|---------|---------------------|
| `FP.DEC2017_L1S1D05.plasmids.fa` | `k119_100265` |
| `bin.bin.9.FP.DEC2017_L1S1D05.fa` | `EC2017L1S05_k119_100265` |

El prefijo en los MAGs (`EC2017L1S05_`) es el identificador de muestra que se añade durante el proceso de binning. El script extrae el patrón `k[kmer]_[número]` del final del ID para hacer el matching entre ambos archivos.

## Archivos incluidos

- `plasmid_mag_analysis.py` — Script principal de análisis
- `plot_seasonal_patterns.py` — Script de visualización estacional

## Uso

### Paso 1: Análisis principal

```bash
# Uso básico — genera TSV + FASTA
python plasmid_mag_analysis.py \
    --plasmids lista_plasmids2.txt \
    --mags lista_mags.txt \
    --out plasmid_mag_analysis/

# Solo TSV (más rápido, sin cargar secuencias en memoria)
python plasmid_mag_analysis.py \
    --plasmids lista_plasmids2.txt \
    --mags lista_mags.txt \
    --out plasmid_mag_analysis/ \
    --no-fasta
```

### Paso 2: Visualización

```bash
python plot_seasonal_patterns.py \
    --input plasmid_mag_analysis/summary_metrics.tsv \
    --out plasmid_mag_analysis/figures/
```

**Requisitos:** `pip install pandas matplotlib numpy`

---

## Archivos de salida

### `plasmid_mag_analysis/plasmids_in_MAGs.fa`  ⭐
FASTA de los contigs plasmídicos **asociados a algún MAG**.
El header original se respeta íntegramente; solo se añade al final una anotación entre corchetes:
```
>k119_100265 evidence=weak_evidence methods=genomad [sample=FP.DEC2017_L1S1D05][MAG_bin=bin.bin.9.FP.DEC2017_L1S1D05.fa][MAG_contig=EC2017L1S05_k119_100265]
ATCGATCGATCG...
```

### `plasmid_mag_analysis/plasmids_NOT_in_MAGs.fa`  ⭐
FASTA de los contigs plasmídicos **NO asociados a ningún MAG**.
Solo se añade `[sample=...]` al final del header original:
```
>k119_999999 evidence=strong_evidence methods=genomad [sample=FP.DEC2017_L1S1D05]
GCTAGCTAGCTA...
```

---

### `plasmid_mag_analysis/summary_metrics.tsv`
Tabla principal con métricas por muestra, ordenada temporalmente. Columnas:

| Columna | Descripción |
|---------|-------------|
| `sample_id` | ID de la muestra |
| `date_label` | Etiqueta temporal (ej: DEC2017) |
| `year` | Año |
| `month_num` | Número de mes (01-12) |
| `month_name` | Nombre del mes |
| `season` | Estación del año (Hemisferio Sur) |
| `sort_key` | Clave de ordenamiento temporal |
| `n_mag_files` | Número de archivos MAG para esa muestra |
| `total_plasmids` | Total de contigs plasmídicos |
| `plasmids_in_mag` | Plásmidos encontrados en algún MAG |
| `plasmids_not_in_mag` | Plásmidos NO en ningún MAG |
| `pct_in_mag` | % de plásmidos en MAGs |
| `pct_not_in_mag` | % de plásmidos NO en MAGs |
| `status` | ok / no_mags / file_not_accessible |

### `plasmid_mag_analysis/plasmids_in_MAGs.tsv`
Listado detallado de cada contig plasmídico encontrado en un MAG:
- `sample_id`, fecha, estación
- `plasmid_contig_base` — ID base del contig (ej: `k119_100265`)
- `plasmid_contig_full` — ID completo del plásmido
- `mag_bin` — nombre del archivo del bin MAG donde se encontró
- `mag_contig_full` — ID completo del contig en el MAG

### `plasmid_mag_analysis/plasmids_NOT_in_MAGs.tsv`
Listado de contigs plasmídicos que NO se encontraron en ningún MAG.

### `plasmid_mag_analysis/samples_without_mags.txt`
Muestras de plásmidos sin archivo MAG correspondiente.

---

## Figuras generadas

- **`seasonal_overview.png`** — 4 paneles: total plásmidos, barras apiladas, % en MAG, boxplot estacional
- **`timeseries_seasonal.png`** — Serie temporal con fondo sombreado por estación
- **`season_statistics.png`** — Estadísticas agregadas por estación (media ± SD)

---

## Estaciones (Hemisferio Sur)

| Mes | Estación |
|-----|---------|
| Diciembre, Enero, Febrero | **Verano** |
| Marzo, Abril, Mayo | **Otoño** |
| Junio, Julio, Agosto | **Invierno** |
| Septiembre, Octubre, Noviembre | **Primavera** |

---

## Notas técnicas

- Las muestras `EC-Bac*` y `BAC4*` se etiquetan como `Lab_Control`
- Las muestras `L151`, `L152`, etc., se etiquetan como `L-series`
- El script es robusto a archivos no encontrados (reporta advertencias sin detenerse)
- Para graficar patrones estacionales, se recomienda usar solo las muestras con `status = ok`
