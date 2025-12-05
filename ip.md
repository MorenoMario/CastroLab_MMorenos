## Acceso servidor utal
```
in= 
ssh mmoreno@172.22.0.21

out= 
ssh mmoreno@190.110.100.21 -p 20223
```
## Datos pato
```

/home/patricio/analisis_comau_separados
(pato: ahi tengo la fase plasmidial y cromosomal)
```
## Datos Vale-bins

```
/home/comau_great_project_saint_trinity/bins

```

## Run concoct 

```
#!/usr/bin/bash

#SBATCH --partition=basic
#SBATCH -J binning_concoct
#SBATCH -e binning_concoct.err
#SBATCH -o binning_concoct.out
#SBATCH --cpus-per-task=40

# Cargar módulos necesarios
module load Miniconda

# Activar el ambiente de concoct
conda activate concoct_env

# Definir rutas
CONTIGS="/home/comau_great_project_saint_trinity/contigs_universe/unique_contigs.fasta"
BAM_DIR="/home/comau_great_project_saint_trinity/bam_files"
OUTPUT_DIR="concoct_output"

# Crear directorio de salida si no existe
mkdir -p ${OUTPUT_DIR}
mkdir -p bins

echo "Iniciando análisis de binning con CONCOCT..."
echo "Fecha y hora de inicio: $(date)"

# Paso 1: Cortar contigs en fragmentos de 10K
echo "Paso 1: Cortando contigs en fragmentos de 10K..."
cut_up_fasta.py ${CONTIGS} -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

if [ $? -eq 0 ]; then
    echo "Contigs cortados exitosamente"
else
    echo "Error al cortar contigs"
    exit 1
fi

# Paso 2: Generar tabla de cobertura
echo "Paso 2: Generando tabla de cobertura..."
concoct_coverage_table.py contigs_10K.bed ${BAM_DIR}/*.bam > coverage_table.tsv

if [ $? -eq 0 ]; then
    echo "Tabla de cobertura generada exitosamente"
else
    echo "Error al generar tabla de cobertura"
    exit 1
fi

# Paso 3: Ejecutar CONCOCT
echo "Paso 3: Ejecutando CONCOCT con 40 threads..."
concoct --threads 40 --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b ${OUTPUT_DIR}/

if [ $? -eq 0 ]; then
    echo "CONCOCT ejecutado exitosamente"
else
    echo "Error al ejecutar CONCOCT"
    exit 1
fi

# Paso 4: Merge de clustering
echo "Paso 4: Mergeando clustering..."
merge_cutup_clustering.py ${OUTPUT_DIR}/clustering_gt1000.csv > ${OUTPUT_DIR}/clustering_merged.csv

if [ $? -eq 0 ]; then
    echo "Clustering mergeado exitosamente"
else
    echo "Error al mergear clustering"
    exit 1
fi

# Paso 5: Extraer bins en archivos FASTA
echo "Paso 5: Extrayendo bins en archivos FASTA..."
extract_fasta_bins.py ${CONTIGS} ${OUTPUT_DIR}/clustering_merged.csv --output_path bins/

if [ $? -eq 0 ]; then
    echo "Bins extraídos exitosamente"
else
    echo "Error al extraer bins"
    exit 1
fi

echo "Análisis completado exitosamente!"
echo "Fecha y hora de finalización: $(date)"
echo "Los bins se encuentran en el directorio: bins/"
```
