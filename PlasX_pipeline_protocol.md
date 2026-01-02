# Pipeline de Predicción de Contigs con PlasX

Este protocolo describe los pasos para la **predicción y clasificación de contigs putativos plasmídicos** utilizando **PlasX** en entornos HPC gestionados por **Slurm**.  
Incluye scripts de exportación, predicción de novo, y análisis integrativo de resultados.

---

## 1. Exportación de anotaciones desde Anvi’o

Archivo: `01_anvio_export_gene_calls_all_n002.slurm`

Este script extrae las anotaciones génicas desde las bases de datos de contigs generadas con Anvi’o.  
Su propósito es obtener las secuencias de genes necesarias para las predicciones de PlasX.

### Comando principal:
```bash
#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=export_gene_calls
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH -o anvio_export_gene_calls.%j.out
#SBATCH -e anvio_export_gene_calls.%j.err

# Exportar las anotaciones génicas de cada contig.db
for db in $(cat list_genomes); do
    anvi-export-gene-calls -c $db -o ${db%.db}_gene_calls.fa
done
```

**Salida esperada:**  
Archivos `.fa` con las secuencias de genes para cada base de datos de contigs.

---

## 2. Predicción de novo con PlasX

Archivo: `02_plasx_denovo_all_test.slurm`

Este paso utiliza el modo *de novo* de PlasX para generar un modelo predictivo a partir de las secuencias génicas exportadas.

### Comando principal:
```bash
#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=plasx_denovo_all
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH -o plasx_denovo_all.%j.out
#SBATCH -e plasx_denovo_all.%j.err

module load plasx

PlasX denovo     --input gene_calls_all.fa     --output plasx_denovo_out     --threads 20
```

**Salida esperada:**  
Modelo de predicción entrenado (`plasx_denovo_out/`) con pesos derivados de las características de los genes.

---

## 3. Predicción de Contigs con PlasX

Archivo: `03_plasx_predict_all.slurm`

Usa el modelo generado en el paso anterior para clasificar los contigs de interés y predecir su probabilidad de ser plasmídicos.

### Comando principal:
```bash
#!/usr/bin/env bash
#SBATCH --partition=basic
#SBATCH --job-name=plasx_predict_all
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH -o plasx_predict_all.%j.out
#SBATCH -e plasx_predict_all.%j.err

module load plasx

PlasX predict     --input contigs_all.fa     --model plasx_denovo_out/model.pkl     --output plasx_predictions_all     --threads 20
```

**Salida esperada:**  
Carpeta `plasx_predictions_all/` con resultados por contig, incluyendo scores de probabilidad plasmídica.

---

## 4. Análisis de resultados

Archivo: `06_analyze_plasmids.py`

Script en Python para consolidar y visualizar los resultados de PlasX.  
Combina las predicciones de múltiples muestras, resume la clasificación de contigs y permite análisis exploratorios adicionales.

### Ejemplo de uso:
```bash
python 06_analyze_plasmids.py --input plasx_predictions_all --output summary_plasmids.tsv
```

**Salida esperada:**  
Archivo tabular `summary_plasmids.tsv` con las categorías de PlasX por contig, incluyendo la probabilidad asignada y estadísticas agregadas.

---

## 5. Notas y recomendaciones

- Verificar que las bases de datos de contigs Anvi’o contengan genes correctamente llamados antes de la exportación.  
- Se recomienda usar nodos con al menos **120–150 GB de RAM** y **≥20 hilos**.  
- Los resultados pueden integrarse posteriormente con anotaciones funcionales (COG, Pfam) o herramientas de clasificación de elementos móviles (p. ej. PlasFlow, geNomad).  
- Mantener una estructura reproducible de carpetas:
  ```
  01_anvio_export/
  02_plasx_denovo/
  03_plasx_predict/
  04_analysis/
  ```

---

**Autor:** Mario Moreno  
**Última actualización:** Enero 2026  
**Licencia:** MIT  
