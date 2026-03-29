#!/usr/bin/env python3
"""
=============================================================
Plasmid-MAG Association Analysis
=============================================================
Determina qué contigs plasmídicos están insertados en los
Metagenome-Assembled Genomes (MAGs) de cada muestra.

Lógica de matching:
- Plásmidos:  >k119_100265 evidence=weak_evidence methods=genomad
- MAGs:       >EC2017L1S05_k119_100265 total_depth_19.96 ...
  El prefijo de muestra se elimina extrayendo el patrón k[kmer]_[num].

Salidas (TSV + FASTA con IDs originales):
  plasmid_mag_analysis/
    ├── summary_metrics.tsv           <- métricas por muestra
    ├── plasmids_in_MAGs.tsv          <- tabla detallada (con MAG bin)
    ├── plasmids_NOT_in_MAGs.tsv      <- tabla detallada
    ├── plasmids_in_MAGs.fa           <- secuencias FASTA (header original + anotación)
    ├── plasmids_NOT_in_MAGs.fa       <- secuencias FASTA (header original)
    └── samples_without_mags.txt      <- muestras sin MAG correspondiente

Uso:
  python plasmid_mag_analysis.py [--plasmids lista_plasmids2.txt]
                                 [--mags lista_mags.txt]
                                 [--out output_dir]
                                 [--no-fasta]        # solo TSV, sin FASTA
=============================================================
"""

import os
import re
import sys
import csv
import argparse
from collections import defaultdict


# ============================================================
# EXTRACCIÓN DE SAMPLE IDs
# ============================================================

def extract_sample_id_from_plasmid_path(filepath):
    """
    Extrae sample_id del path de un archivo de plásmido.

      FP.DEC2017_L1S1D05.plasmids.fa                -> FP.DEC2017_L1S1D05
      CL_FP.BAC4A_ATATCTCG-ACTAAGAT_L00M.plasmids.fa -> CL_FP.BAC4A
      FP.EC-BacJ.plasmids.fa                          -> FP.EC-BacJ
      FP.L151-D5.plasmids.fa                          -> FP.L151-D5
    """
    fname = os.path.basename(filepath)
    name  = fname.replace('.plasmids.fa', '')
    # Remover barcodes Illumina: _[ACGT]+-[ACGT]+_L[0-9]+[letra]
    name  = re.sub(r'_[ACGT]+-[ACGT]+_L\d+\w*$', '', name)
    return name


def extract_sample_id_from_mag_path(filepath):
    """
    Extrae sample_id del path de un archivo MAG.

      bin.bin.9.FP.DEC2017_L1S1D05.fa.FP.DEC2017_L1S1D05.fa -> FP.DEC2017_L1S1D05
      bin.bin.10.CL_FP.BAC4G.fa.CL_FP.BAC4G.fa               -> CL_FP.BAC4G
      bin.10.FP.JUN2017_L4S1D05.fa                            -> FP.JUN2017_L4S1D05
      bin.13.FP.EC-BacJ.fa                                    -> FP.EC-BacJ
    """
    fname = os.path.basename(filepath)
    # Patrón 1: bin.bin.NUM.SAMPLEID.fa.SAMPLEID.fa
    m = re.match(r'^bin\.bin\.\d+\.(.*?)\.fa\.\1\.fa$', fname)
    if m:
        return m.group(1)
    # Patrón 2: bin.bin.NUM.SAMPLEID.fa
    m = re.match(r'^bin\.bin\.\d+\.(.*?)\.fa$', fname)
    if m:
        return m.group(1)
    # Patrón 3: bin.NUM.SAMPLEID.fa
    m = re.match(r'^bin\.\d+\.(.*?)\.fa$', fname)
    if m:
        return m.group(1)
    return None


# ============================================================
# LECTURA FASTA
# ============================================================

def extract_contig_base(contig_id):
    """
    Extrae el ID base (k[kmer]_[num]) eliminando prefijo de muestra.

      k119_100265              -> k119_100265
      EC2017L1S05_k119_100265  -> k119_100265
      BAC4G_k119_55321         -> k119_55321
    """
    m = re.search(r'(k\d+_\d+)$', contig_id)
    return m.group(1) if m else contig_id


def read_fasta_records(filepath):
    """
    Lee registros FASTA completos (header + secuencia).
    Retorna lista de dicts:
      { 'full_id':  str,   # ID sin '>'  (primera palabra del header)
        'base_id':  str,   # ID sin prefijo de muestra
        'header':   str,   # línea de header completa sin '>'
        'sequence': str }  # secuencia concatenada
    """
    records = []
    cur_header = None
    cur_full_id = None
    cur_base_id = None
    cur_seq = []

    def _flush():
        if cur_header is not None:
            records.append({
                'full_id':  cur_full_id,
                'base_id':  cur_base_id,
                'header':   cur_header,
                'sequence': ''.join(cur_seq)
            })

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    _flush()
                    cur_header  = line[1:]                     # sin '>'
                    cur_full_id = cur_header.split()[0]
                    cur_base_id = extract_contig_base(cur_full_id)
                    cur_seq     = []
                elif cur_header is not None:
                    cur_seq.append(line.strip())
        _flush()
    except FileNotFoundError:
        print(f"  ⚠ ARCHIVO NO ENCONTRADO: {filepath}", file=sys.stderr)
    except PermissionError:
        print(f"  ⚠ PERMISO DENEGADO:      {filepath}", file=sys.stderr)

    return records


def read_fasta_headers_only(filepath):
    """
    Lee solo headers FASTA (rápido, sin cargar secuencias).
    Retorna lista de tuplas: (full_id, base_id)
    Usado para indexar MAGs sin cargar sus secuencias en memoria.
    """
    headers = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    full_id = line[1:].split()[0]
                    base_id = extract_contig_base(full_id)
                    headers.append((full_id, base_id))
    except (FileNotFoundError, PermissionError) as e:
        print(f"  ⚠ {e}", file=sys.stderr)
    return headers


def write_fasta_record(fh, header, sequence, annotation='', linewidth=80):
    """
    Escribe un registro FASTA respetando el header original.
    Si se proporciona annotation, se añade al final del header.
    """
    if annotation:
        fh.write(f'>{header} {annotation}\n')
    else:
        fh.write(f'>{header}\n')
    for i in range(0, len(sequence), linewidth):
        fh.write(sequence[i:i + linewidth] + '\n')


# ============================================================
# ANOTACIÓN ESTACIONAL (Hemisferio Sur)
# ============================================================

MONTH_INFO = {
    'JAN': ('January',   '01', 2,  'Summer'),
    'FEB': ('February',  '02', 3,  'Summer'),
    'MAR': ('March',     '03', 4,  'Autumn'),
    'APR': ('April',     '04', 5,  'Autumn'),
    'MAY': ('May',       '05', 6,  'Autumn'),
    'JUN': ('June',      '06', 7,  'Winter'),
    'JUL': ('July',      '07', 8,  'Winter'),
    'AUG': ('August',    '08', 9,  'Winter'),
    'SEP': ('September', '09', 10, 'Spring'),
    'OCT': ('October',   '10', 11, 'Spring'),
    'NOV': ('November',  '11', 12, 'Spring'),
    'DEC': ('December',  '12', 1,  'Summer'),
}

def get_temporal_info(sample_id):
    for abbr, (month_name, month_num, sort_month, season) in MONTH_INFO.items():
        if abbr in sample_id.upper():
            ym = re.search(rf'{abbr}(\d{{4}})', sample_id.upper())
            year = ym.group(1) if ym else 'unknown'
            date_label = f"{abbr}{year}" if year != 'unknown' else abbr
            sort_key   = int(year + month_num) if year != 'unknown' else 999999
            return dict(month_name=month_name, month_num=month_num, year=year,
                        season=season, date_label=date_label, sort_key=sort_key)

    if 'EC-BAC' in sample_id.upper() or 'BAC4' in sample_id.upper():
        return dict(month_name='Lab_Control', month_num='LC', year='lab',
                    season='Lab_Control', date_label='Lab', sort_key=1000000)

    lm = re.search(r'L(\d)(\d)(\d).*?D(\d+)', sample_id)
    if lm:
        return dict(month_name=f'L-series_{sample_id}', month_num='LS', year='L-series',
                    season='L-series',
                    date_label=f'L{lm.group(1)}{lm.group(2)}{lm.group(3)}',
                    sort_key=900000 + int(lm.group(1)) * 10000 + int(lm.group(2)) * 100)

    return dict(month_name='Unknown', month_num='UN', year='unknown',
                season='Unknown', date_label='Unknown', sort_key=999998)


# ============================================================
# ANÁLISIS PRINCIPAL
# ============================================================

def resolve_file_list(path, extensions=('.fa', '.fasta', '.fa.gz')):
    """
    Acepta tanto un archivo de texto con rutas como un directorio.
    - Si es un archivo: lee las rutas línea a línea.
    - Si es un directorio: escanea recursivamente buscando archivos
      con las extensiones indicadas.
    Retorna lista de rutas absolutas existentes.
    """
    path = os.path.abspath(path)

    if os.path.isfile(path):
        with open(path) as f:
            paths = [l.strip() for l in f if l.strip()]
        print(f"  Fuente: archivo de lista  ({path})")
        return paths

    if os.path.isdir(path):
        found = []
        for root, dirs, files in os.walk(path):
            dirs.sort()   # orden reproducible
            for fname in sorted(files):
                if any(fname.endswith(ext) for ext in extensions):
                    found.append(os.path.join(root, fname))
        print(f"  Fuente: directorio        ({path})")
        return found

    raise FileNotFoundError(
        f"'{path}' no es ni un archivo de lista ni un directorio válido.")


def run_analysis(plasmid_list_file, mag_list_file, output_dir, write_fasta=True):
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 65)
    print("  ANÁLISIS DE ASOCIACIÓN PLÁSMIDO-MAG")
    print("=" * 65)

    print("\nResolviendo archivos de plásmidos...")
    plasmid_files = resolve_file_list(plasmid_list_file,
                                      extensions=('.plasmids.fa', '.fa', '.fasta'))
    # Si el directorio contiene .fa genéricos, filtrar solo los de plásmidos
    # (archivos cuyo nombre contiene 'plasmid')
    if os.path.isdir(plasmid_list_file):
        plasmid_files_filtered = [p for p in plasmid_files
                                  if 'plasmid' in os.path.basename(p).lower()]
        if plasmid_files_filtered:
            plasmid_files = plasmid_files_filtered

    print(f"Resolviendo archivos de MAGs...")
    mag_files = resolve_file_list(mag_list_file,
                                  extensions=('.fa', '.fasta'))
    # En directorio de MAGs excluir archivos de plásmidos si los hubiera
    if os.path.isdir(mag_list_file):
        mag_files = [p for p in mag_files
                     if 'plasmid' not in os.path.basename(p).lower()]

    print(f"\nArchivos de plásmidos: {len(plasmid_files)}")
    print(f"Archivos de MAGs:      {len(mag_files)}")

    # Índice sample_id -> [rutas MAG]
    sample_to_mags = defaultdict(list)
    for mpath in mag_files:
        sid = extract_sample_id_from_mag_path(mpath)
        if sid:
            sample_to_mags[sid].append(mpath)
        else:
            print(f"  ⚠ No se extrajo sample_id de MAG: {os.path.basename(mpath)}",
                  file=sys.stderr)

    print(f"Muestras únicas en MAGs: {len(sample_to_mags)}")

    # -------------------------------------------------------
    # Abrir FASTAs de salida globales (todas las muestras)
    # -------------------------------------------------------
    fasta_in_path  = os.path.join(output_dir, 'plasmids_in_MAGs.fa')
    fasta_not_path = os.path.join(output_dir, 'plasmids_NOT_in_MAGs.fa')

    if write_fasta:
        fh_in  = open(fasta_in_path,  'w')
        fh_not = open(fasta_not_path, 'w')
        fh_in.write(f'# plasmids_in_MAGs.fa — contigs plasmídicos asociados a MAGs\n')
        fh_not.write(f'# plasmids_NOT_in_MAGs.fa — contigs plasmídicos NO asociados a MAGs\n')

    # -------------------------------------------------------
    # Procesar cada archivo de plásmidos
    # -------------------------------------------------------
    results_in_mag     = []
    results_not_in_mag = []
    summary_rows       = []
    no_mag_match       = []

    print(f"\n{'Sample':<35} {'Total':>6} {'InMAG':>6} {'%InMAG':>7} {'MAGfiles':>8}")
    print("-" * 65)

    for ppath in plasmid_files:
        psample  = extract_sample_id_from_plasmid_path(ppath)
        temporal = get_temporal_info(psample)

        mag_list_for_sample = sample_to_mags.get(psample, [])
        if not mag_list_for_sample:
            no_mag_match.append(psample)

        # Leer registros FASTA completos del archivo de plásmidos
        plasmid_records = read_fasta_records(ppath)

        if not plasmid_records:
            summary_rows.append({
                'sample_id': psample, **temporal,
                'n_mag_files': len(mag_list_for_sample),
                'total_plasmids': 0, 'plasmids_in_mag': 0,
                'plasmids_not_in_mag': 0, 'pct_in_mag': 0.0,
                'pct_not_in_mag': 0.0, 'status': 'file_not_accessible'
            })
            continue

        # Construir lookup MAG: base_id -> (bin_filename, full_id_en_MAG)
        # Solo se leen headers de MAGs (sin secuencia) para eficiencia
        mag_contig_to_bin = {}
        for mpath in mag_list_for_sample:
            bin_name = os.path.basename(mpath)
            for full_id, base_id in read_fasta_headers_only(mpath):
                if base_id not in mag_contig_to_bin:
                    mag_contig_to_bin[base_id] = (bin_name, full_id)

        mag_base_ids = set(mag_contig_to_bin.keys())

        # Clasificar plásmidos
        n_in_mag = 0
        n_not    = 0

        for rec in plasmid_records:
            full_id  = rec['full_id']
            base_id  = rec['base_id']
            header   = rec['header']    # header original completo
            sequence = rec['sequence']  # secuencia original

            if base_id in mag_base_ids:
                bin_name, mag_full_id = mag_contig_to_bin[base_id]
                results_in_mag.append({
                    'sample_id':           psample,
                    'date_label':          temporal['date_label'],
                    'year':                temporal['year'],
                    'month_num':           temporal['month_num'],
                    'season':              temporal['season'],
                    'plasmid_contig_base': base_id,
                    'plasmid_contig_full': full_id,
                    'mag_bin':             bin_name,
                    'mag_contig_full':     mag_full_id,
                    'plasmid_file':        os.path.basename(ppath),
                })
                n_in_mag += 1

                if write_fasta:
                    # Header original + anotación del MAG al final
                    annotation = (f'[sample={psample}]'
                                  f'[MAG_bin={bin_name}]'
                                  f'[MAG_contig={mag_full_id}]')
                    write_fasta_record(fh_in, header, sequence, annotation)

            else:
                results_not_in_mag.append({
                    'sample_id':           psample,
                    'date_label':          temporal['date_label'],
                    'year':                temporal['year'],
                    'month_num':           temporal['month_num'],
                    'season':              temporal['season'],
                    'plasmid_contig_base': base_id,
                    'plasmid_contig_full': full_id,
                    'plasmid_file':        os.path.basename(ppath),
                })
                n_not += 1

                if write_fasta:
                    # Header original sin modificar
                    annotation = f'[sample={psample}]'
                    write_fasta_record(fh_not, header, sequence, annotation)

        total   = n_in_mag + n_not
        pct_in  = round(100 * n_in_mag / total, 2) if total > 0 else 0.0
        pct_not = round(100 * n_not / total, 2)    if total > 0 else 0.0

        summary_rows.append({
            'sample_id':           psample,
            **temporal,
            'n_mag_files':         len(mag_list_for_sample),
            'total_plasmids':      total,
            'plasmids_in_mag':     n_in_mag,
            'plasmids_not_in_mag': n_not,
            'pct_in_mag':          pct_in,
            'pct_not_in_mag':      pct_not,
            'status':              'ok' if mag_list_for_sample else 'no_mags'
        })

        marker = "" if mag_list_for_sample else " ⚠ no MAGs"
        print(f"  {psample:<33} {total:>6} {n_in_mag:>6} {pct_in:>6.1f}%"
              f" {len(mag_list_for_sample):>8}{marker}")

    # Cerrar FASTAs
    if write_fasta:
        fh_in.close()
        fh_not.close()

    # -------------------------------------------------------
    # Escribir TSVs
    # -------------------------------------------------------
    summary_rows.sort(key=lambda x: x['sort_key'])

    # 1. Métricas resumen
    summary_file = os.path.join(output_dir, 'summary_metrics.tsv')
    _write_tsv(summary_file, summary_rows,
               fields=['sample_id', 'date_label', 'year', 'month_num', 'month_name',
                       'season', 'sort_key', 'n_mag_files', 'total_plasmids',
                       'plasmids_in_mag', 'plasmids_not_in_mag',
                       'pct_in_mag', 'pct_not_in_mag', 'status'])

    # 2. Plásmidos en MAGs
    in_mag_tsv = os.path.join(output_dir, 'plasmids_in_MAGs.tsv')
    _write_tsv(in_mag_tsv, results_in_mag,
               empty_msg='# No se encontraron plásmidos en MAGs\n')

    # 3. Plásmidos NO en MAGs
    not_mag_tsv = os.path.join(output_dir, 'plasmids_NOT_in_MAGs.tsv')
    _write_tsv(not_mag_tsv, results_not_in_mag,
               empty_msg='# Todos los plásmidos están en MAGs\n')

    # 4. Muestras sin MAGs
    if no_mag_match:
        no_match_file = os.path.join(output_dir, 'samples_without_mags.txt')
        with open(no_match_file, 'w') as f:
            f.write("# Plásmidos sin MAG correspondiente\n")
            for s in no_mag_match:
                f.write(s + '\n')
        print(f"\n  ⚠ {len(no_mag_match)} muestras sin MAGs → {no_match_file}")

    # -------------------------------------------------------
    # Resumen global
    # -------------------------------------------------------
    ok_rows  = [r for r in summary_rows if r['status'] == 'ok']
    total_p  = sum(r['total_plasmids']      for r in ok_rows)
    total_in = sum(r['plasmids_in_mag']     for r in ok_rows)
    total_no = sum(r['plasmids_not_in_mag'] for r in ok_rows)
    pct_g    = round(100 * total_in / total_p, 2) if total_p > 0 else 0

    print("\n" + "=" * 65)
    print("  RESUMEN GLOBAL")
    print("=" * 65)
    print(f"  Muestras procesadas (con MAGs): {len(ok_rows)}")
    print(f"  Total contigs plasmídicos:      {total_p:,}")
    print(f"  Plásmidos en MAGs:              {total_in:,}  ({pct_g:.1f}%)")
    print(f"  Plásmidos NO en MAGs:           {total_no:,}  ({100 - pct_g:.1f}%)")
    print(f"\n  Archivos generados en: {output_dir}/")
    print(f"    summary_metrics.tsv")
    print(f"    plasmids_in_MAGs.tsv")
    print(f"    plasmids_NOT_in_MAGs.tsv")
    if write_fasta:
        print(f"    plasmids_in_MAGs.fa       ← FASTA con IDs originales + anotación MAG")
        print(f"    plasmids_NOT_in_MAGs.fa   ← FASTA con IDs originales")
    if no_mag_match:
        print(f"    samples_without_mags.txt")
    print("=" * 65)

    return summary_rows


# ============================================================
# HELPERS
# ============================================================

def _write_tsv(filepath, rows, fields=None, empty_msg='# Sin datos\n'):
    with open(filepath, 'w', newline='') as f:
        if rows:
            if fields is None:
                fields = list(rows[0].keys())
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                                    extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
        else:
            f.write(empty_msg)


# ============================================================
# ENTRY POINT
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description='Análisis Plásmido-MAG',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Con archivos de lista (comportamiento original)
  python plasmid_mag_analysis.py --plasmids lista_plasmids2.txt --mags lista_mags.txt

  # Con directorios directamente
  python plasmid_mag_analysis.py \\
      --plasmids /ruta/a/directorio_plasmidos/ \\
      --mags     /ruta/a/directorio_mags/ \\
      --out      ./resultados/

  # Sin FASTA (más rápido)
  python plasmid_mag_analysis.py --plasmids ... --mags ... --no-fasta
""")
    p.add_argument('--plasmids',
                   help='Archivo de lista de rutas O directorio con archivos .plasmids.fa')
    p.add_argument('--mags',
                   help='Archivo de lista de rutas O directorio con archivos MAG .fa')
    p.add_argument('--out',       default='plasmid_mag_analysis',
                   help='Directorio de salida (default: plasmid_mag_analysis/)')
    p.add_argument('--no-fasta',  action='store_true',
                   help='No generar archivos FASTA (solo TSV, más rápido)')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    run_analysis(args.plasmids, args.mags, args.out,
                 write_fasta=not args.no_fasta)
