# Plasmid consensus pipeline (geNomad + PlasX + PLASMe)

This document is a GitHub-ready, copy/paste pipeline that:
1) builds a **consensus plasmid call** across tools (IDs + per-sample overlaps),  
2) generates **report tables** for downstream plotting, and  
3) produces **publication-ready plots** in R.

> Generated on 2026-01-10.

---

## Inputs (your current paths)

- geNomad outputs:  
  `/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/03_genomad`
- PlasX outputs:  
  `/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/plasx_predict`
- PLASMe reports directory:  
  `/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/01_files_plasme_used`
- Output consensus dir:  
  `/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/out_consensus`
- Output visualization dir:  
  `/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/visualizations`

---

## Step 0 — Minimal environment notes

### Python
- Requires Python 3.
- Both `plasmid_consensus_scores_ids.py` and `plasmid_consensus_visualizer.py` are **stdlib-only** (no third-party Python packages).

### R
The plotting script installs required packages automatically (if missing):
`ggplot2, readr, tidyr, dplyr, VennDiagram, gridExtra, RColorBrewer`.

---

## Step 1 — Compute consensus and export union IDs per sample

Command (exactly as you provided):

```bash
python3 plasmid_consensus_scores_ids.py   --genomad_dir /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/03_genomad   --plasx_dir /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/plasx_predict   --plasme_reports_dir /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/01_files_plasme_used   --outdir /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/out_consensus   --thr_genomad 0.7   --thr_plasx 0.5   --thr_plasme 0   --export_union_ids
```

### Key outputs

- **Global summary (all samples):**  
  `.../out_consensus/plasmid_consensus.global_summary.tsv`

- **Per-sample folders:**  
  `.../out_consensus/<sample>/`  
  plus (because you used `--export_union_ids`) a union list at:  
  `.../out_consensus/<sample>.plasmids_union.txt`

### Global summary columns

The global TSV is written with the following header:

- `sample`
- `thr_genomad, thr_plasx, thr_plasme`
- `genomad_n_ids, plasx_n_ids, plasme_n_ids`
- `G_int_X, G_int_M, X_int_M, G_int_X_int_M`
- `Gpass, Xpass, Mpass`
- `Gpass_int_Xpass, Gpass_int_Mpass, Xpass_int_Mpass, Gpass_int_Xpass_int_Mpass`
- `union_all_tools, only_genomad, only_plasx, only_plasme`
- `genomad_file_used, plasx_file, plasme_report, plasme_delim`

---

## Step 2 — Generate visualization-ready report tables (Python, stdlib-only)

Command (exactly as you provided):

```bash
python3 plasmid_consensus_visualizer.py   --global_summary /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/out_consensus/plasmid_consensus.global_summary.tsv   --outdir /media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/visualizations
```

This writes multiple TSVs to your visualization directory, including:
- `statistics_summary.tsv`
- `per_sample_detailed_report.tsv`
- `tool_comparison.tsv`
- `overlap_matrix.tsv`
- `venn_diagram_data_for_r.tsv`
- `README.txt` (instructions for plotting)

---

## Step 3 — Plot directly from the global summary (R)

You can plot **directly** from the global summary with:

```bash
Rscript plot_plasmid_consensus_direct.R
```

Important:
- `plot_plasmid_consensus_direct.R` uses **hard-coded paths** in a “CONFIGURACIÓN” section:
  - `GLOBAL_SUMMARY_FILE` must point to your `plasmid_consensus.global_summary.tsv`
  - `OUTPUT_DIR` must point to your visualization output directory

With your current setup, those paths already match your commands.

### Plot outputs (written under `OUTPUT_DIR`)

- `01_tool_comparison_barplot.png`
- `02_union_composition_stacked.png`
- `03_detection_heatmap.png`
- `04_detection_distribution_boxplot.png`
- `06_tool_specific_percentages.png`
- `venn_diagrams/<sample>_venn.png` (first 6 samples by default)
- `statistics_summary.tsv`

---

## Repository layout suggestion

A simple layout that stays readable:

```
plasmid-consensus/
  README.md
  scripts/
    plasmid_consensus_scores_ids.py
    plasmid_consensus_visualizer.py
    plot_plasmid_consensus_direct.R
```

---

# Scripts (full source)

## scripts/plasmid_consensus_scores_ids.py

```python
#!/usr/bin/env python3
# plasmid_consensus_scores.py - Enhanced version with union statistics

import argparse
import csv
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Set, Tuple

KEY_RE = re.compile(r"(k\d+_\d+)")

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def diagnose_id_mismatches(G_ids: Set[str], X_ids: Set[str], M_ids: Set[str], 
                          G_raw: Dict[str, str], X_raw: Dict[str, str], M_raw: Dict[str, str],
                          sample: str) -> None:
    """
    Diagnóstico de IDs no pareados entre herramientas.
    G_raw, X_raw, M_raw contienen {extracted_key: original_id}
    """
    eprint(f"\n[ID DIAGNOSTIC] [{sample}]")
    
    # IDs que están en una herramienta pero no en las otras
    only_in_G = G_ids - X_ids - M_ids
    only_in_X = X_ids - G_ids - M_ids
    only_in_M = M_ids - G_ids - X_ids
    
    if only_in_G:
        eprint(f"  IDs solo en Genomad (primeros 5): {sorted(only_in_G)[:5]}")
        if G_raw:
            eprint(f"    Originales: {[G_raw.get(k, 'NA') for k in sorted(only_in_G)[:5]]}")
    
    if only_in_X:
        eprint(f"  IDs solo en PlasX (primeros 5): {sorted(only_in_X)[:5]}")
        if X_raw:
            eprint(f"    Originales: {[X_raw.get(k, 'NA') for k in sorted(only_in_X)[:5]]}")
    
    if only_in_M:
        eprint(f"  IDs solo en PLASMe (primeros 5): {sorted(only_in_M)[:5]}")
        if M_raw:
            eprint(f"    Originales: {[M_raw.get(k, 'NA') for k in sorted(only_in_M)[:5]]}")

def mkdirp(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def extract_key(s: str) -> Optional[str]:
    if not s:
        return None
    m = KEY_RE.search(s)
    return m.group(1) if m else None

def extract_key_with_sample_prefix(s: str, sample: str) -> Optional[str]:
    """
    Extrae el key removiendo el prefijo del sample si existe.
    Ej: 'BAC4A_k141_1000' -> 'k141_1000'
        'k141_1000' -> 'k141_1000'
    """
    if not s:
        return None
    
    # Remover prefijo del sample si existe (case-insensitive)
    s_lower = s.lower()
    sample_prefix = sample.lower() + "_"
    if s_lower.startswith(sample_prefix):
        s = s[len(sample_prefix):]
    
    # Extraer el patrón k\d+_\d+
    m = KEY_RE.search(s)
    return m.group(1) if m else None

def choose_best(cands: List[str]) -> Optional[str]:
    if not cands:
        return None
    cands = sorted(set(cands), key=lambda p: (p.count(os.sep), len(p), p))
    return cands[0]

def find_genomad_aggregated(genomad_dir: str, sample: str) -> Optional[str]:
    pattern = os.path.join(genomad_dir, "**", f"{sample}_final.contigs_aggregated_classification.tsv")
    cands = glob.glob(pattern, recursive=True)
    if not cands:
        pattern2 = os.path.join(genomad_dir, "**", f"{sample}_final.contigs_aggregated_classification*.tsv")
        cands = glob.glob(pattern2, recursive=True)
    return choose_best(cands)

def find_plasx_scores(plasx_dir: str, sample: str) -> Optional[str]:
    cands = glob.glob(os.path.join(plasx_dir, "*-scores.txt"))
    if not cands:
        return None
    s = sample.lower()
    hits = [p for p in cands if s in os.path.basename(p).lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        eprint(f"[WARN] [{sample}] PlasX: múltiples candidatos, usaré: {os.path.basename(choose_best(hits))}")
        return choose_best(hits)
    return None

def find_plasme_report(plasme_reports_dir: str, sample: str) -> Optional[str]:
    cands = glob.glob(os.path.join(plasme_reports_dir, "*.plasme.plasmids.fna_report.csv"))
    if not cands:
        return None
    s = sample.lower()
    hits = [p for p in cands if s in os.path.basename(p).lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        eprint(f"[WARN] [{sample}] PLASMe: múltiples candidatos, usaré: {os.path.basename(choose_best(hits))}")
        return choose_best(hits)
    return None

def infer_samples_from_genomad(genomad_dir: str) -> List[str]:
    cands = glob.glob(os.path.join(genomad_dir, "**", "*_final.contigs_aggregated_classification.tsv"), recursive=True)
    samples: Set[str] = set()
    for p in cands:
        base = os.path.basename(p)
        m = re.match(r"(.+?)_final\.contigs_aggregated_classification\.tsv$", base)
        if m:
            samples.add(m.group(1))
    return sorted(samples)

def write_tsv(path: str, header: List[str], rows: List[List[str]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

def write_id_list(path: str, ids: Set[str]) -> None:
    """Escribe una lista de IDs, uno por línea"""
    with open(path, "w", encoding="utf-8") as f:
        for id_ in sorted(ids):
            f.write(id_ + "\n")

# -------------------------
# LOADERS con columnas fijas
# -------------------------

def load_genomad_plasmid_score_col3(path: str, sample: str = "") -> Tuple[Dict[str, float], Dict[str, str]]:
    """
    geNomad aggregated TSV:
      col1 = seq_name
      col2 = chromosome_score
      col3 = plasmid_score  <- usamos esta
      col4 = virus_score
    
    Devuelve: (scores_dict, raw_ids_dict)
      - scores_dict: {extracted_key: score}
      - raw_ids_dict: {extracted_key: original_seq_name}
    """
    out: Dict[str, float] = {}
    raw_ids: Dict[str, str] = {}
    
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            # salta header si existe
            if i == 0 and ("seq_name" in line.lower() or "plasmid_score" in line.lower()):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            original_id = parts[0].strip()
            k = extract_key_with_sample_prefix(original_id, sample) or original_id
            if not k:
                continue
            try:
                out[k] = float(parts[2])
                raw_ids[k] = original_id
            except Exception:
                continue
    return out, raw_ids

def load_plasx_score_col2(path: str, sample: str = "") -> Tuple[Dict[str, float], Dict[str, str]]:
    """
    PlasX scores.txt (TSV):
      col1 = contig
      col2 = score
    
    Devuelve: (scores_dict, raw_ids_dict)
    """
    out: Dict[str, float] = {}
    raw_ids: Dict[str, str] = {}
    
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            # header típico: "contig\tscore"
            if i == 0 and ("contig" in line.lower() and "score" in line.lower()):
                continue
            # algunas versiones incluyen "..." como línea
            if line.startswith("..."):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                # fallback por whitespace
                parts = line.split()
            if len(parts) < 2:
                continue
            original_id = parts[0].strip()
            k = extract_key_with_sample_prefix(original_id, sample)
            if not k:
                continue
            try:
                out[k] = float(parts[1])
                raw_ids[k] = original_id
            except Exception:
                continue
    return out, raw_ids

def sniff_csv_delim(path: str) -> str:
    # PLASMe suele ser CSV con comas; igualmente soportamos ; o TAB
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        head = f.readline()
        if "\t" in head:
            return "\t"
        if ";" in head and "," not in head:
            return ";"
        return ","

def load_plasme_score_col6(path: str, sample: str = "") -> Tuple[Dict[str, float], str, Dict[str, str]]:
    """
    PLASMe report CSV:
      col1 = contig
      ...
      col6 = score  <- usamos esta
    Devuelve: (scores, delim_usado, raw_ids_dict)
    """
    out: Dict[str, float] = {}
    raw_ids: Dict[str, str] = {}
    delim = sniff_csv_delim(path)

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        for i, row in enumerate(reader):
            if not row:
                continue
            # header típico: contig,length,reference,order,evidence,score,amb_region
            if i == 0 and row[0].strip().lower() == "contig":
                continue
            if len(row) < 6:
                continue
            original_id = (row[0] or "").strip()
            k = extract_key_with_sample_prefix(original_id, sample)
            if not k:
                continue
            try:
                out[k] = float((row[5] or "").strip())
                raw_ids[k] = original_id
            except Exception:
                continue

    return out, delim, raw_ids

def filter_pass(scores: Dict[str, float], thr: float) -> Set[str]:
    return {k for k, v in scores.items() if v >= thr}

def main():
    ap = argparse.ArgumentParser(
        description="Consensus plasmid calls from geNomad(aggregated col3), PlasX(col2), PLASMe(col6). Threshold siempre por herramienta."
    )
    ap.add_argument("--genomad_dir", required=True)
    ap.add_argument("--plasx_dir", required=True)
    ap.add_argument("--plasme_reports_dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", default=None,
                    help="Comma-separated samples OR path to file (one per line). If omitted: infer from geNomad aggregated.")

    # THR obligatorios
    ap.add_argument("--thr_genomad", type=float, required=True)
    ap.add_argument("--thr_plasx", type=float, required=True)
    ap.add_argument("--thr_plasme", type=float, required=True)

    ap.add_argument("--preflight_only", action="store_true")
    ap.add_argument("--export_union_ids", action="store_true",
                    help="Export union of all plasmid IDs to <sample>/<sample>.plasmids_union.txt")
    ap.add_argument("--diagnose_ids", action="store_true",
                    help="Show detailed ID matching diagnostics for troubleshooting")
    args = ap.parse_args()

    mkdirp(args.outdir)

    if args.samples is None:
        samples = infer_samples_from_genomad(args.genomad_dir)
        if not samples:
            eprint("[ERROR] No pude inferir muestras desde geNomad aggregated. Usa --samples.")
            sys.exit(2)
    else:
        if os.path.exists(args.samples):
            with open(args.samples, "r", encoding="utf-8") as f:
                samples = [x.strip() for x in f if x.strip() and not x.startswith("#")]
        else:
            samples = [x.strip() for x in args.samples.split(",") if x.strip()]

    global_header = [
        "sample",
        "thr_genomad", "thr_plasx", "thr_plasme",
        "genomad_n_ids", "plasx_n_ids", "plasme_n_ids",
        "G_int_X", "G_int_M", "X_int_M", "G_int_X_int_M",
        "Gpass", "Xpass", "Mpass",
        "Gpass_int_Xpass", "Gpass_int_Mpass", "Xpass_int_Mpass", "Gpass_int_Xpass_int_Mpass",
        "union_all_tools", "only_genomad", "only_plasx", "only_plasme",
        "genomad_file_used", "plasx_file", "plasme_report",
        "plasme_delim"
    ]
    global_rows: List[List[str]] = []

    for sample in samples:
        out_s = os.path.join(args.outdir, sample)
        mkdirp(out_s)

        g_file = find_genomad_aggregated(args.genomad_dir, sample)
        x_file = find_plasx_scores(args.plasx_dir, sample)
        m_rep = find_plasme_report(args.plasme_reports_dir, sample)

        if g_file is None:
            eprint(f"[WARN] [{sample}] No encontré geNomad aggregated_classification.tsv")
        if x_file is None:
            eprint(f"[WARN] [{sample}] No encontré PlasX -scores.txt")
        if m_rep is None:
            eprint(f"[WARN] [{sample}] No encontré PLASMe report.csv")

        # Cargar con información de IDs originales
        if g_file:
            G, G_raw = load_genomad_plasmid_score_col3(g_file, sample)
        else:
            G, G_raw = {}, {}
            
        if x_file:
            X, X_raw = load_plasx_score_col2(x_file, sample)
        else:
            X, X_raw = {}, {}
            
        if m_rep:
            M, m_delim, M_raw = load_plasme_score_col6(m_rep, sample)
        else:
            M, m_delim, M_raw = {}, "NA", {}

        G_ids = set(G.keys())
        X_ids = set(X.keys())
        M_ids = set(M.keys())

        # Diagnóstico de IDs si está activado
        if args.diagnose_ids:
            diagnose_id_mismatches(G_ids, X_ids, M_ids, G_raw, X_raw, M_raw, sample)

        GX = G_ids & X_ids
        GM = G_ids & M_ids
        XM = X_ids & M_ids
        GXM = G_ids & X_ids & M_ids

        # PASS por herramienta (siempre threshold)
        Gpass = filter_pass(G, args.thr_genomad)
        Xpass = filter_pass(X, args.thr_plasx)
        Mpass = filter_pass(M, args.thr_plasme)

        GpXp = Gpass & Xpass
        GpMp = Gpass & Mpass
        XpMp = Xpass & Mpass
        GpXpMp = Gpass & Xpass & Mpass

        # NUEVA ESTADÍSTICA: UNIÓN Y EXCLUSIVOS
        all_union = Gpass | Xpass | Mpass
        only_G = Gpass - Xpass - Mpass
        only_X = Xpass - Gpass - Mpass
        only_M = Mpass - Gpass - Xpass

        preflight_path = os.path.join(out_s, f"{sample}.id_preflight.tsv")
        pre_rows = [
            ["sample", sample],
            ["thr_genomad", str(args.thr_genomad)],
            ["thr_plasx", str(args.thr_plasx)],
            ["thr_plasme", str(args.thr_plasme)],
            ["genomad_file_used", g_file or "NA"],
            ["plasx_file", x_file or "NA"],
            ["plasme_report", m_rep or "NA"],
            ["plasme_delim", m_delim],
            ["genomad_n_ids", str(len(G_ids))],
            ["plasx_n_ids", str(len(X_ids))],
            ["plasme_n_ids", str(len(M_ids))],
            ["G_int_X", str(len(GX))],
            ["G_int_M", str(len(GM))],
            ["X_int_M", str(len(XM))],
            ["G_int_X_int_M", str(len(GXM))],
            ["Gpass", str(len(Gpass))],
            ["Xpass", str(len(Xpass))],
            ["Mpass", str(len(Mpass))],
            ["Gpass_int_Xpass", str(len(GpXp))],
            ["Gpass_int_Mpass", str(len(GpMp))],
            ["Xpass_int_Mpass", str(len(XpMp))],
            ["Gpass_int_Xpass_int_Mpass", str(len(GpXpMp))],
            ["union_all_tools", str(len(all_union))],
            ["only_genomad", str(len(only_G))],
            ["only_plasx", str(len(only_X))],
            ["only_plasme", str(len(only_M))],
        ]
        write_tsv(preflight_path, ["field", "value"], pre_rows)

        eprint(
            f"[PREFLIGHT] [{sample}] thr(G/X/M)=({args.thr_genomad},{args.thr_plasx},{args.thr_plasme}) | "
            f"G={len(G_ids)} X={len(X_ids)} M={len(M_ids)} | "
            f"Gpass={len(Gpass)} Xpass={len(Xpass)} Mpass={len(Mpass)} | "
            f"UNION={len(all_union)} | Only_G={len(only_G)} Only_X={len(only_X)} Only_M={len(only_M)} | "
            f"Gpass∩Xpass={len(GpXp)} Gpass∩Mpass={len(GpMp)} Xpass∩Mpass={len(XpMp)} | "
            f"Gpass∩Xpass∩Mpass={len(GpXpMp)} | report={preflight_path}"
        )

        global_rows.append([
            sample,
            str(args.thr_genomad), str(args.thr_plasx), str(args.thr_plasme),
            str(len(G_ids)), str(len(X_ids)), str(len(M_ids)),
            str(len(GX)), str(len(GM)), str(len(XM)), str(len(GXM)),
            str(len(Gpass)), str(len(Xpass)), str(len(Mpass)),
            str(len(GpXp)), str(len(GpMp)), str(len(XpMp)), str(len(GpXpMp)),
            str(len(all_union)), str(len(only_G)), str(len(only_X)), str(len(only_M)),
            os.path.basename(g_file) if g_file else "NA",
            os.path.basename(x_file) if x_file else "NA",
            os.path.basename(m_rep) if m_rep else "NA",
            m_delim
        ])

        # Exportar IDs de la unión si se solicita
        if args.export_union_ids and not args.preflight_only:
            union_ids_path = os.path.join(out_s, f"{sample}.plasmids_union.txt")
            write_id_list(union_ids_path, all_union)
            eprint(f"[EXPORT] [{sample}] Union IDs exported to: {union_ids_path}")

        if args.preflight_only:
            continue

        # Si después quieres el consensus.tsv (por contig), se agrega aquí.

    global_path = os.path.join(args.outdir, "plasmid_consensus.global_summary.tsv")
    write_tsv(global_path, global_header, global_rows)
    eprint(f"[OK] Global summary: {global_path}")
    eprint(f"[OK] Outputs por muestra en: {args.outdir}/<sample>/")

if __name__ == "__main__":
    main()

```

## scripts/plasmid_consensus_visualizer.py

```python
#!/usr/bin/env python3
# plasmid_consensus_report_simple.py
# Genera reportes detallados sin dependencias externas (solo stdlib)

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Dict

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def load_global_summary(tsv_path: str) -> List[Dict[str, str]]:
    """Carga el archivo global_summary.tsv"""
    data = []
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            data.append(row)
    return data

def write_tsv(path: str, header: List[str], rows: List[List[str]]) -> None:
    """Escribe un archivo TSV"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(row) + '\n')

def create_statistics_summary(data: List[Dict[str, str]], output_path: str):
    """Crea tabla de estadísticas resumidas"""
    n = len(data)
    
    # Calcular estadísticas
    union_vals = [int(row['union_all_tools']) for row in data]
    gpass_vals = [int(row['Gpass']) for row in data]
    xpass_vals = [int(row['Xpass']) for row in data]
    mpass_vals = [int(row['Mpass']) for row in data]
    gx_vals = [int(row['Gpass_int_Xpass']) for row in data]
    gm_vals = [int(row['Gpass_int_Mpass']) for row in data]
    xm_vals = [int(row['Xpass_int_Mpass']) for row in data]
    gxm_vals = [int(row['Gpass_int_Xpass_int_Mpass']) for row in data]
    only_g_vals = [int(row['only_genomad']) for row in data]
    only_x_vals = [int(row['only_plasx']) for row in data]
    only_m_vals = [int(row['only_plasme']) for row in data]
    
    # Porcentajes
    pct_only_g = [int(row['only_genomad']) / int(row['union_all_tools']) * 100 
                  if int(row['union_all_tools']) > 0 else 0 for row in data]
    pct_only_x = [int(row['only_plasx']) / int(row['union_all_tools']) * 100 
                  if int(row['union_all_tools']) > 0 else 0 for row in data]
    pct_only_m = [int(row['only_plasme']) / int(row['union_all_tools']) * 100 
                  if int(row['union_all_tools']) > 0 else 0 for row in data]
    
    samples_with_triple = sum(1 for val in gxm_vals if val > 0)
    
    stats = [
        ['Metric', 'Value'],
        ['Total samples analyzed', str(n)],
        ['Mean UNION per sample', f"{sum(union_vals)/n:.1f}"],
        ['Median UNION per sample', f"{sorted(union_vals)[n//2]:.1f}"],
        ['Total UNION across all samples', str(sum(union_vals))],
        ['Mean geNomad detection', f"{sum(gpass_vals)/n:.1f}"],
        ['Mean PlasX detection', f"{sum(xpass_vals)/n:.1f}"],
        ['Mean PLASMe detection', f"{sum(mpass_vals)/n:.1f}"],
        ['Mean G∩X overlap', f"{sum(gx_vals)/n:.1f}"],
        ['Mean G∩M overlap', f"{sum(gm_vals)/n:.1f}"],
        ['Mean X∩M overlap', f"{sum(xm_vals)/n:.1f}"],
        ['Mean G∩X∩M overlap', f"{sum(gxm_vals)/n:.1f}"],
        ['Samples with triple overlap > 0', f"{samples_with_triple} ({samples_with_triple/n*100:.1f}%)"],
        ['Mean % only geNomad', f"{sum(pct_only_g)/n:.1f}%"],
        ['Mean % only PlasX', f"{sum(pct_only_x)/n:.1f}%"],
        ['Mean % only PLASMe', f"{sum(pct_only_m)/n:.1f}%"],
    ]
    
    write_tsv(output_path, ['Metric', 'Value'], stats[1:])
    eprint(f"[OK] Statistics summary: {output_path}")
    
    # Imprimir a stderr
    eprint("\n" + "="*70)
    eprint("SUMMARY STATISTICS")
    eprint("="*70)
    for row in stats[1:]:
        eprint(f"{row[0]:.<50} {row[1]:>18}")
    eprint("="*70 + "\n")

def create_per_sample_report(data: List[Dict[str, str]], output_path: str):
    """Crea reporte detallado por muestra"""
    header = [
        'sample', 'union_all_tools', 
        'Gpass', 'Xpass', 'Mpass',
        'only_genomad', 'only_plasx', 'only_plasme',
        'Gpass_int_Xpass', 'Gpass_int_Mpass', 'Xpass_int_Mpass',
        'Gpass_int_Xpass_int_Mpass',
        'pct_only_G', 'pct_only_X', 'pct_only_M'
    ]
    
    rows = []
    for row in data:
        union = int(row['union_all_tools'])
        pct_g = int(row['only_genomad']) / union * 100 if union > 0 else 0
        pct_x = int(row['only_plasx']) / union * 100 if union > 0 else 0
        pct_m = int(row['only_plasme']) / union * 100 if union > 0 else 0
        
        rows.append([
            row['sample'], row['union_all_tools'],
            row['Gpass'], row['Xpass'], row['Mpass'],
            row['only_genomad'], row['only_plasx'], row['only_plasme'],
            row['Gpass_int_Xpass'], row['Gpass_int_Mpass'], row['Xpass_int_Mpass'],
            row['Gpass_int_Xpass_int_Mpass'],
            f"{pct_g:.1f}", f"{pct_x:.1f}", f"{pct_m:.1f}"
        ])
    
    write_tsv(output_path, header, rows)
    eprint(f"[OK] Per-sample report: {output_path}")

def create_venn_data_for_r(data: List[Dict[str, str]], output_path: str):
    """Crea archivo de datos para hacer Venn diagrams en R"""
    header = [
        'sample', 'only_G', 'only_X', 'only_M',
        'GX_not_M', 'GM_not_X', 'XM_not_G', 'GXM'
    ]
    
    rows = []
    for row in data:
        # Calcular valores para Venn
        only_g = int(row['only_genomad'])
        only_x = int(row['only_plasx'])
        only_m = int(row['only_plasme'])
        gxm = int(row['Gpass_int_Xpass_int_Mpass'])
        
        # Intersecciones de pares sin la triple
        gx_not_m = int(row['Gpass_int_Xpass']) - gxm
        gm_not_x = int(row['Gpass_int_Mpass']) - gxm
        xm_not_g = int(row['Xpass_int_Mpass']) - gxm
        
        rows.append([
            row['sample'], str(only_g), str(only_x), str(only_m),
            str(gx_not_m), str(gm_not_x), str(xm_not_g), str(gxm)
        ])
    
    write_tsv(output_path, header, rows)
    eprint(f"[OK] Venn diagram data for R: {output_path}")

def create_tool_comparison(data: List[Dict[str, str]], output_path: str):
    """Crea tabla comparativa de herramientas"""
    header = ['sample', 'geNomad', 'PlasX', 'PLASMe', 'UNION']
    
    rows = []
    for row in data:
        rows.append([
            row['sample'],
            row['Gpass'],
            row['Xpass'],
            row['Mpass'],
            row['union_all_tools']
        ])
    
    write_tsv(output_path, header, rows)
    eprint(f"[OK] Tool comparison: {output_path}")

def create_overlap_matrix(data: List[Dict[str, str]], output_path: str):
    """Crea matriz de overlaps"""
    header = ['sample', 'G∩X', 'G∩M', 'X∩M', 'G∩X∩M']
    
    rows = []
    for row in data:
        rows.append([
            row['sample'],
            row['Gpass_int_Xpass'],
            row['Gpass_int_Mpass'],
            row['Xpass_int_Mpass'],
            row['Gpass_int_Xpass_int_Mpass']
        ])
    
    write_tsv(output_path, header, rows)
    eprint(f"[OK] Overlap matrix: {output_path}")

def create_readme(output_dir: Path):
    """Crea un README explicando los archivos"""
    readme_path = output_dir / "README.txt"
    
    content = """
PLASMID CONSENSUS ANALYSIS REPORT
==================================

This directory contains analysis reports and data files for plasmid detection
across three tools: geNomad, PlasX, and PLASMe.

FILES:
------

1. statistics_summary.tsv
   - Overall statistics across all samples
   - Mean/median values for detection and overlaps
   
2. per_sample_detailed_report.tsv
   - Detailed metrics for each sample
   - Includes percentages of tool-specific detections
   - Ready for plotting in R/Excel
   
3. tool_comparison.tsv
   - Simple comparison of detections per tool
   - Format: sample, geNomad, PlasX, PLASMe, UNION
   - Good for bar plots
   
4. overlap_matrix.tsv
   - Intersection counts between tools
   - Format: sample, G∩X, G∩M, X∩M, G∩X∩M
   
5. venn_diagram_data_for_r.tsv
   - Data formatted for creating Venn diagrams in R
   - Columns: only_G, only_X, only_M, GX_not_M, GM_not_X, XM_not_G, GXM
   
EXAMPLE R CODE FOR VENN DIAGRAMS:
----------------------------------

library(VennDiagram)
library(readr)

data <- read_tsv("venn_diagram_data_for_r.tsv")

# For a specific sample (e.g., first row)
row <- data[1,]

venn.plot <- draw.triple.venn(
  area1 = row$only_G + row$GX_not_M + row$GM_not_X + row$GXM,
  area2 = row$only_X + row$GX_not_M + row$XM_not_G + row$GXM,
  area3 = row$only_M + row$GM_not_X + row$XM_not_G + row$GXM,
  n12 = row$GX_not_M + row$GXM,
  n23 = row$XM_not_G + row$GXM,
  n13 = row$GM_not_X + row$GXM,
  n123 = row$GXM,
  category = c("geNomad", "PlasX", "PLASMe"),
  fill = c("red", "blue", "green"),
  alpha = 0.3
)

EXAMPLE R CODE FOR BAR PLOTS:
------------------------------

library(ggplot2)
library(readr)
library(tidyr)

data <- read_tsv("tool_comparison.tsv")
data_long <- pivot_longer(data, cols = c(geNomad, PlasX, PLASMe, UNION),
                          names_to = "Tool", values_to = "Count")

ggplot(data_long, aes(x = sample, y = Count, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

"""
    
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    eprint(f"[OK] README created: {readme_path}")

def main():
    ap = argparse.ArgumentParser(
        description="Generate reports from plasmid consensus analysis (no external dependencies)"
    )
    ap.add_argument("--global_summary", required=True,
                    help="Path to plasmid_consensus.global_summary.tsv")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for reports")
    
    args = ap.parse_args()
    
    # Crear directorio de salida
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Cargar datos
    eprint(f"[INFO] Loading global summary from: {args.global_summary}")
    data = load_global_summary(args.global_summary)
    eprint(f"[INFO] Loaded {len(data)} samples")
    
    # Crear reportes
    eprint("\n[INFO] Creating reports...")
    create_statistics_summary(data, str(outdir / "statistics_summary.tsv"))
    create_per_sample_report(data, str(outdir / "per_sample_detailed_report.tsv"))
    create_tool_comparison(data, str(outdir / "tool_comparison.tsv"))
    create_overlap_matrix(data, str(outdir / "overlap_matrix.tsv"))
    create_venn_data_for_r(data, str(outdir / "venn_diagram_data_for_r.tsv"))
    create_readme(outdir)
    
    eprint(f"\n[OK] All reports saved to: {outdir}")
    eprint(f"[OK] Read README.txt for instructions on plotting in R")

if __name__ == "__main__":
    main()

```

## scripts/plot_plasmid_consensus_direct.R

```r
#!/usr/bin/env Rscript
# plot_plasmid_consensus_direct.R
# Script para visualizar resultados del análisis de consenso de plasmidos
# Lee DIRECTAMENTE el archivo plasmid_consensus.global_summary.tsv

# Instalar paquetes si no están instalados
packages <- c("ggplot2", "readr", "tidyr", "dplyr", "VennDiagram", "gridExtra", "RColorBrewer")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

# Cargar librerías
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(gridExtra)
library(RColorBrewer)

# =====================================================================
# CONFIGURACIÓN - CAMBIA ESTOS PATHS
# =====================================================================
GLOBAL_SUMMARY_FILE <- "/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/out_consensus/plasmid_consensus.global_summary.tsv"
OUTPUT_DIR <- "/media/mmorenos/disk3/012_eduCastroLab/02_plasmidos/04_plasmid_consensus_scores/visualizations"

# Crear directorio de salida
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# CARGAR DATOS
# =====================================================================
cat("Cargando datos desde:", GLOBAL_SUMMARY_FILE, "\n")

if(!file.exists(GLOBAL_SUMMARY_FILE)) {
  stop("ERROR: El archivo no existe: ", GLOBAL_SUMMARY_FILE)
}

df <- read_tsv(GLOBAL_SUMMARY_FILE, show_col_types = FALSE)
cat(sprintf("  ✓ %d muestras cargadas\n", nrow(df)))

# =====================================================================
# PREPARAR DATOS DERIVADOS
# =====================================================================

# Calcular porcentajes
df <- df %>%
  mutate(
    pct_only_G = ifelse(union_all_tools > 0, only_genomad / union_all_tools * 100, 0),
    pct_only_X = ifelse(union_all_tools > 0, only_plasx / union_all_tools * 100, 0),
    pct_only_M = ifelse(union_all_tools > 0, only_plasme / union_all_tools * 100, 0)
  )

# Preparar datos para Venn
venn_data <- df %>%
  mutate(
    only_G = only_genomad,
    only_X = only_plasx,
    only_M = only_plasme,
    GX_not_M = Gpass_int_Xpass - Gpass_int_Xpass_int_Mpass,
    GM_not_X = Gpass_int_Mpass - Gpass_int_Xpass_int_Mpass,
    XM_not_G = Xpass_int_Mpass - Gpass_int_Xpass_int_Mpass,
    GXM = Gpass_int_Xpass_int_Mpass
  )

# =====================================================================
# 1. GRÁFICO DE BARRAS: COMPARACIÓN DE HERRAMIENTAS
# =====================================================================
cat("\n1. Creando gráfico de comparación de herramientas...\n")

# Preparar datos en formato largo
tool_comp_long <- df %>%
  select(sample, Gpass, Xpass, Mpass, union_all_tools) %>%
  rename(geNomad = Gpass, PlasX = Xpass, PLASMe = Mpass, UNION = union_all_tools) %>%
  pivot_longer(cols = c(geNomad, PlasX, PLASMe, UNION),
               names_to = "Tool", values_to = "Count")

# Ordenar factores
tool_comp_long$Tool <- factor(tool_comp_long$Tool, 
                              levels = c("geNomad", "PlasX", "PLASMe", "UNION"))

# Limitar a primeras 20 muestras si hay muchas
n_samples <- min(20, length(unique(df$sample)))
samples_to_plot <- unique(df$sample)[1:n_samples]
tool_comp_subset <- tool_comp_long %>% filter(sample %in% samples_to_plot)

p1 <- ggplot(tool_comp_subset, aes(x = sample, y = Count, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("geNomad" = "#FF9999", 
                               "PlasX" = "#99CCFF", 
                               "PLASMe" = "#99FF99",
                               "UNION" = "#FFCC00")) +
  labs(title = "Plasmid Detection by Tool",
       subtitle = sprintf("First %d samples (Thresholds: G=%.1f, X=%.1f, M=%.1f)", 
                         n_samples, df$thr_genomad[1], df$thr_plasx[1], df$thr_plasme[1]),
       x = "Sample", y = "Number of Plasmids") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "01_tool_comparison_barplot.png"), 
       p1, width = 12, height = 6, dpi = 300)
cat("  ✓ Guardado: 01_tool_comparison_barplot.png\n")

# =====================================================================
# 2. GRÁFICO DE BARRAS APILADAS: COMPOSICIÓN DEL UNION
# =====================================================================
cat("\n2. Creando gráfico de composición del UNION...\n")

union_comp_long <- df %>%
  select(sample, only_genomad, only_plasx, only_plasme, 
         Gpass_int_Xpass, Gpass_int_Mpass, Xpass_int_Mpass,
         Gpass_int_Xpass_int_Mpass) %>%
  pivot_longer(cols = -sample, names_to = "Category", values_to = "Count")

# Renombrar categorías
union_comp_long$Category <- factor(union_comp_long$Category,
  levels = c("only_genomad", "only_plasx", "only_plasme",
             "Gpass_int_Xpass", "Gpass_int_Mpass", "Xpass_int_Mpass",
             "Gpass_int_Xpass_int_Mpass"),
  labels = c("Only geNomad", "Only PlasX", "Only PLASMe",
             "G ∩ X", "G ∩ M", "X ∩ M", "G ∩ X ∩ M"))

# Subset de muestras
union_comp_subset <- union_comp_long %>% filter(sample %in% samples_to_plot)

p2 <- ggplot(union_comp_subset, aes(x = sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "UNION Composition: Tool-Specific vs Shared Plasmids",
       subtitle = sprintf("First %d samples", n_samples),
       x = "Sample", y = "Number of Plasmids") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "02_union_composition_stacked.png"), 
       p2, width = 14, height = 6, dpi = 300)
cat("  ✓ Guardado: 02_union_composition_stacked.png\n")

# =====================================================================
# 3. HEATMAP: DETECCIÓN POR MUESTRA
# =====================================================================
cat("\n3. Creando heatmap de detección...\n")

heatmap_data <- df %>%
  select(sample, Gpass, Xpass, Mpass, union_all_tools,
         Gpass_int_Xpass, Gpass_int_Mpass, Xpass_int_Mpass,
         Gpass_int_Xpass_int_Mpass) %>%
  pivot_longer(cols = -sample, names_to = "Metric", values_to = "Count")

# Renombrar métricas
heatmap_data$Metric <- factor(heatmap_data$Metric,
  levels = c("Gpass", "Xpass", "Mpass", "union_all_tools",
             "Gpass_int_Xpass", "Gpass_int_Mpass", "Xpass_int_Mpass",
             "Gpass_int_Xpass_int_Mpass"),
  labels = c("geNomad", "PlasX", "PLASMe", "UNION",
             "G ∩ X", "G ∩ M", "X ∩ M", "G ∩ X ∩ M"))

p3 <- ggplot(heatmap_data, aes(x = Metric, y = sample, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red", name = "Plasmids") +
  labs(title = "Plasmid Detection Heatmap",
       x = "Detection Category", y = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "03_detection_heatmap.png"), 
       p3, width = 10, height = max(8, nrow(df) * 0.2), dpi = 300)
cat("  ✓ Guardado: 03_detection_heatmap.png\n")

# =====================================================================
# 4. BOXPLOTS: DISTRIBUCIÓN DE DETECCIONES
# =====================================================================
cat("\n4. Creando boxplots de distribución...\n")

box_data <- df %>%
  select(Gpass, Xpass, Mpass, union_all_tools) %>%
  rename(geNomad = Gpass, PlasX = Xpass, PLASMe = Mpass, UNION = union_all_tools) %>%
  pivot_longer(cols = everything(), names_to = "Tool", values_to = "Count")

box_data$Tool <- factor(box_data$Tool, 
                        levels = c("geNomad", "PlasX", "PLASMe", "UNION"))

p4 <- ggplot(box_data, aes(x = Tool, y = Count, fill = Tool)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("geNomad" = "#FF9999", 
                               "PlasX" = "#99CCFF", 
                               "PLASMe" = "#99FF99",
                               "UNION" = "#FFCC00")) +
  labs(title = "Distribution of Plasmid Detections Across Samples",
       x = "Tool", y = "Number of Plasmids") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "04_detection_distribution_boxplot.png"), 
       p4, width = 8, height = 6, dpi = 300)
cat("  ✓ Guardado: 04_detection_distribution_boxplot.png\n")

# =====================================================================
# 5. DIAGRAMAS DE VENN (primeras 6 muestras)
# =====================================================================
cat("\n5. Creando diagramas de Venn...\n")

venn_dir <- file.path(OUTPUT_DIR, "venn_diagrams")
dir.create(venn_dir, showWarnings = FALSE)

# Función para crear Venn diagram
create_venn <- function(row_data, sample_name, output_file) {
  
  # Extraer valores
  only_G <- as.numeric(row_data$only_G)
  only_X <- as.numeric(row_data$only_X)
  only_M <- as.numeric(row_data$only_M)
  GX_not_M <- as.numeric(row_data$GX_not_M)
  GM_not_X <- as.numeric(row_data$GM_not_X)
  XM_not_G <- as.numeric(row_data$XM_not_G)
  GXM <- as.numeric(row_data$GXM)
  
  # Calcular áreas totales
  area1 <- only_G + GX_not_M + GM_not_X + GXM  # geNomad
  area2 <- only_X + GX_not_M + XM_not_G + GXM  # PlasX
  area3 <- only_M + GM_not_X + XM_not_G + GXM  # PLASMe
  
  n12 <- GX_not_M + GXM
  n23 <- XM_not_G + GXM
  n13 <- GM_not_X + GXM
  n123 <- GXM
  
  union_total <- as.numeric(row_data$union_all_tools)
  
  # Crear diagrama
  venn.plot <- draw.triple.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    n12 = n12,
    n23 = n23,
    n13 = n13,
    n123 = n123,
    category = c("geNomad", "PlasX", "PLASMe"),
    fill = c("#FF9999", "#99CCFF", "#99FF99"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
  )
  
  # Guardar
  png(output_file, width = 800, height = 800)
  grid.draw(venn.plot)
  title(main = sprintf("Sample: %s\nTotal UNION: %d plasmids", 
                       sample_name, union_total),
        cex.main = 1.5)
  dev.off()
}

# Crear Venns para primeras 6 muestras
n_venns <- min(6, nrow(venn_data))
for(i in 1:n_venns) {
  sample_name <- venn_data$sample[i]
  output_file <- file.path(venn_dir, sprintf("%s_venn.png", sample_name))
  
  tryCatch({
    create_venn(venn_data[i,], sample_name, output_file)
    cat(sprintf("  ✓ Venn guardado: %s_venn.png\n", sample_name))
  }, error = function(e) {
    cat(sprintf("  ✗ Error en %s: %s\n", sample_name, e$message))
  })
}

# =====================================================================
# 6. GRÁFICO DE PROPORCIÓN: % SOLO DETECTADO POR CADA HERRAMIENTA
# =====================================================================
cat("\n6. Creando gráfico de proporciones...\n")

prop_data <- df %>%
  select(sample, pct_only_G, pct_only_X, pct_only_M) %>%
  pivot_longer(cols = -sample, names_to = "Tool", values_to = "Percentage")

prop_data$Tool <- factor(prop_data$Tool,
  levels = c("pct_only_G", "pct_only_X", "pct_only_M"),
  labels = c("Only geNomad", "Only PlasX", "Only PLASMe"))

prop_subset <- prop_data %>% filter(sample %in% samples_to_plot)

p6 <- ggplot(prop_subset, aes(x = sample, y = Percentage, fill = Tool)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Only geNomad" = "#FF9999", 
                               "Only PlasX" = "#99CCFF", 
                               "Only PLASMe" = "#99FF99")) +
  labs(title = "Percentage of Tool-Specific Plasmids",
       subtitle = "Plasmids detected only by one tool (not shared)",
       x = "Sample", y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUTPUT_DIR, "06_tool_specific_percentages.png"), 
       p6, width = 12, height = 6, dpi = 300)
cat("  ✓ Guardado: 06_tool_specific_percentages.png\n")

# =====================================================================
# 7. ESTADÍSTICAS RESUMIDAS
# =====================================================================
cat("\n7. Calculando estadísticas resumidas...\n")

stats_summary <- data.frame(
  Metric = c(
    "Total samples analyzed",
    "Mean UNION per sample",
    "Median UNION per sample",
    "Total UNION across all samples",
    "Mean geNomad detection",
    "Mean PlasX detection",
    "Mean PLASMe detection",
    "Mean G∩X overlap",
    "Mean G∩M overlap",
    "Mean X∩M overlap",
    "Mean G∩X∩M overlap",
    "Samples with triple overlap > 0",
    "Mean % only geNomad",
    "Mean % only PlasX",
    "Mean % only PLASMe"
  ),
  Value = c(
    nrow(df),
    round(mean(df$union_all_tools), 1),
    median(df$union_all_tools),
    sum(df$union_all_tools),
    round(mean(df$Gpass), 1),
    round(mean(df$Xpass), 1),
    round(mean(df$Mpass), 1),
    round(mean(df$Gpass_int_Xpass), 1),
    round(mean(df$Gpass_int_Mpass), 1),
    round(mean(df$Xpass_int_Mpass), 1),
    round(mean(df$Gpass_int_Xpass_int_Mpass), 1),
    sum(df$Gpass_int_Xpass_int_Mpass > 0),
    round(mean(df$pct_only_G), 1),
    round(mean(df$pct_only_X), 1),
    round(mean(df$pct_only_M), 1)
  )
)

write_tsv(stats_summary, file.path(OUTPUT_DIR, "statistics_summary.tsv"))
cat("  ✓ Guardado: statistics_summary.tsv\n")

# =====================================================================
# RESUMEN FINAL
# =====================================================================
cat("\n" %+% strrep("=", 70) %+% "\n")
cat("GRAFICACIÓN COMPLETADA\n")
cat(strrep("=", 70) %+% "\n")
cat(sprintf("Directorio de salida: %s\n", OUTPUT_DIR))
cat("\nArchivos creados:\n")
cat("  1. 01_tool_comparison_barplot.png\n")
cat("  2. 02_union_composition_stacked.png\n")
cat("  3. 03_detection_heatmap.png\n")
cat("  4. 04_detection_distribution_boxplot.png\n")
cat(sprintf("  5. venn_diagrams/ (%d diagramas)\n", n_venns))
cat("  6. 06_tool_specific_percentages.png\n")
cat("  7. statistics_summary.tsv\n")
cat(strrep("=", 70) %+% "\n")

# Imprimir estadísticas en consola
cat("\nESTADÍSTICAS RESUMIDAS:\n")
cat(strrep("=", 70) %+% "\n")
for(i in 1:nrow(stats_summary)) {
  cat(sprintf("%-45s %15s\n", stats_summary$Metric[i], stats_summary$Value[i]))
}
cat(strrep("=", 70) %+% "\n\n")

```
