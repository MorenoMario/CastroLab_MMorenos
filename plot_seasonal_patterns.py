#!/usr/bin/env python3
"""
=============================================================
Visualización de Patrones Estacionales - Plásmidos en MAGs
=============================================================
Lee el archivo summary_metrics.tsv generado por plasmid_mag_analysis.py
y genera gráficos de patrones estacionales.

Uso:
  python plot_seasonal_patterns.py [--input plasmid_mag_analysis/summary_metrics.tsv]
                                   [--out figures/]
=============================================================
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # sin GUI
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.gridspec import GridSpec


# Paleta estacional (Hemisferio Sur)
SEASON_COLORS = {
    'Summer':      '#e74c3c',   # rojo
    'Autumn':      '#e67e22',   # naranja
    'Winter':      '#3498db',   # azul
    'Spring':      '#2ecc71',   # verde
    'Lab_Control': '#95a5a6',   # gris
    'L-series':    '#9b59b6',   # violeta
    'Unknown':     '#bdc3c7',   # gris claro
}

SEASON_ORDER = ['Summer', 'Autumn', 'Winter', 'Spring', 'L-series', 'Lab_Control', 'Unknown']


def load_data(summary_file):
    df = pd.read_csv(summary_file, sep='\t')

    # Filtrar muestras con datos y MAGs correspondientes
    df = df[df['status'] == 'ok'].copy()
    df = df[df['total_plasmids'] > 0].copy()

    # Ordenar por sort_key
    df = df.sort_values('sort_key').reset_index(drop=True)

    # Color por estación
    df['color'] = df['season'].map(SEASON_COLORS).fillna('#bdc3c7')

    return df


def plot_overview(df, out_dir):
    """Panel principal: 4 subplots con métricas estacionales."""
    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, :])    # ancho completo
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, :])    # ancho completo

    x = np.arange(len(df))
    colors = df['color'].tolist()

    # ----------------------------------------------------------
    # 1. Total de plásmidos por muestra (barras)
    # ----------------------------------------------------------
    bars1 = ax1.bar(x, df['total_plasmids'], color=colors, edgecolor='white', linewidth=0.5)
    ax1.set_title('Total de contigs plasmídicos por muestra', fontsize=13, fontweight='bold')
    ax1.set_ylabel('Número de plásmidos', fontsize=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels(df['sample_id'], rotation=90, fontsize=6)
    ax1.tick_params(axis='x', pad=1)
    ax1.grid(axis='y', alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ----------------------------------------------------------
    # 2. Plásmidos en MAGs vs No en MAGs (barras apiladas)
    # ----------------------------------------------------------
    ax2.bar(x, df['plasmids_in_mag'],     color='#2ecc71', label='En MAG',    alpha=0.85)
    ax2.bar(x, df['plasmids_not_in_mag'], bottom=df['plasmids_in_mag'],
            color='#e74c3c', label='No en MAG', alpha=0.85)
    ax2.set_title('Plásmidos: en MAG vs no en MAG', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Número de plásmidos', fontsize=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(df['sample_id'], rotation=90, fontsize=6)
    ax2.legend(fontsize=9)
    ax2.grid(axis='y', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ----------------------------------------------------------
    # 3. % de plásmidos en MAGs por muestra
    # ----------------------------------------------------------
    ax3.bar(x, df['pct_in_mag'], color=colors, edgecolor='white', linewidth=0.5)
    ax3.axhline(df['pct_in_mag'].mean(), color='black', linestyle='--',
                linewidth=1.2, label=f"Media: {df['pct_in_mag'].mean():.1f}%")
    ax3.set_title('% plásmidos asociados a MAGs', fontsize=12, fontweight='bold')
    ax3.set_ylabel('% en MAG', fontsize=10)
    ax3.set_ylim(0, 105)
    ax3.set_xticks(x)
    ax3.set_xticklabels(df['sample_id'], rotation=90, fontsize=6)
    ax3.legend(fontsize=9)
    ax3.grid(axis='y', alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # ----------------------------------------------------------
    # 4. Boxplot por estación: % en MAG
    # ----------------------------------------------------------
    season_df = df[df['season'].isin(['Summer', 'Autumn', 'Winter', 'Spring'])]
    seasons_present = [s for s in SEASON_ORDER if s in season_df['season'].unique()]

    if len(seasons_present) > 0:
        data_by_season = [season_df[season_df['season'] == s]['pct_in_mag'].values
                          for s in seasons_present]
        bp = ax4.boxplot(data_by_season,
                         patch_artist=True,
                         medianprops=dict(color='black', linewidth=2),
                         whiskerprops=dict(linewidth=1.2),
                         capprops=dict(linewidth=1.2),
                         flierprops=dict(marker='o', markersize=5, alpha=0.6))
        for patch, season in zip(bp['boxes'], seasons_present):
            patch.set_facecolor(SEASON_COLORS.get(season, '#bdc3c7'))
            patch.set_alpha(0.75)

        # Añadir puntos individuales (jitter)
        for i, (season, patch) in enumerate(zip(seasons_present, bp['boxes']), 1):
            vals = season_df[season_df['season'] == season]['pct_in_mag'].values
            jitter = np.random.uniform(-0.2, 0.2, size=len(vals))
            ax4.scatter([i + j for j in jitter], vals,
                        color=SEASON_COLORS.get(season, '#bdc3c7'),
                        edgecolors='grey', linewidths=0.5,
                        s=40, alpha=0.7, zorder=3)

        ax4.set_xticklabels(seasons_present, fontsize=11)
        ax4.set_title('% plásmidos en MAGs por estación (Hemisferio Sur)',
                      fontsize=12, fontweight='bold')
        ax4.set_ylabel('% en MAG', fontsize=10)
        ax4.set_ylim(-5, 110)
        ax4.grid(axis='y', alpha=0.3)
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
    else:
        ax4.text(0.5, 0.5, 'Datos estacionales insuficientes',
                 ha='center', va='center', fontsize=12, transform=ax4.transAxes)

    # Leyenda estacional
    legend_patches = [mpatches.Patch(color=SEASON_COLORS[s], label=s)
                      for s in SEASON_ORDER if s in df['season'].unique()]
    fig.legend(handles=legend_patches, title='Estación', loc='lower center',
               ncol=len(legend_patches), fontsize=9, bbox_to_anchor=(0.5, -0.01),
               frameon=True, fancybox=True)

    fig.suptitle('Análisis de Asociación Plásmido-MAG: Patrones Estacionales',
                 fontsize=15, fontweight='bold', y=1.01)

    out_path = os.path.join(out_dir, 'seasonal_overview.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Guardado: {out_path}")


def plot_seasonal_timeseries(df, out_dir):
    """Serie de tiempo: plásmidos totales y % en MAG."""
    # Solo muestras temporales (no lab, no L-series)
    temp_df = df[df['season'].isin(['Summer', 'Autumn', 'Winter', 'Spring'])].copy()

    if len(temp_df) < 2:
        print("  ⚠ Insuficientes muestras temporales para serie de tiempo.")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    fig.subplots_adjust(hspace=0.15)

    x = np.arange(len(temp_df))
    colors = temp_df['color'].tolist()

    # Total plásmidos
    ax1.fill_between(x, temp_df['total_plasmids'], alpha=0.15, color='steelblue')
    ax1.plot(x, temp_df['total_plasmids'], 'o-', color='steelblue',
             linewidth=1.5, markersize=7, zorder=5)
    # Colorear puntos por estación
    for i, (xi, yi, c) in enumerate(zip(x, temp_df['total_plasmids'], colors)):
        ax1.scatter(xi, yi, color=c, s=80, zorder=6, edgecolors='white', linewidth=0.8)

    ax1.set_ylabel('Total plásmidos', fontsize=11)
    ax1.set_title('Serie Temporal: Plásmidos Totales y % en MAGs', fontsize=13, fontweight='bold')
    ax1.grid(alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # % en MAGs
    ax2.fill_between(x, temp_df['pct_in_mag'], alpha=0.15, color='green')
    ax2.plot(x, temp_df['pct_in_mag'], 's-', color='green',
             linewidth=1.5, markersize=7, zorder=5)
    for i, (xi, yi, c) in enumerate(zip(x, temp_df['pct_in_mag'], colors)):
        ax2.scatter(xi, yi, color=c, s=80, zorder=6, edgecolors='white', linewidth=0.8)

    ax2.axhline(temp_df['pct_in_mag'].mean(), color='darkgreen', linestyle='--',
                linewidth=1.2, alpha=0.8, label=f"Media: {temp_df['pct_in_mag'].mean():.1f}%")
    ax2.set_ylabel('% en MAG', fontsize=11)
    ax2.set_ylim(0, 105)
    ax2.legend(fontsize=10)
    ax2.grid(alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax2.set_xticks(x)
    ax2.set_xticklabels(temp_df['sample_id'], rotation=90, fontsize=7)

    # Sombrear por estación
    current_season = None
    start_x = 0
    for i, row in enumerate(temp_df.itertuples()):
        if row.season != current_season:
            if current_season is not None:
                for ax in [ax1, ax2]:
                    ax.axvspan(start_x - 0.5, i - 0.5,
                               alpha=0.06,
                               color=SEASON_COLORS.get(current_season, '#bdc3c7'))
            current_season = row.season
            start_x = i
    # Último segmento
    if current_season:
        for ax in [ax1, ax2]:
            ax.axvspan(start_x - 0.5, len(temp_df) - 0.5,
                       alpha=0.06,
                       color=SEASON_COLORS.get(current_season, '#bdc3c7'))

    # Leyenda estacional
    legend_patches = [mpatches.Patch(color=SEASON_COLORS[s], label=s, alpha=0.7)
                      for s in SEASON_ORDER if s in temp_df['season'].unique()]
    fig.legend(handles=legend_patches, title='Estación', loc='upper right',
               fontsize=9, frameon=True, fancybox=True)

    out_path = os.path.join(out_dir, 'timeseries_seasonal.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Guardado: {out_path}")


def plot_season_stats(df, out_dir):
    """Estadísticas agregadas por estación."""
    season_df = df[df['season'].isin(['Summer', 'Autumn', 'Winter', 'Spring'])].copy()
    if len(season_df) == 0:
        return

    stats = season_df.groupby('season').agg(
        n_samples=('sample_id', 'count'),
        mean_total=('total_plasmids', 'mean'),
        std_total=('total_plasmids', 'std'),
        mean_in_mag=('plasmids_in_mag', 'mean'),
        mean_pct=('pct_in_mag', 'mean'),
        std_pct=('pct_in_mag', 'std'),
        sum_total=('total_plasmids', 'sum'),
        sum_in_mag=('plasmids_in_mag', 'sum'),
    ).reindex([s for s in SEASON_ORDER if s in season_df['season'].unique()])
    stats['pct_global'] = 100 * stats['sum_in_mag'] / stats['sum_total']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    seasons = stats.index.tolist()
    x = np.arange(len(seasons))
    season_colors = [SEASON_COLORS.get(s, '#bdc3c7') for s in seasons]

    # Plásmidos totales por estación (media ± SD)
    axes[0].bar(x, stats['mean_total'], color=season_colors, alpha=0.8,
                edgecolor='grey', linewidth=0.5)
    axes[0].errorbar(x, stats['mean_total'], yerr=stats['std_total'].fillna(0),
                     fmt='none', color='black', capsize=5, linewidth=1.5)
    axes[0].set_title('Plásmidos totales\n(media ± SD)', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Nº plásmidos', fontsize=10)
    axes[0].set_xticks(x); axes[0].set_xticklabels(seasons, fontsize=10)
    axes[0].grid(axis='y', alpha=0.3)
    axes[0].spines['top'].set_visible(False); axes[0].spines['right'].set_visible(False)

    # % en MAG por estación (media ± SD)
    axes[1].bar(x, stats['mean_pct'], color=season_colors, alpha=0.8,
                edgecolor='grey', linewidth=0.5)
    axes[1].errorbar(x, stats['mean_pct'], yerr=stats['std_pct'].fillna(0),
                     fmt='none', color='black', capsize=5, linewidth=1.5)
    axes[1].axhline(season_df['pct_in_mag'].mean(), color='black', linestyle='--',
                    linewidth=1.2, label='Media global')
    axes[1].set_title('% plásmidos en MAG\n(media ± SD)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('% en MAG', fontsize=10)
    axes[1].set_ylim(0, 110)
    axes[1].set_xticks(x); axes[1].set_xticklabels(seasons, fontsize=10)
    axes[1].legend(fontsize=9)
    axes[1].grid(axis='y', alpha=0.3)
    axes[1].spines['top'].set_visible(False); axes[1].spines['right'].set_visible(False)

    # % global de plásmidos en MAG (conteo real, no media de %)
    axes[2].bar(x, stats['pct_global'], color=season_colors, alpha=0.8,
                edgecolor='grey', linewidth=0.5)
    axes[2].set_title('% global plásmidos en MAG\n(suma total / suma en MAG)', fontsize=11, fontweight='bold')
    axes[2].set_ylabel('% global', fontsize=10)
    axes[2].set_ylim(0, 110)
    axes[2].set_xticks(x); axes[2].set_xticklabels(seasons, fontsize=10)
    axes[2].grid(axis='y', alpha=0.3)
    axes[2].spines['top'].set_visible(False); axes[2].spines['right'].set_visible(False)

    # Añadir n= encima de cada barra
    for ax, col in zip(axes, ['mean_total', 'mean_pct', 'pct_global']):
        for i, (s, v) in enumerate(zip(seasons, stats[col])):
            n = stats.loc[s, 'n_samples']
            ax.text(i, v + ax.get_ylim()[1] * 0.02,
                    f'n={n}', ha='center', va='bottom', fontsize=9, color='grey')

    fig.suptitle('Estadísticas por Estación', fontsize=14, fontweight='bold')
    out_path = os.path.join(out_dir, 'season_statistics.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Guardado: {out_path}")


def main():
    parser = argparse.ArgumentParser(description='Visualización estacional plásmido-MAG')
    parser.add_argument('--input', default='plasmid_mag_analysis/summary_metrics.tsv',
                        help='Archivo summary_metrics.tsv')
    parser.add_argument('--out', default='plasmid_mag_analysis/figures',
                        help='Directorio de figuras (default: plasmid_mag_analysis/figures/)')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: No se encuentra {args.input}")
        print("  Primero ejecuta: python plasmid_mag_analysis.py")
        sys.exit(1)

    os.makedirs(args.out, exist_ok=True)

    print(f"\nCargando datos desde: {args.input}")
    df = load_data(args.input)
    print(f"  Muestras con datos: {len(df)}")
    print(f"  Estaciones encontradas: {df['season'].unique().tolist()}\n")

    print("Generando figuras:")
    plot_overview(df, args.out)
    plot_seasonal_timeseries(df, args.out)
    plot_season_stats(df, args.out)

    print(f"\nFiguras guardadas en: {args.out}/")
    print("  seasonal_overview.png   - Panel completo de métricas")
    print("  timeseries_seasonal.png - Serie temporal")
    print("  season_statistics.png   - Estadísticas por estación")


if __name__ == '__main__':
    main()
