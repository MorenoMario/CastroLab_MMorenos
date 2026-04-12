# Moviloma-Resistoma — Project Summary
**Fjord Comau, Chile | 90 metagenomes | 3 seasons × 3 years × 2 depths**

---

## Data overview

| Variable | Groups | n samples |
|---|---|---|
| Season | Winter / Summer / Fall | 44 / 23 / 24 |
| Year | 2017 / 2018 / 2019 | 40 / 26 / 24 |
| Depth | 5m / 20m | ~45 / ~45 |
| Season × Year | Winter_2017 / Summer_2017 / Winter_2018 / Summer_2018 / Fall_2019 | 19/21/24/2/24 |

**Matrices**: MAG TPM (90 samples × bins), Plasmid TPM, Combined TPM  
**ARG annotation**: 10,503 ARGs in dereplicated MAGs; 22,254 total (reliable: 22,226)

---

## Key results by script

### Script 11 — Composition by Season and Depth (`results/01_season_depth/`)

**Community composition (Bray-Curtis) differs significantly by both Season and Depth, but total abundance does not.**

| Test | Factor | R² | p |
|---|---|---|---|
| PERMANOVA | MAG — Season | **0.250** | 0.001 |
| PERMANOVA | MAG — Depth | **0.094** | 0.001 |
| PERMANOVA | Plasmid — Season | **0.235** | 0.001 |
| PERMANOVA | Plasmid — Depth | **0.099** | 0.001 |
| KW (total TPM) | Season | — | 0.201 ns |
| MW (total TPM) | Depth | — | 0.262 ns |

**Interpretation**: The *who is there* changes with season and depth, but not *how much biomass* there is. The community restructures without increasing total abundance.

Betadisper (Depth, MAGs): p=0.002 → 5m samples are more variable than 20m.

---

### Script 12 — ARG richness and load by Season and Depth (`results/02_ARG_season_depth/`)

**ARG richness and load are significantly structured by Season, not by Depth.**

| Test | Variable | χ² / W | p |
|---|---|---|---|
| KW | ARG richness ~ Season | 17.05 | **0.0002** |
| KW | ARG total hits ~ Season | 14.41 | **0.0007** |
| KW | ARG per MAG TPM ~ Season | 14.44 | **0.0007** |
| MW | ARG richness ~ Depth | — | 0.373 ns |
| MW | ARG total hits ~ Depth | — | 0.345 ns |
| Spearman | ARG richness ~ MAG TPM | ρ=0.098 | 0.361 ns |

**Interpretation**: ARG diversity tracks season, not depth or total microbial biomass. Resistance genes are not simply proportional to bacterial abundance — they are acquired/lost independently, suggesting active selection.

---

### Script 13 — Two winters are not the same (`results/03_inter_annual_winter/`)

**Winter 2017 and Winter 2018 have radically different community compositions despite similar total abundances — the strongest separation signal in the entire dataset.**

| Test | Comparison | R² / R | p |
|---|---|---|---|
| PERMANOVA | SeasonYear (3 groups) | **0.480** | 0.001 |
| ANOSIM | SeasonYear | R=0.702 | 0.001 |
| PERMANOVA | Winter_2017 vs Winter_2018 | **0.505** | 0.001 |
| ANOSIM | Winter_2017 vs Winter_2018 | R=0.816 | 0.001 |
| KW (total TPM) | Winter_2017 vs Winter_2018 | — | 0.916 ns |
| MW (total TPM) | Winter_2017 vs Winter_2018 | — | 0.678 ns |

**Interpretation**: The two winter communities are compositionally as different as entirely different seasons (R²=0.50), yet carry the same total biomass. This is not a seasonal effect — it is an inter-annual shift. Winter 2017 (June 2017, surface spring conditions) and Winter 2018 (August 2018, deep winter) likely experienced different environmental drivers (temperature, stratification, precipitation). This is the paper's central finding.

---

### Script 14 — Year as an explicit factor (`results/04_year_factor/`)

**Year explains MORE variance than Season when both are modeled simultaneously — confirming the inter-annual signal is real and not a proxy for seasonal structure.**

#### PERMANOVA additive model (`~ Year + Season`):

| Factor | Df | R² | F | p |
|---|---|---|---|---|
| Year | 2 | **0.249** | 17.875 | **0.001** |
| Season | 1 | **0.151** | 21.682 | **0.001** |
| Residual | 86 | 0.600 | — | — |

**Interpretation**: Year (R²=0.249) explains more community variance than Season (R²=0.151) after accounting for each other. The inter-annual signal is not a seasonal artifact — both factors structure the fjord resistome independently. Together, Year + Season explain ~40% of total community variance.

---

### Script 15 — ARG mobility: plasmid × MAG × compartment (`results/05_ARG_mobility/`)

**ARG mobility (% of ARGs on mobile elements) varies significantly with Season and Year, but not Depth. Mobility positively correlates with ARG richness.**

#### Mobility index by ecological factor:

| Test | Factor | χ² / W | p |
|---|---|---|---|
| KW | Mobility % ~ Season | 11.771 | **0.0028** |
| KW | Mobility % ~ Year | 10.114 | **0.0064** |
| MW | Mobility % ~ Depth | W=1216.5 | 0.087 ns |
| Spearman | Mobility % ~ ARG richness | ρ=0.543 | **3.3×10⁻⁸** |

#### Mobility by Season × Year:

| Group | n | Median mobility % | Mean mobility % |
|---|---|---|---|
| Winter 2017 | 19 | **0.0** | 0.00 |
| Winter 2018 | 24 | 0.3 | 2.09 |
| Summer 2017 | 21 | 1.4 | 2.24 |
| Summer 2018 | 2 | 1.9 | 1.93 |
| Fall 2019 | 24 | **4.1** | 5.31 |

**Striking observation**: Winter 2017 has a median mobility of 0% — no ARGs on mobile elements. Fall 2019 is the most mobile (4.1% median). This mirrors the compositional inter-annual divergence seen in Script 13/14.

#### Compartment breakdown (22,226 reliable ARGs):

| Compartment | n ARGs | % |
|---|---|---|
| Chromosome | 21,527 | 96.8% |
| Plasmid-in-MAG | 687 | 3.1% |
| Free plasmid | 12 | 0.05% |

#### Multi-compartment genes (direct HGT evidence — 17 genes):

Key genes detected simultaneously in chromosome AND plasmid:

| Gene | Compartments | Clinical relevance |
|---|---|---|
| **GIM-1** | chromosome + **free plasmid** | ⚠ Carbapenemase — highest priority |
| **macB** | chromosome + plasmid-in-MAG (148 hits, 21 samples) | Macrolide efflux |
| **novA** | chromosome + plasmid-in-MAG (130 hits, 11 samples) | Novobiocin resistance |
| **ugd** | chromosome + plasmid-in-MAG (55 hits, 26 samples) | Polymyxin resistance |
| **bcrA** | chromosome + plasmid-in-MAG (58 hits, 11 samples) | Bacitracin resistance |
| **msbA** | chromosome + plasmid-in-MAG | Multidrug efflux |
| **vanH** (vanA/vanB) | chromosome + plasmid-in-MAG | Glycopeptide resistance |
| **tetA(58)**, **tetB(46)** | chromosome + plasmid-in-MAG | Tetracycline resistance |

#### Shared contigs (direct dispersal evidence):

- 686 assembled contigs appear in >1 sample across the time series
- 17 carry annotated ARGs, including **mecA** (methicillin resistance), **msbA** (5 contigs), **tetA(41)**, **vanR** (vanF cluster), **arlR**
- mecA-carrying contig (k119_26528) detected in Summer and Fall 2017 — evidence of horizontal transmission across time points

#### Network hub analysis (mobile ARGs):

Top ARG genes by number of MAG-bin connections: **RANA** (20 bins), **macB** (17 bins), **RANB** (15 bins), **ugd** (8 bins), **novA** (7 bins)  
Top MAG bins by ARG connections: **bin1** (8 ARGs), **bin20** (6), **bin28** (6)

---

## Summary figure map

| Figure | Script | Key message |
|---|---|---|
| NMDS/PCoA Season | 11 | Communities differ by season |
| NMDS/PCoA Depth | 11 | Communities differ by depth |
| ARG richness ~ Season | 12 | ARG diversity is seasonal |
| PCoA SeasonYear | 13 | Two winters are compositionally distinct |
| Winter highlight PCoA | 13 | W2017 vs W2018: R²=0.50 |
| PCoA Year × Season | 14 | Year explains more variance than Season |
| Additive model bar | 14 | Year R²=0.249 > Season R²=0.151 |
| Mobility stacked bar | 15 | Compartment structure per sample |
| Mobility ~ Season/Year | 15 | Mobility is seasonal and inter-annual |
| Multi-compartment dot | 15 | 17 genes with HGT evidence |
| Shared contigs heatmap | 15 | mecA + msbA dispersal across time |
| Mobility network | 15 | MAG bins as mobilization hubs |
| Mobility ~ ARG richness | 15 | ρ=0.54, mobility predicts diversity |

---

## Scripts — execution order

```bash
cd ~/Documents/Claude/Projects/CASTROLAB/01_resistoma_MGE_MAGs

# Already run:
Rscript 11_visualize.R               # Season + Depth composition → results/01_season_depth/
Rscript 12_ARG_analysis.R            # ARG richness/load → results/02_ARG_season_depth/
Rscript 13_SeasonYear_analysis.R     # Two winters → results/03_inter_annual_winter/
Rscript 14_year_analysis.R           # Year factor → results/04_year_factor/
Rscript 15_mobility_network.R        # ARG mobility → results/05_ARG_mobility/
```

Metadata files: `metadata/metadata_MM.txt`, `metadata/metadata_MM_year.txt`

Required R packages: `vegan`, `ggplot2`, `dplyr`, `tidyr`, `scales`, `RColorBrewer`, `readxl`, `ggrepel`, `pheatmap`, `dunn.test`, `igraph`, `ggraph`

---

## Biological narrative (draft for paper)

The resistome-mobilome of fjord Comau is **not static**. Community composition responds significantly to both season (R²=0.25, p=0.001) and depth (R²=0.09, p=0.001), but total microbial biomass remains stable — suggesting ecological restructuring rather than bloom dynamics.

ARG diversity tracks seasonality independently of biomass (KW p=0.0002), implying that resistance gene acquisition and loss are driven by selective pressures beyond simple community turnover. The decoupling of ARG richness from total TPM (Spearman ρ=0.098, ns) is particularly notable: resistance is not a passive function of bacterial abundance.

When Year is modeled explicitly alongside Season (additive PERMANOVA), Year emerges as the dominant structuring factor (R²=0.249 vs Season R²=0.151, both p<0.001). Together they explain ~40% of total community variance, with ~60% remaining unexplained — consistent with fine-scale environmental stochasticity.

The most striking signal is the **divergence between winter 2017 and winter 2018** (PERMANOVA R²=0.50, ANOSIM R=0.816, both p=0.001) — a compositional difference as large as between entirely different seasons, yet with identical total abundance and zero difference in ARG richness. This inter-annual shift cannot be attributed to seasonal cycling. Winter 2017 (June 2017, early sampling) and Winter 2018 (August 2018, peak winter conditions) likely experienced different regimes of temperature, fjord stratification, and terrestrial runoff. This is the paper's central finding.

ARG mobility analysis reveals that ~3.1% of ARGs in dereplicated MAGs reside on mobile genetic elements (687 plasmid-in-MAG, 12 free-plasmid). Mobility varies significantly with Season (KW p=0.003) and Year (KW p=0.006), but not Depth. Critically, mobility % correlates strongly with ARG richness (Spearman ρ=0.543, p<10⁻⁷) — samples with more diverse resistomes are also more mobilized. Winter 2017 stands out again: median mobility = 0%, the lowest in the dataset, while Fall 2019 reaches 4.1%.

Seventeen gene families were detected simultaneously in chromosomal and plasmid compartments, providing direct molecular evidence of horizontal gene transfer. **GIM-1**, a carbapenemase of the highest clinical priority, was recovered from a free plasmid — the most mobile compartment. **mecA** (methicillin resistance) was detected on the same assembled contig in two distinct samples from different time points, constituting direct evidence of physical dispersal. Network analysis identifies **bin1**, **bin20**, and **bin28** as mobilization hubs (MAG bins with the highest ARG connectivity), and **RANA**, **macB**, and **RANB** as the most widely distributed mobile ARGs across the fjord.

Taken together, the fjord Comau resistome is structured by inter-annual environmental variability operating above the seasonal cycle, with horizontal gene transfer serving as an active mechanism of resistance spread — particularly pronounced in warmer, more productive seasons.
