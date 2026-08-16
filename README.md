## Introduction

This repository holds code and data associated to the manuscript:

[**Environmental Niches in Dual-mycorrhizal and Non-mycorrhizal Australian Forest Trees**]()

**Authors:**
Luke Florence<sup>1</sup>, John W. Morgan<sup>1</sup>, Peter A. Vesk<sup>2</sup>, Jen L. Wood<sup>3</sup>, Camille Truong<sup>4,5</sup>

**Affiliations:**
1. Department of Ecological, Plant and Animal Sciences, La Trobe University, Bundoora, Victoria, Australia. 
2. School of Agriculture, Food and Ecosystem Sciences, University of Melbourne, Parkville, Victoria, Australia. 
3. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, Victoria, Australia. 
4. Royal Botanic Gardens Victoria, Melbourne, Victoria, Australia.
5. School of BioSciences, University of Melbourne, Parkville, Victoria, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au) 

## Repository contents

* `./code/` — All R scripts required to reproduce the analysis: raw-data acquisition, data preparation, statistical/spatial modelling, and figure/table generation.
* `./generated_data/` — The processed datasets needed to reproduce every figure in the main text without acquiring or reprocessing raw data (see "Reproducing the analysis" below).
* `./envs/` — Conda environment set-up scripts for the species distribution modelling step, which is computationally intensive (approximately two days on 12 cores).
* `./covariates.txt` — A mapping of covariate names to their sources, descriptions and access routes.
* `./plant_data.txt` — A mapping of plant database names to their sources, descriptions and access routes.

### Raw data is not included

This repository and its Figshare deposit do not include the raw data (occurrence records, taxonomic checklists, or covariate rasters, including georeferenced predictors such as the national forest-cover raster and bioregion boundary layers). `data/` is excluded from version control and is not part of the Figshare deposit due to size and in some cases terms-of-use restrictions:

1. **Size.** The raw inputs total around 84 GB (versus around 130 MB for `generated_data/`). This exceeds what is practical to version-control or host as a supplementary deposit.
2. **Redistribution terms.** Some raw sources are openly licensed and could in principle be redistributed (e.g. the Harmonised Australian Vegetation Plot dataset and Forests of Australia (2023) are both CC BY 4.0). Others including the GBIF occurrence download and BGCI's GlobalTreeSearch are aggregated or terms-of-use-restricted third-party products that we have not confirmed we can redistribute in bulk. Rather than apply different handling per source, no raw data is provided. Licences and access routes for every covariate are listed in `covariates.txt` and plant database in `./plant_data.txt`.

## Reproducing the analysis

There are two ways to use this repository, depending on whether you need to regenerate the processed datasets or only primary figures.

### Path A — primnary figures from the generated data

Every script listed as `TRUE` in the "Generated data only" column of the Outputs table below reads exclusively from `generated_data/` (and, for context/colour scales, `covariates.txt`). Clone the repository and run the relevant script directly.

### Path B — the full pipeline from raw data

Reproducing `generated_data/` itself, or anything marked `FALSE` below, means acquiring the raw data first. Run the following in order; each stage's output feeds the next.

1. **`01.Harmonise_tree_data.R`** — Downloads and harmonises the World Flora Online backbone, FungalRoot, GlobalTreeSearch, HAVPlot, the Australian Plant Census (APC) and GBIF, then builds the tree-occurrence dataset and `generated_data/global_tree_mycorrhizal_types.txt`. Three sources have no stable direct-download URL and require a manual step (see each section's inline instructions):
   - The APC export.
   - GlobalTreeSearch's Australia-filtered species list (`global_tree_search_trees_1_9_australia.csv`) — the global list downloads automatically, but the Australia-only export has to be produced from BGCI's web tool. This file is the sole source of the `native_status` field used throughout the rest of the pipeline.
   - The GBIF occurrence download itself (`data/gbif/0071622-*.csv`) and its GBIF-backbone-matched species list (`data/gbif/normalized.csv`) — the script cites the download's DOI but does not fetch it automatically.
2. **`build_phylogenetic_tree.R`** — Builds `generated_data/phylo_tree_mycorrhizal_types.tre` from (1)'s output. Required by `05a`, `05b`, `07d` and `Figure_1.R`.
3. **`02a.Upscale_predictor_rasters.R`** — Downloads and upscales the forest-cover raster and climate/soil covariate rasters.
4. **`02b.Prepare_prediction_grids.R`** — Builds the 10 km prediction grid and site-level presence data. Produces Figs. S2–S4 and Table S1.
5. **`03a.Niche_cache_enmeval.R` → `03b.Niche_est_enmeval.R` → `03c.Niche_run_enmeval.R`** — Species distribution modelling with ENMeval. This is the computationally intensive step referenced above (~2 days on 12 cores; see `envs/`).
6. **`03d.Niche_est_summary.R`** — Summarises the modelling output into `generated_data/niche_estimates.txt`.
7. **`04a.Relative_richness_analysis.R`, `04b.Relative_richness_latitude.R`, `04c.Absolute_richness_analysis.R`, `04d.Absolute_richness_latitude.R`** — Richness analyses. Produce `generated_data/figure_2*`, `figure_3.RData` and `figure_S8*`, plus Figs. S5, S6, S9–S11 and Table S2.
8. **`05a.Species_niche_breadth_differences.R`, `05b.Species_niche_position_differences.R`** — Niche breadth and position analyses. Produce `generated_data/figure_4.RData` and `figure_5.RData`, plus Tables S3–S7.
9. **`06a.Sensitivity_relative_richness.R`** — Sensitivity analysis run directly on raw presence data (Figs. S14–S16).
10. **`07a.Generate_empirical_dataset.R` → `07b.Prepare_empirical_prediction_grids.R` → `07c.Relative_richness_empirical.R` → `07d.Niche_breadth_differences_empirical.R`** — The empirical pathway, using measured rather than genus-inferred mycorrhizal types. Produces Table S8 and Dataset S2.
11. **`08_Mid-domain_null_models.R`** — Fits the observed latitude-breadth models and runs the 1,000-simulation mid-domain-effect null model, in parallel across mycorrhizal types. Only needs `generated_data/niche_estimates.txt` (stage 6's output), so it can run any time after that. Produces `generated_data/figure_6.RData`; `Figure_6.R` reads that directly and does not run any modelling itself.

## Outputs
To reproduce the outputs of this project, see the following scripts. "Generated data only" indicates whether the script reads exclusively from `generated_data/` (`TRUE`), or additionally requires raw data acquired via Path B above (`FALSE`).

| Output | Script | Generated data only |
|--------|--------|:---:|
| Fig. 1 | [code/Figure_1.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_1.R) | TRUE |
| Fig. 2 | [code/Figure_2.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_2.R) | TRUE |
| Fig. 3 | [code/Figure_3.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_3.R) | TRUE |
| Fig. 4 | [code/Figure_4.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_4.R) | TRUE |
| Fig. 5 | [code/Figure_5.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_5.R) | TRUE |
| Fig. 6 | [code/Figure_6.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_6.R) | TRUE |
| Fig. S1 | [code/Figure_S1.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S1.R) | FALSE |
| Fig. S2 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) | FALSE |
| Fig. S3 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) | FALSE |
| Fig. S4 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) | FALSE |
| Fig. S5 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) | FALSE |
| Fig. S6 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) | FALSE |
| Fig. S7 | [code/Figure_S7.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S7.R) | TRUE |
| Fig. S8 | [code/Figure_S8.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S8.R) | TRUE |
| Fig. S9 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) | FALSE |
| Fig. S10 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) | FALSE |
| Fig. S11 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) | FALSE |
| Fig. S12 | [code/Figure_S12.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S12.R) | FALSE |
| Fig. S13 | [code/Figure_S13.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S13.R) | TRUE |
| Fig. S14 | [code/06a.Sensitivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensitivity_relative_richness.R) | FALSE |
| Fig. S15 | [code/06a.Sensitivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensitivity_relative_richness.R) | FALSE |
| Fig. S16 | [code/06a.Sensitivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensitivity_relative_richness.R) | FALSE |
| Fig. S17 | [code/Figure_S17.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S17.R) | TRUE |
| Fig. S18 | [code/Figure_S18.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S18.R) | TRUE |
| Fig. S19 | [code/06b.Sensitivity_niche_breadth.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06b.Sensitivity_niche_breadth.R) | TRUE |
| Fig. S20 | [code/Figure_S20.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S20.R) | FALSE |
| Table S1 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) | FALSE |
| Table S2 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) | FALSE |
| Table S3 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) | FALSE |
| Table S4 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) | FALSE |
| Table S5 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) | FALSE |
| Table S6 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) | FALSE |
| Table S7 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) | FALSE |
| Table S8 | [code/07a.Generate_empirical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_empirical_dataset.R) | FALSE |
| Dataset S1 | [code/01.Harmonise_tree_data.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/01.Harmonise_tree_data.R) | FALSE |
| Dataset S2 | [code/07a.Generate_empirical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_empirical_dataset.R) | FALSE |

## Data dictionary for generated_data/

Every file in `generated_data/`, what it is, and the script that produces it. Column lists are given for the `.txt` files; `.RData` files bundle several R objects, so their key objects and columns are summarised rather than exhaustively listed (several carry hundreds of `PEM1...PEMn` phylogenetic eigenvector map columns, one per predictor).

**Base taxonomic and phylogenetic data**

| File | Produced by | Contents |
|---|---|---|
| `global_tree_mycorrhizal_types.txt` | `01.Harmonise_tree_data.R` | One row per tree species (57,081 rows). Columns: `family`, `genus`, `scientific_name`, `native_status` (native/non-native, from GlobalTreeSearch's Australia-filtered export), `mycorrhizal_type` (from FungalRoot, genus-level). The source of truth for native status used throughout the pipeline. |
| `phylo_tree_mycorrhizal_types.tre` | `build_phylogenetic_tree.R` | Newick-format phylogenetic tree, 3,802 tips, one per native Australian tree species with a resolved mycorrhizal type. |
| `niche_estimates.txt` | `03d.Niche_est_summary.R` | One row per modelled species (2,334 rows). Columns: `species`, `family`, `genus`, `mycorrhizal_type`, `AOO`/`EOO` (area/extent of occurrence), `env_B2`/`env_B2_corrected` (environmental niche breadth, raw and range-corrected), `geo_B2`/`geo_B2_corrected` (geographic niche breadth), `ex_dent`, `RC1_position`/`RC2_position`/`RC3_position` (position on the three rotated climate/soil components), `lat_position`, `lat_range`, `biome`, `climate_zone`. |

**Figure 2**

| File | Produced by | Contents |
|---|---|---|
| `figure_2a.tif` | `04a.Relative_richness_analysis.R` | Raster, 4 layers (AM, EcM, EcM-AM, NM), 373 x 492 cells. Predicted relative richness surfaces. |
| `figure_2b.RData` | `04b.Relative_richness_latitude.R` | `latitude_gradient_data` (47,120 x 3: `latitude`, `mycorrhizal_type`, `relative_richness` -- the raw latitude-richness observations), `latitude_gradient_marginal_effects` (400 x 5: model-predicted `relative_richness` with `lower`/`upper` CIs across the latitude gradient), `latitude_gradient_coeficients` (4 x 4: the fitted slope annotation per mycorrhizal type). |
| `figure_2c.RData` | `04a.Relative_richness_analysis.R` | `data_figure_2c` (290 x 6: `RC1_bin`, `RC2_bin` and predicted richness per mycorrhizal type across the RC1-RC2 environmental grid), `data_figure_2c_limits` (1 x 6: the RC1/RC2 axis ranges used to build that grid). |

**Figure 3**

| File | Produced by | Contents |
|---|---|---|
| `figure_3.RData` | `04a.Relative_richness_analysis.R` | `relative_richness_data` (47,120 x 5: `RC1`, `RC2`, `RC3`, `mycorrhizal_type`, `response` -- observed relative richness against each rotated component), `relative_richness_marginal_effects` (1,200 x 6: fitted response curves with CIs), `relative_richness_coeficients` (12 x 5: fitted slope annotations). |

**Figure 4 and 5 (niche position and breadth)**

| File | Produced by | Contents |
|---|---|---|
| `figure_4.RData` | `05b.Species_niche_position_differences.R` | Niche-position analysis (2,334 species, split into tropical/non-tropical subsets). Key objects: `data_position*` (species-level `RC1`-`RC3` position by `mycorrhizal_type`/`biome`), `pem_position*` (phylogenetic eigenvector map predictors, one `PEM` column per eigenvector), `phylo_prox_position*` (phylogenetic proximity matrices), `RC1_results`/`RC2_results`/`RC3_results` (fitted model results per component), plus model parameter lists. |
| `figure_5.RData` | `05a.Species_niche_breadth_differences.R` | Niche-breadth analysis, same structure as `figure_4.RData`: `data_breadth*` (species-level `env_breadth` by `mycorrhizal_type`/`biome`), `pem_breadth*`, `phylo_prox_breadth*`, `pairwise_results_breadth` (18 x 5: pairwise comparisons between mycorrhizal types), `env_breadth_results` (fitted model results). |
| `figure_6.RData` | `08_Mid-domain_null_models.R` | `data` (2,334 x 7: the observed latitude/`env_breadth` points), `pred_combined` (13,604 x 7: fitted latitude-breadth model with CIs), `mde_combined` (13,604 x 5: the mid-domain-effect null model's predicted breadth with CIs), `prop_within_mde` and `r2_combined` (4 x 2 each: summary statistics per mycorrhizal type). |

**Figure S8 (absolute richness)**

| File | Produced by | Contents |
|---|---|---|
| `figure_S8a.tif` | `04c.Absolute_richness_analysis.R` | Raster, 4 layers (AM, EcM, EcM-AM, NM), 373 x 492 cells. Predicted absolute richness surfaces. |
| `figure_S8b.RData` | `04d.Absolute_richness_latitude.R` | Same structure as `figure_2b.RData` but for absolute richness: `latitude_gradient_data_abs` (47,116 x 3, with `richness` in place of `relative_richness`), `latitude_gradient_marginal_effects_abs`, `latitude_gradient_coeficients_abs`. |
| `figure_S8c.RData` | `04c.Absolute_richness_analysis.R` | Same structure as `figure_2c.RData` but for absolute richness: `data_figure_c` (290 x 6), `data_figure_c_limits` (1 x 6). |

**Empirical pathway** (measured rather than genus-inferred mycorrhizal types)

| File | Produced by | Contents |
|---|---|---|
| `empirical_measurements.txt` | Not regenerated by this pipeline | 1,436 rows, 22 columns: `fungal_root_id`, `family`, `genus`, `species`, `species_verbatim`, `mycorrhizal_type`, `replicates`, `habitat`, `host_age`, `host_age_group`, `am_evaluated`/`am_criteria`/`am_structures`, `ecm_evaluated`/`ecm_structures`, `erm_structure`, `nm_structures`, `latitude`, `longitude`, `comment`, `source`, `link`. A curated compilation of literature-sourced mycorrhizal measurements. Originally produced by `get_fungal_root_sources.R`, which is kept locally but is no longer tracked or part of this pipeline (see that script's note in the repository) -- treat this file as a fixed input rather than something reproducible from the other scripts here. |
| `empirical_mycorrhizal_type.txt` | `07a.Generate_empirical_dataset.R` | 290 rows: `family`, `genus`, `species`, `mycorrhizal_type` -- the curated measurements above, resolved to one type per species. |
| `empirical_genus_consensus.txt` | `07a.Generate_empirical_dataset.R` | 168 rows: `family`, `genus`, per-genus counts by type (`AM`, `EcM`, `Dual`, `NM`), `consensus_mycorrhizal_type`, `consensus_type_proportion`, plus `fungal_root_mycorrhizal_type` and `match` for comparison against FungalRoot's genus-level call. |
| `empirical_niche_estimates.txt` | `07a.Generate_empirical_dataset.R` | 266 rows, same structure as `niche_estimates.txt` above, restricted to species with a measured (not genus-inferred) mycorrhizal type. |
| `figure_2a_empirical.tif` | `07c.Relative_richness_empirical.R` | Same structure as `figure_2a.tif`, built from the empirical species set. |
| `figure_2c_empirical.RData` | `07c.Relative_richness_empirical.R` | Same structure as `figure_2c.RData`, built from the empirical species set. |
| `figure_3_empirical.RData` | `07c.Relative_richness_empirical.R` | Same structure as `figure_3.RData` (39,040 x 5 for `relative_richness_data`, reflecting the smaller empirical species set), built from the empirical species set. |
| `figure_5_empirical.RData` | `07d.Niche_breadth_differences_empirical.R` | Same structure as `figure_5.RData` (264 species rather than 2,334), built from the empirical species set. |

## Dependencies
This project is conducted using R version 4.3.3 (2024-02-29) -- "Angel Food Cake" and the following packages:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) version 1.17.8
* [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html) version 2.0.0
* [CoordinateCleaner](https://cran.r-project.org/web/packages/CoordinateCleaner/index.html) version 3.0.1
* [WorldFlora](https://cran.r-project.org/web/packages/WorldFlora/index.html) version 1.14.5
* [tidyverse](https://www.tidyverse.org/blog/2023/03/tidyverse-2-0-0/) version 2.0.0
* [terra](https://rspatial.github.io/terra/index.html) version 1.8.54
* [rnaturalearth](https://cran.r-project.org/web/packages/rnaturalearth/index.html) version 1.0.1
* [sf](https://r-spatial.r-universe.dev/sf) version 1.0.21
* [psych](https://cran.r-project.org/web/packages/psych/index.html) version 2.4.6.26
* [viridis](https://cran.r-project.org/web/packages/viridis/index.html) version 0.6.5
* [fst](https://cran.r-project.org/web/packages/fst/index.html) version 0.9.8
* [blockCV](https://github.com/rvalavi/blockCV) version 3.1.6
* [ENMeval](https://jamiemkass.github.io/ENMeval/) version 2.0.4
* [dsmextra](https://densitymodelling.github.io/dsmextra/) version 1.1.5
* [ecospat](https://ecospat.r-universe.dev/ecospat) version 4.1.2
* [maxnet](https://mrmaxent.r-universe.dev/maxnet) version 0.1.4
* [lhs](https://bertcarnell.github.io/lhs/) version 1.2.0
* [INLA](https://www.r-inla.org) version 24.2.9
* [fmesher](https://cran.r-project.org/web/packages/fmesher/index.html) version 0.2.0
* [boot](https://cran.r-project.org/web/packages/boot/index.html) version 1.3.30
* [emmeans](https://cran.r-project.org/web/packages/emmeans/index.html) version 1.10.1
* [ape](https://cran.r-project.org/web/packages/ape/index.html) version 5.8
* [adephylo](https://cran.r-project.org/web/packages/adephylo/index.html) version 1.1.16
* [MPSEM](https://cran.r-project.org/web/packages/MPSEM/index.html) version 0.6.1
* [V.PhyloMaker2](https://github.com/jinyizju/V.PhyloMaker2) version 0.1.0
* [parameters](https://easystats.github.io/parameters/) version 0.24.2
* [ggtext](https://cran.r-project.org/web/packages/ggtext/index.html) version 0.1.2
* [ggtree](https://guangchuangyu.github.io/software/ggtree/) version 3.10.1
* [ggtreeExtra](https://github.com/YuLab-SMU/ggtreeExtra) version 1.19.990
* [ggnewscale](https://eliocamp.github.io/ggnewscale/) version
* [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html) version 1.1.3
* [patchwork](https://github.com/thomasp85/patchwork) version 1.2.0
* [ggpubr](https://github.com/kassambara/ggpubr) version 0.6.0
* [betareg](https://cran.r-project.org/web/packages/betareg/index.html) version 3.2.2
* [performance](https://easystats.github.io/performance/) version 0.13.0
* [ggeffects](https://strengejacke.github.io/ggeffects/) version 1.5.2
* [ggridges](https://cran.r-project.org/web/packages/ggridges/index.html) version 0.5.6
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html) version 1.5.2
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) version 1.0.17
* [future.apply](https://future.apply.futureverse.org/) version 1.11.3
* [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html) version 1.6.5
* [RANN](https://cran.r-project.org/web/packages/RANN/index.html) version 2.6.2
* [gstat](https://cran.r-project.org/web/packages/gstat/index.html) version 2.1.1
* [sp](https://cran.r-project.org/web/packages/sp/index.html) version 2.2.0
* [phyloMEM](https://github.com/LukeLikesDirt/phyloMEM) version 0.0.0.9000 -- not on CRAN, install with `devtools::install_github("LukeLikesDirt/phyloMEM")`
