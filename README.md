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

1. **Size.** The raw inputs total over 100 GB (versus around 130 MB for `generated_data/`), exceeding what is practical to version-control or host as a supplementary deposit.
2. **Redistribution terms.** Some raw sources are openly licensed and could in principle be redistributed (e.g. the Harmonised Australian Vegetation Plot dataset and Forests of Australia (2023) are both CC BY 4.0). Others — including the GBIF occurrence download and BGCI's GlobalTreeSearch — are aggregated or terms-of-use-restricted third-party products that we have not confirmed we can redistribute in bulk. Rather than apply different handling per source, no raw data is provided; `01.Harmonise_tree_data.R` and the other raw-data scripts listed below acquire it directly from source. Licences and access routes for every covariate are listed in `covariates.txt` and plant database in `./plant_data.txt`.

## Reproducing the analysis

There are two ways to use this repository, depending on whether you need to regenerate the processed datasets or only primary figures.

### Path A — primnary figures from the generated data

Every script listed as `TRUE` in the "Generated data only" column of the Outputs table below reads exclusively from `generated_data/` (and, for context/colour scales, `covariates.txt`). Clone the repository and run the relevant script directly.

### Path B — the full pipeline from raw data

Reproducing `generated_data/` itself, or anything marked `FALSE` below, means acquiring the raw data first. Run the following in order; each stage's output feeds the next.

1. **`01.Harmonise_tree_data.R`** — Downloads and harmonises the World Flora Online backbone, FungalRoot, GlobalTreeSearch, HAVPlot and GBIF, then builds the tree-occurrence dataset and `generated_data/global_tree_mycorrhizal_types.txt`. The Australian Plant Census (APC) step requires a manual export (see the script's inline instructions) because the APC has no stable direct-download URL.
2. **`build_phylogenetic_tree.R`** — Builds `generated_data/phylo_tree_mycorrhizal_types.tre` from (1)'s output. Required by `05a`, `05b`, `07d` and `Figure_1.R`.
3. **`02a.Upscale_predictor_rasters.R`** — Downloads and upscales the forest-cover raster and climate/soil covariate rasters.
4. **`02b.Prepare_prediction_grids.R`** — Builds the 10 km prediction grid and site-level presence data. Produces Figs. S2–S4 and Table S1.
5. **`03a.Niche_cache_enmeval.R` → `03b.Niche_est_enmeval.R` → `03c.Niche_run_enmeval.R`** — Species distribution modelling with ENMeval. This is the computationally intensive step referenced above (~2 days on 12 cores; see `envs/`).
6. **`03d.Niche_est_summary.R`** — Summarises the modelling output into `generated_data/niche_estimates.txt`.
7. **`04a.Relative_richness_analysis.R`, `04b.Relative_richness_latitude.R`, `04c.Absolute_richness_analysis.R`, `04d.Absolute_richness_latitude.R`** — Richness analyses. Produce `generated_data/figure_2*`, `figure_3.RData` and `figure_S8*`, plus Figs. S5, S6, S9–S11 and Table S2.
8. **`05a.Species_niche_breadth_differences.R`, `05b.Species_niche_position_differences.R`** — Niche breadth and position analyses. Produce `generated_data/figure_4.RData` and `figure_5.RData`, plus Tables S3–S7.
9. **`06a.Sensetivity_relative_richness.R`** — Sensitivity analysis run directly on raw presence data (Figs. S14–S16).
10. **`07a.Generate_emperical_dataset.R` → `07b.Prepare_emperical_prediction_grids.R` → `07c.Relative_richness_emperical.R` → `07d.Niche_breadth_differences_emperical.R`** — The empirical pathway, using measured rather than genus-inferred mycorrhizal types. Produces Table S8 and Dataset S2.

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
| Fig. S14 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) | FALSE |
| Fig. S15 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) | FALSE |
| Fig. S16 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) | FALSE |
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
| Table S8 | [code/07a.Generate_emperical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_emperical_dataset.R) | FALSE |
| Dataset S1 | [code/01.Harmonise_tree_data.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/01.Harmonise_tree_data.R) | FALSE |
| Dataset S2 | [code/07a.Generate_emperical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_emperical_dataset.R) | FALSE |

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
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) version 2.3
* [gtable](https://r-lib.r-universe.dev/gtable) version 0.3.6
* [patchwork](https://github.com/thomasp85/patchwork) version 1.2.0
* [ggpubr](https://github.com/kassambara/ggpubr) version 0.6.0
