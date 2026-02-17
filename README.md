## Introduction

This repository holds code and data associated to the manuscript:

[**Uncovering Environmental Niches in Dual-mycorrhizal and Non-mycorrhizal Forest Trees**]()

**Authors:**
Luke Florence<sup>1</sup>, John W. Morgan<sup>1</sup>, Peter A. Vesk<sup>2</sup>, Jen L. Wood<sup>3</sup>, Camille Truong<sup>4</sup>

**Affiliations:**
1. Department of Environment and Genetics, La Trobe University, Bundoora, Victoria, Australia. 
2. School of Agriculture, Food and Ecosystem Sciences, University of Melbourne, Parkville, Victoria, Australia. 
3. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, Victoria, Australia. 
4. Royal Botanic Gardens Victoria, Melbourne, Victoria, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au) 

## Repository contents

* `./code/` — Contains all R scripts required to reproduce the analysis, including data sourcing, preparation, execution of statistical workflows, and figure generation
* `./generated_data/` — Includes all datasets required to reproduce the primary and figures and results presented in the manuscript
* `./envs/` — Contains the Conda environment set-up scripts used for modelling environmental breadth. This analysis was computationally intensive and takes around 2 days on 12 cores.
* `./covariates.txt` — A mapping of covariate names to their sources and descriptions.
* **Note** — Raw tree occurrence data retrieval and curation can be reproduced in script `./code/01.Harmonise_tree_data.R`

## Outputs
To reproduce the outputs of this project, see the following scripts:

| Output | Script |
|--------|--------|
| Fig. 1 | [code/Figure_1.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_1.R) |
| Fig. 2 | [code/Figure_2.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_2.R) |
| Fig. 3 | [code/Figure_3.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_3.R) |
| Fig. 4 | [code/Figure_4.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_4.R) |
| Fig. 5 | [code/Figure_5.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_5.R) |
| Fig. 6 | [code/Figure_6.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_6.R) |
| Fig. S1 | [code/Figure_S1.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S1.R) |
| Fig. S2 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) |
| Fig. S3 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) |
| Fig. S4 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) |
| Fig. S5 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) |
| Fig. S6 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) |
| Fig. S7 | [code/Figure_S7.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S7.R) |
| Fig. S8 | [code/Figure_S8.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S8.R) |
| Fig. S9 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) |
| Fig. S10 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) |
| Fig. S11 | [code/04c.Absolute_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04c.Absolute_richness_analysis.R) |
| Fig. S12 | [code/Figure_S12.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S12.R) |
| Fig. S13 | [code/Figure_S13.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S13.R) |
| Fig. S14 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) |
| Fig. S15 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) |
| Fig. S16 | [code/06a.Sensetivity_relative_richness.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06a.Sensetivity_relative_richness.R) |
| Fig. S17 | [code/Figure_S17.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S17.R) |
| Fig. S18 | [code/Figure_S18.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S18.R) |
| Fig. S19 | [code/06b.Sensitivity_niche_breadth.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/06b.Sensitivity_niche_breadth.R) |
| Fig. S20 | [code/Figure_S20.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/Figure_S20.R) |
| Table S1 | [code/02b.Prepare_prediction_grids.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/02b.Prepare_prediction_grids.R) |
| Table S2 | [code/04a.Relative_richness_analysis.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/04a.Relative_richness_analysis.R) |
| Table S3 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) |
| Table S4 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) |
| Table S5 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) |
| Table S6 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) |
| Table S7 | [code/05a.Species_niche_breadth_differences.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/05a.Species_niche_breadth_differences.R) |
| Table S8 | [code/07a.Generate_emperical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_emperical_dataset.R) |
| Dataset S1 | [code/01.Harmonise_tree_data.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/01.Harmonise_tree_data.R) |
| Dataset S2 | [code/07a.Generate_emperical_dataset.R](https://github.com/LukeLikesDirt/MycorrhizaAusForests/blob/main/code/07a.Generate_emperical_dataset.R) |

## Dependencies
This project is conducted using R version 4.3.3 (2024-02-29) -- "Angel Food Cake" and the following packages:

* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) version 1.17.8
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