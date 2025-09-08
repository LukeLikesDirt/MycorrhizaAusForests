## Introduction

This repository holds code and data associated to the manuscript:

[**Rethinking Mycorrhizal Ecology: Niche Patterns of Dual-mycorrhizal and Non-Mycorrhizal Trees**]()

**Authors:**
Luke Florence<sup>1</sup>, John W. Morgan<sup>1</sup>, Peter A. Vesk<sup>2</sup>, Jen L. Wood<sup>3</sup> Camille Truong<sup>4</sup>

**Affiliations:**
1. Department of Environment and Genetics, La Trobe University, Bundoora, Victoria, Australia. 
2. School of Agriculture, Food and Ecosystem Sciences, University of Melbourne, Parkville, Victoria, Australia. 
3. Department of Department of Microbiology, Anatomy, Physiology and Pharmacology, La Trobe University, Bundoora, Victoria, Australia. 
4. Royal Botanic Gardens Victoria, Melbourne, Victoria, Australia.

Corresponding author: Luke Florence (L.Florence@latrobe.edu.au) 

## Repository contents

* `code/` — Contains all R scripts used to reproduce the analysis, including data sourcing, preparation, execution of statistical workflows, and figure generation
* `data/` — Will store all sourced and prepared data required for datasets required statistical workflows
* `code/Figure*.R` - Individual scripts for generating each figure presented in the manuscript
* `output/generated_data/` — Includes all datasets required to reproduce the primary figures and results
* `envs/` — Contains the Conda environment set-up scripts used for modelling environmental breadth. This analysis was executed on a high-performance computing cluster due to computational intensity

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


