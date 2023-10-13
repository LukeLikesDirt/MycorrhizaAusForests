## Introduction

This repository is dedicated to the first chapter of my PhD at La Trobe University—a project aimed at mapping patterns of root-trait diversity and endemism in Australian forests.

This project encompasses eight main tasks:

1. **Bioinformatics Analysis:** Processing and analysing fungal (ITS) metabarcoding data from the [Australian Microbiome Initiative (AMI)](https://www.australianmicrobiome.com/).

2. **Bacterial Metabarcoding:** Similar to the ITS bioinformatics analysis, but this time with bacterial (16S) metabarcoding data from AMI.

3. **Data Harmonisation:** Cleaning and harmonising plant species names from various databases, including the [Biomass Plot Library](https://portal.tern.org.au/metadata/23218#:~:text=The%20Biomass%20Plot%20Library%20is,private%20companies%20and%20other%20agencies.), [Harmonised Australian Vegetation Plot (HAVPlot)](https://data.csiro.au/collection/csiro:54461?_st=browse&_str=6&_si=2&browseType=kw&browseValue=vegetation), [FungalRoot](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18207), and [NodDB](https://onlinelibrary.wiley.com/doi/10.1111/jvs.12627).

4. **Root Trait Prediction:** Estimating proportions of root traits (arbuscular mycorrhizal, ectomycorrhizal, and nitrogen-fixing) across Australian forests based on the Biomass Plot Library dataset. These predictions will serve as covariates in subsequent statistical analyses.

5. **Covariate Extraction:** Extracting environmental and vegetation covariates for each AMI and HAVPlot site to be used in subsequent statistical analyses.

6. **Diversity Drivers:** Determining the dominant drivers of root trait diversity for both microbes and plants.

7. **Vulnerability Assessment:** Evaluating the vulnerability of root traits to global change drivers for both microbes and plants.

8. **Endemism Mapping:** Creating maps that showcase patterns of root trait endemism for both microbes and plants.

This repository holds the code and tools necessary to carry out these tasks and contribute to the understanding of root trait diversity in Australian forests. The repository is designed to be entirely reproducible, ensuring that the research can be independently verified and built upon.