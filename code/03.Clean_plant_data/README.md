# Harmonise plant databases according to World Flora Online

## Overview:

[World Flora Online](http://www.worldfloraonline.org/) is an online flora of all known plants, used to match
a list of plant names against a static copy of the World Flora Online Taxonomic Backbone. 

Here I use the WorldFlora package to harmonise the plant names in the [Biomass Plot Library](https://portal.tern.org.au/metadata/23218#:~:text=The%20Biomass%20Plot%20Library%20is,private%20companies%20and%20other%20agencies.), [Harmonised Australian Vegetation Plot (HAVPlot)](https://data.csiro.au/collection/csiro:54461?_st=browse&_str=6&_si=2&browseType=kw&browseValue=vegetation), [FungalRoot](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18207), [Root nodulation (NodDB)](https://onlinelibrary.wiley.com/doi/10.1111/jvs.12627), and [GlobalTreeSearch](https://tools.bgci.org/global_tree_search.php) databases accoding to the World Flora Online taxonomy.

The World Flora Online provides plant name synonyms in its output. The Biomass Plot Library and HAVPlot will subsequently be used for modeling the distribution and diversity of Australian forest tree species. Therefore, the harmonised outputs for species names in these databases will be reduced to a single species name. The FungalRoot, NodDB, and GlobalTreeSearch databases, however, will be used to identify the mycorrhizal type, root nodulation potential, and whether a plant species is classified as a tree, respectively. As such, all synonyms for species names in these databases will be retained to facilitate the correct annotation of mycorrhizal type, root nodulation potential, and tree status in the Biomass Plot Library and HAVPlot databases.