##############################################################################*
##### Harmonise data using the World Flora database ##########################
##############################################################################*

require(WorldFlora)
require(stringr)
require(tidyverse)

# World Flora taxonomic backbone dataset: Downloaded from www.worldfloraonline.org/downloadData
WFO = read.csv ('data/WorldFloraOnline/classification.csv', sep = '')
# My list of tree species
treelist = read.csv('data/BiomassPlotLib/biolib_treelist.csv') %>%
  rename(taxon = species) %>%
  # remove 'subspecies', 'variety', and 'authority' annotations (keep genus and spices names only)
  mutate(taxon = word(taxon, 1, 2)) %>%
  mutate(genus_old = word(taxon, 1)) %>%
  mutate(species_old = word(taxon, 2)) %>%
  # remove observations without generic annotations
  mutate(taxon = na_if(taxon,'')) %>%
  filter(!taxon %in% NA) %>%
  as_tibble()
glimpse(treelist)

spec.list = treelist %>%
  select(spec.name = taxon) %>%
  unique()

spec.check = WFO.match(spec.data = spec.list, WFO.data = WFO, counter = 1, verbose = TRUE)

gen.list = treelist %>%
  select(Genus = genus_old) %>%
  unique()

gen.check = WFO.match(gen.list$Genus[1:291], WFO.file = WFO, counter = 1, verbose = T)




require(Taxonstand)
data(bryophytes)

spec.check = WFO.match(spec.list[1:1031,], WFO.data = WFO, 
                       spec.name = 'spec.name', counter = 1, verbose = TRUE)

w1 <- WFO.match(bryophytes[1:20, ], WFO.file=WFO, spec.name="Full.name", counter=1)
w1


