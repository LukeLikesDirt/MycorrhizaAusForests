
require(tidyverse, lib.loc = '/data/group/frankslab/project/LFlorence/MycorrhizaAusForests/envs/R-packages')

setwd("/data/group/frankslab/project/LFlorence/MycorrhizaAusForests")

otu <- "data/AusMicrobiome/ITS"

OTUs = read.csv("data/AusMicrobiome/ITS/08.Clustered/OTUs.txt", header = T, sep = " ")

