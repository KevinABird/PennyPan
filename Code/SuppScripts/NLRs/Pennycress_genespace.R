library(GENESPACE)
#R version 4.4.0 (2024-04-24)
#GENESPACE_1.3.1
Sys.getenv("PATH")
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "~minimap2", sep = ":"))
minimap2=".~/minimap2/minimap2"

path2mcscanx ="~/MCScanX"

setwd("~/Pennycress")

# working directories 
wd <- "~/Pennycress/workingDirectory_with_TAIR"

##### Run GENESPACE across T. arvense accessions and other Brassica relatives #####

# -- parse the annotations to fastas with headers that match a gene bed file - add or remove genome IDs
parsedPaths <- parse_annotations(
  rawGenomeRepo = "~/Pennycress/genomes",
  genomeDirs = c("TAIR11", "TarvensevarAK34W", "TarvensevarAmes32873", "TarvensevarLorettoMN", "TarvensevarMN106", "TarvensevarMN134", "TarvensevarPI650286", "TarvensevarTibet33", "A_lyrata", "B_rapa"),
  genomeIDs = c("TAIR11", "TarvensevarAK34W", "TarvensevarAmes32873", "TarvensevarLorettoMN", "TarvensevarMN106", "TarvensevarMN134", "TarvensevarPI650286", "TarvensevarTibet33", "A_lyrata", "B_rapa"),
  presets = "phytozome",
  genespaceWd = wd)

parsedPaths <- parse_annotations(
  rawGenomeRepo = "~/Pennycress/genomes",
  genomeDirs = c( "TarvensevarAK34W", "TarvensevarAmes32873", "TarvensevarLorettoMN", "TarvensevarMN106", "TarvensevarMN134", "TarvensevarPI650286", "TarvensevarTibet33"),
  genomeIDs = c( "TarvensevarAK34W", "TarvensevarAmes32873", "TarvensevarLorettoMN", "TarvensevarMN106", "TarvensevarMN134", "TarvensevarPI650286", "TarvensevarTibet33"),
  presets = "phytozome",
  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar, overwrite = F)
save(out, file="Pennycress_genespace_var_genomes_with_TAIR11_20241024.RDA")

