#############################################
##### Set directories
#############################################
if(Sys.info()['sysname'] == "Windows"){
  # paths Win
  wd <- "C:/Users/Manuel/Nextcloud/Masterarbeit/"
  db <- "C:/Users/Manuel/Dropbox/Apps/Overleaf/Masterarbeit"
}else if(Sys.info()['sysname'] == "Linux"){
  # paths Lin
  wd <- "/home/manuel/Nextcloud/Masterarbeit/"
  db <- "/home/manuel/Dropbox/Apps/Overleaf/Masterarbeit"
}else{
  print("Error: OS not identified.")
}

dir_fig <- file.path(wd, "fig")

#############################################
##### Load packages
#############################################
packages <- c("sp", "sf", "dplyr", "reshape2", "betapart", "vegan",
              "ggplot2", "WorldFlora", "xtable", "stringi")
for(i in 1:NROW(packages)){
  if(!require(packages[i], character.only = TRUE)){
    install.packages(packages[i])
    library(packages[i], character.only = TRUE)
  }
}

# ArcGIS connection
if(!require("arcgisbinding", character.only = TRUE)){
  install.packages("arcgisbinding", repos = "https://r.esri.com", type = "win.binary")
  library("arcgisbinding", character.only = TRUE)
}

#############################################
##### Set style options
#############################################
theme_set(theme_bw())
cols <- c(
  rgb(0, 150, 130, alpha = 254.99, maxColorValue = 255), #kit colour
  rgb(70, 100, 170, alpha = 254.99, maxColorValue = 255), #kit blue
  rgb(223, 155, 27, alpha = 254.99, maxColorValue = 255), #kit orange
  rgb(140, 182, 60, alpha = 254.99, maxColorValue = 255), #kit Mai green
  rgb(162, 34, 35, alpha = 254.99, maxColorValue = 255), #kit red
  rgb(163, 16, 124, alpha = 254.99, maxColorValue = 255), #kit violet
  rgb(167, 130, 46, alpha = 254.99, maxColorValue = 255) #kit braun
)
colz <- c(
  rgb(0, 150, 130, alpha = 255*0.5, maxColorValue = 255), #kit colour
  rgb(70, 100, 170, alpha = 255*0.5, maxColorValue = 255), #kit blue
  rgb(223, 155, 27, alpha = 255*0.5, maxColorValue = 255), #kit orange
  rgb(140, 182, 60, alpha = 255*0.5, maxColorValue = 255), #kit Mai green
  rgb(162, 34, 35, alpha = 255*0.5, maxColorValue = 255), #kit red
  rgb(163, 16, 124, alpha = 255*0.5, maxColorValue = 255), #kit violet
  rgb(167, 130, 46, alpha = 255*0.5, maxColorValue = 255) #kit braun
)

#############################################
##### Load data
#############################################
# ArcGIS Online data
arc.check_product()

Blocks <- arc.data2sf(
  arc.select(
    arc.open(
      "https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/LELEBlockOutline/FeatureServer/0"
    )
  )
)[, c("OBJECTID", "LELE_Block", "Block", "Plot", "Treatment")]
#Blocks <- st_transform(Blocks, crs = "+proj=longlat +datum=WGS84")

Trees <- arc.data2sf(
  arc.select(
    arc.open(
      "https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/Trees_from_scratch/FeatureServer/0"
    )
  )
)[, c("OBJECTID", "Species", "TreeDiameter1", "SpeciesID")]
Trees <- st_transform(Trees, crs = st_crs(Blocks))
Trees <- Trees[Trees$Species != "Other",]

# Shapefiles
require("rgdal")
Lapalala <- readOGR(file.path(wd, "gis", "QGIS", "Shapefiles", "Lapalala", "LapalalaBorders.shp"))

#############################################
##### Analyse species composition in plots
#############################################
TreesWithPlots <- st_join(x = Trees, y = Blocks)
AbMatr <- data.frame(TreesWithPlots) %>%
  group_by(Block, Plot) %>%
  summarise(Species) %>%
  count(Species)

AbundanceMatrix <- dcast(data = AbMatr, Block + Plot ~ Species, value.var = "n")
AbundanceMatrix[is.na(AbundanceMatrix)] <- 0
row.names(AbundanceMatrix) <- paste("B", AbundanceMatrix$Block, "_", AbundanceMatrix$Plot, sep = "")

if(length(which(AbundanceMatrix$Block == 0 | AbundanceMatrix$Plot == 0)) > 0){
  AbundanceMatrix <- AbundanceMatrix[-which(AbundanceMatrix$Block == 0 | AbundanceMatrix$Plot == 0),]
}
if(length(which(colnames(AbundanceMatrix) %in% c("Block", "Plot"))) > 1){
  AbundanceMatrix <- AbundanceMatrix[, -c(which(colnames(AbundanceMatrix) %in% c("Block", "Plot")))]
}
if(length(which(colSums(AbundanceMatrix) == 0)) > 0){
  AbundanceMatrix <- AbundanceMatrix[, -which(colSums(AbundanceMatrix) == 0)]
}

AbundanceMatrix <- AbundanceMatrix[which(rowSums(AbundanceMatrix) > 40),]

# Bray-distance
bray_matr <- beta.pair.abund(AbundanceMatrix)$beta.bray

paMatrix <- AbundanceMatrix
paMatrix[paMatrix >= 1] <- 1
sor_matr <- beta.pair(paMatrix)$beta.sor
sim_matr <- beta.pair(paMatrix)$beta.sim
sne_matr <- beta.pair(paMatrix)$beta.sne

# beta diversity on abundance
range(bray_matr)

# beta diversity on pa
range(sor_matr)

# alpha diversity
range(rowSums(paMatrix))

#############################################
##### Print species list
#############################################
options(timeout = max(300000, getOption("timeout")))
WFO.download(WFO.url = "http://104.198.143.165/files/WFO_Backbone/_WFOCompleteBackbone/WFO_Backbone.zip",
             save.dir = tempdir(), WFO.remember = TRUE)
WFO.remember(WFO.file = file.path(tempdir(), "WFO_Backbone.zip"), WFO.data = "WFO.data", WFO.pos = 1)
dir.create(file.path(tempdir(), "WFO_Backbone"), showWarnings = FALSE)
setwd(file.path(tempdir(), "WFO_Backbone"))
unzipped <- unzip(file.path(tempdir(), "WFO_Backbone.zip"))
WFO.file.RK <- file.path(tempdir(), "WFO_Backbone", "classification.txt")
WFO.data1 <- data.table::fread(WFO.file.RK, encoding = "UTF-8")

all_species <- unique(TreesWithPlots$Species)
all_species[which(all_species == "Euclea linearis")] <- "Euclea crispa"

df_all_spec <- data.frame(Family = rep(NA, length(all_species)),
                          Genus = rep(NA, length(all_species)),
                          Epithet = rep(NA, length(all_species)),
                          Authority = rep(NA, length(all_species)),
                          Source = rep(NA, length(all_species)))

get_infos <- function(Species = NULL){
  WFO.browse(Species, WFO.data = WFO.data1, accepted.only = TRUE)
  spec_info <- WFO.match(Species, WFO.data = WFO.data1)
  return(c(spec_info$family[1], spec_info$genus[1], spec_info$specificEpithet[1],
           spec_info$scientificNameAuthorship[1],
           spec_info$namePublishedIn[1]))
}

for(i in 1:length(all_species)){
  df_all_spec[i, ] <- get_infos(all_species[i])
}

df_all_spec$Genus <- paste0("openitshape", df_all_spec$Genus, "closeitshape")
df_all_spec$Epithet <- paste0("openitshape", df_all_spec$Epithet, "closeitshape")

df_all_spec[which(df_all_spec$Epithet == "openitshapecrispacloseitshape")[1],
            c("Epithet", "Authority", "Source")] <- c(
              "openitshapecrispacloseitshape subsp. openitshapelineariscloseitshape",
              "(Zeyh. ex Hiern) F.White",
              "Bull. Jard. Bot. Natl. Belg.")

# order by species name
df_all_spec <- df_all_spec[with(df_all_spec,
                                order(df_all_spec$Genus, df_all_spec$Epithet)), ]

print(xtable(df_all_spec,
             caption = "List of all tree species found in the research plots. Taxon names and information were semiautomatically matched and manuelly checked following World Flora Online (Kindt, 2021; WFO, 2022).",
             label = "tab:Secies_list"),
      table.placement = "pbht",
      caption.placement = "top",
      include.rownames = FALSE,
      booktabs = TRUE,
      file = file.path(db, "tab", "Species_List.tex"),
      add.to.row = list(list(nrow(df_all_spec)),  
                        '\\bottomrule\n\\multicolumn{5}{l}{Kindt, Roeland (2020). \\enquote{WorldFlora: An R Package for Exact and Fuzzy Matching of Plant Names against
the World Flora Online Taxonomic Backbone Data}.}\\\\\n \\multicolumn{5}{l}{\\hspace{1em} In: openitshapeApplications in Plant Sciencescloseitshape 8. 9, e11388. issn: 2168-0450. doi: \\href{https://doi.org/10.1002/aps3.11388}{10.1002/aps3.11388}.}\\\\\n \\multicolumn{5}{l}{WFO (2022): World Flora Online. Published on the Internet; \\url{http://www.worldfloraonline.org}. Accessed on: 06 Feb 2022.}\\\\'))
# adjust table size
text_file <- readLines(file.path(db, "tab", "Species_List.tex"))
text_file <- c(text_file[1:6], "\\resizebox{\\textwidth}{!}{", text_file[7:length(text_file)])
text_file <- c(text_file[1:length(text_file)-1], "}", text_file[length(text_file)])
text_file <- gsub("1947 \\[Sep 1947\\]", "(1947)", text_file, perl = TRUE)
text_file <- gsub("caption", "caption[Species list]", text_file)
text_file <- gsub("ä", '\\\\"a', text_file)
text_file <- gsub("ü", '\\\\"u', text_file)
text_file <- gsub("é", "\\\\'e", text_file)
text_file <- stri_replace_last(text_file, fixed = "\\\\ \\bottomrule", replacement = "\\\\")
text_file <- gsub("openitshape", "\\\\textit\\{", text_file)
text_file <- gsub("closeitshape", "\\}", text_file)
writeLines(text_file, file.path(db, "tab", "Species_List.tex"))

#############################################
##### Diversity analyses
#############################################
plot(hclust(bray_matr))
plot(hclust(sor_matr))

dca <- decorana(bray_matr)
summary(dca)

RDA <- rda(AbundanceMatrix)
pdf(file.path(wd, "fig", "RDA_trees.pdf"), width = 7, height = 7)
biplot(RDA, display = c("species", "sites"), type = c("text", "text"), col = cols[c(2, 1)])
#ordihull(RDA, groups = "species")
dev.off()

RDA <- rda(paMatrix)
biplot(RDA, display = c("species", "sites"), type = c("text", "text"))

CCA <- cca(AbundanceMatrix)
plot(CCA, type = "n", dis = "sp", col = cols[1])
points(CCA, dis = "sp")
decorana(AbundanceMatrix)
decorana(paMatrix)

plot(as.matrix(bray_matr)[, 1] ~ seq(1:length(as.matrix(bray_matr)[, 1])),
     xlab = "Plot number", ylab = expression(beta[bray]))
plot(as.matrix(sor_matr)[, 1] ~ seq(1:length(as.matrix(sor_matr)[, 1])),
     xlab = "Plot number", ylab = expression(beta[SOR]),
     col = cols[1])
points(as.matrix(sim_matr)[, 1] ~ seq(1:length(as.matrix(sim_matr)[, 1])),
       col = cols[2])
points(as.matrix(sne_matr)[, 1] ~ seq(1:length(as.matrix(sne_matr)[, 1])),
       col = cols[3])

#############################################
##### Vegetation Sampling
#############################################
Herb_Species <- Veg_Species[Veg_Species$layer == "grass",]
PlotGlobalIDs <- unique(Herb_Species$parentglobalid)
get_quadrat_ID <- function(globalID){
  block <- Veg_Plots[Veg_Plots$globalid == globalID, "blockNumber"]$blockNumber
  plot <- Veg_Plots[Veg_Plots$globalid == globalID, "plotNumber"]$plotNumber
  quadrat <- Veg_Plots[Veg_Plots$globalid == globalID, "quadratName"]$quadratName
  quadratID <- paste0("B", as.character(block), as.character(plot), quadrat,
                      collapse = "")
  as.character(quadratID)
  return(quadratID)
}
get_quadrat_info <- function(globalID){
  coverTrees <- Veg_Plots[Veg_Plots$globalid == globalID, "coverTrees"]$coverTrees
  coverShrubs <- Veg_Plots[Veg_Plots$globalid == globalID, "coverShrubs"]$coverShrubs
  coverRocks <- Veg_Plots[Veg_Plots$globalid == globalID, "coverRocks"]$coverRocks
  return(c(coverTrees, coverShrubs, coverRocks))
}

Herb_Species$quadratID <- vapply(X = Herb_Species$parentglobalid,
                              FUN = get_quadrat_ID, FUN.VALUE = character(1))
Plot_infos <- data.frame(Tree_cover = rep(NA, nrow(Herb_Species)), Shrub_cover = NA, Rock_cover = NA)
Herb_Species <- Herb_Species[nchar(Herb_Species$quadratID) < 10,]

SpeciesList <- unique(Herb_Species$epithet5)
QuadratList <- unique(Herb_Species$quadratID)

VegAbundanceMatrix <- matrix(nrow = length(QuadratList),
                             ncol = length(SpeciesList))
row.names(VegAbundanceMatrix) <- QuadratList
colnames(VegAbundanceMatrix) <- SpeciesList

"
writeAbundance <- function(quadrat, taxon, abundance){
  VegAbundanceMatrix[QuadratList == quadrat,
                     SpeciesList == taxon] <<- abundance
}

mapply(writeAbundance,
       quadrat = Herb_Species$quadratID,
       taxon = Herb_Species$epithet5,
       abundance = Herb_Species$cover)
"
for(i in 1:nrow(Herb_Species)){
  VegAbundanceMatrix[which(QuadratList == Herb_Species$quadratID[i]),
                     which(SpeciesList == Herb_Species$epithet5[i])] <- Herb_Species$cover[i]
}
VegAbundanceMatrix[is.na(VegAbundanceMatrix)] <- 0

RDA_herbs <- rda(VegAbundanceMatrix)
pdf(file.path(wd, "fig", "RDA_herbs.pdf"), width = 7, height = 7)
biplot(RDA_herbs, display = c("species", "sites"), type = c("text", "text"),
       col = cols[c(2, 1)])
#ordihull(RDA, groups = "species")
dev.off()

veg_bray_matr <- beta.pair.abund(VegAbundanceMatrix)$beta.bray
plot(hclust(veg_bray_matr))

dbRDA <- dbrda(veg_bray_matr)
#############################################
##### Landcover classification
#############################################
library("raster")
library("cluster")
imgs <- list.files(file.path(wd, "gis", "QGIS", "Sentinel2"), pattern = "\\.tif$", full.names = TRUE)[c(1, 3)]
stk <- stack(imgs)
names(stk) <- paste(substr(names(stk), 13, 15), substr(names(stk), 24, 31), substr(names(stk), 84, 85), sep = "")
Sent2 <- crop(stk, extent(Lapalala))
rm(stk)

idx <- 1:ncell(Sent2)
set.seed(99)
K <- 9
clust <- cluster::clara(na.omit(scale(as.data.frame(Sent2))), k = K)
classified <- Sent2[[1]]
classified[idx] <- clust$clustering
x11()
plot(classified, breaks = seq(1, K, 1), col = cols)
writeRaster(classified, filename = paste(wd, "gis/QGIS/Classification/Classified.tif", sep = ""), format = "GTiff", overwrite = TRUE)
