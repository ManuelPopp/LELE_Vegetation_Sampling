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
arc.check_product()

#############################################
##### Set directories
#############################################
if(Sys.info()['sysname'] == "Windows"){
  # paths Win
  wd <- "D:/Dateien/Studium_KIT/Master_GOEK/Masterarbeit"
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
  rgb(167, 130, 46, alpha = 254.99, maxColorValue = 255), #kit brown
  rgb(252, 229, 0, alpha = 254.99, maxColorValue = 255), #kit yellow
  rgb(25, 161, 224, alpha = 254.99, maxColorValue = 255)# kit cyan
)
colz <- c(
  rgb(0, 150, 130, alpha = 255*0.5, maxColorValue = 255), #kit colour
  rgb(70, 100, 170, alpha = 255*0.5, maxColorValue = 255), #kit blue
  rgb(223, 155, 27, alpha = 255*0.5, maxColorValue = 255), #kit orange
  rgb(140, 182, 60, alpha = 255*0.5, maxColorValue = 255), #kit Mai green
  rgb(162, 34, 35, alpha = 255*0.5, maxColorValue = 255), #kit red
  rgb(163, 16, 124, alpha = 255*0.5, maxColorValue = 255), #kit violet
  rgb(167, 130, 46, alpha = 255*0.5, maxColorValue = 255) #kit brown
)

#############################################
##### Load data
#############################################
# ArcGIS Online data
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
      #"https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/Trees_from_scratch/FeatureServer/0"
      "https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/Tree_mapping_Nov_2021/FeatureServer/0"
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
all_species <- all_species[all_species != "Other"]
all_species[which(all_species == "Euclea linearis")] <- "Euclea crispa"

df_all_spec <- data.frame(Family = rep(NA, length(all_species)),
                          Genus = rep(NA, length(all_species)),
                          Epithet = rep(NA, length(all_species)),
                          Authority = rep(NA, length(all_species)),
                          Source = rep(NA, length(all_species)))

get_infos <- function(Species = NULL){
    #WFO.browse(Species, WFO.data = WFO.data1, accepted.only = TRUE)
    spec_info <- WFO.match(Species, WFO.data = WFO.data1)
    return(c(spec_info$family[1], spec_info$genus[1], spec_info$specificEpithet[1],
             spec_info$scientificNameAuthorship[1],
             spec_info$namePublishedIn[1]))
}

for(i in 1:length(all_species)){
  print(paste("Species:", all_species[i]))
  df_all_spec[i, ] <- get_infos(all_species[i])
  print(paste("Found:", paste(df_all_spec[i, c(2, 3)], collapse = " ")))
}

df_all_spec$Genus <- paste0("openitshape", df_all_spec$Genus, "closeitshape")
df_all_spec$Epithet <- paste0("openitshape", df_all_spec$Epithet, "closeitshape")

df_all_spec[which(df_all_spec$Epithet == "openitshapecrispacloseitshape")[1],
            c("Epithet", "Authority", "Source")] <- c(
              "openitshapecrispacloseitshape subsp. openitshapelineariscloseitshape",
              "(Zeyh. ex Hiern) F.White",
              "Bull. Jard. Bot. Natl. Belg.")

df_all_spec[which(df_all_spec$Genus == "openitshapeBrachylaenacloseitshape" &
                    df_all_spec$Epithet == "openitshapediscolorcloseitshape")[1],
            c("Epithet")] <- c(
              "openitshapediscolorcloseitshape subsp. openitshaperotundatacloseitshape")

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
the World Flora Online Taxonomic Backbone Data}.}\\\\\n \\multicolumn{5}{l}{\\hspace{1em} In: openitshapeApplications in Plant Sciencescloseitshape 8. 9, e11388. issn: 2168-0450. doi: \\href{https://doi.org/10.1002/aps3.11388}{10.1002/aps3.11388}.}\\\\\n \\multicolumn{5}{l}{WFO (2022): World Flora Online. Published on the Internet; \\url{http://www.worldfloraonline.org}. Accessed on: 08 Apr 2022.}\\\\'))
# adjust table size
text_file <- readLines(file.path(db, "tab", "Species_List.tex"))
text_file <- c(text_file[1:6], "\\resizebox{\\textwidth}{!}{", text_file[7:length(text_file)])
text_file <- c(text_file[1:length(text_file)-1], "}", text_file[length(text_file)])

# remove duplicate date
text_file <- gsub("1947 \\[Sep 1947\\]", "(1947)", text_file, perl = TRUE)

# edit caption
text_file <- gsub("caption", "caption[Species list]", text_file)

# adjust umlauts
text_file <- gsub("ä", '\\\\"a', text_file)
text_file <- gsub("ü", '\\\\"u', text_file)
text_file <- gsub("é", "\\\\'e", text_file)

# add bottom line and italics
text_file <- stri_replace_last(text_file, fixed = "\\\\ \\bottomrule", replacement = "\\\\")
text_file <- gsub("openitshape", "\\\\textit\\{", text_file)
text_file <- gsub("closeitshape", "\\}", text_file)

# line breaks at some long citations
text_file <- gsub("Rast. 8", "\\\\\\\\ & & & & Rast. 8", text_file)
text_file <- gsub("Akad. Wiss.", "Akad. Wiss. \\\\\\\\ & & & & ", text_file)
text_file <- gsub("149: 150", "149: \\\\\\\\ & & & & 150", text_file)

# use raw inputencoding
text_file[1] <- "\\UseRawInputEncoding"
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
Veg_Plots <- arc.data2sf(
  arc.select(
    arc.open(
      "https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/service_d48d561dfb0942418475096d297205c3/FeatureServer/0"
    )
  )
)

Veg_Species <- arc.select(
    arc.open(
      "https://services9.arcgis.com/ogIE9wkhkJFrskxj/arcgis/rest/services/service_d48d561dfb0942418475096d297205c3/FeatureServer/1"
    )
  )[, c("objectid", "parentglobalid", "layer", "cover", "epithet5", "iNaturalistObs")]

# functions
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

##############################
# Quadrats complete statistics
PlotGlobalIDs <- unique(Veg_Species$parentglobalid)
Veg_Species$quadratID <- vapply(X = Veg_Species$parentglobalid,
                                 FUN = get_quadrat_ID, FUN.VALUE = character(1))
Veg_Species[is.na(Veg_Species$epithet5), "epithet5"] <- "Unknown_spec."
Veg_Species$epithet5[Veg_Species$epithet5 == "Unknown_spec."] <- Veg_Species$iNaturalistObs[Veg_Species$epithet5 == "Unknown_spec."]

# mean and sd species richness per quadrat
Species_per_quadrat <- Veg_Species %>%
  group_by(quadratID) %>%
  summarize(Spec_count = n_distinct(epithet5))
Species_per_quadrat <- as.data.frame(Species_per_quadrat)
mean(Species_per_quadrat$Spec_count)
sd(Species_per_quadrat$Spec_count)

# abundance matrix
AbMatr <- data.frame(Veg_Species) %>%
  group_by(quadratID, epithet5) %>%
  summarise(mean_cover = mean(cover))

AbundanceMatrix <- dcast(data = AbMatr, quadratID ~ epithet5, value.var = "mean_cover")
AbundanceMatrix[is.na(AbundanceMatrix)] <- 0
row.names(AbundanceMatrix) <- AbundanceMatrix$quadratID
AbundanceMatrix <- AbundanceMatrix[, -which(names(AbundanceMatrix) == "quadratID")]
which(rowSums(AbundanceMatrix) == 0)
which(colSums(AbundanceMatrix) == 0)

# PCA on Hellinger transformed abundance data (just for fun)
AbundanceMatrix_hell <- decostand(AbundanceMatrix, "hellinger")
PCA <- rda(AbundanceMatrix_hell)
biplot(PCA, scaling = "species", display = c("sites", "species"),
       type = c("text", "text"), col = cols)

require("factoextra")
pca <- prcomp(x = AbundanceMatrix_hell, scale = TRUE)
fviz_eig(pca)
fviz_pca_ind(pca,
             col.ind = "cos2",
             repel = TRUE)
fviz_pca_var(pca,
             col.var = "contrib")

fviz_pca_biplot(pca,
                # Individuals
                geom.ind = "point",
                fill.ind = substr(row.names(AbundanceMatrix), 1, 2), col.ind = "black",
                pointshape = 21,
                pointsize = 2,
                palette = cols,
                addEllipses = TRUE,
                # Variables
                alpha.var = "contrib",
                col.var = "contrib",
                gradient.cols = c("lightblue", "black"),#cols[c(1, 2, 3)],
                legend.title = list(fill = "Type", color = "Contr.",
                                    alpha = "Contr."))

# Herb layer
Herb_Species <- Veg_Species[Veg_Species$layer == "grass",]
PlotGlobalIDs <- unique(Herb_Species$parentglobalid)

Herb_Species$quadratID <- vapply(X = Herb_Species$parentglobalid,
                                 FUN = get_quadrat_ID, FUN.VALUE = character(1))

#
Plot_infos <- data.frame(Tree_cover = rep(NA, nrow(Herb_Species)), Shrub_cover = NA, Rock_cover = NA)
Herb_Species <- Herb_Species[nchar(Herb_Species$quadratID) < 10,]

SpeciesList <- unique(Herb_Species$epithet5)
QuadratList <- unique(Herb_Species$quadratID)

AbMatrHerb <- data.frame(Herb_Species) %>%
  group_by(quadratID, epithet5) %>%
  summarise(mean_cover = mean(cover))

AbundanceMatrixHerb <- dcast(data = AbMatrHerb, quadratID ~ epithet5, value.var = "mean_cover")
AbundanceMatrixHerb[is.na(AbundanceMatrixHerb)] <- 0
row.names(AbundanceMatrixHerb) <- AbundanceMatrixHerb$quadratID
AbundanceMatrixHerb <- AbundanceMatrixHerb[, -which(names(AbundanceMatrixHerb) == "quadratID")]
