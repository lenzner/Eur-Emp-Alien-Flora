#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Analysis script 01                                                      #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#


#================================================================================#
# Packages ----
library("rgdal")
library("vegan")
library("betapart")
library("fuzzySim")
library("zetadiv")
library("dplyr")
library("sf")
library("lwgeom")
library("rgeos")
library("geosphere")
library("tidyverse")
library("tictoc")
library("janitor")
library("viridis")
library("netrankr")
library("tidygraph")
library("assertthat")
library("purrr")
library("ggplot2")
library("GGally")
library("ggraph")
library("ggmap")
library("ggpmisc")
library("igraph")
library("gridExtra")
library("ggrepel")
library("mapproj")


#================================================================================#
# Project paths ----
repo <- paste0('~') # Set path to the downloaded directory
code.path <- paste0(repo, '/R_Scripts') # path for directory with code
data.path <- paste0(repo, '/Datasets') # path for directory with analysis datasets
fig.path <- paste0(repo, '/Results/Figures') # path for directory with figures
shape.path <- paste0(repo, "/Datasets/Shapefile") # path for directory with region shapefile
res.path <- paste0(repo, '/Results') # path for directory with analysis datasets
fun.path <- paste0(repo, '/Functions') # path for directory with custom functions

#================================================================================#
# Datasets ----

### Shapefile
shape <- readOGR(dsn = shape.path, layer = "Shapefile")

# GloNAF species data ----
glonaf <- read.csv(paste0(data.path,"/GloNAF.csv"), sep = ",", header = T, stringsAsFactors = F)
head(glonaf)

# Predictors dataset
pred <- read.csv(paste0(data.path,"/Predictors.csv"), sep = ",", header = T, stringsAsFactors = F)
head(pred)

# Empire data
empires <- read.csv(paste0(data.path,"/Empires.csv"), sep = ",", header = T, stringsAsFactors = F)
head(empires)

# Regions object to run null model code
regions <- read.csv(paste0(data.path,"/Regions.csv"), sep = ",", header = T, stringsAsFactors = F)
head(regions)


#================================================================================#
# Null model ----

### Load Null Model functions from source ----
source(paste0(fun.path,"/Null_Model_Functions.R")) # Loads functions random.empire() and random.empire.seq()

### Function description

#### random.empire()
##### The null model selects random regions globally. The regions selection is constrained by the number of regions in a geospatial region in the original empire and by the number of mainland and island regions in the geospatial region. Besically the null model builds a random empire with the same spatial structure (i.e. number of mainland and island regions per geospatial region as in the observed empire)

#### ranodm.empire.seq() --> Null model used in the manuscript
##### The regions selection is constrained by the number of regions in a geospatial region in the original empire and by the number of mainland and island regions in the geospatial region. The first region per geospatial region is selected randomly, all subsequent regions will be selected sequentially based on the shortest distance to the previous regions. 


### Calculate centroid coordinates for each region
trueCentroids = gCentroid(shape,byid=TRUE) # calculate mass centroids from polygons
coord <- as.data.frame(cbind(shape@data[,c("OBJIDsic")],coordinates(trueCentroids)))
colnames(coord)[1] <- "OBJIDsic"
colnames(coord)[2] <- "LON"
colnames(coord)[3] <- "LAT"
coord <- coord[match(pred$OBJIDsic, coord$OBJIDsic),] # Order coord similar to pred

### Calculate pairwise (geographic) distance between regions
m <- distm(coord[,c("LON","LAT")], coord[,c("LON","LAT")], fun = distGeo) # create a distance matrix
diag(m) <- NA # replace the diagonal with NA
colnames(m) <- coord$OBJIDsic # make column names for the distance matrix
rownames(m) <- coord$OBJIDsic # make column names for the distance matrix


# Build datasets for each empire
emp.GBR <- empires %>%
  filter(colonizer  == "Great_Britain" | OBJIDsic == 162) %>%
  select(OBJIDsic, Region_unique, MainIsl, colonizer) %>%
  distinct()

emp.PRT <- empires %>%
  filter(colonizer  == "Portugal" | OBJIDsic == 1361) %>%
  select(OBJIDsic, Region_unique, MainIsl, colonizer) %>%
  distinct()

emp.ESP <- empires %>%
  filter(colonizer  == "Spain" | OBJIDsic == 1360) %>%
  select(OBJIDsic, Region_unique, MainIsl, colonizer) %>%
  distinct()

emp.NED <- empires %>%
  filter(colonizer  == "Netherlands" | OBJIDsic == 1029) %>%
  select(OBJIDsic, Region_unique, MainIsl, colonizer) %>%
  distinct()



# Run the null model
iter = 10 # decide on a number of random empires to calculate

#### Fully random empire
emp.ran.GBR <- random.empire("Great_Britain", iter)
emp.ran.PRT <- random.empire("Portugal", iter)
emp.ran.ESP <- random.empire("Spain", iter)
emp.ran.NED <- random.empire("Netherlands", iter)

#### Random empire with sequential colonization structure --> Null model used in the manuscript
emp.ran.seq.GBR <- random.empire.seq("Great_Britain", iter)
emp.ran.seq.PRT <- random.empire.seq("Portugal", iter)
emp.ran.seq.ESP <- random.empire.seq("Spain", iter)
emp.ran.seq.NED <- random.empire.seq("Netherlands", iter)


#================================================================================#
# Build presence-absence matrices for all GloNAF species and for each empire ----
glonaf.pa <- splist2presabs(glonaf, sites.col = "OBJIDsic", sp.col = "Species", keep.n = F)
rownames(glonaf.pa) <- glonaf.pa$OBJIDsic # set region id (= OBJIDsic) as rownames
glonaf.pa <- glonaf.pa[,-1] # remove column with region id (= OBJIDsic)

# Build presence absence matrix for each empire ----
pa.GBR <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.GBR$OBJIDsic),]
pa.PRT <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.PRT$OBJIDsic),]
pa.ESP <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.ESP$OBJIDsic),]
pa.NED <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.NED$OBJIDsic),]



#================================================================================#
# Zeta diversity ----
### Observed zeta diversity ----
orders <- 40 # Define how many zeta order should be calculated

#### British Empire
pa.obs <- pa.GBR
ptm <- proc.time()
zeta.dec.ex.obs.GBR <- Zeta.decline.ex(pa.obs, orders = 1:orders)
proc.time() - ptm

date <- Sys.Date()
save(zeta.dec.ex.obs.GBR, file=paste0(res.path, "/", date, "_zeta_obs_GBR_order40.Rdata"))

#### Spanish Empire
pa.obs <- pa.ESP
ptm <- proc.time()
zeta.dec.ex.obs.ESP <- Zeta.decline.ex(pa.obs, orders = 1:orders)
proc.time() - ptm

date <- Sys.Date()
save(zeta.dec.ex.obs.ESP, file=paste0(res.path, "/", date, "_zeta_obs_ESP_order40.Rdata"))

#### Portuguese Empire
pa.obs <- pa.PRT
ptm <- proc.time()
zeta.dec.ex.obs.PRT <- Zeta.decline.ex(pa.obs, orders = 1:orders)
proc.time() - ptm

date <- Sys.Date()
save(zeta.dec.ex.obs.PRT, file=paste0(res.path, "/", date, "_zeta_obs_PRT_order40.Rdata"))

#### Netherlands
pa.obs <- pa.NED
ptm <- proc.time()
zeta.dec.ex.obs.NED <- Zeta.decline.ex(pa.obs, orders = 1:orders)
proc.time() - ptm

date <- Sys.Date()
save(zeta.dec.ex.obs.NED, file=paste0(res.path, "/", date, "_zeta_obs_NED_order40.Rdata"))



#================================================================================#
### Random zeta diversity ----
orders = 40 # Define how many zeta order should be calculated
iter = 10 # Define number of iterations of how many random empires should be used


#### British Empire
emp.ran.seq <- emp.ran.seq.GBR 
zeta.ran.seq.GBR <- list()

for(i in 1:iter){
  
  ptm <- proc.time()
  pa <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.ran.seq[[i]]$OBJIDsic),]
  zeta.random <- Zeta.decline.ex(pa, orders = 1:orders)
  
  zeta.ran.seq.GBR[[i]] <- zeta.random
  end <- proc.time() - ptm
  
  print(i)
  print(end)
  
}


#### Spanish Empire
emp.ran.seq <- emp.ran.seq.ESP
zeta.ran.seq.ESP <- list()

for(i in 1:iter){
  
  ptm <- proc.time()
  pa <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.ran.seq[[i]]$OBJIDsic),]
  zeta.random <- Zeta.decline.ex(pa, orders = 1:orders)
  
  zeta.ran.seq.ESP[[i]] <- zeta.random
  end <- proc.time() - ptm
  
  print(i)
  print(end)
  
}


#### Portuguese Empire
emp.ran.seq <- emp.ran.seq.PRT
zeta.ran.seq.PRT <- list()

for(i in 1:iter){
  
  ptm <- proc.time()
  pa <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.ran.seq[[i]]$OBJIDsic),]
  zeta.random <- Zeta.decline.ex(pa, orders = 1:orders)
  
  zeta.ran.seq.PRT[[i]] <- zeta.random
  end <- proc.time() - ptm
  
  print(i)
  print(end)
  
}


#### Dutch Empire
emp.ran.seq <- emp.ran.seq.NED
zeta.ran.seq.NED <- list()

for(i in 1:iter){
  
  ptm <- proc.time()
  pa <- glonaf.pa[which(rownames(glonaf.pa) %in% emp.ran.seq[[i]]$OBJIDsic),]
  zeta.random <- Zeta.decline.ex(pa, orders = 1:orders)
  
  zeta.ran.seq.NED[[i]] <- zeta.random
  end <- proc.time() - ptm
  
  print(i)
  print(end)
  
}


date <- Sys.Date()
save(zeta.ran.seq.GBR, file=paste0(res.path, "/", date, "_zeta_ran_GBR_ran.", iter, "_order.", orders, ".Rdata"))

save(zeta.ran.seq.ESP, file=paste0(res.path, "/", date, "_zeta_ran_ESP_ran.", iter, "_order.", orders, ".Rdata"))

save(zeta.ran.seq.PRT, file=paste0(res.path, "/", date, "_zeta_ran_PRT_ran.", iter, "_order.", orders, ".Rdata"))

save(zeta.ran.seq.NED, file=paste0(res.path, "/", date, "_zeta_ran_NED_ran.", iter, "_order.", orders, ".Rdata"))




#================================================================================#
### Export zeta ratio results for the observed empires ----
zeta.rat.obs <- cbind(zeta.dec.ex.obs.GBR$ratio,zeta.dec.ex.obs.ESP$ratio,zeta.dec.ex.obs.PRT$ratio,zeta.dec.ex.obs.NED$ratio)

date <- Sys.Date()
write.table(zeta.rat.obs, file=paste0(res.path, "/", date, "_zeta_ratio_obs.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")

### Export zeta and zeta ratio results for the random empires

#### British Empire
##### Zeta order
zeta.ran.GBR <- zeta.ran.seq.GBR[[1]]$zeta.val[1:5]
for(i in 2:iter){
  x <- rbind(x,zeta.ran.seq.GBR[[i]]$zeta.val[1:5])
}
apply(zeta.ran.GBR, 2,mean)
apply(zeta.ran.GBR, 2,sd)


#### Zeta ratio
zeta.rat.ran.GBR <- zeta.ran.seq.GBR[[1]]$ratio
for(j in 2:iter){
  y <- rbind(y,zeta.ran.seq.GBR[[j]]$ratio)
}
zeta.rat.ran.GBR.m <- round(apply(zeta.rat.ran.GBR, 2,mean),2)
zeta.rat.ran.GBR.s <- round(apply(zeta.rat.ran.GBR, 2,sd),2)
zeta.rat.ran.GBR <- cbind(zeta.rat.ran.GBR.m,zeta.rat.ran.GBR.s, paste0(zeta.rat.ran.GBR.m," (", '\u00B1 ', zeta.rat.ran.GBR.s, ")"))


#### Spanish Empire
##### Zeta order
zeta.ran.ESP <- zeta.ran.seq.ESP[[1]]$zeta.val[1:5]
for(i in 2:iter){
  x <- rbind(x,zeta.ran.seq.ESP[[i]]$zeta.val[1:5])
}
apply(zeta.ran.ESP, 2,mean)
apply(zeta.ran.ESP, 2,sd)


#### Zeta ratio
zeta.rat.ran.ESP <- zeta.ran.seq.ESP[[1]]$ratio
for(j in 2:iter){
  y <- rbind(y,zeta.ran.seq.ESP[[j]]$ratio)
}
zeta.rat.ran.ESP.m <- round(apply(zeta.rat.ran.ESP, 2,mean),2)
zeta.rat.ran.ESP.s <- round(apply(zeta.rat.ran.ESP, 2,sd),2)
zeta.rat.ran.ESP <- cbind(zeta.rat.ran.ESP.m,zeta.rat.ran.ESP.s, paste0(zeta.rat.ran.ESP.m," (", '\u00B1 ', zeta.rat.ran.ESP.s, ")"))


#### Portuguese Empire
##### Zeta order
zeta.ran.PRT <- zeta.ran.seq.PRT[[1]]$zeta.val[1:5]
for(i in 2:iter){
  x <- rbind(x,zeta.ran.seq.PRT[[i]]$zeta.val[1:5])
}
apply(zeta.ran.PRT, 2,mean)
apply(zeta.ran.PRT, 2,sd)


#### Zeta ratio
zeta.rat.ran.PRT <- zeta.ran.seq.PRT[[1]]$ratio
for(j in 2:iter){
  y <- rbind(y,zeta.ran.seq.PRT[[j]]$ratio)
}
zeta.rat.ran.PRT.m <- round(apply(zeta.rat.ran.PRT, 2,mean),2)
zeta.rat.ran.PRT.s <- round(apply(zeta.rat.ran.PRT, 2,sd),2)
zeta.rat.ran.PRT <- cbind(zeta.rat.ran.PRT.m,zeta.rat.ran.PRT.s, paste0(zeta.rat.ran.PRT.m," (", '\u00B1 ', zeta.rat.ran.PRT.s, ")"))


#### Dutch Empire
##### Zeta order
zeta.ran.NED <- zeta.ran.seq.NED[[1]]$zeta.val[1:5]
for(i in 2:iter){
  x <- rbind(x,zeta.ran.seq.NED[[i]]$zeta.val[1:5])
}
apply(zeta.ran.NED, 2,mean)
apply(zeta.ran.NED, 2,sd)


#### Zeta ratio
zeta.rat.ran.NED <- zeta.ran.seq.NED[[1]]$ratio
for(j in 2:iter){
  y <- rbind(y,zeta.ran.seq.NED[[j]]$ratio)
}
zeta.rat.ran.NED.m <- round(apply(zeta.rat.ran.NED, 2,mean),2)
zeta.rat.ran.NED.s <- round(apply(zeta.rat.ran.NED, 2,sd),2)
zeta.rat.ran.NED <- cbind(zeta.rat.ran.NED.m,zeta.rat.ran.NED.s, paste0(zeta.rat.ran.NED.m," (", '\u00B1 ', zeta.rat.ran.NED.s, ")"))




date <- Sys.Date()
write.table(zeta.ran.GBR, file=paste0(res.path, "/", date,"_zeta_value_GBR.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")
write.table(zeta.rat.ran.GBR, file=paste0(res.path, "/", date,"_zeta_ratio_GBR.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")

write.table(zeta.ran.ESP, file=paste0(res.path, "/", date,"_zeta_value_ESP.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")
write.table(zeta.rat.ran.ESP, file=paste0(res.path, "/", date,"_zeta_ratio_ESP.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")

write.table(zeta.ran.PRT, file=paste0(res.path, "/", date,"_zeta_value_PRT.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")
write.table(zeta.rat.ran.PRT, file=paste0(res.path, "/", date,"_zeta_ratio_PRT.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")

write.table(zeta.ran.NED, file=paste0(res.path, "/", date,"_zeta_value_NED.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")
write.table(zeta.rat.ran.NED, file=paste0(res.path, "/", date,"_zeta_ratio_NED.csv"), sep =";", row.names = F, fileEncoding = "UTF-8")


#================================================================================#
#### Visualize zeta observed and random ----
##### Load Null Model functions from source ----
source(paste0(fun.path,"/Zeta_value_ratio_visualization.R")) # Loads function to visualize the zeta values and zeta ratios for each empire for their observed and random runs

Plot.zeta.emp("Great_Britain")
Plot.zeta.emp("Spain")
Plot.zeta.emp("Portugal")
Plot.zeta.emp("Netherlands")


#================================================================================#
#================================================================================#
#================================================================================#
# Zeta dissimilarity modelling ----
### Prepare predictor dataset for each empire

#### British Empire
pred.calc.GBR <- pred %>%
  select(OBJIDsic, LON, LAT, island, area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.GB, crop.area.prop) %>%
  mutate_if(is.integer, as.numeric) # change all integer to numeric

pred.calc.GBR <- pred.calc.GBR[match(rownames(pa.GBR), pred.calc.GBR$OBJIDsic),]

pred.calc.GBR2 <- pred.calc.GBR %>%
  select(area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.GB, crop.area.prop)

pred.calc.GBR2.s <- pred.calc.GBR2 %>%
  mutate(
    area = scale(log(area)), 
    dem.sd = scale(log(dem.sd+0.1)),  # removed to answer reviewer comment
    temp.mean = scale(temp.mean), 
    aridity.mean = scale(log(aridity.mean+0.1)), 
    pop.dens.2010.mean = scale(log(pop.dens.2010.mean+0.1)), 
    gdp.pc.2011.mean = scale(log(gdp.pc.2011.mean+0.1)), 
    time.emp.GB = scale(time.emp.GB+0.1), 
    crop.area.prop = scale(log(crop.area.prop+0.1)))

# remove attributes from dataframe
pred.calc.GBR2.s[] = lapply(pred.calc.GBR2.s, function(x) { attributes(x) <- NULL; x })
pred.calc.GBR2[] = lapply(pred.calc.GBR2, function(x) { attributes(x) <- NULL; x })
pred.calc.GBR[] = lapply(pred.calc.GBR, function(x) { attributes(x) <- NULL; x })


#### Spanish Empire
pred.calc.ESP <- pred %>%
  select(OBJIDsic, LON, LAT, island, area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.ESP, crop.area.prop) %>%
  mutate_if(is.integer, as.numeric) # change all integer to numeric
pred.calc.ESP <- pred.calc.ESP[match(rownames(pa.ESP), pred.calc.ESP$OBJIDsic),]

pred.calc.ESP2 <- pred.calc.ESP %>%
  select(area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.ESP, crop.area.prop)

pred.calc.ESP2.s <- pred.calc.ESP2 %>%
  mutate(
    area = scale(log(area)), 
    dem.sd = scale(log(dem.sd+0.1)), 
    temp.mean = scale(temp.mean), 
    aridity.mean = scale(log(aridity.mean+0.1)), 
    pop.dens.2010.mean = scale(log(pop.dens.2010.mean+0.1)), 
    gdp.pc.2011.mean = scale(log(gdp.pc.2011.mean+0.1)), 
    time.emp.ESP = scale(time.emp.ESP+0.1), 
    crop.area.prop = scale(log(crop.area.prop+0.1)))

# remove attributes from dataframe
pred.calc.ESP2.s[] = lapply(pred.calc.ESP2.s, function(x) { attributes(x) <- NULL; x })
pred.calc.ESP2[] = lapply(pred.calc.ESP2, function(x) { attributes(x) <- NULL; x })
pred.calc.ESP[] = lapply(pred.calc.ESP, function(x) { attributes(x) <- NULL; x })


#### Portuguese Empire
pred.calc.PRT <- pred %>%
  select(OBJIDsic, LON, LAT, island, area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.PRT, crop.area.prop) %>%
  mutate_if(is.integer, as.numeric) # change all integer to numeric

pred.calc.PRT <- pred.calc.PRT[match(rownames(pa.PRT), pred.calc.PRT$OBJIDsic),]

pred.calc.PRT2 <- pred.calc.PRT %>%
  select(area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.PRT, crop.area.prop)

pred.calc.PRT2.s <- pred.calc.PRT2 %>%
  mutate(
    area = scale(log(area)), 
    dem.sd = scale(log(dem.sd+0.1)), 
    temp.mean = scale(temp.mean), 
    aridity.mean = scale(log(aridity.mean+0.1)), 
    pop.dens.2010.mean = scale(log(pop.dens.2010.mean+0.1)), 
    gdp.pc.2011.mean = scale(log(gdp.pc.2011.mean+0.1)), 
    time.emp.PRT = scale(time.emp.PRT+0.1), 
    crop.area.prop = scale(log(crop.area.prop+0.1)))

# remove attributes from dataframe
pred.calc.PRT2.s[] = lapply(pred.calc.PRT2.s, function(x) { attributes(x) <- NULL; x })
pred.calc.PRT2[] = lapply(pred.calc.PRT2, function(x) { attributes(x) <- NULL; x })
pred.calc.PRT[] = lapply(pred.calc.PRT, function(x) { attributes(x) <- NULL; x })


#### Dutch Empire
pred.calc.NED <- pred %>%
  select(OBJIDsic, LON, LAT, island, area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.NL, crop.area.prop) %>%
  mutate_if(is.integer, as.numeric) # change all integer to numeric

pred.calc.NED <- pred.calc.NED[match(rownames(pa.NED), pred.calc.NED$OBJIDsic),]

pred.calc.NED2 <- pred.calc.NED %>%
  select(area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, time.emp.NL, crop.area.prop)

pred.calc.NED2.s <- pred.calc.NED2 %>%
  mutate(
    area = scale(log(area)), 
    dem.sd = scale(log(dem.sd+0.1)), 
    temp.mean = scale(temp.mean), 
    aridity.mean = scale(log(aridity.mean+0.1)), 
    pop.dens.2010.mean = scale(log(pop.dens.2010.mean+0.1)), 
    gdp.pc.2011.mean = scale(log(gdp.pc.2011.mean+0.1)), 
    time.emp.NL = scale(time.emp.NL+0.1), 
    crop.area.prop = scale(log(crop.area.prop+0.1)))

# remove attributes from dataframe
pred.calc.NED2.s[] = lapply(pred.calc.NED2.s, function(x) { attributes(x) <- NULL; x })
pred.calc.NED2[] = lapply(pred.calc.NED2, function(x) { attributes(x) <- NULL; x })
pred.calc.NED[] = lapply(pred.calc.NED, function(x) { attributes(x) <- NULL; x })



### Check predictor variable correlations for each empire

pairs.GBR <- pred.calc.GBR2.s %>%
  rename(Area = area,
         Habitat.heterogeneity = dem.sd,
         Mean.annual.temperature = temp.mean,
         Aridity.index = aridity.mean,
         Population.density = pop.dens.2010.mean,
         GDPpc = gdp.pc.2011.mean,
         Occupation.time = time.emp.GB,
         Agricultural.land = crop.area.prop)

pairs.ESP <- pred.calc.ESP2.s %>%
  rename(Area = area,
         Habitat.heterogeneity = dem.sd,
         Mean.annual.temperature = temp.mean,
         Aridity.index = aridity.mean,
         Population.density = pop.dens.2010.mean,
         GDPpc = gdp.pc.2011.mean,
         Occupation.time = time.emp.ESP,
         Agricultural.land = crop.area.prop)

pairs.PRT <- pred.calc.PRT2.s %>%
  rename(Area = area,
         Habitat.heterogeneity = dem.sd,
         Mean.annual.temperature = temp.mean,
         Aridity.index = aridity.mean,
         Population.density = pop.dens.2010.mean,
         GDPpc = gdp.pc.2011.mean,
         Occupation.time = time.emp.PRT,
         Agricultural.land = crop.area.prop)


pairs.NED <- pred.calc.NED2.s %>%
  rename(Area = area,
         Habitat.heterogeneity = dem.sd,
         Mean.annual.temperature = temp.mean,
         Aridity.index = aridity.mean,
         Population.density = pop.dens.2010.mean,
         GDPpc = gdp.pc.2011.mean,
         Occupation.time = time.emp.NL,
         Agricultural.land = crop.area.prop)


p.pairs.GBR = ggpairs(pairs.GBR)
p.pairs.ESP = ggpairs(pairs.ESP)
p.pairs.PRT = ggpairs(pairs.PRT)
p.pairs.NED = ggpairs(pairs.NED)

date <- Sys.Date()
png(file=paste0(fig.path, "/", date, "_Predictor_correlation_GBR.png"), width=1500, height=700)
  p.pairs.GBR
dev.off()

png(file=paste0(fig.path, "/", date, "_Predictor_correlation_ESP.png"), width=1500, height=700)
  p.pairs.ESP
dev.off()

png(file=paste0(fig.path, "/", date, "_Predictor_correlation_PRT.png"), width=1500, height=700)
  p.pairs.PRT
dev.off()

png(file=paste0(fig.path, "/", date, "_Predictor_correlation_NED.png"), width=1500, height=700)
  p.pairs.NED
dev.off()



#================================================================================#
# MS-GDM Model ----
### British Empire

tic()
mod.GBR.ispl.rep <- list()

for(j in 1:6){ # Run 6 replicates
  
  print(j)
  
  mod.GBR.ispl <- list()
  
  for(i in 2:5){ # Run for zeta orders 2 to 5
    
    model <- Zeta.msgdm(pa.GBR, pred.calc.GBR2.s, xy = pred.calc.GBR[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "log"))
    mod.GBR.ispl[[i-1]] <- model
    print(i)
  }
  
  mod.GBR.ispl.rep[[j]] <- mod.GBR.ispl
  
  
}
toc() 

###########################
# Run model without habitat heterogeneity to test reviewer comment
pred.calc.GBR2.s <- pred.calc.GBR2.s %>% select(!dem.sd)


tic()

mod.GBR.ispl.rep <- list()

for(j in 1:6){ # Run 6 replicates
  
  print(j)
  
  mod.GBR.ispl <- list()
  
  for(i in 2:5){ # Run for zeta orders 2 to 5
    
    model <- Zeta.msgdm(pa.GBR, pred.calc.GBR2.s, xy = pred.calc.GBR[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "logit"))
    mod.GBR.ispl[[i-1]] <- model
    print(i)
    
  }
  
  mod.GBR.ispl.rep[[j]] <- mod.GBR.ispl
  
}

toc() 

### Spanish Empire

tic()

mod.ESP.ispl.rep <- list()

for(j in 1:6){ # Run 6 replicates
  
  print(j)
  
  mod.ESP.ispl <- list()
  
  for(i in 2:5){ # Run for zeta orders 2 to 5
    
    model <- Zeta.msgdm(pa.ESP, pred.calc.ESP2.s, xy = pred.calc.ESP[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "logit"))
    mod.ESP.ispl[[i-1]] <- model
    print(i)
    
  }
  
  mod.ESP.ispl.rep[[j]] <- mod.ESP.ispl
  
}

toc() 


### Portuguese Empire

tic()

mod.PRT.ispl.rep <- list()

for(j in 1:6){ # Run 6 replicates
  
  print(j)
  
  mod.PRT.ispl <- list()
  
  for(i in 2:5){ # Run for zeta orders 2 to 5
    
    model <- Zeta.msgdm(pa.PRT, pred.calc.PRT2.s, xy = pred.calc.PRT[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "logit"))
    mod.PRT.ispl[[i-1]] <- model
    print(i)
    
  }
  
  mod.PRT.ispl.rep[[j]] <- mod.PRT.ispl
  
}

toc() 


### Dutch Empire
tic()

mod.NED.ispl.rep <- list()

for(j in 1:6){
  
  print(j)
  
  mod.NED.ispl <- list()
  
  for(i in 2:5){ # Run for zeta orders 2 to 5
    
    model <- Zeta.msgdm(pa.NED, pred.calc.NED2.s, xy = pred.calc.NED[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "logit"))
    mod.NED.ispl[[i-1]] <- model
    print(i)
    
  }
  
  mod.NED.ispl.rep[[j]] <- mod.NED.ispl
  
}

toc() 

###########################
# Run model without habitat heterogeneity to test reviewer comment
pred.calc.NED2.s <- pred.calc.NED2.s %>% select(!dem.sd)

tic()

mod.NED.ispl.rep <- list()

for(j in 1:6){
  
  print(j)
  
  mod.NED.ispl <- list()
  for(i in 2:5){
    model <- Zeta.msgdm.test(pa.NED, pred.calc.NED2.s, xy = pred.calc.NED[,c("LAT","LON")], sam = 8000, order = i, reg.type = "ispline", rescale = T, rescale.pred = T, normalize = "Simpson", family = binomial(link = "logit"))
    mod.NED.ispl[[i-1]] <- model
    print(i)
  }
  
  mod.NED.ispl.rep[[j]] <- mod.NED.ispl
  
  
}
toc() 

###########################



date <- Sys.Date()
save(mod.GBR.ispl.rep, file=paste0(res.path, "/", date, "_MSGDM_ispline_GBR_8000.rep.6.Rdata"))

save(mod.ESP.ispl.rep, file=paste0(res.path, "/", date, "_MSGDM_ispline_ESP_8000.rep.6.Rdata"))

save(mod.PRT.ispl.rep, file=paste0(res.path, "/", date, "_MSGDM_ispline_PRT_8000.rep.6.Rdata"))

save(mod.NED.ispl.rep, file=paste0(res.path, "/", date, "_MSGDM_ispline_NED_8000.rep.6.Rdata"))

#####SAVE MSGDM MODEL RUNS!!!!!







# variance partitioning
reps <- 6
zetas <- 4

### GBR
var.GBR <- data.frame(matrix(NA, nrow = reps, ncol = zetas))
colnames(var.GBR) <- paste0("zeta", 1:4)

for(i in 1:reps){
  
  for(j in 1:zetas){
    
    var.GBR[i, j] <- Zeta.varpart(mod.GBR.ispl.rep[[i]][[j]])[1,]
    
  }
  
}
var.GBR

colMeans(var.GBR)

### ESP
var.ESP <- data.frame(matrix(NA, nrow = reps, ncol = zetas))
colnames(var.ESP) <- paste0("zeta", 1:4)

for(i in 1:reps){
  
  for(j in 1:zetas){
    
    var.ESP[i, j] <- Zeta.varpart(mod.ESP.ispl.rep[[i]][[j]])[1,]
    
  }
  
}
var.ESP

### PRT
var.PRT <- data.frame(matrix(NA, nrow = reps, ncol = zetas))
colnames(var.PRT) <- paste0("zeta", 1:4)

for(i in 1:reps){
  
  for(j in 1:zetas){
    
    var.PRT[i, j] <- Zeta.varpart(mod.PRT.ispl.rep[[i]][[j]])[1,]
    
  }
  
}
var.PRT

### NED
var.NED <- data.frame(matrix(NA, nrow = reps, ncol = zetas))
colnames(var.NED) <- paste0("zeta", 1:4)

for(i in 1:reps){
  
  for(j in 1:zetas){
    
    var.NED[i, j] <- Zeta.varpart(mod.NED.ispl.rep[[i]][[j]])[1,]
    
  }
  
}
var.NED

colMeans(var.NED)


# Summary of model runs
var.GBR.sum <- var.GBR %>%
  summarise(mean.GBR = colMeans(var.GBR), sd.GBR = apply(var.GBR,2,sd))
var.ESP.sum <- var.ESP %>%
  summarise(mean.ESP = colMeans(var.ESP), sd.ESP = apply(var.ESP,2,sd))
var.PRT.sum <- var.PRT %>%
  summarise(mean.PRT = colMeans(var.PRT), sd.PRT = apply(var.PRT,2,sd))
var.NED.sum <- var.NED %>%
  summarise(mean.NED = colMeans(var.NED), sd.NED = apply(var.NED,2,sd))

var.sum <- cbind(var.GBR.sum, var.ESP.sum, var.PRT.sum, var.NED.sum)


#x11()
date <- Sys.Date()
name <- paste0(res.path,"/Figures/",date,"_Explained_variance.png")
png(name, width = 842, height = 842,res=72)
  par(mar = c(6, 6, 6, 6))  
  plot(var.sum[,"mean.GBR"] ~ c(2:5), ylim = c(0,1), xaxt = "n", xlab = "Zeta order", ylab = "Explained variance", col = "#0066CC", pch = 19, type = "b", las = 2, cex = 2, cex.lab = 2, cex.axis = 1.5)
  axis(1, labels = c(2:5), at = 2:5, cex.axis = 1.5)
  points(var.sum[,"mean.ESP"] ~ c(2:5), col = "#CC3300", pch = 15, type = "b", cex = 2)
  points(var.sum[,"mean.PRT"] ~ c(2:5), col = "#660066", pch = 16, type = "b", cex = 2)
  points(var.sum[,"mean.NED"] ~ c(2:5), col = "#FF9933", pch = 17, type = "b", cex = 2)
  legend("topright", c("British Empire", "Spanish Empire", "Portuguese Empire", "Dutch Empire"), lty = 1, pch = c(19,15:17), col = c("#0066CC", "#CC3300", "#660066", "#FF9933"), cex = 2)
dev.off()
  
#

mod.ESP.ispl.rep[[1]]

#






#####
## Most widespread species per empire
spp.per.reg.GBR <- data.frame(colSums(pa.GBR))
spp.per.reg.GBR <- cbind(rownames(spp.per.reg.GBR), spp.per.reg.GBR)
colnames(spp.per.reg.GBR) <- c("species", "n_regions")
spp.per.reg.GBR <- spp.per.reg.GBR[order(spp.per.reg.GBR$n_regions, decreasing = T),]

spp.per.reg.ESP <- data.frame(colSums(pa.ESP))
spp.per.reg.ESP <- cbind(rownames(spp.per.reg.ESP), spp.per.reg.ESP)
colnames(spp.per.reg.ESP) <- c("species", "n_regions")
spp.per.reg.ESP <- spp.per.reg.ESP[order(spp.per.reg.ESP$n_regions, decreasing = T),]

spp.per.reg.PRT <- data.frame(colSums(pa.PRT))
spp.per.reg.PRT <- cbind(rownames(spp.per.reg.PRT), spp.per.reg.PRT)
colnames(spp.per.reg.PRT) <- c("species", "n_regions")
spp.per.reg.PRT <- spp.per.reg.PRT[order(spp.per.reg.PRT$n_regions, decreasing = T),]

spp.per.reg.NED <- data.frame(colSums(pa.NED))
spp.per.reg.NED <- cbind(rownames(spp.per.reg.NED), spp.per.reg.NED)
colnames(spp.per.reg.NED) <- c("species", "n_regions")
spp.per.reg.NED <- spp.per.reg.NED[order(spp.per.reg.NED$n_regions, decreasing = T),]


spp.top.20 <- cbind(spp.per.reg.GBR[1:20,], spp.per.reg.ESP[1:20,], spp.per.reg.PRT[1:20,], spp.per.reg.NED[1:20,])
rownames(spp.top.20) <- NULL
colnames(spp.top.20) <- c("GBR.spp", "GBR.reg", "ESP.spp", "ESP.reg", "PRT.spp", "PRT.reg", "NED.spp", "NED.reg")

## date <- Sys.Date()
## write.table(spp.top.20, file=paste0(res.path, "/", date, "_Top20_widespread_spp_per_region.csv"), row.names = F, sep = ";")
#
















############################################################################################
## Richness in dependence of residence time ----

head(glonaf)
head(pred)

spp.nat.SR <- data.frame(table(glonaf[,"OBJIDsic"]))
colnames(spp.nat.SR) <- c("OBJIDsic", "SR.nat")
spp.nat.SR[,1] <- as.numeric(as.character(spp.nat.SR[,1]))


pred.new <- pred %>%
  mutate(time.emp.all = rowSums(.[47:50])) %>%
  left_join(spp.nat.SR, by = "OBJIDsic") %>%
  left_join(empires[,c("OBJIDsic", "UN_geospat")], by = "OBJIDsic") %>%
  select(area, dem.sd, temp.mean, aridity.mean, pop.dens.2010.mean, gdp.pc.2011.mean, crop.area.prop, time.emp.all, SR.nat, UN_geospat)

head(pred.new)

pred.new <- pred.new %>%
  mutate(
    area.s = scale(log(area)), 
    dem.sd.s = scale(log(dem.sd+0.1)), 
    temp.mean.s = scale(temp.mean), 
    aridity.mean.s = scale(log(aridity.mean+0.1)), 
    pop.dens.2010.mean.s = scale(log(pop.dens.2010.mean+0.1)), 
    gdp.pc.2011.mean.s = scale(log(gdp.pc.2011.mean+0.1)), 
    crop.area.prop.s = scale(log(crop.area.prop+0.1)),
    time.emp.all.s = scale(log(time.emp.all+0.1)),
    SR.nat.s = scale(log(SR.nat+0.1)))
  



# remove attributes from dataframe
pred.new[] = lapply(pred.new, function(x) { attributes(x) <- NULL; x })

pred.new


# Univariate model
mod <- lm(SR.nat.s ~ time.emp.all.s, pred.new)
summary(mod)

x11()
par(mfrow=c(2,2))
plot(mod)

x11()
plot(SR.nat.s ~ time.emp.all.s, pred.new, las = 1, xlab = "Time in Empire [log]", ylab = "Number of naturalized species [log]", col = "grey")
abline(mod, col = "deepskyblue2", lwd = 2)

mod2 <- glm(SR.nat ~ area.s + dem.sd.s + temp.mean.s + aridity.mean.s + pop.dens.2010.mean.s + crop.area.prop.s + time.emp.all.s, pred.new, family = "poisson")
summary(mod2)

# Generalized linear mixed model
library(lme4)
mod3 <- glmer(SR.nat ~ area.s + dem.sd.s + temp.mean.s + aridity.mean.s + pop.dens.2010.mean.s + crop.area.prop.s + time.emp.all.s + (1|UN_geospat), pred.new, family = "poisson")
summary(mod3)


library("sjPlot")
x11()
plot_model(mod3, type = "pred")

x11()
par(mfrow=c(2,2))
plot(mod2)

#

























