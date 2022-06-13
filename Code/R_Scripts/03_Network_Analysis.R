#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Analysis script 03 (Network analysis)                                   #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#

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
# Load relevant function for network visualization ----
source(paste0(fun.path,"/Network.stats.R"))
source(paste0(fun.path,"/Plot.network2.R"))


#================================================================================#
# Plot networks with statistics based on full empire and reduced links ----

stat.thresh = 0 # Define pairwise similarity threshold to calculate network centrality and clusters
net.sim.thresh = 0.2   # Define pairwise similarity threshold used to plot the network

net.stat.GBR <- Network.stats(ruler = "Great_Britain", sim.thresh.stats = stat.thresh)
net.stat.ESP <- Network.stats(ruler = "Spain", sim.thresh.stats = stat.thresh)
net.stat.PRT <- Network.stats(ruler = "Portugal", sim.thresh.stats = stat.thresh)
net.stat.NED <- Network.stats(ruler = "Netherlands", sim.thresh.stats = stat.thresh)


#================================================================================#
# Rename top 5 regions per cluster for plotting ----

### British Empire
t5.GBR.names.new <- c("Port Curtis (AUS)","Fujian (CHN)","Odisha (IND)","Kennedy South (AUS)","Kennedy North (AUS)","Darling Downs (AUS)","Moreton (AUS)","Wide Bay (AUS)","Burnett (AUS)","Mpumalanga (ZAF)","Tamil Nadu (IND)","Himachal Pradesh (IND)","Maharashtra (IND)","Uttarakhand (IND)","Karnataka (IND)")

t5.GBR <- net.stat.GBR %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>% 
  bind_cols(Region = t5.GBR.names.new)
  
for(i in 1:dim(t5.GBR)[1]){
  
  x <- net.stat.GBR %>% filter(OBJIDsic == t5.GBR$OBJIDsic[i]) %>%
    mutate(Region_unique = t5.GBR$Region[i])
  
  net.stat.GBR[net.stat.GBR$OBJIDsic == t5.GBR$OBJIDsic[i],] <- x
  
}

### Spanish Empire
t5.ESP.names.new <- c("Guerrero (MEX)","Colima (MEX)","Sinaloa (MEX)","Cundinamarca (COL)","Nariño (COL)","Distrito Federal (MEX)","Tlaxcala (MEX)","Baja California Norte (MEX)","Santiago (CHL)","Valparaiso (CHL)","Queretaro (MEX)","Jalisco (MEX)","Puebla (MEX)","Michoacan (MEX)","Mexico State (MEX)")

t5.ESP <- net.stat.ESP %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>% 
  bind_cols(Region = t5.ESP.names.new)

for(i in 1:dim(t5.ESP)[1]){
  
  x <- net.stat.ESP %>% filter(OBJIDsic == t5.ESP$OBJIDsic[i]) %>%
    mutate(Region_unique = t5.ESP$Region[i])
  
  net.stat.ESP[net.stat.ESP$OBJIDsic == t5.ESP$OBJIDsic[i],] <- x
  
}

### Portuguese Empire
t5.PRT.names.new <- c("Manica (MOZ)","Maputo (MOZ)","Sofala (MOZ)","Tete (MOZ)","Niassa (MOZ)","Karnataka (IND)","Tamil Nadu (IND)","Kerala (IND)","Maharashtra (IND)","East Timor (TLS)","Paraiba (BRA)","Pernambuco (BRA)","Bahia (BRA)","Ceara (BRA)","Mato Grosso do Sul (BRA)","Zhejiang (CHN)","Flores (PRT)","Sao Jorge (PRT)","St Helena (GBR)","Graciosa (PRT)")

t5.PRT <- net.stat.PRT %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>% 
  bind_cols(Region = t5.PRT.names.new)

for(i in 1:dim(t5.PRT)[1]){
  
  x <- net.stat.PRT %>% filter(OBJIDsic == t5.PRT$OBJIDsic[i]) %>%
    mutate(Region_unique = t5.PRT$Region[i])
  
  net.stat.PRT[net.stat.PRT$OBJIDsic == t5.PRT$OBJIDsic[i],] <- x
  
}

### Dutch Empire
t5.NED.names.new <- c("Paraiba (BRA)","Pernambuco (BRA)","Ceara (BRA)","Rio Grande do Norte (BRA)","St John (USA)","Mpumalanga (ZAF)","Limpopo (ZAF)","KwaZulu Natal (ZAF)","Gauteng (ZAF)","North West (ZAF)","Java (IND)","Penninsular Malaysia (MYS)","Sulawesi (IDN)","Cocos Keeling Islands (AUS)","Sumatra (IDN)","Kerala (IND)","Tamil Nadu (IND)","Andhra Pradesh (IND)","West Bengal (IND)","Gujarat (IND)")

t5.NED <- net.stat.NED %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>% 
  bind_cols(Region = t5.NED.names.new)

for(i in 1:dim(t5.NED)[1]){
  
  x <- net.stat.NED %>% filter(OBJIDsic == t5.NED$OBJIDsic[i]) %>%
    mutate(Region_unique = t5.NED$Region[i])
  
  net.stat.NED[net.stat.NED$OBJIDsic == t5.NED$OBJIDsic[i],] <- x
  
}

#================================================================================#
# Plot networks for all empires ----
emp0.GBR <- Plot.network.sim.empire2(ruler = "Great_Britain", sim.thresh.plot = net.sim.thresh, net.stats = net.stat.GBR, sim.thresh.stats = stat.thresh, weight = F, top.reg = 5, clust = T)
emp0.ESP <- Plot.network.sim.empire2(ruler = "Spain", sim.thresh.plot = net.sim.thresh, net.stats = net.stat.ESP, sim.thresh.stats = stat.thresh, weight = F, top.reg = 5, clust = T)
emp0.PRT <- Plot.network.sim.empire2(ruler = "Portugal", sim.thresh.plot = net.sim.thresh, net.stats = net.stat.PRT, sim.thresh.stats = stat.thresh, weight = F, top.reg = 5, clust = T)
emp0.NED <- Plot.network.sim.empire2(ruler = "Netherlands", sim.thresh.plot = net.sim.thresh, net.stats = net.stat.NED, sim.thresh.stats = stat.thresh, weight = F, top.reg = 5, clust = T)


date <- Sys.Date()
png(file=paste0(fig.path, "/", date, "_Networks.png"), width=1740, height=980)

grid.arrange(emp0.GBR, emp0.ESP, emp0.PRT, emp0.NED, nrow = 2)

dev.off()


#===================================================================#
# Extract 5 most important nodes per cluster ----

# Get region with highest centrality per cluster module
top.reg.GBR <- net.stat.GBR %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

top.reg.ESP <- net.stat.ESP %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

top.reg.PRT <- net.stat.PRT %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

top.reg.NED <- net.stat.NED %>%
  group_by(cluster.mod) %>%
  top_n(5, cent.eig) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

top.reg <- rbind(top.reg.GBR, top.reg.ESP, top.reg.PRT, top.reg.NED)


date <- Sys.Date()
write.table(top.reg, paste0(res.path,"/",Sys.Date(),"_Network_table_top_regions.csv"), sep = ";", row.names = F)


#===================================================================#
# Extract full network stats ----

# Get region with highest centrality per cluster module
net.stat.GBR.exp <- net.stat.GBR %>%
  group_by(cluster.mod) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

net.stat.ESP.exp <- net.stat.ESP %>%
  group_by(cluster.mod) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

net.stat.PRT.exp <- net.stat.PRT %>%
  group_by(cluster.mod) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)

net.stat.NED.exp <- net.stat.NED %>%
  group_by(cluster.mod) %>%
  arrange(cluster.mod, desc(cent.eig)) %>%
  select(Region_unique, cluster.mod, cent.eig) %>%
  mutate(cent.eig = round(cent.eig,2)) %>%
  rename(Region = Region_unique, Modularity = cluster.mod, Centrality = cent.eig)


date <- Sys.Date()
 write.table(net.stat.GBR.exp, paste0(res.path,"/",Sys.Date(),"_Network_stats_GBR.csv"), sep = ";", row.names = F)
 write.table(net.stat.ESP.exp, paste0(res.path,"/",Sys.Date(),"_Network_stats_ESP.csv"), sep = ";", row.names = F)
 write.table(net.stat.PRT.exp, paste0(res.path,"/",Sys.Date(),"_Network_stats_PRT.csv"), sep = ";", row.names = F)
 write.table(net.stat.NED.exp, paste0(res.path,"/",Sys.Date(),"_Network_stats_NED.csv"), sep = ";", row.names = F)
 
 
