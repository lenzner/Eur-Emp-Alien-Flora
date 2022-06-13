#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Analysis script 02  (MS-GDM visualization)                              #
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
# Return the mean I-splines for each replicate per empire ----

### Load function to calculate mean I-splines
source(paste0(fun.path,"/Return.ispline.2.R"))

### Calculate mean I-splines for British Empire
mod.ispline.GBR <- list()
# j = number of zeta orders
# i = number of replicates

for(j in 1:4){
  mod.ispline.GBR[[j]] <- list()
  for(i in 1:6){
    mod.ispline.GBR[[j]][[i]] <- Return.ispline.2(msgdm=mod.GBR.ispl.rep[[i]][[j]],data.env=pred.calc.GBR2.s,distance=TRUE)
  }
}


### Calculate mean I-splines for Spanish Empire
mod.ispline.ESP <- list()
# j = number of zeta orders
# i = number of replicates

for(j in 1:4){
    mod.ispline.ESP[[j]] <- list()
    for(i in 1:6){
      mod.ispline.ESP[[j]][[i]] <- Return.ispline.2(msgdm=mod.ESP.ispl.rep[[i]][[j]],data.env=pred.calc.ESP2.s,distance=TRUE)
    }
  }
  

### Calculate mean I-splines for Portuguese Empire
mod.ispline.PRT <- list()
# j = number of zeta orders
# i = number of replicates

for(j in 1:4){
  mod.ispline.PRT[[j]] <- list()
  for(i in 1:6){
    mod.ispline.PRT[[j]][[i]] <- Return.ispline.2(msgdm=mod.PRT.ispl.rep[[i]][[j]],data.env=pred.calc.PRT2.s,distance=TRUE)
  }
}


### Calculate mean I-splines for Dutch Empire
mod.ispline.NED <- list()
# j = number of zeta orders
# i = number of replicates

for(j in 1:4){
  mod.ispline.NED[[j]] <- list()
  for(i in 1:6){
    mod.ispline.NED[[j]][[i]] <- Return.ispline.2(msgdm=mod.NED.ispl.rep[[i]][[j]],data.env=pred.calc.NED2.s,distance=TRUE)
  }
}



#================================================================================#
# Plot I-splines for each empire for the individual replicates and the mean trend ----

### Load function to plots individual and mean I-splines
source(paste0(fun.path,"/Plot.mean.ispline.2.R"))


### British Empire
date <- Sys.Date()
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_GBR.png"), width=1740, height=980)

  par(mfrow=c(4,9))
  for(j in 1:4){
    for(spl in 1:(ncol(pred.calc.GBR2.s)+1)){##the +1 is for distance
      #x11()
      Plot.ispline.2(mod.ispline.GBR[[j]],pred.calc.GBR2.s,distance=TRUE,biotic=FALSE,num.spline=spl,num.rep=6)
    }
  }

dev.off()
  


### Spanish Empire
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_ESP.png"), width=1740, height=980)

par(mfrow=c(4,9))
for(j in 1:4){
  for(spl in 1:(ncol(pred.calc.ESP2.s)+1)){##the +1 is for distance
    #x11()
    Plot.ispline.2(mod.ispline.ESP[[j]],pred.calc.ESP2.s,distance=TRUE,biotic=FALSE,num.spline=spl,num.rep=6)
  }
}

dev.off()

### Portuguese Empire
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_PRT.png"), width=1740, height=980)

par(mfrow=c(4,9))
for(j in 1:4){
  for(spl in 1:(ncol(pred.calc.PRT2.s)+1)){##the +1 is for distance
    #x11()
    Plot.ispline.2(mod.ispline.PRT[[j]],pred.calc.PRT2.s,distance=TRUE,biotic=FALSE,num.spline=spl,num.rep=6)
  }
}

dev.off()

### Dutch Empire
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_NED.png"), width=1740, height=980)

par(mfrow=c(4,9))
for(j in 1:4){
  for(spl in 1:(ncol(pred.calc.NED2.s)+1)){##the +1 is for distance
    #x11()
    Plot.ispline.2(mod.ispline.NED[[j]],pred.calc.NED2.s,distance=TRUE,biotic=FALSE,num.spline=spl,num.rep=6)
  }
}

dev.off()



#================================================================================#
# Plot mean I-splines for all empires ----

### Load function to plots mean I-splines
source(paste0(fun.path,"/Mean_Splines_Plotting_Dataframe.R"))

#### Plot mean splines of all empires across zeta orders ----
Isplines_to_plot <- Plot_mean_isplines(models = list(mod.ispline.GBR, mod.ispline.ESP, mod.ispline.PRT, mod.ispline.NED), num.zeta = 4, pred.num = 9)


Isplines_to_plot <- Plot_mean_isplines(models = list(mod.ispline.GBR, mod.ispline.NED), num.zeta = 4, pred.num = 8)


Splines_max_val <- Isplines_to_plot %>%
  group_by(empire, zeta, predictor) %>%
  summarize(max = max(spline)) %>%
  mutate(col = as.numeric(as.factor(predictor)))


Splines_max_val2 <- Isplines_to_plot %>%
  group_by(empire, zeta, predictor) %>%
  filter(zeta %in% c("Zeta 2", "Zeta 5")) %>%
  summarize(max = max(spline)) %>%
  mutate(col = as.numeric(as.factor(predictor)))


# Prepare and export max values (amplitudes)
pred.amp <- Splines_max_val2 %>%
  select(empire, predictor, zeta, max) %>% 
  spread(zeta, max) %>%
  mutate(`Zeta 2` = round(`Zeta 2`,2),
         `Zeta 5` = round(`Zeta 5`,2))

date <- Sys.Date()
write.table(pred.amp, paste0(res.path,"/", date,"_Predictor_amplitude.csv"), sep = ";", row.names = F)

# Display all Zeta orders
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_all_zeta_orders.png"), width=1740, height=980)

Isplines_to_plot %>%
  ggplot(aes(x = pred.range, y = spline, color = predictor)) +
  ylim(0, 3.5) +
  ylab("Mean I-splines") +
  xlab("Rescaled predictor range") +
  facet_grid(zeta ~ empire) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        # Change axis label size
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Change facet label size
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  geom_point(data=Splines_max_val, aes(x = 1, y = max), size = 3, shape = Splines_max_val$col) +
  guides(colour = guide_legend(override.aes = list(shape = c(1:9))))

dev.off()


# Display Zeta 2 and Zeta 5 only
png(file=paste0(fig.path, "/", date, "_Mean_Ispline_zeta_orders_2_and_5.png"), width=1740, height=980)

Isplines_to_plot %>%
  filter(zeta %in% c("Zeta 2", "Zeta 5")) %>%
  ggplot(aes(x = pred.range, y = spline, color = predictor)) +
  ylim(0, 4) +
  ylab("Mean I-splines \n (turnover)") +
  xlab("Rescaled predictor range") +
  facet_grid(zeta ~ empire) +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        # Change axis label size
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        # Change facet label size
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  geom_point(data=Splines_max_val2, aes(x = 1, y = max), size = 3, shape = Splines_max_val2$col) +
  guides(colour = guide_legend(override.aes = list(shape = c(1:9))))

dev.off()


#================================================================================#
# Calculate variable importance based on height of the I-Spline ----

### Load function to calculate variable importance
source(paste0(fun.path,"/Variable_Importance_MSGDM.R"))

### British Empire
preds <- colnames(mod.ispline.GBR[[1]][[1]]$Ispline) # predictor variables
zetas <- 4 # number of zeta orders of interest

var.part.preds.GBR <- data.frame(cbind(preds, matrix(NA, ncol = 4, nrow = length(preds))))
colnames(var.part.preds.GBR)[2:5] <- paste0("zeta", 2:5)

for(i in 1:zetas){
  for(j in 1:length(preds)){
    var.part.preds.GBR[j,i+1] <- round(max(mean.Ispline.val(mod.ispline.GBR[[i]],pred.calc.GBR2.s,distance=TRUE,biotic=FALSE,num.spline=j,num.rep=6)),2)
  }
}

### Spanish Empire
preds <- colnames(mod.ispline.ESP[[1]][[1]]$Ispline) # predictor variables
zetas <- 4 # number of zeta orders of interest

var.part.preds.ESP <- data.frame(cbind(preds, matrix(NA, ncol = 4, nrow = length(preds))))
colnames(var.part.preds.ESP)[2:5] <- paste0("zeta", 2:5)

for(i in 1:zetas){
  for(j in 1:length(preds)){
    var.part.preds.ESP[j,i+1] <- round(max(mean.Ispline.val(mod.ispline.ESP[[i]],pred.calc.ESP2.s,distance=TRUE,biotic=FALSE,num.spline=j,num.rep=6)),2)
  }
}

### Portuguese Empire
preds <- colnames(mod.ispline.PRT[[1]][[1]]$Ispline) # predictor variables
zetas <- 4 # number of zeta orders of interest

var.part.preds.PRT <- data.frame(cbind(preds, matrix(NA, ncol = 4, nrow = length(preds))))
colnames(var.part.preds.PRT)[2:5] <- paste0("zeta", 2:5)

for(i in 1:zetas){
  for(j in 1:length(preds)){
    var.part.preds.PRT[j,i+1] <- round(max(mean.Ispline.val(mod.ispline.PRT[[i]],pred.calc.PRT2.s,distance=TRUE,biotic=FALSE,num.spline=j,num.rep=6)),2)
  }
}

### Dutch Empire
preds <- colnames(mod.ispline.NED[[1]][[1]]$Ispline) # predictor variables
zetas <- 4 # number of zeta orders of interest

var.part.preds.NED <- data.frame(cbind(preds, matrix(NA, ncol = 4, nrow = length(preds))))
colnames(var.part.preds.NED)[2:5] <- paste0("zeta", 2:5)

for(i in 1:zetas){
  for(j in 1:length(preds)){
    var.part.preds.NED[j,i+1] <- round(max(mean.Ispline.val(mod.ispline.NED[[i]],pred.calc.NED2.s,distance=TRUE,biotic=FALSE,num.spline=j,num.rep=6)),2)
  }
}

# Change variable names
var.part.preds.GBR$preds <- c("Area", "Habitat Heterogeneity", "Temperature", "Aridity", "Population Density", "per capita GDP", "Time in Empire", "Cropland Area", "Distance")
var.part.preds.ESP$preds <- c("Area", "Habitat Heterogeneity", "Temperature", "Aridity", "Population Density", "per capita GDP", "Time in Empire", "Cropland Area", "Distance")
var.part.preds.PRT$preds <- c("Area", "Habitat Heterogeneity", "Temperature", "Aridity", "Population Density", "per capita GDP", "Time in Empire", "Cropland Area", "Distance")
var.part.preds.NED$preds <- c("Area", "Habitat Heterogeneity", "Temperature", "Aridity", "Population Density", "per capita GDP", "Time in Empire", "Cropland Area", "Distance")

date <- Sys.Date()
write.table(var.part.preds.GBR, file=paste0(res.path, "/", date, "_MSGDM_variable_importance_GBR.csv"), sep = ";", row.names = F)
write.table(var.part.preds.ESP, file=paste0(res.path, "/", date, "_MSGDM_variable_importance_ESP.csv"), sep = ";", row.names = F)
write.table(var.part.preds.PRT, file=paste0(res.path, "/", date, "_MSGDM_variable_importance_PRT.csv"), sep = ";", row.names = F)
write.table(var.part.preds.NED, file=paste0(res.path, "/", date, "_MSGDM_variable_importance_NED.csv"), sep = ";", row.names = F)
#





