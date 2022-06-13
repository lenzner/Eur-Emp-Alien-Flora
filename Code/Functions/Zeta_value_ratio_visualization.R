#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#


# Visualize zeta value and zeta ratio for the observed and random empires

Plot.zeta.emp <- function(ruler){
  
  if(ruler == "Great_Britain") {
    zeta.dec.obs.emp <- zeta.dec.ex.obs.GBR
    zeta.ran.seq.emp <- zeta.ran.seq.GBR
    COL = "#0066CC"}
  if(ruler == "Spain") {
    zeta.dec.obs.emp <- zeta.dec.ex.obs.ESP
    zeta.ran.seq.emp <- zeta.ran.seq.ESP
    COL = "#CC3300"}
  if(ruler == "Portugal") {
    zeta.dec.obs.emp <- zeta.dec.ex.obs.PRT
    zeta.ran.seq.emp <- zeta.ran.seq.PRT
    COL = "#660066"}
  if(ruler == "Netherlands") {
    zeta.dec.obs.emp <- zeta.dec.ex.obs.NED
    zeta.ran.seq.emp <- zeta.ran.seq.NED
    COL = "#FF9933"}
  
  
  
  zeta.emp <- data.frame(zeta.ran.seq.emp[[1]]$zeta.val)
  colnames(zeta.emp) <- "zeta.val.r1"
  
  for(i in 2:length(zeta.ran.seq.emp)){
    zeta.emp[,i] <- zeta.ran.seq.emp[[i]]$zeta.val
    colnames(zeta.emp)[i] <- paste0("zeta.val.","r",i)
  }
  
  zeta.rat.emp <- data.frame(zeta.ran.seq.emp[[1]]$ratio)
  colnames(zeta.rat.emp) <- "zeta.rat.r1"
  
  for(i in 2:length(zeta.ran.seq.emp)){
    zeta.rat.emp[,i] <- zeta.ran.seq.emp[[i]]$ratio
    colnames(zeta.rat.emp)[i] <- paste0("zeta.rat.","r",i)
  }
  
  x11(width = 11, height = 6)
  par(mfrow=c(1,2))
  plot(zeta.dec.obs.emp$zeta.val, type = "b", col = COL, ylim = c(-100,500), xlim = c(0,10), las = 1, main = paste0("Zeta diversity decline - ", ruler), ylab = "Zeta value", xlab = "Zeta order", lwd = 3)
  points(zeta.dec.obs.emp$zeta.val-zeta.dec.obs.emp$zeta.val.sd, type = "l", col = COL, lty = 2)
  points(zeta.dec.obs.emp$zeta.val+zeta.dec.obs.emp$zeta.val.sd, type = "l", col = COL, lty = 2)
  legend("topright", c("observed", "random", "sd"), pch = c(1,1,26), col = c(COL, "darkgrey",COL), lty = c(1,2,2), ncol = 2)
  
  for(i in 1:dim(zeta.emp)[2]){
    points(zeta.emp[,i], type = "b", col = "darkgrey")
  }
  #points(zeta.ran[[i]]$zeta.val+zeta.ran[[i]]$zeta.val.sd, type = "l", col = "darkgrey", lty = 2)
  #points(zeta.ran[[i]]$zeta.val-zeta.ran[[i]]$zeta.val.sd, type = "l", col = "darkgrey", lty = 2)
  
  
  plot(zeta.dec.obs.emp$ratio, type = "b", col = COL, ylim = c(0,1), xlim = c(0,40), las = 1, main = paste0("Ratio of zeta diversity decline - ", ruler), ylab = "Zeta ratio - retention rate", xlab = "Zeta order", lwd = 3)
  legend("topright", c("observed", "random"), pch = 1, col = c(COL, "darkgrey"))
  for(i in 1:dim(zeta.rat.emp)[2]){
    points(zeta.rat.emp[,i], type = "b", col = "darkgrey")
  }
  
}