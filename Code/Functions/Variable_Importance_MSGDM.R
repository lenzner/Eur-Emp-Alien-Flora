#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#



mean.Ispline.val <- function (isplines = NULL, msgdm, data.env, distance = FALSE, 
                                 biotic = 0, pch = NULL, lty = 1, legend = TRUE, lwd = 1, num.spline=0,col="grey",num.rep=0) 
{
  
  
  # isplines = mod.ispline[[1]]
  # #msgdm=pred.calc.GBR2
  # data.env=pred.calc.GBR2
  # distance = FALSE
  # biotic = 0
  # pch = NULL
  # lty = NULL
  # legend = TRUE
  # lwd = 1
  cex = 1
  # num.quantiles = 11
  # num.spline=10
  # col="grey"
  # num.rep=5
  
  env.resc <- matrix(NA,nrow(isplines[[1]]$env),num.rep)
  Isplines.pred <- matrix(NA,nrow(isplines[[1]]$env),num.rep)
  #distance <- matrix(NA,nrow(isplines[[1]]$env),num.rep)
  for(i in 1:num.rep){
    env.resc[,i] <- isplines[[i]]$env[,num.spline]
    Isplines.pred[,i] <- isplines[[i]]$Ispline[,num.spline]
    #distance[,i] <- isplines[[i]]$distance
  }
  
  #my.order <- isplines[[1]]$my.order
  
  
  
  #if (distance == FALSE) {
  #graphics::plot(sort(env.resc[, 1]), Isplines.pred[order(env.resc[,1]), 1], type = "l", ylim = c(0, max(Isplines.pred)),main = names(isplines[[1]]$env)[num.spline], xlab = "Rescaled range", ylab = "I-splines",cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex = cex,lwd = lwd,col=col)
  
  #### EDIT BERND: main = "" to main = names(isplines[[1]]$env)[num.spline]
  
  #for(i in 2:num.rep){
  #  graphics::lines(sort(env.resc[, i]), Isplines.pred[order(env.resc[,i]), i],lwd = lwd,col=col)
  #}
  
  ####EDIT Bernd
  mean.spline <- rowMeans(Isplines.pred)
 
  return(mean.spline) 
}
