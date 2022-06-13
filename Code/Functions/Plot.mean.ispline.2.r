#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#


Plot.mean.ispline.2 <- function (isplines = NULL, msgdm, data.env, distance = FALSE, 
                            biotic = 0, pch = NULL, lty = 1, legend = TRUE, lwd = 1, num.spline=0,col="grey",num.rep=0, mean.col = "deepskyblue2", pch.p = 2) 
{
  
  cex = 1
  
  env.resc <- matrix(NA,nrow(isplines[[1]]$env),num.rep)
  Isplines.pred <- matrix(NA,nrow(isplines[[1]]$env),num.rep)
  
  for(i in 1:num.rep){
    env.resc[,i] <- isplines[[i]]$env[,num.spline]
    Isplines.pred[,i] <- isplines[[i]]$Ispline[,num.spline]
    
  }

  mean.spline <- rowMeans(Isplines.pred)

  #### draw lines only
  graphics::lines(env.resc[, 1], mean.spline, col = mean.col, lwd = 2, type = "l", ylim = c(0, max(Isplines.pred)),main = names(isplines[[1]]$env)[num.spline], xlab = "Rescaled range", ylab = "I-splines",cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, cex = cex, las = 1)
  
  graphics::points(quantile(env.resc[, 1], probs = seq(.1, .9, by = .1)), quantile(mean.spline, probs = seq(.1, .9, by = .1)), pch = pch.p, col = mean.col, cex = 2)
  
   
  invisible(isplines)
}
