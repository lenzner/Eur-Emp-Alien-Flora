#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#



Return.ispline.2 <- function (msgdm, data.env, distance = FALSE, biotic = 0)
{
  
  
  # msgdm <- mod.GBR.ispl.rep[[1]][[1]]
  # data.env <- pred.calc.GBR2
  # distance = TRUE
  # biotic = 0
  
  
  if (biotic == FALSE) 
    biotic <- 0
  if (biotic == TRUE) 
    biotic <- 1
  my.order <- msgdm$my.order
  #order.ispline = msgdm$order.ispline
  #kn.ispline = msgdm$kn.ispline
  order.ispline = 2
  kn.ispline = 1
  num.splines <- order.ispline + kn.ispline
  data.env.num <- as.data.frame(data.env[, (which(sapply(data.env, 
                                                         class) == "numeric"))])
  names(data.env.num) <- names(data.env)[(which(sapply(data.env, 
                                                       class) == "numeric"))]
  Fa <- as.data.frame(data.env[, which(sapply(data.env, class) == 
                                         "factor")])
  names(Fa) <- names(data.env)[which(sapply(data.env, class) == 
                                       "factor")]
  if (distance == FALSE & biotic == 0) {
    XX <- as.data.frame(matrix(rep(seq(0, 1, 0.01), ncol(data.env)), 
                               101, ncol(data.env)))
    names(XX) <- names(data.env)
  }else if (distance == TRUE & biotic == 0) {
    XX <- as.data.frame(matrix(rep(seq(0, 1, 0.01), ncol(data.env) + 
                                     1), 101, ncol(data.env) + 1))
    names(XX) <- c(names(data.env), "Distance")
  }else if (distance == FALSE & biotic > 0) {
    XX <- as.data.frame(matrix(rep(seq(0, 1, 0.01), ncol(data.env) + 
                                     biotic), 101, ncol(data.env) + biotic))
    names(XX) <- c(names(data.env), paste("Biotic", 1:biotic, 
                                          sep = ""))
  }else {
    XX <- as.data.frame(matrix(rep(seq(0, 1, 0.01), ncol(data.env) + 
                                     biotic + 1), 101, ncol(data.env) + biotic + 1))
    names(XX) <- c(names(data.env), paste("Biotic", 1:biotic, 
                                          sep = ""), "Distance")
  }
  env.ispline <- Ispline(data.env.num, rescale = 1, order.ispline = order.ispline, 
                         kn.ispline = kn.ispline)$splines
  subsamp <- 0
  if (length(msgdm$val) < nrow(data.env)) {
    ind.sel <- sample(nrow(data.env), length(msgdm$val))
    env.ispline <- env.ispline[ind.sel, ]
    data.env <- data.env[ind.sel, ]
    subsamp <- 1
  }
  if (distance == TRUE) {
    d.ind <- c(sample(1:length(msgdm$val), min(nrow(data.env) - 2, length(msgdm$val) - 2)), which.max(msgdm$distance),which.min(msgdm$distance))
    d <- msgdm$distance[d.ind]
    d.spline <- msgdm$predictors[d.ind, (ncol(msgdm$predictors) - 
                                           (num.splines - 1)):ncol(msgdm$predictors)]
    if (biotic > 0) {
      for (b in 1:biotic) {
        bio.ind <- c(sample(1:length(msgdm$val), min(nrow(data.env) - 
                                                       2, length(msgdm$val) - 2)), which.max(msgdm$biotic[, 
                                                                                                          b]), which.min(msgdm$biotic[, b]))
        if (b == 1) {
          bio <- msgdm$biotic[bio.ind, b]
          bio <- matrix(bio, length(bio), 1)
          bio.spline <- msgdm$predictors[bio.ind, (ncol(msgdm$predictors) - 
                                                     ((biotic - (b - 2)) * 3 - 1)):(ncol(msgdm$predictors) - 
                                                                                      ((biotic - (b - 1)) * 3))]
        }else {
          bio <- cbind(bio, msgdm$biotic[bio.ind, b])
          bio.spline <- cbind(bio.spline, msgdm$predictors[bio.ind, 
                                                           (ncol(msgdm$predictors) - ((biotic - (b - 
                                                                                                   2)) * 3 - 1)):(ncol(msgdm$predictors) - 
                                                                                                                    ((biotic - (b - 1)) * 3))])
        }
      }
      if (length(which(sapply(data.env, class) == "factor")) == 
          0) {
        X.ispline <- cbind(env.ispline, bio.spline, 
                           d.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  biotic + 1, nrow(data.env))
        for (i in 1:(ncol(data.env) + biotic + 1)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                          X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                     num.splines)])
        }
      }else {
        X.ispline <- cbind(env.ispline, matrix(rep(1, 
                                                   ncol(Fa) * nrow(env.ispline)), nrow(env.ispline), 
                                               ncol(Fa)), bio.spline, d.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  biotic + 1, nrow(data.env))
        for (i in 1:ncol(data.env.num)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                          X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                     num.splines)])
        }
        ii <- 0
        for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                            ncol(Fa))) {
          ii <- ii + 1
          Isplines.pred[i, ] <- -stats::coef(msgdm$model)[(1 + 
                                                             ncol(data.env.num) * num.splines) + ii] * 
            X.ispline[, ncol(data.env.num) * num.splines + 
                        ii]
        }
        for (b in 1:biotic) {
          Isplines.pred[ncol(data.env) + b, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                                             (ncol(data.env.num) + b - 1) * num.splines + 
                                                                                             ncol(Fa)):((ncol(data.env.num) + b) * num.splines + 
                                                                                                          ncol(Fa) + 1)], nrow(data.env), num.splines, 
                                                                byrow = TRUE) * X.ispline[, (1 + (ncol(data.env.num) + 
                                                                                                    b - 1) * num.splines + ncol(Fa)):((ncol(data.env.num) + 
                                                                                                                                         b) * num.splines + ncol(Fa))])
        }
        Isplines.pred[ncol(data.env) + biotic + 1, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                                                    (ncol(data.env.num) + biotic - 1) * num.splines + 
                                                                                                    ncol(Fa) + num.splines):((ncol(data.env.num) + 
                                                                                                                                biotic) * num.splines + ncol(Fa) + 1 + num.splines)], 
                                                                       nrow(data.env), num.splines, byrow = TRUE) * 
                                                                  X.ispline[, (1 + (ncol(data.env.num) + biotic - 
                                                                                      1) * num.splines + ncol(Fa) + num.splines):((ncol(data.env.num) + 
                                                                                                                                     biotic) * num.splines + ncol(Fa) + num.splines)])
      }
    }else {
      if (length(which(sapply(data.env, class) == "factor")) == 
          0) {
        X.ispline <- cbind(env.ispline, d.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  1, min(nrow(data.env), length(msgdm$val)))
        for (i in 1:(ncol(data.env.num) + 1)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], min(nrow(data.env), length(msgdm$val)), 
                                               num.splines, byrow = TRUE) * X.ispline[, 
                                                                                      (1 + (i - 1) * num.splines):(i * num.splines)])
        }
      }else {
        X.ispline <- cbind(env.ispline, matrix(rep(1, 
                                                   ncol(Fa) * nrow(env.ispline)), nrow(env.ispline), 
                                               ncol(Fa)), d.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  1, nrow(data.env))
        for (i in 1:(ncol(data.env.num) + 1)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                          X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                     num.splines)])
        }
        ii <- 0
        for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                            ncol(Fa))) {
          ii <- ii + 1
          Isplines.pred[i, ] <- -stats::coef(msgdm$model)[(1 + 
                                                             ncol(data.env.num) * num.splines) + ii] * 
            X.ispline[, ncol(data.env.num) * num.splines + 
                        ii]
        }
        Isplines.pred[ncol(data.env) + 1, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                                           (ncol(data.env.num)) * num.splines + ncol(Fa)):((ncol(data.env.num) + 
                                                                                                                                              1) * num.splines + ncol(Fa) + 1)], nrow(data.env), 
                                                              num.splines, byrow = TRUE) * X.ispline[, (1 + 
                                                                                                          ncol(data.env.num) * num.splines + ncol(Fa)):((ncol(data.env.num) + 
                                                                                                                                                           1) * num.splines + ncol(Fa))])
      }
    }
  }else {
    if (biotic > 0) {
      for (b in 1:biotic) {
        bio.ind <- c(sample(1:length(msgdm$val), min(nrow(data.env) - 
                                                       2, length(msgdm$val) - 2)), which.max(msgdm$biotic[, 
                                                                                                          b]), which.min(msgdm$biotic[, b]))
        if (b == 1) {
          bio <- msgdm$biotic[bio.ind, b]
          bio <- matrix(bio, length(bio), 1)
          bio.spline <- msgdm$predictors[bio.ind, (ncol(msgdm$predictors) - 
                                                     ((biotic - (b - 1)) * 3 - 1)):(ncol(msgdm$predictors) - 
                                                                                      ((biotic - (b)) * 3))]
        }else {
          bio <- cbind(bio, msgdm$biotic[bio.ind, b])
          bio.spline <- cbind(bio.spline, msgdm$predictors[bio.ind, 
                                                           (ncol(msgdm$predictors) - ((biotic - (b - 
                                                                                                   1)) * 3 - 1)):(ncol(msgdm$predictors) - 
                                                                                                                    ((biotic - (b)) * 3))])
        }
      }
      if (length(which(sapply(data.env, class) == "factor")) == 
          0) {
        X.ispline <- cbind(env.ispline, bio.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  biotic, nrow(data.env))
        for (i in 1:(ncol(data.env) + biotic)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                          X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                     num.splines)])
        }
      }else {
        X.ispline <- cbind(env.ispline, matrix(rep(1, 
                                                   ncol(Fa) * nrow(env.ispline)), nrow(env.ispline), 
                                               ncol(Fa)), bio.spline)
        Isplines.pred <- matrix(NA, ncol(data.env) + 
                                  biotic, nrow(data.env))
        for (i in 1:ncol(data.env.num)) {
          Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                            (i - 1) * num.splines):(i * num.splines + 
                                                                                                      1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                          X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                     num.splines)])
        }
        ii <- 0
        for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                            ncol(Fa))) {
          ii <- ii + 1
          Isplines.pred[i, ] <- -stats::coef(msgdm$model)[(1 + 
                                                             ncol(data.env.num) * num.splines) + ii] * 
            X.ispline[, ncol(data.env.num) * num.splines + 
                        ii]
        }
        for (b in 1:biotic) {
          Isplines.pred[ncol(data.env) + b, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                                             (ncol(data.env.num) + b - 1) * num.splines + 
                                                                                             ncol(Fa)):((ncol(data.env.num) + b) * num.splines + 
                                                                                                          ncol(Fa) + 1)], nrow(data.env), num.splines, 
                                                                byrow = TRUE) * X.ispline[, (1 + (ncol(data.env.num) + 
                                                                                                    b - 1) * num.splines + ncol(Fa)):((ncol(data.env.num) + 
                                                                                                                                         b) * num.splines + ncol(Fa))])
        }
      }
    }else {
      if (length(which(sapply(data.env, class) == "factor")) == 
          0) {
        X.ispline <- env.ispline
      }else {
        X.ispline <- cbind(env.ispline, matrix(rep(1, 
                                                   ncol(Fa) * nrow(env.ispline)), nrow(env.ispline), 
                                               ncol(Fa)))
      }
      Isplines.pred <- matrix(NA, ncol(data.env), nrow(data.env))
      for (i in 1:ncol(data.env.num)) {
        Isplines.pred[i, ] <- rowSums(matrix(-stats::coef(msgdm$model)[(2 + 
                                                                          (i - 1) * num.splines):(i * num.splines + 
                                                                                                    1)], nrow(data.env), num.splines, byrow = TRUE) * 
                                        X.ispline[, (1 + (i - 1) * num.splines):(i * 
                                                                                   num.splines)])
      }
      if (length(which(sapply(data.env, class) == "factor")) > 
          0) {
        ii <- 0
        for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                            ncol(Fa))) {
          ii <- ii + 1
          Isplines.pred[i, ] <- -stats::coef(msgdm$model)[(1 + 
                                                             ncol(data.env.num) * num.splines) + ii] * 
            X.ispline[, ncol(data.env.num) * num.splines + 
                        ii]
        }
      }
    }
  }
  if (subsamp == 1) {
    env.resc <- data.env.num[ind.sel, ]
  }else {
    env.resc <- data.env.num
  }
  for (i in 1:ncol(env.resc)) {
    env.resc[, i] <- (env.resc[i] - min(env.resc[i]))/(max(env.resc[i]) - 
                                                         min(env.resc[i]))
  }
  splines.out <- list()
  env.resc.out <- env.resc
  Isplines.pred.out <- Isplines.pred
  for (i in 1:ncol(data.env.num)) {
    env.resc.out[, i] <- sort(env.resc[, i])
    Isplines.pred.out[i, ] <- Isplines.pred[i, order(env.resc[, 
                                                              i])]
  }
  splines.out$env <- data.frame(env.resc.out)
  if (ncol(Fa) > 0) {
    Fa.out <- matrix(0, ncol(Fa), ncol(Isplines.pred))
    ii <- 0
    for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                        ncol(Fa))) {
      ii <- ii + 1
      Fa.out[ii, ] <- seq(0, 1, 1/(ncol(Isplines.pred) - 
                                     1))
      Isplines.pred.out[i, ] <- seq(0, max(Isplines.pred[i, 
      ]), max(Isplines.pred[i, ])/(ncol(Isplines.pred) - 
                                     1))
    }
    splines.out$env <- cbind(splines.out$env, data.frame(t(Fa.out)))
    for (i in (ncol(data.env.num) + 1):(ncol(data.env.num) + 
                                        ncol(Fa))) {
      names(splines.out$env)[i] <- paste("Fa", i - ncol(data.env.num), 
                                         sep = "")
    }
  }
  if (biotic > 0) {
    for (b in 1:biotic) {
      i <- i + 1
      bio.out <- sort(bio[, b])
      Isplines.pred.out[i, ] <- Isplines.pred[i, order(bio[, 
                                                           b])]
      splines.out$env <- cbind(splines.out$env, bio.out)
      names(splines.out$env)[i] <- paste("biotic", b, 
                                         sep = "")
    }
  }
  if (distance == TRUE) {
    i <- i + 1
    d.out <- sort(d/max(d))
    Isplines.pred.out[i, ] <- Isplines.pred[i, order(d)]
    splines.out$env <- cbind(splines.out$env, d.out)
    names(splines.out$env)[i] <- "distance"
  }
  splines.out$Ispline <- data.frame(t(Isplines.pred.out))
  names(splines.out$Ispline) <- names(splines.out$env)
  splines.out$env.num <- data.env.num
  splines.out$Fa <- Fa
  splines.out$distance <- distance
  splines.out$biotic <- biotic
  splines.out$my.order <- my.order
  return(splines.out)
}
