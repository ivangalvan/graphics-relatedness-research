#
#
#    March 15, 2017  
#
#    Graphics for relatedness research     ######
#
#    Galvan-Femenia I., Graffelman J., Barcelo-i-Vidal C.
#
#
#
#    CPU: Linux 3.16.0-4-amd64 #1 SMP Debian 3.16.36-1+deb8u2 (2016-10-19) x86_64 GNU/Linux
#    Model name: Intel(R) Xeon(R) CPU E5-2609 v2 @ 2.50GHz
#    Memory: 32.00 GB
#    Cores: 8
#
#
#    Data from Rosenberg lab at Standford University: https://rosenberglab.stanford.edu/diversity.html
#
#    multiplot function from Peter Haschke: http://www.peterhaschke.com/r/2013/04/24/MultiPlot.html
#
#    stat_chull function from A. Kassambara : https://github.com/kassambara/ggpubr/blob/master/R/stat_chull.R
#





onerowperAll <- function (x, id, colmarkers) 

  {
  if (missing(colmarkers)) 
    colmarkers <- 1:ncol(x)
  else colmarkers <- colmarkers
  x <- dplyr::select(x, colmarkers)
  x <- split(x, f = id)
  return(do.call("rbind", lapply(x, function(aux) {
    aux <- apply(aux, MARGIN = 2, FUN = function(a) {
      if (any(is.na(a))) {
        NA
      } else paste(a, collapse = "/")
    })
  })))
}


ibs <- function (x, y) 
{
  m <- cbind(x, y)
  apply(m, MARGIN = 1, FUN = function(md) {
    max((md[1] == md[3]) + (md[2] == md[4]), (md[1] == md[4]) + 
          (md[2] == md[3]))
  })
}


allelesharing <- function (X, colmarkers, splits = 4) 

  {
  require(compiler)
  ibs <- cmpfun(ibs)
  if (missing(colmarkers)) 
    colmarkers <- 1:col(X)
  else colmarkers <- colmarkers
  X <- X[, colmarkers]
  xrnames <- rownames(X)
  if (is.null(xrnames)) 
    xrnames <- 1:nrow(X)
  xcnames <- colnames(X)
  if (is.null(xcnames)) 
    xcnames <- colmarkers
  n <- nrow(X)
  rdsfiles <- paste("tmp-", 1:splits, ".rds", sep = "")
  tsplits <- round(seq(1, ncol(X), length.out = splits + 1), 
                   0)
  for (i in 1:length(rdsfiles)) {
    if (i != length(rdsfiles)) {
      saveRDS(X[, tsplits[i]:(tsplits[i + 1] - 1)], file = rdsfiles[i])
    }
    else saveRDS(X[, tsplits[i]:tsplits[i + 1]], file = rdsfiles[i])
  }
  rm(X)
  gc()
  xc <- combn(1:n, m = 2)
  for (r in rdsfiles) {
    X <- readRDS(r)
    sX <- lapply(1:ncol(X), function(i) {
      s <- do.call("rbind", strsplit(X[, i], split = "/"))
    })
    rm(X)
    sX <- do.call("cbind", lapply(sX, function(s) {
      ibs(s[xc[1, ], ], s[xc[2, ], ])
    }))
    saveRDS(sX, file = r)
    rm(sX)
    gc()
  }
  x <- do.call("cbind", lapply(1:length(rdsfiles), function(i) {
    readRDS(rdsfiles[[i]])
  }))
  rownames(x) <- paste(xrnames[xc[1, ]], xrnames[xc[2, ]], 
                       sep = "-")
  colnames(x) <- xcnames
  unlink(rdsfiles, recursive = TRUE)
  return(x)
}



percentages <- function (x, labels = c("P0", "P1", "P2")) 
{
  res <- apply(x, MARGIN = 1, FUN = function(xrow) {
    xrow <- factor(xrow, levels = c(0, 1, 2), labels = labels)
    prop.table(table(xrow))
  })
  return(as.data.frame(t(res)))
}



mean_sd <- function(x){
  
  
  n <- dim(x)[1]
  
  mean_ibs <- rep(NA, n)
  sd_ibs <- rep(NA, n)
  
  for (i in 1:n) {
    mean_ibs[i] <- mean(x[i, ], na.rm = T)
    sd_ibs[i] <- sd(x[i, ], na.rm = T)
  }
  
  return(data.frame(mean=mean_ibs,sd=sd_ibs))
  
  
}


children <- function (X, father = NULL, mother = NULL) 
{
  m <- ncol(X)
  fatherRow <- X[father, ]
  motherRow <- X[mother, ]
  fatherAll <- do.call(rbind, strsplit(fatherRow, "/"))
  motherAll <- do.call(rbind, strsplit(motherRow, "/"))
  fSample <- sample(1:2, size = m, replace = TRUE, prob = c(0.5,0.5))
  mSample <- sample(1:2, size = m, replace = TRUE, prob = c(0.5,0.5))
  child <- sapply(1:m, function(i) {
    paste(fatherAll[i, fSample[i]], motherAll[i, mSample[i]], 
          sep = "/")
  })
  X <- rbind(X, child)
  rownames(X) <- paste("ID", 1:nrow(X), sep = "")
  return(X)
}





StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y"))


stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



probIBD <- function (X, pairs, colmarkers = NULL, ibs = NULL) {
  number_pairs = length(pairs)
  mylist = list()
  m <- length(colmarkers)
  n <- dim(X)[1]
  X <- X[, colmarkers]
  for (h in 1:number_pairs) {
    pair = NULL
    pair = pairs[h]
    if (require(genetics) == F) {
      install.packages("genetics")
      require(genetics)
    }
    else {
      require(genetics)
    }
    id1 <- identifyID(pair, n)[2]
    id2 <- identifyID(pair, n)[3]
    s <- ibs[pair, ]
    ML <- matrix(NA, m, 3)
    for (k in 1:m) {
      p <- summary(genotype(X[, k]))$allele.freq
      g1 <- strsplit(X[id1, k], "/")[[1]]
      g2 <- strsplit(X[id2, k], "/")[[1]]
      u <- unique(c(g1, g2))
      if (s[k] == 2 && is.na(s[k]) == FALSE) {
        if (length(u) == 2) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          ML[k, 1] <- 4 * (pi^2) * (pj^2)
          ML[k, 2] <- pi * pj * (pi + pj)
          ML[k, 3] <- 2 * pi * pj
        }
        else if (length(u) == 1) {
          pi <- p[u[1], 2]
          ML[k, 1] <- (pi^4)
          ML[k, 2] <- (pi^3)
          ML[k, 3] <- (pi^2)
        }
      }
      else if (s[k] == 1 && is.na(s[k]) == FALSE) {
        if (length(u) == 3) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          pk <- p[u[3], 2]
          ML[k, 1] <- 4 * (pi^2) * pj * pk
          ML[k, 2] <- pi * pj * pk
          ML[k, 3] <- 0
        }
        else if (length(u) == 2) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          ML[k, 1] <- 2 * (pi^3) * pj
          ML[k, 2] <- (pi^2) * pj
          ML[k, 3] <- 0
        }
      }
      else if (s[k] == 0 && is.na(s[k]) == FALSE) {
        if (length(u) == 4) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          pk <- p[u[3], 2]
          pl <- p[u[4], 2]
          ML[k, 1] <- 4 * pi * pj * pk * pl
          ML[k, 2] <- 0
          ML[k, 3] <- 0
        }
        else if (length(u) == 3) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          pk <- p[u[3], 2]
          ML[k, 1] <- 2 * (pi^2) * pj * pk
          ML[k, 2] <- 0
          ML[k, 3] <- 0
        }
        else if (length(u) == 2) {
          pi <- p[u[1], 2]
          pj <- p[u[2], 2]
          ML[k, 1] <- (pi^2) * (pj^2)
          ML[k, 2] <- 0
          ML[k, 3] <- 0
        }
      }
    }
    ml <- ML[!is.na(ML[, 1]), ]
    ml <- as.matrix(ml)
    mylist[[h]] = ml
  }
  return(mylist)
}



MaximumLikelihood <- function (X, ...) 
{
  if (require(Rsolnp) == F) {
    install.packages("Rsolnp")
    require(Rsolnp)
  }
  else {
    require(Rsolnp)
  }
  eqn1 = function(y) {
    z1 = y[1] + y[2] + y[3]
    return(z1)
  }
  eqn2 = function(y) {
    z1 = y[2]^2 - 4 * y[1] * y[3]
    return(z1)
  }
  fn1 = function(y) {
    -sum(log(as.matrix(X) %*% y))
  }
  out <- solnp(c(0.25, 0.5, 0.25), fun = fn1, eqfun = eqn1, 
               eqB = c(1), ineqfun = eqn2, ineqLB = c(0), ineqUB = c(1), 
               LB = c(0, 0, 0), UB = c(1, 1, 1), ...)
  return(out$pars)
}



identifyID <- function (x, n) 
{
  id <- matrix(NA, nrow = length(x), ncol = 3)
  R <- matrix(NA, n, n)
  C <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      R[i, j] <- i
      C[i, j] <- j
    }
  }
  first <- C[lower.tri(C, diag = F)]
  second <- R[lower.tri(R, diag = F)]
  id[, 1] <- x
  id[, 2] <- first[x]
  id[, 3] <- second[x]
  colnames(id) <- c("pair", "ID1", "ID2")
  rownames(id) <- rep("", length(x))
  return(id)
}



cotterman <- function (X, pairs, colmarkers = 1:377, ibs = NULL) 
{
  IBDprob = probIBD(X, pairs, colmarkers, ibs)
  cot_list = lapply(IBDprob, MaximumLikelihood)
  return_matrix = matrix(unlist(cot_list), ncol = 3, byrow = TRUE)
  colnames(return_matrix) = c("k0", "k1", "k2")
  rownames(return_matrix) = pairs
  return(data.frame(return_matrix))
}




ilr_maximum <- function(x){
  
  a <- sqrt(2)
  b <- sqrt(6)
  
  z1 <- x[1]
  z2 <- x[2]
  
  d0 <- prob.rel[,1]
  d1 <- prob.rel[,2]
  d2 <- prob.rel[,3]
  
  -sum(log(d0*exp(a*z1) + d1*exp(0.5*(a*z1-b*z2)) + d2) - log( 1 + exp(a*z1) + exp(0.5*(a*z1-b*z2))) )
}



