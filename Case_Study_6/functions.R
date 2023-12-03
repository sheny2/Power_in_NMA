## For binary outcomes, the input data have four columns:
##   sid (study ID), tid (treatment ID), r (event number), and n (sample size)
##
## For continuous outcomes, the input data have five columns:
##   sid (study ID), tid (treatment ID), y (mean effect), s (standard error), and n (sample size)

#####################################################################################
## Direct evidence for binary outcomes
dir.binary <- function(data){
  N <- max(data$sid)
  K <- max(data$tid)
  correction <- FALSE
  if(any(data$r == 0) | any(data$n - data$r == 0)){
    data$r.c <- data$r + 0.5
    data$n.c <- data$n + 1
    correction <- TRUE
  }
  studies <- matrix(0, K, K)
  sz <- matrix(0, K, K)
  prec <- matrix(0, K, K)
  for(i in 1:N){
    temp <- data[data$sid == i,]
    t.temp <- temp$tid
    for(j in 1:length(t.temp)){
      for(k in 1:length(t.temp)){
        if(j != k){
          studies[t.temp[j], t.temp[k]] <- studies[t.temp[j], t.temp[k]] + 1
          sz[t.temp[j], t.temp[k]] <- sz[t.temp[j], t.temp[k]] + temp$n[j] + temp$n[k]
          if(correction){
            prec[t.temp[j], t.temp[k]] <- prec[t.temp[j], t.temp[k]] + 1/(1/(temp$r.c[j]) + 1/(temp$n.c[j] - temp$r.c[j]) + 1/(temp$r.c[k]) + 1/(temp$n.c[k] - temp$r.c[k]))
          }else{
            prec[t.temp[j], t.temp[k]] <- prec[t.temp[j], t.temp[k]] + 1/(1/(temp$r[j]) + 1/(temp$n[j] - temp$r[j]) + 1/(temp$r[k]) + 1/(temp$n[k] - temp$r[k]))
          }
        }
      }
    }
  }
  out <- list(studies = studies, sz = sz, prec = prec)
  return(out)
}

#####################################################################################
## Direct evidence for continuous outcomes
dir.cont <- function(data){
  N <- max(data$sid)
  K <- max(data$tid)
  studies <- matrix(0, K, K)
  sz <- matrix(0, K, K)
  prec <- matrix(0, K, K)
  for(i in 1:N){
    temp <- data[data$sid == i,]
    t.temp <- temp$tid
    for(j in 1:length(t.temp)){
      for(k in 1:length(t.temp)){
        if(j != k){
          studies[t.temp[j], t.temp[k]] <- studies[t.temp[j], t.temp[k]] + 1
          sz[t.temp[j], t.temp[k]] <- sz[t.temp[j], t.temp[k]] + temp$n[j] + temp$n[k]
          prec[t.temp[j], t.temp[k]] <- prec[t.temp[j], t.temp[k]] + 1/(1/(temp$n[j]) + 1/(temp$n[k]))
        }
      }
    }
  }
  out <- list(studies = studies, sz = sz, prec = prec)
  return(out)
}

#####################################################################################
## Overall evidence for binary outcomes
eff.binary <- function(data){
  studies <- dir.binary(data)$studies
  sz <- dir.binary(data)$sz
  prec <- dir.binary(data)$prec
  conn <- (studies > 0)
  K <- dim(conn)[1]
  gra <- graph_from_adjacency_matrix(conn, mode = "undirected")

  eff.studies <- matrix(NA, K, K)
  eff.sz <- matrix(NA, K, K)
  eff.prec <- matrix(NA, K, K)

  for(i in 1:(K - 1)){
    for(j in (i + 1):K){
      paths <- all_simple_paths(gra, from = i, to = j)
      n.paths <- length(paths)
      if(n.paths == 1){
        path.temp <- paths[[1]]
        eff.studies.temp <- 0
        eff.sz.temp <- 0
        eff.prec.temp <- 0
        for(m in 1:(length(path.temp) - 1)){
          eff.studies.temp <- eff.studies.temp + 1/studies[path.temp[m], path.temp[m + 1]]
          eff.sz.temp <- eff.sz.temp + 1/sz[path.temp[m], path.temp[m + 1]]
          eff.prec.temp <- eff.prec.temp + 1/prec[path.temp[m], path.temp[m + 1]]
        }
        eff.studies[i, j] <- eff.studies[j, i] <- 1/eff.studies.temp
        eff.sz[i, j] <- eff.sz[j, i] <- 1/eff.sz.temp
        eff.prec[i, j] <- eff.prec[j, i] <- 1/eff.prec.temp
      }
      if(n.paths > 1){
        S.eff.studies <- matrix(NA, n.paths, n.paths)
        S.eff.sz <- matrix(NA, n.paths, n.paths)
        S.eff.prec <- matrix(NA, n.paths, n.paths)
        for(k in 1:n.paths){
          path.temp <- as.numeric(paths[[k]])
          eff.studies.temp <- 0
          eff.sz.temp <- 0
          eff.prec.temp <- 0
          for(m in 1:(length(path.temp) - 1)){
            eff.studies.temp <- eff.studies.temp + 1/studies[path.temp[m], path.temp[m + 1]]
            eff.sz.temp <- eff.sz.temp + 1/sz[path.temp[m], path.temp[m + 1]]
            eff.prec.temp <- eff.prec.temp + 1/prec[path.temp[m], path.temp[m + 1]]
          }
          S.eff.studies[k, k] <- eff.studies.temp
          S.eff.sz[k, k] <- eff.sz.temp
          S.eff.prec[k, k] <- eff.prec.temp
        }
        for(k in 1:(n.paths - 1)){
          for(m in (k + 1):n.paths){
            path1 <- as.numeric(paths[[k]])
            path2 <- as.numeric(paths[[m]])
            commp <- common.path(path1, path2)
            common <- commp$common
            pm <- commp$pm
            if(is.null(common)){
              S.eff.studies[k, m] <- S.eff.studies[m, k] <- 0
              S.eff.sz[k, m] <- S.eff.sz[m, k] <- 0
              S.eff.prec[k, m] <- S.eff.prec[m, k] <- 0
            }else{
              eff.studies.temp <- 0
              eff.sz.temp <- 0
              eff.prec.temp <- 0
              for(p in 1:length(common)){
                eff.studies.temp <- eff.studies.temp + pm[[p]]/studies[common[[p]][1], common[[p]][2]]
                eff.sz.temp <- eff.sz.temp + pm[[p]]/sz[common[[p]][1], common[[p]][2]]
                eff.prec.temp <- eff.prec.temp + pm[[p]]/prec[common[[p]][1], common[[p]][2]]
              }
              S.eff.studies[k, m] <- S.eff.studies[m, k] <- eff.studies.temp
              S.eff.sz[k, m] <- S.eff.sz[m, k] <- eff.sz.temp
              S.eff.prec[k, m] <- S.eff.prec[m, k] <- eff.prec.temp
            }
          }
        }

        eigen.studies <- eigen(S.eff.studies)
        eigen.studies.val <- eigen.studies$values
        if(all(eigen.studies.val > 10^(-6))){
          eff.studies[i, j] <- eff.studies[j, i] <- sum(solve(S.eff.studies))
        }else{
          pos.eigen <- sum(eigen.studies.val > 10^(-6))
          Q.pos <- eigen.studies$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.studies.val[1:pos.eigen])
          eff.studies[i, j] <- eff.studies[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }

        eigen.sz <- eigen(S.eff.sz)
        eigen.sz.val <- eigen.sz$values
        if(all(eigen.sz.val > 10^(-6))){
          eff.sz[i, j] <- eff.sz[j, i] <- sum(solve(S.eff.sz))
        }else{
          pos.eigen <- sum(eigen.sz.val > 10^(-6))
          Q.pos <- eigen.sz$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.sz.val[1:pos.eigen])
          eff.sz[i, j] <- eff.sz[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }

        eigen.prec <- eigen(S.eff.prec)
        eigen.prec.val <- eigen.prec$values
        if(all(eigen.prec.val > 10^(-6))){
          eff.prec[i, j] <- eff.prec[j, i] <- sum(solve(S.eff.prec))
        }else{
          pos.eigen <- sum(eigen.prec.val > 10^(-6))
          Q.pos <- eigen.prec$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.prec.val[1:pos.eigen])
          eff.prec[i, j] <- eff.prec[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }
      }
    }
  }
  out <- list(eff.studies = eff.studies, eff.sz = eff.sz, eff.prec = eff.prec)
  return(out)
}

#####################################################################################
## Overall evidence for continuous outcomes
eff.cont <- function(data){
  studies <- dir.cont(data)$studies
  sz <- dir.cont(data)$sz
  prec <- dir.cont(data)$prec
  conn <- (studies > 0)
  K <- dim(conn)[1]
  gra <- graph_from_adjacency_matrix(conn, mode = "undirected")

  eff.studies <- matrix(NA, K, K)
  eff.sz <- matrix(NA, K, K)
  eff.prec <- matrix(NA, K, K)

  for(i in 1:(K - 1)){
    for(j in (i + 1):K){
      paths <- all_simple_paths(gra, from = i, to = j)
      n.paths <- length(paths)
      if(n.paths == 1){
        path.temp <- paths[[1]]
        eff.studies.temp <- 0
        eff.sz.temp <- 0
        eff.prec.temp <- 0
        for(m in 1:(length(path.temp) - 1)){
          eff.studies.temp <- eff.studies.temp + 1/studies[path.temp[m], path.temp[m + 1]]
          eff.sz.temp <- eff.sz.temp + 1/sz[path.temp[m], path.temp[m + 1]]
          eff.prec.temp <- eff.prec.temp + 1/prec[path.temp[m], path.temp[m + 1]]
        }
        eff.studies[i, j] <- eff.studies[j, i] <- 1/eff.studies.temp
        eff.sz[i, j] <- eff.sz[j, i] <- 1/eff.sz.temp
        eff.prec[i, j] <- eff.prec[j, i] <- 1/eff.prec.temp
      }
      if(n.paths > 1){
        S.eff.studies <- matrix(NA, n.paths, n.paths)
        S.eff.sz <- matrix(NA, n.paths, n.paths)
        S.eff.prec <- matrix(NA, n.paths, n.paths)
        for(k in 1:n.paths){
          path.temp <- as.numeric(paths[[k]])
          eff.studies.temp <- 0
          eff.sz.temp <- 0
          eff.prec.temp <- 0
          for(m in 1:(length(path.temp) - 1)){
            eff.studies.temp <- eff.studies.temp + 1/studies[path.temp[m], path.temp[m + 1]]
            eff.sz.temp <- eff.sz.temp + 1/sz[path.temp[m], path.temp[m + 1]]
            eff.prec.temp <- eff.prec.temp + 1/prec[path.temp[m], path.temp[m + 1]]
          }
          S.eff.studies[k, k] <- eff.studies.temp
          S.eff.sz[k, k] <- eff.sz.temp
          S.eff.prec[k, k] <- eff.prec.temp
        }
        for(k in 1:(n.paths - 1)){
          for(m in (k + 1):n.paths){
            path1 <- as.numeric(paths[[k]])
            path2 <- as.numeric(paths[[m]])
            commp <- common.path(path1, path2)
            common <- commp$common
            pm <- commp$pm
            if(is.null(common)){
              S.eff.studies[k, m] <- S.eff.studies[m, k] <- 0
              S.eff.sz[k, m] <- S.eff.sz[m, k] <- 0
              S.eff.prec[k, m] <- S.eff.prec[m, k] <- 0
            }else{
              eff.studies.temp <- 0
              eff.sz.temp <- 0
              eff.prec.temp <- 0
              for(p in 1:length(common)){
                eff.studies.temp <- eff.studies.temp + pm[[p]]/studies[common[[p]][1], common[[p]][2]]
                eff.sz.temp <- eff.sz.temp + pm[[p]]/sz[common[[p]][1], common[[p]][2]]
                eff.prec.temp <- eff.prec.temp + pm[[p]]/prec[common[[p]][1], common[[p]][2]]
              }
              S.eff.studies[k, m] <- S.eff.studies[m, k] <- eff.studies.temp
              S.eff.sz[k, m] <- S.eff.sz[m, k] <- eff.sz.temp
              S.eff.prec[k, m] <- S.eff.prec[m, k] <- eff.prec.temp
            }
          }
        }

        eigen.studies <- eigen(S.eff.studies)
        eigen.studies.val <- eigen.studies$values
        if(all(eigen.studies.val > 10^(-6))){
          eff.studies[i, j] <- eff.studies[j, i] <- sum(solve(S.eff.studies))
        }else{
          pos.eigen <- sum(eigen.studies.val > 10^(-6))
          Q.pos <- eigen.studies$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.studies.val[1:pos.eigen])
          eff.studies[i, j] <- eff.studies[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }

        eigen.sz <- eigen(S.eff.sz)
        eigen.sz.val <- eigen.sz$values
        if(all(eigen.sz.val > 10^(-6))){
          eff.sz[i, j] <- eff.sz[j, i] <- sum(solve(S.eff.sz))
        }else{
          pos.eigen <- sum(eigen.sz.val > 10^(-6))
          Q.pos <- eigen.sz$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.sz.val[1:pos.eigen])
          eff.sz[i, j] <- eff.sz[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }

        eigen.prec <- eigen(S.eff.prec)
        eigen.prec.val <- eigen.prec$values
        if(all(eigen.prec.val > 10^(-6))){
          eff.prec[i, j] <- eff.prec[j, i] <- sum(solve(S.eff.prec))
        }else{
          pos.eigen <- sum(eigen.prec.val > 10^(-6))
          Q.pos <- eigen.prec$vectors[,1:pos.eigen]
          Lambda.pos.inv <- diag(1/eigen.prec.val[1:pos.eigen])
          eff.prec[i, j] <- eff.prec[j, i] <- sum(Q.pos %*% Lambda.pos.inv %*% t(Q.pos))
        }
      }
    }
  }
  out <- list(eff.studies = eff.studies, eff.sz = eff.sz, eff.prec = eff.prec)
  return(out)
}

#####################################################################################
## Finding common treatment comparisons shared by two sources of evidence
common.path <- function(path1, path2){
  if(length(path1) > length(path2)){
    temp <- path1
    path1 <- path2
    path2 <- temp
  }
  len1 <- length(path1) - 1
  len2 <- length(path2) - 1
  common <- NULL
  pm <- NULL
  k <- 1
  for(i in 1:len1){
    from1 <- path1[i]
    to1 <- path1[i + 1]
    if(is.element(from1, path2)){
      idx <- which(path2 == from1)
      if(idx <= len2){
        if(path2[idx + 1] == to1){
          common[[k]] <- sort(c(from1, to1))
          pm[[k]] <- 1
          k <- k + 1
        }
      }
      if(idx > 1){
        if(path2[idx - 1] == to1){
          common[[k]] <- sort(c(from1, to1))
          pm[[k]] <- -1
          k <- k + 1
        }
      }
    }
  }
  out <- list(common = common, pm = pm)
  return(out)
}

#####################################################################################
## Plotting treatment networks with edge width proportional to given weights
## wt.mat: a K x K matrix; each element is the weight of the corresponding comparison
## wt.min: the minimum weight that corresponds to the thinnest edge
## other arguments adjust sizes and colors of nodes, edges, and title
networkplot <- function(wt.mat, wt.min, title = "", cex.main = 1, adjust.title = 1,
  trtname, adjust.thick = 5, weight.node = FALSE, adjust.node.size = 10,
  node.col = "orange", edge.col = "black", text.cex = 1,
  adjust.figsizex = 1.1, adjust.figsizey = 1.1){
  if((dim(wt.mat)[1] != dim(wt.mat)[2]) | any(wt.mat - t(wt.mat) > 10^(-6), na.rm = TRUE)){
    stop("wt.mat is not a symmetric matrix.")
  }
  diag(wt.mat) <- 0
  K <- dim(wt.mat)[1]
  if(K <= 2) stop("there are less than 3 treatments, no need for network plot.")
  if(missing(trtname)){
    trtname <- 1:K
  }

  polar <- pi/2 - 2*pi/K*(0:(K - 1))
  x <- cos(polar)
  y <- sin(polar)
  plot(x, y, axes = FALSE, xlab = "", ylab = "", cex = 0.1,
    xlim = c(-adjust.figsizex, adjust.figsizex),
    ylim = c(-adjust.figsizey, adjust.figsizey))
  title(main = title, cex.main = cex.main, font.main = 1, line = adjust.title)

  wt <- c(wt.mat)
  wt.unique <- unique(wt[wt > 0])
  if(missing(wt.min)){
    wt.min <- min(wt.unique)
  }
  wt.max <- max(wt.unique)
  if(wt.min < wt.max){
    wt[wt > 0] <- round(1 + adjust.thick*(wt[wt > 0] - wt.min)/(wt.max - wt.min))
  }else{
    wt[wt > 0] <- 2
  }

  wt <- matrix(wt, K, K)

  for(t1 in 2:K){
    for(t2 in 1:(t1 - 1)){
      if(t1 != t2 & wt[t1, t2] > 0){
        lines(x = x[c(t1, t2)], y = y[c(t1, t2)], lwd = wt[t1, t2], col = edge.col)
      }
    }
  }

  if(weight.node){
    wt <- wt.mat
    wt <- colSums(wt)
    node.sizes <- 3 + (wt - min(wt))/(max(wt) - min(wt))*adjust.node.size
    points(x, y, pch = 20, cex = node.sizes, col = node.col)
  }else{
    points(x, y, pch = 20, cex = 3, col = node.col)
  }

  sides <- numeric(K)
  eps <- 10^(-4)
  for(t in 1:K){
    if((polar[t] <= pi/2 & polar[t] > pi/4) | (polar[t] < -5*pi/4 & polar[t] >= -3*pi/2)){
      sides[t] <- 3
    }
    if(polar[t] <= pi/4 & polar[t] >= -pi/4){
      sides[t] <- 4
    }
    if(polar[t] < -pi/4 & polar[t] > -3*pi/4){
      sides[t] <- 1
    }
    if(polar[t] <= -3*pi/4 & polar[t] >= -5*pi/4){
      sides[t] <- 2
    }
  }
  for(t in 1:K){
    if(weight.node){
      text(x = x[t], y = y[t], labels = trtname[t], cex = text.cex)
    }else{
      text(x = x[t], y = y[t], labels = trtname[t], pos = sides[t], cex = text.cex)
    }
  }
}