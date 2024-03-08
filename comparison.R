#############
# Faithful graphical representations of local independence
# comparison of algorithms
# source code is below

library(FIAR)
library(ggplot2)
library(igraph)

set.seed(595959595)

### parameters

n.v <- c(5,7,9) # number of coordinate processes
M <- 100 # number of repetitions
N <- 100 # length of observed time series
p.v <- c(0.001,0.01,0.05,0.10,0.20) # significance level
sigma <- .5 # noise level

# algorithms
# CA (ca)
# dSGS (dS)
# parent screening (pa)
# causal min (cm)

ll <- list()

t0 <- proc.time()[3]
for (n in n.v){
  
  nodes <- letters[1:n]  

    ll[[length(ll) + 1]] <- simfunc(n, M, N, p.v, sigma)

  print(n)
}
t1 <- proc.time()[3]

#save(ll, file = "comparison")

res <- data.frame(n = integer(), M = integer(), N = integer(),
                  sigma = double(), G0edges = integer(),
                  extra = integer(), missing = integer(), nTests = integer(),
                  p = double(), type = character())

for (i in 1:length(ll)){
  
  param <- ll[[i]][[1]]
  
  for (j in 2:length(ll[[i]])){
    
    gg <- ll[[i]][[j]]
    G0 <- gg[["G0"]]
    G0edges <- sum(G0)
    
    g1 <- gg[["outCA"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                    gg[["outCA"]]$nTests, NA))
    res[nrow(res),ncol(res)] <- "CA"
    res[nrow(res),ncol(res)-1] <- gg[["outCA"]]$p
    
    g1 <- gg[["outdS"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                    NA, NA))
    res[nrow(res),ncol(res)] <- "dSGS"
    res[nrow(res),ncol(res)-1] <- gg[["outdS"]]$p
    
    g1 <- gg[["outPA"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                    NA, NA))
    res[nrow(res),ncol(res)] <- "PA"
    res[nrow(res),ncol(res)-1] <- gg[["outPA"]]$p
    
    g1 <- gg[["outCM"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                    NA, NA))
    res[nrow(res),ncol(res)] <- "CM"
    res[nrow(res),ncol(res)-1] <- gg[["outCM"]]$p
    
  }
  
}

res$diff <- res$extra + res$missing

aa <- aggregate(res$diff, list(n = res$n, p = res$p, type = res$type), FUN = mean)
bb <- aggregate(res$extra, list(n = res$n, p = res$p, type = res$type), FUN = mean)
cc <- aggregate(res$missing, list(n = res$n, p = res$p, type = res$type), FUN = mean)

aa$metric <- "difference"
bb$metric <- "surplus"
cc$metric <- "missing"

dd <- rbind(aa, bb, cc)
dd$type[dd$type == "dSGS"] <- "dS"
dd$type[dd$type == "PA"] <- "CS"

ggplot(dd, aes(x = type, y = x, col = metric, group = metric, shape = metric)) + 
  geom_point(size = 2, position = position_dodge(width=0.3)) + 
  facet_grid(n ~ p, scales = "free_y",
             labeller = labeller(
    n = c(`5` = "n = 5", `7` = "n = 7", `9` = "n = 9"),
    p = c(`0.001` = "p = 0.001", `0.01` = "p = 0.01", `0.05` = "p = 0.05", `0.1` = "p = 0.1", `0.2` = "p = 0.2"),
    type = c(`dSGS` = "dS")
  )) +
  labs(color = NULL) + labs(shape = NULL) +
  xlab(NULL) + ylab(NULL) +
  theme(text=element_text(size=16))

ggsave("comparison.pdf", height = 5, width = 10)

# number of tests

ee <- res[res$type == "CA",]

aggregate(ee$nTests, list(ee$n), min)
aggregate(ee$nTests, list(ee$n), max)


########################
######################
# partial observation

rm(list = ls())

set.seed(1919191919)

n.v <- c(5) # number of coordinate processes
M <- 50 # number of repetitions
N <- 100 # length of observed time series
p.v <- c(0.001,0.01,0.05,0.10,0.20) # significance level
sigma <- .5 # noise level

nodes <- letters[1:n]

# algorithms
# CA (ca)
# dSGS (dS)
# parent screening (pa)
# causal min (cm)
# dFCI (dF)

ll <- list()

t0 <- proc.time()[3]
for (n in n.v){
  
  nodes <- letters[1:n]  
  
  ll[[length(ll) + 1]] <- simfuncPartial(n, M, N, p.v, sigma)
  
  print(n)
}
t1 <- proc.time()[3]

#save(ll, file = "comparisonPartial")

res <- data.frame(n = integer(), M = integer(), N = integer(),
                  sigma = double(), G0edges = integer(),
                  extra = integer(), missing = integer(), nTests = integer(),
                  p = double(), type = character())

for (i in 1:length(ll)){
  
  param <- ll[[i]][[1]]
  
  for (j in 2:length(ll[[i]])){
    
    gg <- ll[[i]][[j]]
    G0 <- gg[["G0"]]
    G0edges <- sum(G0)
    
    g1 <- gg[["outCA"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                              gg[["outCA"]]$nTests, NA))
    res[nrow(res),ncol(res)] <- "CA"
    res[nrow(res),ncol(res)-1] <- gg[["outCA"]]$p
    
    g1 <- gg[["outdS"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                              NA, NA))
    res[nrow(res),ncol(res)] <- "dSGS"
    res[nrow(res),ncol(res)-1] <- gg[["outdS"]]$p
    
    g1 <- gg[["outPA"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                              NA, NA))
    res[nrow(res),ncol(res)] <- "PA"
    res[nrow(res),ncol(res)-1] <- gg[["outPA"]]$p
    
    g1 <- gg[["outCM"]]$G - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                              NA, NA))
    res[nrow(res),ncol(res)] <- "CM"
    res[nrow(res),ncol(res)-1] <- gg[["outCM"]]$p
    
    g1 <- gg[["outdFCI"]]$G$M1 - G0
    res[nrow(res)+1,-ncol(res)] <- c(param, c(G0edges, sum(g1 > 0), sum(g1 < 0),
                                              gg[["outdFCI"]]$notests, NA))
    res[nrow(res),ncol(res)] <- "dF"
    res[nrow(res),ncol(res)-1] <- gg[["outdFCI"]]$p
    
  }
  
}

res$diff <- res$extra + res$missing

aa <- aggregate(res$diff, list(n = res$n, p = res$p, type = res$type), FUN = mean)
bb <- aggregate(res$extra, list(n = res$n, p = res$p, type = res$type), FUN = mean)
cc <- aggregate(res$missing, list(n = res$n, p = res$p, type = res$type), FUN = mean)

aa$metric <- "difference"
bb$metric <- "surplus"
cc$metric <- "missing"

dd <- rbind(aa, bb, cc)
dd$type[dd$type == "dSGS"] <- "dS"
dd$type[dd$type == "PA"] <- "CS"

ggplot(dd, aes(x = type, y = x, col = metric, group = metric, shape = metric)) + 
  geom_point(size = 2, position = position_dodge(width=0.3)) + 
  facet_grid(n ~ p, scales = "free_y",
             labeller = labeller(
               n = c(`5` = "n = 5", `7` = "n = 7", `9` = "n = 9"),
               p = c(`0.001` = "p = 0.001", `0.01` = "p = 0.01", `0.05` = "p = 0.05", `0.1` = "p = 0.1", `0.2` = "p = 0.2"),
               type = c(`dSGS` = "dS")
             )) +
  labs(color = NULL) + labs(shape = NULL) +
  xlab(NULL) + ylab(NULL) +
  theme(text=element_text(size=16))

#ggsave("comparisonPartial.pdf", height = 2.5, width = 10)

# number of tests

ee <- res[res$type == "CA",]

aggregate(ee$nTests, list(ee$n), min)
aggregate(ee$nTests, list(ee$n), max)


#################
# source


pcParent <- function(nodes, testFunc, thres){
  
  nNodes <- length(nodes)
  S <- list(M1 = matrix(1, nrow = nNodes, ncol = nNodes, dimnames = list(nodes, nodes)),
            M2 = matrix(1, nrow = nNodes, ncol = nNodes, dimnames = list(nodes, nodes)))

  sep.m <- S$M1
  sep.m[] <- NA
  
  tests.l <- list()
  n <- 0
  
  while (n <= max(colSums(S$M1))){
    for (i in nodes){
      for (j in setdiff(nodes, i)){
        if (S$M1[i,j] == 0) next
        Osub <- setdiff(nodes[which(S$M1[,j] == 1)], c(i,j))
        if (length(Osub) < n) next
        for (C in combn(Osub, n, simplify=FALSE)){
          sep.bool <- FALSE
          sep.bool <- testFunc(i, j, setdiff(union(C, j), i)) > thres
          tests.l[[length(tests.l) + 1]] <- paste(paste(i, j, sep = "-"), paste(setdiff(union(C, j), i), collapse = ""), sep = "-")
          if (sep.bool){
            S$M1[i,j] <- S$M2[i,j] <- S$M2[j,i] <- 0
            sep.m[i,j] <- paste(setdiff(union(C, j), i), collapse = "")
            break
          }
        }
      }
    }
    n <- n + 1
  }
  
  
  list(S = S, sep.m = sep.m, tests = tests.l)
}

###


learndFCI <- function(nodes, testFunc, thres) {
  
  nNodes <- length(nodes)
  G <- S <- list(M1 = matrix(1, nrow = nNodes, ncol = nNodes, dimnames = list(nodes, nodes)),
                M2 = matrix(1, nrow = nNodes, ncol = nNodes, dimnames = list(nodes, nodes)))
  
  sep.m <- G$M1
  

  
  ###
  # PC step
  
    pc.l <- pcParent(nodes, testFunc, thres)
    sep.m <- pc.l$sep.m
    tests.l <- pc.l[[3]]
    S <- pc.l$S

  ###
  # FCI step
  # In the oracle case, S is a supergraph of the true graph, G0,
  # then D_G0(alpha,beta) is a subset of D_S(alpha,beta)

    for (i in nodes){
      for (j in setdiff(nodes,i)){
        if (S$M1[i,j] == 0) next
        
        C0 <- findD(S, i, j)
        if (length(C0) == 0) next # this is handled by PC-like step already
        
        for (m in length(C0):1){
          for (C in combn(C0,m,simplify=FALSE)){
            sep.bool <- FALSE
            sep.bool <- testFunc(i, j, C) > thres
            tests.l[[length(tests.l) + 1]] <- paste(paste(i, j, sep = "-"), paste(C, collapse = ""), sep = "-")
            if (sep.bool){
              S$M1[i,j] <- S$M2[i,j] <- S$M2[j,i] <- 0
              sep.m[i,j] <- paste(C, collapse = "")
              break
            }
            if (sep.bool)
              break
          }
        }
      }
    }

G <- S
  
    
    ###
    # define G^c to mark those certainly in N
    certainG <- G
    certainG$M1[] <- certainG$M2[] <- 0
    
    ###
    # loop over each edge in G
    
    
    # order by "connectivity"
    
    tmp <- G
    diag(tmp$M1) <- diag(tmp$M2) <- 0
    ord <- rowSums(tmp$M1) + colSums(tmp$M1) + rowSums(tmp$M2)
    
    for (alpha in nodes){
      for (beta in nodes){
        
        # check all known dependences, if removing e violates one of them, then leave e in
        Gdir <- Gbidir <- G
        Gdir$M1[alpha,beta] <- 0
        Gbidir$M2[alpha,beta] <- Gbidir$M2[beta,alpha] <- 0
        
        dirSep <- FALSE
        bidirSep <- FALSE
        
        if (length(tests.l) > 0){
          u.l <- unique(lapply(tests.l, function(x) x[[1]]))
          for (li in 1:length(u.l)){
            
            te <- u.l[[li]]
            tmp <- strsplit(te, split = "-")[[1]]
            if (length(tmp) == 2) tmp <- c(tmp, "")
            mDir <- muSep(Gdir, tmp[[1]], tmp[[2]], 
                          strsplit(tmp[[3]], "")[[1]])
            mBidir <- muSep(Gbidir, tmp[[1]], tmp[[2]], 
                            strsplit(tmp[[3]], "")[[1]])
            m0 <- testFunc(tmp[[1]], tmp[[2]], 
                        strsplit(tmp[[3]], "")[[1]]) > thres
            if (!m0 && mDir) certainG$M1[alpha,beta] <- 1
            if (!m0 && mBidir) certainG$M2[alpha,beta] <- certainG$M2[beta,alpha] <- 1
            
          }
        }
        
        
        # consider those that are possibly in N, remove all of them, and check for Markov equivalence
        # repeat this after each PP iteration
        
        Gminus <- G
        Gminus$M1[certainG$M1 == 0] <- 0
        Gminus$M2[certainG$M2 == 0] <- 0
        
        if (!areME(G, Gminus)){
          
          
          # check that these are still uncertain    
          
          breakBool <- FALSE
          
          #potential sibling 
          
          if (G$M2[alpha,beta] == 1 &
              certainG$M2[alpha,beta] == 0){ # testing the third condition
            # checking that there's inseparability both ways
            
            # this should add to certainG$M2 too, but needs to do it both ways, then
            
            for (delta in nodes){ #setdiff(nodes, c(alpha,beta))){
              for (C in all.subsets(setdiff(nodes, delta))){
                if (!(alpha %in% C)) next
                if (!muSep(G, delta, alpha, C)){
                  if (delta == alpha) next
                  if (delta == beta) next
                  if (testFunc(delta, alpha, C) < thres && testFunc(delta, beta, C) > thres) {
                    G$M2[alpha, beta] <- G$M2[beta, alpha] <- 0
                    breakBool <- TRUE
                    tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                    tests.l[[length(tests.l) + 1]] <- paste(paste(delta, beta, sep = "-"), paste(C, collapse = ""), sep = "-")
                    break
                  }
                  tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                  tests.l[[length(tests.l) + 1]] <- paste(paste(delta, beta, sep = "-"), paste(C, collapse = ""), sep = "-")
                }
                if (breakBool) break
              }
              if (breakBool) break
            }
            
            
          }
          
          breakBool <- FALSE
          
          #############
          #potential parent
          
            if (G$M1[alpha,beta] == 1 & 
                certainG$M1[alpha,beta] == 0){ # testing the third condition
              
              for (delta in nodes){
                for (C in all.subsets(setdiff(nodes, delta))){
                  if (alpha %in% C) next
                  
                  #PP2
                  if (!muSep(G, delta, alpha, C)){
                    if (delta == alpha) next
                    if (delta == beta) next
                    if (testFunc(delta, alpha, C) < thres && testFunc(delta, beta, C) > thres) { 
                      G$M1[alpha, beta] <- 0
                      breakBool <- TRUE
                      tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                      tests.l[[length(tests.l) + 1]] <- paste(paste(delta, beta, sep = "-"), paste(C, collapse = ""), sep = "-")
                      break
                    }
                    tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                    tests.l[[length(tests.l) + 1]] <- paste(paste(delta, beta, sep = "-"), paste(C, collapse = ""), sep = "-")
                  }
                  
                  #PP3
                  if (!(alpha %in% C) && (beta %in% C)){
                    for (epsilon in nodes){ # setdiff(nodes, c(alpha, beta))){
                      if (!muSep(G, delta, beta, setdiff(union(C, epsilon), delta)) &&
                          !muSep(G, delta, epsilon, C)){
                        #!testFunc(G0, delta, beta, setdiff(C, epsilon)) && 
                        if (delta == beta) next
                        if (alpha == epsilon) next
                        if (delta == epsilon) next
                        if (testFunc(delta, beta, C) < thres && 
                            testFunc(alpha, epsilon, C) < thres && testFunc(delta, epsilon, C) > thres){
                          G$M1[alpha, beta] <- 0
                          tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                          tests.l[[length(tests.l) + 1]] <- paste(paste(alpha, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                          tests.l[[length(tests.l) + 1]] <- paste(paste(delta, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                          breakBool <- TRUE
                          break
                        }
                        tests.l[[length(tests.l) + 1]] <- paste(paste(delta, alpha, sep = "-"), paste(C, collapse = ""), sep = "-")
                        tests.l[[length(tests.l) + 1]] <- paste(paste(alpha, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                        tests.l[[length(tests.l) + 1]] <- paste(paste(delta, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                      }
                      if (breakBool) break
                    }
                  }
                  
                  #PP4
                  if (!(beta %in% C)){ # note that delta is not used here - not very efficient
                    for (epsilon in nodes){ # setdiff(nodes, c(alpha, beta))){
                      if (!muSep(G, alpha, epsilon, C)){
                        if (alpha == epsilon) next
                        if (beta == epsilon) next
                        if (testFunc(alpha, epsilon, C) < thres && testFunc(beta, epsilon, C) > thres) { 
                          G$M1[alpha, beta] <- 0
                          breakBool <- TRUE
                          tests.l[[length(tests.l) + 1]] <- paste(paste(alpha, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                          tests.l[[length(tests.l) + 1]] <- paste(paste(beta, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                          break
                        }
                        tests.l[[length(tests.l) + 1]] <- paste(paste(alpha, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                        tests.l[[length(tests.l) + 1]] <- paste(paste(beta, epsilon, sep = "-"), paste(C, collapse = ""), sep = "-")
                      }
                      if (breakBool) break
                    }
                  }
                  
                }
                if (breakBool) break
              }
              
              # updating certainG, alpha -> beta is in N
              certainG$M1[alpha,beta] <- 1
              
            }
          
          # alpha, beta loops end:
        }
        
      }
    }
  
  
  ######
  # make output
  
  # count distinct tests
  notests <- length(unique(lapply(tests.l, function(x) x[[1]])))
  
  
  list(G = G, tests = tests.l, notests = notests, p = thres)
}


###################################
# function that returns a list of all subsets of the input set, C

all.subsets <- function(C){
  
  set.l <- list()
  
  for (k in 0:length(C)){
    for (Z in combn(C,k,simplify=FALSE)){
      set.l[[length(set.l) + 1]] <- Z
    }
  }
  
  set.l
}


###################################
# determine if two graphs are Markov equivalent
# if O is provided, then only equality on O is checked (e.g. for checking marginalizations)

areME <- function(G1, G2, verbose = FALSE, O = NULL){
  
  if (!isTRUE(all.equal(names(G1), c("M1", "M2"))))
    stop("G1 is not a DMG object.")
  if (!isTRUE(all.equal(names(G2), c("M1", "M2"))))
    stop("G2 is not a DMG object.")
  if (!isTRUE(all.equal(rownames(G1$M1), rownames(G2$M1))))
    warning("Node names provided differ between the two graphs.")
  
  nodes <- union(rownames(G1$M1), rownames(G2$M1))
  mu1 <- mu2 <- c()
  res.l <- list()
  
  if (is.null(O)){
    checkNodes <- nodes
  } else {
    checkNodes <- O
  }
  
  
  for (alpha in checkNodes) {
    for (beta in checkNodes) {
      C <- setdiff(checkNodes, alpha)
      tmp1 <- muSep(G1, alpha, beta, c())
      tmp2 <- muSep(G2, alpha, beta, c())
      mu1 <- c(mu1, tmp1)
      mu2 <- c(mu2, tmp2)
      if (tmp1 != tmp2)
        res.l[[length(res.l) + 1]] <- paste(alpha, "-", beta, "-", paste("", collapse = ""), "-", tmp1, tmp2, collapse = " ")
      if (length(C) > 0){
        for (k in 1:length(C)){
          for (Z in combn(C, k, simplify=FALSE)){
            tmp1 <- muSep(G1, alpha, beta, Z)
            tmp2 <- muSep(G2, alpha, beta, Z)
            mu1 <- c(mu1, tmp1)
            mu2 <- c(mu2, tmp2)
            if (tmp1 != tmp2)
              res.l[[length(res.l) + 1]] <- paste(alpha, "-", beta, "-", paste(Z, collapse = ""), "-", tmp1, tmp2, collapse = " ")
          }
        }
      }
      
    }
  }
  
  if (verbose) {
    return(res.l)
  } else {
    isTRUE(all.equal(mu1, mu2))
  }
}



###################################
# unDirSep
# M is directed adjacency matrix

unDirSep <- function(M, A, B, C){
  if (length(C) > 0) {
    tmp <- which(rownames(M) %in% C)
    sepM <- M[-tmp, -tmp, drop = FALSE]
  } else {
    sepM <- M
  }
  
  !max(ancM(sepM)[A, B])
}

testUnDirSep <- function(nrep, nNodes = 7){
  
  res.v <- rep(NA, nrep)
  
  for (i in 1:nrep){
    
    G <- createDMG(nNodes)
    M <- G$M2
    
    ss <- sample(rownames(M), sample(2:nNodes, 1))
    cut1 <- sample(1:(length(ss)-1), 1)
    cut2 <- cut1 + sample(1:(length(ss)-cut1), 1)
    A <- ss[1:cut1]
    B <- ss[(cut1+1):cut2]
    C <- setdiff(ss, union(A,B))
    
    
    
    uds <- unDirSep(M, A, B, C)
    
    # a hack when C is empty as RGBL::separates does not accept this
    if (length(C) == 0) {
      
      names <- rownames(M)
      M <- cbind(M, 0)
      M <- rbind(M, 0)
      rownames(M) <- colnames(M) <- c(names, "zz")
      g <- graph_from_adjacency_matrix(M)
      sep <- separates(A, B, "zz", igraph.to.graphNEL(g))
      
    } else {
      g <- graph_from_adjacency_matrix(M)
      sep <- separates(A, B, C, igraph.to.graphNEL(g))
    }
    
    
    
    
    res.v[i] <-  uds == sep
    if (!res.v[i]) return(list(G = G, A = A, B = B, C = C,
                               uds, sep))
  }
  
  res.v
  
}

###################################
# muSep
# mSep: use m-separation (instead of mu-separation)

muSep <- function(G, A, B, C, ..., mSep = FALSE){
  
  M1 <- G$M1
  M2 <- G$M2
  Bp <- paste(B, "p", sep = "")
  
  A <- setdiff(A,C)
  
  if (length(A) == 0)
    return(TRUE)
  
  if (!mSep) {
    bH <- bHist(G, B)
  } else {
    bH <- G
    Bp <- B
    
    B <- setdiff(B,C)
    
    if (length(B) == 0)
      return(TRUE)
  }
  
  # find G(B)_An(A\cup Bp \cup C)
  An <- ancM(bH$M1)
  isAn <- rowSums(An[,c(A,Bp,C),drop = FALSE]) > 0
  M1An <- bH$M1[isAn, isAn, drop = FALSE]
  M2An <- bH$M2[isAn, isAn, drop = FALSE]
  
  # find the augmented graph
  g <- augGraph(list(M1 = M1An, M2 = M2An))
  
  unDirSep(g, A, Bp, C)
  
}


###################################
# create B-history version of a DMG
#

bHist <- function(G, B){
  
  nB <- length(B)
  if (nB == 0) stop("Empty B-set provided.")
  
  n <- nrow(G$M1)
  M1 <- M2 <- matrix(0, nrow = n + nB, ncol = n + nB)
  M1[1:n, 1:n] <- G$M1
  M2[1:n, 1:n] <- G$M2
  M1[1:n, (n+1):(n+nB)] <- G$M1[,B]
  M2[1:n, (n+1):(n+nB)] <- G$M2[,B]
  M2[(n+1):(n+nB), 1:n] <- G$M2[B,]
  
  rownames(M1) <- colnames(M1) <- rownames(M2) <- colnames(M2) <- c(rownames(G$M1), paste(B, "p", sep = ""))
  list(M1 = M1, M2 = M2)
}


###################################
# trimS
#

### use W-structures to trim down S

trimS <- function(nodes, S, sep.m){
  
  # iterate over separable pairs, alpha->beta
  # look for W-structures, (alpha,gamma,beta)
  
  N <- S
  
  for (alpha in nodes){
    for (beta in nodes){
      if (alpha == beta) next
      if (is.na(sep.m[alpha,beta])) next
      sep.set <- strsplit(sep.m[alpha,beta], split = "")[[1]]
      for (gamma in setdiff(nodes, beta)){
        # in S (separability matrix) there's never a bidirected without a directed each way, therefore this is enough
        if (S$M1[alpha,gamma] == 0) next
        # we should keep evaluating the above in S, not in the updated S, or we might miss some
        if (!(gamma %in% sep.set)) N$M1[gamma,beta] <- 0
        if (gamma %in% sep.set) N$M2[gamma,beta] <- N$M2[beta,gamma] <- 0
      }
    }
  }
  
  N
}

#################
# function for matrix exponentiation (the function "truncates" entry-wise at 1 to avoid overflow)
# m: matrix
# k: power

mexp <- function(m, k){
  
  if (k == 0) return(diag(nrow(m)))
  count <- 1
  m0 <- m
  while (count < k){
    m0 <- m0%*%m
    m0 <- 1*(m0 > 0)
    count <- count + 1
  }
  m0
} 



###################################
# find ancestors
# M is a directed adjacency matrix (directed part of a DMG)
# output: an "ancestor" matrix, M: M(i,j) = 1 iff i is an ancestor of j

ancM <- function(M) {
  
  Man <- M
  Man[] <- 0
  
  if (sum(M) > 0) {
    for (k in 1:sum(M)) {
      Man <- Man + mexp(M,k)
    }
  }
  
  Man <- 1*(Man > 0)
  diag(Man) <- 1
  Man 
}


###
##############
# function to find D(alpha,beta)-set from graph
# directedly: use definition with directedly collider connected (TRUE), 
# or just collider connected

findD <- function(G, alpha, beta, directedly = TRUE){
  
  nodes <- rownames(G$M1)
  if (directedly) G$M1[beta,] <- 0
  v1 <- augGraph(G)[,beta]
  v2 <- apply(ancM(G$M1)[,c(alpha,beta)], 1, max)
  C0 <- nodes[v1 + v2 == 2]
  setdiff(union(C0, beta), alpha)
  
}

###################################
# function that returns the moral graph of a DMG,
# that is, the symmetric adjacency of an undirected graph that has an edge 
# between alpha and beta if and only if there's a collider path of length <=2 connecting them
# (alpha = beta allowed, that is, a single edge e: alpha *-> beta, can constitute a collider path:
# alpha, e, beta, e, alpha)

moralGraph <- function(G, noWarning = FALSE){
  
  M1 <- G$M1
  M2 <- G$M2
  if (is.null(rownames(G$M1))) stop("Need named nodes.")
  if (!noWarning)
    if (!isTRUE(all.equal(sum(abs(M2)), 0))) warning("Should not be used on (non-LIG) MIGs")
  nodes <- rownames(G$M1)
  
  # find (asymmetric) adjacency matrix
  M0 <- 1*(M1 + M2 > 0)
  # find (symmetric) adjacency matrix
  M <- 1*(M0 + t(M0) > 0)
  
  # go through each column of M0, two non-zero entries in a column means a collider
  for (alpha in nodes){
    for (r1 in nodes){
      for (r2 in nodes){
        if (M0[r1,alpha] == 1 & M0[r2,alpha] == 1) M[r1,r2] <- M[r2,r1] <- 1
      }
    }
  }
  
  M
}


###################################
# create augmented graph (normally of a B-history graph)
# (self-edges are included in the moral graph which is used in augGraph)

augGraph <- function(G){
  
  M1 <- G$M1
  M2 <- G$M2
  if (is.null(rownames(M1))) stop("Need named nodes.")
  V <- rownames(M1)
  
  # start from the moralGraph
  moralG <- moralGraph(G, noWarning = TRUE)
  
  # then add the collider connected nodes (collision path of length > 2)
  for (i in 1:length(V)){
    for (j in i:length(V)){
      alpha <- V[i]
      beta <- V[j]
      if (moralG[alpha,beta] > 0) next
      M1tmp <- M1
      # remove the other directed arrows, can never be in a collider path of length > 2 (length = 2 is already in the moral graph)
      M1tmp[-which(V %in% c(alpha,beta)), ] <- 0
      M0tmp <- 1*(M1tmp + t(M1tmp) + M2 > 0)
      Mconn <- M0tmp
      Mconn[] <- 0
      for (k in 1:sum(M0tmp)) {
        Mconn <- Mconn + mexp(M0tmp,k)
      }
      if (Mconn[alpha,beta] > 0) moralG[alpha,beta] <- moralG[beta,alpha] <- 1
    }
  }
  
  moralG
}

###

simfuncPartial <- function(n, M, N, p.v, sigma){
  
  ll <- list()
  ll[["param"]] <- c(n = n, M = M, N = N, sigma = sigma)

  
for (i in 1:M){
  
  # generate graph (connected)
  notcon <- TRUE
  nn <- sample(0:n, size = 1) # number of unobserved coordinate processes
  while(notcon){
    pr <- runif(1)/2
    G0 <- matrix(sample(0:1, size = (n+nn)^2, replace = TRUE, prob = c(1-pr, pr)), nrow = n+nn)
    diag(G0) <- 1 # always has all directed loops
    
    g <- graph_from_adjacency_matrix(G0)
    notcon <- !is.connected(g)
  }
  
  nodes <- letters[1:(n+nn)]
  rownames(G0) <- colnames(G0) <- nodes
  
  # generate data and test results
  X <- genVAR(N = N, G = G0, sigma = sigma)
  GG <- margDMG(list(M1 = G0, M2 = matrix(0, nrow = n + nn, ncol = n + nn, dimnames = list(nodes, nodes))), nodes[1:n]) # marginalization
  G0 <- findMaxBrute(GG)$M1
  X <- X[,1:n] # marginalization
  nodes <- nodes[1:n]
  at <- allTests(X, k = n, order = 2)
  as <- allSubsets(n, k = n)
  
  
  testFunc <- function(A, B, C){
    V <- nodes
    li <- V %in% C
    l <- which.min(colSums(abs(t(as) - li)))
    at[A,B,l]
  }

  
  for (p in p.v){

    lli <- list()  
    lli[["G0"]] <- G0
    
  # learning algorithms
  
  # CA (ca)
  outCA <- learnCA(nodes, testFunc = testFunc, p)
  lli[["outCA"]] <- outCA
  
  # dSGS (dS)
  outdS <- learnSGS(nodes, at, p)
  lli[["outdS"]] <- outdS
  
  # parent screening (pa)
  outPA <- learnParent(nodes, testFunc, p)
  lli[["outPA"]] <- outPA
  
  # causal min (cm)
  outCM <- learnMin(nodes, testFunc, p)
  lli[["outCM"]] <- outCM
  
  # dFCI (dF)
  outdFCI <- learndFCI(nodes, testFunc, p)
  lli[["outdFCI"]] <- outdFCI
  
  ll[[length(ll) + 1]] <- lli
  
  }
  
  
}
  
  ll
}

######################
# function that can marginalize a DMG
# G, a DMG on nodes V
# O, an observed subset of V

margDMG <- function(G, O){
  
  V <- rownames(G$M1)
  M <- setdiff(V, O)
  
  Gold <- G
  Gold$M1[] <- Gold$M2[] <- 0
  
  while (!isTRUE(all.equal(G, Gold))){
    Gold <- G
    
    for (m in M){
      
      for (alpha in setdiff(V, m)) {
        for (beta in setdiff(V, m)) {
          if (G$M1[alpha, m] == 1 && G$M1[m, beta] == 1)
            G$M1[alpha, beta] <- 1
          if (G$M2[alpha, m] == 1 && G$M1[m, beta] == 1)
            G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
          if (G$M1[m, alpha] == 1 && G$M1[m, beta] == 1)
            G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
        }
      }
    }
    
  }
  
  list(M1 = (G$M1)[O,O, drop = FALSE], M2 = (G$M2)[O,O, drop = FALSE])
  
}

###

simfunc <- function(n, M, N, p.v, sigma){
  
  ll <- list()
  ll[["param"]] <- c(n = n, M = M, N = N, sigma = sigma)
  
  
  for (i in 1:M){
    
    
    
    # generate graph (connected)
    notcon <- TRUE
    while(notcon){
      pr <- runif(1)/2
      G0 <- matrix(sample(0:1, size = n^2, replace = TRUE, prob = c(1-pr, pr)), nrow = n)
      diag(G0) <- 1 # always has all directed loops
      
      g <- graph_from_adjacency_matrix(G0)
      notcon <- !is.connected(g)
    }
    
    nodes <- letters[1:n]
    rownames(G0) <- colnames(G0) <- nodes
    
    # generate data and test results
    X <- genVAR(N = N, G = G0, sigma = sigma)
    at <- allTests(X, k = n, order = 2)
    as <- allSubsets(n, k = n)
    
    
    testFunc <- function(A, B, C){
      V <- nodes
      li <- V %in% C
      l <- which.min(colSums(abs(t(as) - li)))
      at[A,B,l]
    }
    
    for (p in p.v){
      
      lli <- list()  
      lli[["G0"]] <- G0
      
      # learning algorithms
      
      # CA (ca)
      outCA <- learnCA(nodes, testFunc = testFunc, p)
      lli[["outCA"]] <- outCA
      
      # dSGS (dS)
      outdS <- learnSGS(nodes, at, p)
      lli[["outdS"]] <- outdS
      
      # parent screening (pa)
      outPA <- learnParent(nodes, testFunc, p)
      lli[["outPA"]] <- outPA
      
      # causal min (cm)
      outCM <- learnMin(nodes, testFunc, p)
      lli[["outCM"]] <- outCM
      
      ll[[length(ll) + 1]] <- lli
      
    }
    
    
  }
  
  ll
}

###############
# parent learning
# uses at most two tests for each pair (alpha,beta)

learnParent <- function(nodes, testFunc, thres){
  
  n <- length(nodes)
  G <- matrix(1, nrow = n, ncol = n)
  rownames(G) <- colnames(G) <- nodes

  # Subalgorithm 1
  for (alpha in nodes){
    for (beta in setdiff(nodes, alpha)){ # beta \neq alpha
      
      if (testFunc(alpha, beta, beta) > thres){
        G[alpha, beta] <- 0
      }
      
    }
  }
  

  
 
  # Subalgorithm 2
  for (alpha in nodes){ 
    for (beta in setdiff(nodes, alpha)){
      if(G[alpha, beta] == 0) next
      
      paBeta <- nodes[G[ ,beta] == 1]
      if (testFunc(alpha, beta, setdiff(paBeta, alpha)) > thres){
        G[alpha, beta] <- 0
      }
      
    }
  }
  
  list(G = G, p = thres)
}



#########
# causal minimality learning (condition on V\setminus alpha)
# uses only a single test for each pair (alpha,beta)

learnMin <- function(nodes, testFunc, thres){
  
  n <- length(nodes)
  
  G <- matrix(1, nrow = n, ncol = n)
  rownames(G) <- colnames(G) <- nodes
  
  for (alpha in nodes){
    for (beta in setdiff(nodes, alpha)){
      
      if (testFunc(alpha, beta, setdiff(nodes,alpha)) > thres)
        G[alpha,beta] <- 0
      
    }
  }
  
  
  list(G = G, p = thres)
}




#########
# dSGS algorithm
# uses all possible tests

learnSGS <- function(nodes, at, thres){


  n <- length(nodes)
  ii <- 0
  
  G <- matrix(1, nrow = n, ncol = n)
  rownames(G) <- colnames(G) <- nodes
  
  for (alpha in nodes){
    for (beta in setdiff(nodes, alpha)){
      
      if (max(na.omit(at[alpha,beta,])) > thres)
        G[alpha,beta] <- 0
      
    }
  }
  
  
  list(G = G, p = thres)
}



#########
# Meek's CA algorithm

learnCA <- function(nodes, testFunc, thres){
  
  n <- length(nodes)
  ii <- 0
  
  G <- matrix(1, nrow = n, ncol = n)
  rownames(G) <- colnames(G) <- nodes
  
  no <- 0 # counter (number of tests)
  
  while(ii <= max(colSums(G))){
    for (alpha in nodes){
      for (beta in setdiff(nodes, alpha)){
        if (G[alpha,beta] == 0) next
        Osub <- setdiff(nodes[which(G[,beta] == 1)], alpha)
        if (length(Osub) < ii) next
        for (C in combn(Osub, ii, simplify=FALSE)){
          if(testFunc(alpha, beta, C) > thres){
            G[alpha,beta] <- 0
            no <- no + 1
            break
          } else {
            no <- no + 1
          }
        }
      }
    }
    ii <- ii + 1
  }
  
  list(G = G, nTests = no, p = thres)
}





################
# generate VAR(1) from graph G (N observed time points, burn-in of N)
# n is number of coordinate processes (inferred from graph size)

genVAR <- function(N, G, sigma = 1){
  
  n <- nrow(G)
  unstable <- TRUE
  while(unstable){
    A <- matrix(sample(c(-1,1), n^2, replace = TRUE)*runif(n^2,0,1), nrow = n)
    A[t(G) == 0] <- 0
    unstable <- max(abs(eigen(A)$values)) > 1
  }
  
  X <- matrix(0, nrow = 2*N, ncol = n)
  
  for (i in 2:(2*N)){
    X[i,] <- A%*%X[i-1,] + rnorm(n, sd = sigma)
  }
  
  colnames(X) <- colnames(G)
  X[(N+1):(2*N),]
}

################
# test Granger causality
# A,B,C should be colnames of X
# B should be singleton
# requires the FIAR package

condGrangerWrap <- function(X, A, B, C, order = 1){
  
  XX <- cbind(X[,A],X[,B],X[,C])
  condGranger(XX, length(A), length(B), order = order)
  
}

################
# subset list (base set is letters[1:n], all subsets are returned if is.null(k), otherwise subsets of size at most k)

allSubsets <- function(n, k){
  
  ll <- list()
  for (i in 1:n) ll[[length(ll) + 1]] <- c(0,1)
  ee <- expand.grid(ll)
  
  colnames(ee) <- letters[1:n]
  
  if (!is.null(k)){
    ee <- ee[rowSums(ee) < k+1,]
  }
  
  ee
  
}

################
# generate all Granger test results from data (X) (of C size at most k)
# Columns in X should have names

allTests <- function(X, k, order = 1){
  
  n <- ncol(X)
  V <- colnames(X)
  as <- allSubsets(n, k)
  
  res <- array(NA, dim = c(n, n, nrow(as)),
               dimnames = list(V, V, 1:nrow(as)))
  
  for (i in 1:n){
    for (j in 1:n){
      for (l in 1:nrow(as)){
        res[i,j,l] <- condGrangerWrap(X,V[i],V[j],V[as[l,] == 1], order = order)$prob
      }
    }
  }
  
  res
}

######
# find maximal Markov equivalent DMG from known DMG, G
# G: DMG

# brute force:
findMaxBrute <- function(G){
  
  nodes <- rownames(G$M1)
  
  for (alpha in nodes){
    for (beta in nodes){
      if (G$M1[alpha, beta] == 0){
        Gtmp <- G
        Gtmp$M1[alpha, beta] <- 1
        aM <- areME(G, Gtmp)
        if (aM) {
          G$M1[alpha, beta] <- 1
        }
      }
      if (G$M2[alpha, beta] == 0){
        Gtmp <- G
        Gtmp$M2[alpha, beta] <- Gtmp$M2[beta, alpha] <- 1
        aM <- areME(G, Gtmp)
        if (aM) {
          G$M2[alpha, beta] <- G$M2[beta, alpha] <- 1
        }
      }
    }
  }
  
  G
  
}