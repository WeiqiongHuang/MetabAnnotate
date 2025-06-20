ML.gamma <- function(data,k.latent,n.metabolite=3,var.threshold = 1e-3){
  IDs <- data$metabolites
  t <- table(IDs$SUB_PATHWAY)
  subp.big <- names(t)[t>=n.metabolite]
  IDs.big <- IDs[IDs$SUB_PATHWAY %in% subp.big,]
  # var.g <- sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){var(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent][abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3])})-sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){mean( sapply(which(IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i])[abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3], function(m){data$var[m,k.latent]}) )})
  if(is.list(data$var$U)){
    var.g <- sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){var(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent][abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3])})-sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){mean( sapply(which(IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i])[abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3], function(xxx){data$var$U[[xxx]][k.latent,k.latent]})*data$var$V[k.latent] )})
  }else{
    var.g <- sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){var(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent][abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3])})-sapply(1:length(unique(IDs.big$SUB_PATHWAY)), function(i){mean( sapply(which(IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i])[abs((data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]-mean(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))/sd(data$L[IDs$SUB_PATHWAY==unique(IDs.big$SUB_PATHWAY)[i],k.latent]))<3], function(xxx){data$var$U[xxx]})*data$var$V[k.latent] )})
  }
  var.g <- var.g[var.g>var.threshold]
  x <- univariateML::mlgamma(var.g)
  # plot(density(var.g));lines(c(1:2000)/1000,y = dgamma(x = c(1:2000)/1000,shape = 0.5143  ,rate = 1.6972  ),col="red")
  # x[1]
  result <- list(x[1],x[2]);names(result) <- c("shape","rate")
  return(result)
}

listOfPathway <- function(data,supname,subname,L.hat,var.g,k){
  G <- length(unique(data[,subname]))
  data.list <- list()
  g <- 1
  for (i in 1:length(unique(data[,supname]))) {
    sup.name <- unique(data[,supname])[i]
    subindices <- which(data[,supname]==sup.name)

    Nj <- length( unique(data[,subname][subindices]) )#num of subpathways in superpathway
    # mode_j <- estimate_mode(L.hat[which(data$SUPER_PATHWAY==sup.name),k])

    subpathways <- list()
    for (j in 1:Nj) {
      sub.name <- sort(unique(data[,subname][subindices]))[j]
      indices <- which(data$SUB_PATHWAY==sub.name & data$SUPER_PATHWAY==sup.name)
      theta <-  rbeta(n = 1,shape1 = 1,shape2 = 1)
      #mu <- new.para(w)
      if(length(indices)==1){
        nu <- 0
      }else{
        nu <- rep( sqrt(var( L.hat[indices,k] )),length(indices) )
      }
      if(is.list(var.g$U)){
        var.m = sapply(indices, function(xxx){ var.g$U[[xxx]][k,k] })*var.g$V[k]
      }else{
        var.m = var.g$U[indices]*var.g$V[k]
      }
      list1 <- list(sub.name,
                    sup.name,
                    g,
                    indices,
                    L.hat[indices,k],
                    var.m,#variance of L.hat[m,k] given L[m,k]
                    L.hat[indices,k],
                    0,
                    rep(2,length(indices)),
                    length(indices),#Mg
                    rep(0,length(indices)),#z=1:outlier
                    rep(0,length(indices)),#cluster
                    rep(0,length(indices)),#t
                    rep(0,length(indices)),#k
                    theta,1,1,
                    rep(0,length(indices)),#mu
                    nu,
                    vector())#nu

      names(list1) <- c("Subpathway","superPathway","g","indices","Lgm.hat",
                        "var.gm","Lgm","mode_j","indices_0",
                        "Mg","z","c","t","k","theta","a","b","mu","sigma_g","is.spike")
      subpathways[[j]] <- list1
      g <- g + 1
    }
    names(subpathways) <- sort(unique(data[,subname][subindices]))
    data.list[[i]] <- subpathways
  }
  names(data.list) <- sort(unique(data[,supname]))
  return(data.list)
}

InitDishes <- function(){
  dishes <- matrix(ncol = 3);colnames(dishes) <- c("k","mu","sigma");dishes <- dishes[-1,]
  return(dishes)
}

InitRestaurants <- function(pathways){
  J <- length(pathways)
  rest.list <- list()
  for (j in 1:J) {
    matrix_j <- matrix(ncol = 4)
    colnames(matrix_j) <- c("table","dish","mu","sigma")
    matrix_j <- matrix_j[-1,]
    rest.list[[j]] <- matrix_j
  }
  names(rest.list) <- names(pathways)
  return(rest.list)
}

InitClusters <- function(pathways){
  J <- length(pathways)
  rest.list <- list()
  for (j in 1:J) {
    rest_j <- list()
    for (g in 1:length(pathways[[j]])) {
      matrix_jg <- matrix(ncol = 5)
      colnames(matrix_jg) <- c("cluster","table","dish","mu","sigma")
      matrix_jg <- matrix_jg[-1,]
      rest_j[[g]] <- matrix_jg
    }
    names(rest_j) <- names(pathways[[j]])
    rest.list[[j]] <- rest_j
  }
  names(rest.list) <- names(pathways)
  return(rest.list)
}

newTable <- function(restaurants,j,g){
  if(length(restaurants[[j]][[g]][,1])>0){
    return(nextIndix(restaurants[[j]][[g]][,1]))
  }else{
    return(1)
  }
}

pi.post <- function(M_k,alpha,gamma,Stirling.num){
  if(ncol(M_k)==0){
    return(1)
  }else{
    pi0 <- rep(1/ncol(M_k),ncol(M_k))
    r <- matrix(data = 0,nrow = nrow(M_k),ncol = ncol(M_k))
    for (iter in 1:50) {
      for (j in 1:nrow(M_k)) {
        for (k in 1:ncol(M_k)) {
          #j <- 1;k <- 1
          if(M_k[j,k]>0){
            prob <- Stirling.num[M_k[j,k],1:M_k[j,k]]*(alpha*pi0[k])^c(1:M_k[j,k])
            r[j,k] <- sample(x = c(1:M_k[j,k]),size = 1,prob = prob)
          }else{
            r[j,k] <- 0
          }
        }
      }
      pi0 <- rgamma(n = ncol(M_k)+1,shape = c(colSums(r),gamma),rate = 1);pi0 <- pi0/sum(pi0)
    }
    return(pi0)
  }
}

p.spike.spike <- function(L,var.m){
  return( prod(dnorm(x=L,mean=0,sd=sqrt(var.m))) )
}

p.spike.slab <- function(L,var.m,taus){
  x <- sapply(1:length(taus), function(l){prod(dnorm(x=L,mean=0,sd=sqrt(taus[l]^2+var.m)))})
  return(x)
}

p.slab.slab <- function(L,var.m,taus,w){
  x <- sapply(1:length(taus), function(l){prod(dnorm(x=L,mean=0,sd=sqrt(taus[l]^2+var.m+w^2)))})
  return(x)
}

p.L.given.theta <- function(L,mu,sigma,var.m){
  x <- prod(dnorm(x=L,mean=mu,sd=sqrt(sigma^2+var.m)))
  return(x)
}

p.L <- function(L,var.m,p.00,taus,p.taus,p.0,w,p.w){
  if(length(L)==1){
    x = sum( sapply(1:length(taus),function(l){ sapply(1:length(w),function(b){ p.w[b]*p.taus[l]*prod(dnorm(x=L,mean=0,sd=sqrt(taus[l]^2+var.m+w[b]^2))) }) }) )
  }else{
    x = sum( sapply(1:length(taus),function(l){ sapply(1:length(w),function(b){ p.w[b]*p.taus[l]*mvtnorm::dmvnorm(x = as.vector(L),mean = rep(0,length(L)),sigma = diag(taus[l]^2+as.vector(var.m))+w[b]^2) }) }) )
  }
  x <- (1-p.0)*x + p.0*sum(sapply(1:length(taus), function(h){p.taus[h]*prod(dnorm(x=L,mean=0,sd=sqrt(taus[h]^2+var.m)))}))
  x <- (1-p.00)*x + p.00*prod(dnorm(x=L,mean=0,sd=sqrt(var.m)))
  return(x)
}

is.outlier <- function(pathways,j,g,m,p.outlier,alpha1_mu,alpha2_mu,clusters,restaurants,dishes,M_k,N_t,N_c,taus,w,p.w,p.taus,p.00,p.0,b.l,b.u){
  if(p.outlier==0){
    return(0)
  }else if(pathways[[j]][[g]]$Mg==1){
    return(0)
  }else{
    #j <- 8;g <- 7;m <- 50
    p.m.outlier <- p.outlier/(b.u-b.l)*( pnorm((b.u-pathways[[j]][[g]]$Lgm.hat[m])/sqrt(pathways[[j]][[g]]$var.gm[m])) - pnorm((b.l-pathways[[j]][[g]]$Lgm.hat[m])/sqrt(pathways[[j]][[g]]$var.gm[m])))
    L = pathways[[j]][[g]]$Lgm.hat[m];var.m = pathways[[j]][[g]]$var.gm[m]
    if(nrow(clusters[[j]][[g]])==0){
      if(nrow(restaurants[[j]])==0){
        if(nrow(dishes)==0){
          x <- p.L(L = L,var.m = var.m,p.00 = p.00,taus = taus,p.taus = p.taus,p.0 = p.0,w = w,p.w = p.w)
        }else{
          x <- ( alpha1_mu / (sum(M_k)+alpha1_mu) )*p.L(L = L,var.m = var.m,p.00 = p.00,taus = taus,p.taus = p.taus,p.0 = p.0,w = w,p.w = p.w)
          x <- x + sum(sapply(1:nrow(dishes), function(k){sum(M_k[,k])/(sum(M_k)+alpha1_mu)*p.L.given.theta(L = L,mu = dishes[k,2],sigma = dishes[k,3],var.m = var.m)}))
        }
      }else{
        x <- ( alpha1_mu / (sum(M_k)+alpha1_mu) )*p.L(L = L,var.m = var.m,p.00 = p.00,taus = taus,p.taus = p.taus,p.0 = p.0,w = w,p.w = p.w)
        x <- x + sum(sapply(1:nrow(dishes), function(k){sum(M_k[,k])/(sum(M_k)+alpha1_mu)*p.L.given.theta(L = L,mu = dishes[k,2],sigma = dishes[k,3],var.m = var.m)}))
        likelis <- sapply(1:nrow(restaurants[[j]]), function(t){p.L.given.theta(L = L,mu = restaurants[[j]][t,3],sigma = restaurants[[j]][t,4],var.m = var.m)})
        n.jt <- sapply(1:nrow(restaurants[[j]]), function(t){sum(sapply(clusters[[j]], function(g){sum(g[,2]==restaurants[[j]][t,1])}))})
        x <- x*alpha2_mu/( sum(n.jt) + alpha2_mu ) + sum(likelis*n.jt/(sum(n.jt)+alpha2_mu))
      }
      p.m.HDP <- (1-p.outlier)*x
    }else{
      p.m.HDP <- (1-p.outlier)*dnorm(x = L,mean = clusters[[j]][[g]][1,4],sd = sqrt(clusters[[j]][[g]][1,5]^2+var.m))
    }
    return(rbinom(n = 1,size = 1,prob = p.m.outlier/(p.m.outlier+p.m.HDP)))
  }
}

#count the number of clusters except for c_m sitting at each table in j
n.jt.except.c <- function(restaurants,clusters,j,g,c.m){
  if(nrow(restaurants[[j]])==0){
    return(0)
  }else{
    n.jt <- sapply(X = restaurants[[j]][,1],FUN = function(x){sum(sapply(X = clusters[[j]],FUN = function(tc,t){sum(tc[,2]==t)},t=x))})
    if( c.m%in% clusters[[j]][[g]][,1] ){
      t.c <- clusters[[j]][[g]][clusters[[j]][[g]][,1]==c.m,2]
      n.jt <- n.jt - (restaurants[[j]][,1]==t.c)
    }
    return(n.jt)
  }
}

sample_t_c <- function(pathways,j,g,c.m,alpha1_mu,alpha2_mu,restaurants,dishes,clusters,M_k,taus,p.taus,p.00,p.0,w,p.w){
  if(nrow(restaurants[[j]])==0){#if there's no existing tables, return 1
    return( 1 )
  }else{
    L = pathways[[j]][[g]]$Lgm.hat[pathways[[j]][[g]]$c==c.m];var.m = pathways[[j]][[g]]$var.gm[pathways[[j]][[g]]$c==c.m]
    likelis <- sapply(1:nrow(restaurants[[j]]), function(t)(p.L.given.theta(L = L,mu = restaurants[[j]][t,3],sigma = restaurants[[j]][t,4],var.m = var.m)))
    n.jt.c <- n.jt.except.c(restaurants = restaurants,clusters = clusters,j = j,g = g,c.m = c.m)
    p.t <- likelis*n.jt.c
    x <- ( alpha1_mu / (sum(M_k)+alpha1_mu) )*p.L(L = L,var.m = var.m,p.00 = p.00,taus = taus,p.taus = p.taus,p.0 = p.0,w = w,p.w = p.w)
    x <- x + sum(sapply(1:nrow(dishes), function(k){ sum(M_k[,k]) / (sum(M_k)+alpha1_mu) *p.L.given.theta(L = L,mu = dishes[k,2],sigma = dishes[k,3],var.m = var.m)}))
    p.t <- c(p.t,x*alpha2_mu)
    return(sample(x = c(restaurants[[j]][,1],nextIndix(restaurants[[j]][,1])),size = 1,prob = p.t))
    #plot(c(likelis,x))
    #plot(p.t)
  }
}

kForNewTable4 <- function(j,dishes,L,taus,p.taus,var.m,M_k,alpha1_mu,alpha2_mu,w,p.w,p.00,p.0){
  if(length(dishes[,1])>0){
    probs <- sapply(1:nrow(dishes), function(k){p.L.given.theta(L = L,mu = dishes[k,2],sigma = dishes[k,3],var.m = var.m)})
    m.k <- sapply(1:nrow(dishes), function(k){sum(M_k[,k])/(sum(M_k)+alpha1_mu)})
    probs <- probs*m.k
    probs <- c(probs,  alpha1_mu/(sum(M_k)+alpha1_mu)*p.L(L = L,var.m = var.m,p.00 = p.00,taus = taus,p.taus = p.taus,p.0 = p.0,w = w,p.w = p.w)  )
    k.t.new <- sample(x = c( dishes[,1], nextIndix(dishes[,1])),size = 1,prob = probs)
  }else{
    k.t.new <-  1
  }
  return(k.t.new)
}

new.parameter <- function(L,var.m,w,p.w,p.00,p.0,taus,p.taus){
  log.post.00 <- log(p.00)+sum(log( dnorm(x = L,mean = 0,sd = sqrt(var.m)) ))
  log.post.0h <- log(1-p.00)+log(p.0)+log(p.taus)+sapply(1:length(taus), function(h){ sum(dnorm(x = L,mean = 0,sd = sqrt(taus[h]^2+var.m),log = T)) })
  log.post.bh = log(1-p.00)+log(1-p.0)+sapply(1:length(taus),function(h){ sapply(1:length(w),function(b){ log(p.w[b])+log(p.taus[h])+ ifelse(length(L)==1,dnorm(x = L,mean = 0,sd = sqrt(taus[h]^2+var.m+w[b]^2),log = T),mvtnorm::dmvnorm(x = as.vector(L),mean = rep(0,length(L)),sigma = diag(taus[h]^2+as.vector(var.m))+w[b]^2,log = T)) }) })
  log.max <- max(c( log.post.00, max(log.post.0h),max(log.post.bh) ))
  p.spike <- exp(log.post.00-log.max)/( exp(log.post.00-log.max) + sum(exp(log.post.0h-log.max)) + sum(exp(log.post.bh-log.max))  )
  z.00 <- rbinom(n = 1,size = 1,prob = p.spike)
  if(z.00==1){
    new.mu <- 0;new.sigma <- 0
  }else{
    log.max <- max( max(log.post.0h),max(log.post.bh) )
    p.mu.0 <- sum(exp(log.post.0h-log.max))/(sum(exp(log.post.0h-log.max)) + sum(exp( log.post.bh -log.max))  )
    z.0 <- rbinom(n = 1,size = 1,prob = p.mu.0)
    if(z.0==1){
      new.mu <- 0
      new.sigma <- sample(x = taus,size = 1,prob = exp(log.post.0h-max(log.post.0h))/ sum(exp(log.post.0h-max(log.post.0h))) )
    }else{
      xxx = log.post.bh-max(log.post.bh)
      mu.b <- sample(x = 1:length(w),size = 1,prob = rowSums(exp(xxx)))
      xxx = log.post.bh[mu.b,];xxx = xxx-max(xxx)
      new.sigma <- sample(x = taus,size = 1,prob = exp(xxx) )
      post.var <- 1/( 1/w[mu.b]^2+ sum( 1/(new.sigma^2+var.m) ) )
      post.mean <- sum(L/(new.sigma^2+var.m))*post.var
      new.mu <- rnorm(n = 1,mean = post.mean,sd = sqrt(post.var))
    }
  }
  result <- list(new.mu,new.sigma);names(result) <-c("mu","sigma")
  return(result)
}

nextIndix <- function(x){return( ifelse(length(x)==0,1,min(setdiff(c(1:(max(x)+1)),x))))}
#number of tables serving dish k
updateM_k <- function(dishes,M_k,j,restaurants){
  k.new <- setdiff(dishes[,1],colnames(M_k))
  colnames.M_k <- c(colnames(M_k),k.new)
  M_k <- cbind(M_k,matrix(data = 0,nrow = nrow(M_k),ncol = length(k.new)))
  colnames(M_k) <- colnames.M_k
  M_k[j,] <- sapply(c(1:ncol(M_k)), function(k){sum(restaurants[[j]][,2] == colnames(M_k)[k])})
  return(M_k)
}
#number of clusters sitting at table t
updateN_t <- function(N_t,j,restaurants,clusters){
  t.new <- setdiff(restaurants[[j]][,1],colnames(N_t[[j]]))
  colnames.Ntj <- c(colnames(N_t[[j]]),t.new)
  N_t[[j]] <- cbind(N_t[[j]],matrix(data = 0,nrow = nrow(N_t[[j]]),ncol = length(t.new)))
  colnames(N_t[[j]]) <- colnames.Ntj
  for (t in 1:ncol(N_t[[j]])) {
    N_t[[j]][,t] <- sapply(clusters[[j]], function(g){sum(g[,2]==colnames(N_t[[j]])[t])})
  }
  return(N_t)
}
#which cluster metabolite m belongs to
updateN_c <- function(N_c,clusters,pathways,j,g){
  new.c <- setdiff(clusters[[j]][[g]][,1],colnames(N_c[[j]][[g]]))
  colnames.Nc.jg <- c(colnames(N_c[[j]][[g]]),new.c)
  N_c[[j]][[g]] <- cbind(N_c[[j]][[g]],matrix(data = 0,nrow = nrow(N_c[[j]][[g]]),ncol = length(new.c)))
  colnames(N_c[[j]][[g]]) <- colnames.Nc.jg
  for (c in 1:ncol(N_c[[j]][[g]])) {
    N_c[[j]][[g]][,c] <- pathways[[j]][[g]]$c==colnames(N_c[[j]][[g]])[c]
  }
  return(N_c)
}

updataTauCounts2 <- function(taus,clusters,pathways){
  counts <-rep(0,length(taus))
  for (j in 1:J) {
    for (g in 1:length(clusters[[j]])) {
      #j <- 1;g <- 1
      for (c in clusters[[j]][[g]][,1]) {
        if( sum(pathways[[j]][[g]]$c==c)>1 ){
          tau <- clusters[[j]][[g]][which(clusters[[j]][[g]][,1]==c),5]
          counts[taus==tau] <- counts[taus==tau]+1
        }
      }
    }
  }
  return(counts)
}

update.abcd <- function(clusters,J){
  a <- 1 + sum(sapply(1:J, function(j){sum(sapply(1:length(clusters[[j]]), function(g){sum(clusters[[j]][[g]][,4]==0)}))}))
  b <- 1 + sum(sapply(1:J, function(j){sum(sapply(1:length(clusters[[j]]), function(g){sum(clusters[[j]][[g]][,4]!=0)}))}))
  c <- 1 + sum(sapply(1:J, function(j){sum(sapply(1:length(clusters[[j]]), function(g){sum(clusters[[j]][[g]][,5]==0)}))}))
  d <- 1 + sum(sapply(1:J, function(j){sum(sapply(1:length(clusters[[j]]), function(g){sum(clusters[[j]][[g]][,5]!=0)}))}))
  result <- list(a,b,c,d);names(result) <- c("a","b","c","d")
  return(result)
}

select.all.metabolites <- function(k_para,pathways){
  Lgm <- vector();var_m <- vector()
  J <- length(pathways)
  for (j in 1:J) {
    for (g in 1:length(pathways[[j]])) {
      Lgm <- c(Lgm,pathways[[j]][[g]]$Lgm[pathways[[j]][[g]]$k==k_para])
      var_m <- c(var_m,pathways[[j]][[g]]$var.gm[pathways[[j]][[g]]$k==k_para])
    }
  }
  m.list <- list(Lgm,var_m)
  names(m.list) <- c("L.hat","var.m")
  return(m.list)
}

estIG.mom <- function(x){
  alpha <- mean(x)^2/var(x)+2
  beta <- (alpha-1)*mean(x)
  result <- list(alpha,beta);names(result) <- c("alpha","beta")
  return(result)
}
