#' Title
#'
#' @param IDs A table listing metabolites along with their associated sub-pathways and super-pathways. Entries with missing values (NAs) will be excluded. The table must include columns named 'SUPER_PATHWAY' and 'SUB_PATHWAY'
#' @param Y A raw abundance matrix of dimensions p × n (metabolites × samples), allowing missing values (NAs).
#' @param X Covariates to be controlled (n × q).
#' @param B Standardized summary level data from a association study, p × s (metabolites × snps).
#' @param N Sample size of the association study.
#' @param K Number of latent factors to be estimated.
#' @param COMD.column The column in IDs whose values correspond to the row names of Y or B.
#' @param impute Whether impute missing values in Y.Defaults to False.
#' @param imputation.methods If impute = TRUE, specifies the imputation method to use. Available options are "MSNIMBLE" and "Mean".
#' @param Miss.Mech The missing data mechanism for Y. It will be inferred by MSNIMBLE if not explicitly specified
#' @param max.miss Maximum missing proportion for MSNIMBLE to consider.
#'
#' @return
#' * `data` A list of three elements: L, the loading matrix of K factors whose rows correspond to annotated metabolites;
#' metabolites, rows of IDs corresponding to annotated metabolites;
#' Var, a list of two vectors including the information to calculate the observation variance of L.
#'
#'* `data.unknown` A list of three elements: L, the loading matrix of K factors whose rows correspond to unknown metabolites;
#' metabolites, rows of IDs corresponding to unknown metabolites;
#' Var, a list of two vectors including the information to calculate the observation variance of L.
#' @export
#'
#' @examples
preProcess = function(IDs,Y=NULL,X=NULL,B=NULL,N=NULL,K=NULL,COMD.column,impute = F,imputation.methods = c("MSNIMBLE","Mean"),Miss.Mech = NULL,max.miss=NULL){
  if(!is.null(Y) & !is.null(B)){
    print("Only one of Y or B is needed")
  }else if(is.null(Y) & is.null(B)){
    print("You must provide Y or L")
  }else if(!is.null(Y)){
    if(is.null(X)){
      X = matrix(data = rep(1,ncol(Y)),ncol = 1)
    }
    if(impute){
      if(imputation.methods=="MSNIMBLE"){
        if(is.null(Miss.Mech)){
          Miss.Mech = EstimateMissing(Y = Y,Cov = X,K = K,max.missing.consider = max.miss)
        }
        MSNIMBLE <- FactorAnalysis(Y = Y,Cov = X,K = K,Miss.Mech = Miss.Mech,max.miss = max.miss)
        NAs <- unique(which(is.na(MSNIMBLE$L),arr.ind = T)[,1])
        unknown.sub <- which(is.na(IDs$SUB_PATHWAY))
        L.hat <- MSNIMBLE$L;rownames(L.hat) <- rownames(Y);L.hat <- L.hat[-NAs,]
        Var.m <-  MSNIMBLE$Var.L[-NAs]
        p <- dim(L.hat)[1]
        IDs_0 <- IDs[which(IDs[,COMD.column]%in%rownames(L.hat)),]
        IDs_known <- sort_by_sup_sub(data = IDs_0,supname = "SUPER_PATHWAY","SUB_PATHWAY")
        D <- p/colSums(L.hat^2)
        L.hat <- L.hat%*%diag(sqrt(D))
        knowns <- match(IDs_known[,COMD.column],rownames(L.hat))
        Var <- list(Var.m[knowns],D);names(Var) <- c("U","V")#
        data <- list(matrix(L.hat[knowns,],ncol = K),Var,IDs_known);names(data) <- c("L","var","metabolites")
        IDs <- data$metabolites
        Var <- list(Var.m[-knowns],D);names(Var) <- c("U","V")#
        data.unknown <- list(matrix(L.hat[-knowns,],ncol = K),Var,IDs_0[-knowns,]);names(data.unknown) <- c("L","var","metabolites")
        result = list(data,data.unknown);names(result) = c("data","data.unknown")
        return(result)
      }else if(imputation.methods=="Mean"){
        NAs <- unique(which(is.na(Y),arr.ind = T)[,1])
        if(length(NAs)>0){ Y[which(is.na(Y),arr.ind = T)] = mean(Y,na.rm=T);IDs = IDs[IDs[,COMD.column]%in%rownames(Y),] }
        unknown.sub <- which(is.na(IDs$SUB_PATHWAY))
        Y = Y - Y%*%X%*%solve(t(X)%*%X)%*%t(X)
        Y.svd = svd(Y,nu = K,nv = K)
        if(K==1){
          L.hat <- Y.svd$u%*%Y.svd$d[1]/sqrt(ncol(Y))
          Var.m <-  rowMeans((Y-Y.svd$u%*%Y.svd$d[1]%*%t(Y.svd$v))^2)/ncol(Y)
          p <- nrow(L.hat);D <- p/colSums(L.hat^2)
          L.hat <- L.hat%*%sqrt(D)
        }else{
          L.hat <- Y.svd$u%*%diag(Y.svd$d[1:K])/sqrt(ncol(Y))
          Var.m <-  rowMeans((Y-Y.svd$u%*%diag(Y.svd$d[1:K])%*%t(Y.svd$v))^2)/ncol(Y)
          p <- nrow(L.hat);D <- p/colSums(L.hat^2)
          L.hat <- L.hat%*%diag(sqrt(D))
        }
        rownames(L.hat) <- rownames(Y)
        IDs_known <- sort_by_sup_sub(data = IDs,supname = "SUPER_PATHWAY","SUB_PATHWAY")
        knowns <- match(IDs_known[,COMD.column],rownames(L.hat))
        Var <- list(Var.m[knowns],D);names(Var) <- c("U","V")#
        data <- list(matrix(L.hat[knowns,],ncol = K),Var,IDs_known);names(data) <- c("L","var","metabolites")

        Var <- list(Var.m[-knowns],D);names(Var) <- c("U","V")#
        data.unknown <- list(matrix(L.hat[-knowns,],ncol = K),Var,IDs[-knowns,]);names(data.unknown) <- c("L","var","metabolites")
        result = list(data,data.unknown);names(result) = c("data","data.unknown")
        return(result)
      }

    }else{
      NAs <- unique(which(is.na(Y),arr.ind = T)[,1])
      if(length(NAs)>0){ Y = Y[-NAs,];IDs = IDs[IDs[,COMD.column]%in%rownames(Y),] }
      unknown.sub <- which(is.na(IDs$SUB_PATHWAY))
      Y = Y - Y%*%X%*%solve(t(X)%*%X)%*%t(X)
      Y.svd = svd(Y,nu = K,nv = K)
      if(K==1){
        L.hat <- Y.svd$u%*%Y.svd$d[1]/sqrt(ncol(Y))
        Var.m <-  rowMeans((Y-Y.svd$u%*%Y.svd$d[1]%*%t(Y.svd$v))^2)/ncol(Y)
        p <- nrow(L.hat);D <- p/colSums(L.hat^2)
        L.hat <- L.hat%*%sqrt(D)
      }else{
        L.hat <- Y.svd$u%*%diag(Y.svd$d[1:K])/sqrt(ncol(Y))
        Var.m <-  rowMeans((Y-Y.svd$u%*%diag(Y.svd$d[1:K])%*%t(Y.svd$v))^2)/ncol(Y)
        p <- nrow(L.hat);D <- p/colSums(L.hat^2)
        L.hat <- L.hat%*%diag(sqrt(D))
      }
      rownames(L.hat) <- rownames(Y)
      IDs_known <- sort_by_sup_sub(data = IDs,supname = "SUPER_PATHWAY","SUB_PATHWAY")
      knowns <- match(IDs_known[,COMD.column],rownames(L.hat))
      Var <- list(Var.m[knowns],D);names(Var) <- c("U","V")#
      data <- list(matrix(L.hat[knowns,],ncol = K),Var,IDs_known);names(data) <- c("L","var","metabolites")

      Var <- list(Var.m[-knowns],D);names(Var) <- c("U","V")#
      data.unknown <- list(matrix(L.hat[-knowns,],ncol = K),Var,IDs[-knowns,]);names(data.unknown) <- c("L","var","metabolites")
      result = list(data,data.unknown);names(result) = c("data","data.unknown")
      return(result)
    }
  }else if(!is.null(B)){
    B.svd <- svd(x = B, nu =K, nv=K)
    M = nrow(B);S = ncol(B)
    L.hat <- B.svd$u[,1:K]*sqrt(M)
    Sigma <- vector()
    for (m in 1:nrow(B)) {
      Sigma[m] <- (sum(B[m,]^2)-sum( (t(B[m,])%*%B.svd$v[,1:K])%*%t((t(B[m,])%*%B.svd$v[,1:K]))  ) )/(S-1)
    }
    IDs_known <- sort_by_sup_sub(data = IDs,supname = "SUPER_PATHWAY","SUB_PATHWAY")
    knowns <- match(IDs_known[,COMD.column],rownames(B))

    #knowns
    L.known <- L.hat[knowns,]
    Sigma.known <- Sigma[knowns]
    Var <- list(Sigma.known,M*(1+S/N)/B.svd$d[1:K]^2);names(Var) <- c("U","V")
    data <- list(L.known,Var,IDs_known);names(data) <- c("L","var","metabolites")

    #unknowns
    unknowns = match(IDs[is.na(IDs$SUPER_PATHWAY),COMD.column],rownames(B))
    L.unknown <- L.hat[unknowns,]
    Sigma.unknown <- Sigma[unknowns]
    Var <- list(Sigma.unknown,M*(1+S/N)/B.svd$d[1:K]^2);names(Var) <- c("U","V")
    IDs.unknown <- IDs[unknowns,]
    data.unknown <- list(L.unknown,Var,IDs.unknown);names(data.unknown) <- c("L","var","metabolites")
    result = list(data,data.unknown);names(result) = c("data","data.unknown")
    return(result)
  }
}

estAlpha <- function(IDs){
  superpathways <- unique(IDs$SUPER_PATHWAY)
  J <- length(superpathways)
  m <- sapply(1:J, function(j){sum(IDs$SUPER_PATHWAY==superpathways[j])})
  B <- sapply(1:J, function(j){length(unique(IDs$SUB_PATHWAY[IDs$SUPER_PATHWAY==superpathways[j]]))})
  alpha.grid <- c(1:500)/100
  log.p <- sapply(alpha.grid, function(x){sum(sapply(1:J, function(j){B[j]*log(x)-sum(log(x+c(1:m[j])-1))}))})
  log.p <- log.p-max(log.p)
  p <- exp(log.p);p <- p/sum(p)
  alpha.meta <- alpha.grid[which.max(p)]
  return(alpha.meta)
}

Input_Gibbs <- function(Gibbs_Samples){
  K = length(Gibbs_Samples)
  MU <- list();SIGMA <- list();DISH <- list();DISH.assign <- list();ALPHA.sup <- list();ALPHA.sub <- list();
  P.0 <- list();P.00 <- list();P.out <- list()
  W <- list();P.W <- list();Taus <- list();P.Taus <- list()
  PCH <- list()
  B.l <- vector();B.u <- vector()
  for (k.latent in 1:K) {
    b.l <- ((range(Gibbs_Samples[[k.latent]]$data$L[,k.latent])-mean(Gibbs_Samples[[k.latent]]$data$L[,k.latent]))*1.25+mean(Gibbs_Samples[[k.latent]]$data$L[,k.latent]))[1]
    b.u <- ((range(Gibbs_Samples[[k.latent]]$data$L[,k.latent])-mean(Gibbs_Samples[[k.latent]]$data$L[,k.latent]))*1.25+mean(Gibbs_Samples[[k.latent]]$data$L[,k.latent]))[2]
    B.l[k.latent] <- b.l;B.u[k.latent] <- b.u
    MU[[k.latent]] <- Gibbs_Samples[[k.latent]]$Mu.g
    SIGMA[[k.latent]] <- Gibbs_Samples[[k.latent]]$Sigma.g
    DISH[[k.latent]] <- Gibbs_Samples[[k.latent]]$Dishes
    DISH.assign[[k.latent]] <- Gibbs_Samples[[k.latent]]$Dish.assign
    ALPHA.sup[[k.latent]] <-Gibbs_Samples[[k.latent]]$Alphas.sup
    ALPHA.sub[[k.latent]] <-Gibbs_Samples[[k.latent]]$Alphas.sub
    PCH[[k.latent]] <- Gibbs_Samples[[k.latent]]$PCHs
    P.0[[k.latent]] <-Gibbs_Samples[[k.latent]]$p.0s
    P.00[[k.latent]] <-Gibbs_Samples[[k.latent]]$p.00s
    P.out[[k.latent]] <- Gibbs_Samples[[k.latent]]$p.out
    W[[k.latent]] <- Gibbs_Samples[[k.latent]]$w
    P.W[[k.latent]] <- Gibbs_Samples[[k.latent]]$p.w
    Taus[[k.latent]] <- Gibbs_Samples[[k.latent]]$taus
    P.Taus[[k.latent]] <- Gibbs_Samples[[k.latent]]$p.taus
  }
  Result = list(K,data,MU,SIGMA,DISH,DISH.assign,ALPHA.sup,ALPHA.sub,PCH,Gibbs_Samples[[k.latent]]$record,Gibbs_Samples[[k.latent]]$Iter,p.L,P.0,P.00,P.out,W,P.W,Taus,P.Taus,B.l,B.u,Gibbs_Samples[[k.latent]]$p.outlier.alpha,Gibbs_Samples[[k.latent]]$p.outlier.beta)
  names(Result) = c("K","data","MU","SIGMA","DISH","DISH.assign","ALPHA.sup","ALPHA.sub","PCH","record","Iter","p.L","P.0","P.00","P.out","W","P.W","Taus","P.Taus","B.l","B.u","p.outlier.alpha","p.outlier.beta")
  return(Result)
}

Input_naive <- function(data,K){
  MU <- list();SIGMA <- list();DISH <- list();DISH.assign <- list();ALPHA.sup <- list();ALPHA.sub <- list();
  P.0 <- list();P.00 <- list();P.out <- list()
  W <- list();P.W <- list();Taus <- list();P.Taus <- list()
  PCH <- list()
  B.l <- vector();B.u <- vector()
  for (k.latent in 1:K) {
    b.l <- ((range(data$L[,k.latent])-mean(data$L[,k.latent]))*1.25+mean(data$L[,k.latent]))[1]
    b.u <- ((range(data$L[,k.latent])-mean(data$L[,k.latent]))*1.25+mean(data$L[,k.latent]))[2]
    B.l[k.latent] <- b.l;B.u[k.latent] <- b.u
    P.out[[k.latent]] <- 0.01
    x <- ML.gamma(data = data,k.latent = k.latent,var.threshold = 0.05)
    MU[[k.latent]] <- matrix(data = rep(sapply(unique(data$metabolites$SUB_PATHWAY), function(xxx){mean(data$L[data$metabolites$SUB_PATHWAY==xxx,k.latent])}),100),nrow = 100,byrow = T)
    SIGMA[[k.latent]] <- matrix(data = rep(sapply(unique(data$metabolites$SUB_PATHWAY), function(xxx){sd(data$L[data$metabolites$SUB_PATHWAY==xxx,k.latent])}),100),nrow = 100,byrow = T)
    for (i in 1:ncol(SIGMA[[k.latent]])) {
      if(is.na(SIGMA[[k.latent]][1,i])){
        SIGMA[[k.latent]][,i] <- rgamma(n = 100,shape = x$shape,rate = x$rate)
      }
    }
  }
  Result = list(K,data,MU,SIGMA,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,P.out,NULL,NULL,NULL,NULL,B.l,B.u,NULL,NULL)
  names(Result) = c("K","data","MU","SIGMA","DISH","DISH.assign","ALPHA.sup","ALPHA.sub","PCH","record","Iter","p.L","P.0","P.00","P.out","W","P.W","Taus","P.Taus","B.l","B.u","p.outlier.alpha","p.outlier.beta")
  return(Result)
}

estimate_mode <- function(x) {
  if(length(x)==1){
    return(x)
  }else{
    d <- density(x)
    return(d$x[which.max(d$y)])
  }
}

log.p2p <- function(x){
  x <- x-max(x)
  x <- exp(x)
  x <- x/sum(x)
  return(x)
}

sort_by_sup_sub <- function(data,supname,subname){
  sorted.data <- data.frame()
  for (sup.name in sort(unique(data[,supname]))) {
    subindices <- which(data[,supname]==sup.name)
    for (sub.name in sort(unique(data[,subname][subindices]))) {
      sorted.data <- rbind(sorted.data,data[which( data[,subname]==sub.name & data[,supname]==sup.name),])
    }
  }
  return(sorted.data)
}

matrix.up <- function(x){
  if(is.matrix(x)){
    x[col(x) > row(x)]
  }else{
    x
  }
}

logGamma <- function(alpha,n=0){
  #calculate log of Gamma(alpha+n), when n is a large integer
  if(n<2){
    x <- log(gamma(alpha+n))
  }else{
    x <- log(gamma(alpha))+sum(log( c(1:(n-1))+alpha ))
  }
  return(x)
}

find_sup_for_sub <- function(IDs,sub){
  return( unique( IDs$SUPER_PATHWAY[IDs$SUB_PATHWAY==sub] ) )
}

pred.sub <- function(unknowns,data,m,factors,alpha.meta,Iter,record,MU,SIGMA,ALPHA.sub,ALPHA.sup,DISH,DISH.assign,P.out,P.00,P.0,Taus,P.Taus,W,P.W,B.u,B.l,alpha.out=1/5,beta.out=999/5,PCH,plt=F,new.pathway=FALSE,plot.true = F){
  M <- nrow(data$metabolites)
  superpathways <- unique(data$metabolites$SUPER_PATHWAY)
  J <- length(superpathways)
  subpathways <- unique(data$metabolites$SUB_PATHWAY)
  input <- unknowns$L[m,factors]
  if(is.list(unknowns$var$U)){
    post.sub <- sapply(1:length(subpathways), function(b){sum(sapply(factors, function(k){ log(mean( P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+(1-P.out[[k]])*dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) )) }))})
    +log(sapply(subpathways, function(x){sum(data$metabolites$SUB_PATHWAY==x)}))
  }else{
    post.sub <- sapply(1:length(subpathways), function(b){sum(sapply(factors, function(k){ log(mean( P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+(1-P.out[[k]])*dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]]*unknowns$var$V[k])) )) }))})
    +log(sapply(subpathways, function(x){sum(data$metabolites$SUB_PATHWAY==x)}))
  }
  if(sum(is.nan(post.sub))>0){
    post.sub[is.nan(post.sub)] <- -Inf
  }
  if(sum(is.na(post.sub))>0){
    post.sub[is.na(post.sub)] <- -Inf
  }
  if(isTRUE(new.pathway)){
    G <- Iter-record+1
    post.new <- vector()
    for (j in 1:J) {
      superpathway <- superpathways[j]
      subpathways.j <- unique(data$metabolites$SUB_PATHWAY[data$metabolites$SUPER_PATHWAY==superpathway])
      if(is.list(unknowns$var$U)){
        post.new[j] <- sum(sapply(factors, function(k){ log(mean(P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+  (1-P.out[[k]])*(rowSums(sapply(which(unique(data$metabolites$SUB_PATHWAY)%in%subpathways.j), function(b){ dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]][k,k]*unknowns$var$V[k]))/(ALPHA.sub[[k]]+length(subpathways.j)) }))+  ALPHA.sub[[k]]/(ALPHA.sub[[k]]+length(subpathways.j))*sapply(1:G, function(g){(sum(colSums(DISH.assign[[k]][[g]])*dnorm(x = input[k],mean = DISH[[k]][[g]][,2],sd = sqrt(DISH[[k]][[g]][,3]^2+unknowns$var$U[[m]][k,k]*unknowns$var$V[k])))+ALPHA.sup[[k]][g]*p.L(L = input[k],var.m = unknowns$var$U[[m]][k,k]*unknowns$var$V[k],p.00 = P.00[[k]][g],taus =Taus[[k]],p.taus = P.Taus[[k]],p.0 = mean(P.0[[k]]),w = W[[k]],p.w = P.W[[k]]))/(sum(DISH.assign[[k]][[g]])+ALPHA.sup[[k]][g])})))) }))+log(alpha.meta)+log(length(subpathways.j))-log(length(subpathways))
      }else{
        post.new[j] <- sum(sapply(factors, function(k){ log(mean(P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+  (1-P.out[[k]])*(rowSums(sapply(which(unique(data$metabolites$SUB_PATHWAY)%in%subpathways.j), function(b){ dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]]*unknowns$var$V[k]))/(ALPHA.sub[[k]]+length(subpathways.j)) }))+  ALPHA.sub[[k]]/(ALPHA.sub[[k]]+length(subpathways.j))*sapply(1:G, function(g){(sum(colSums(DISH.assign[[k]][[g]])*dnorm(x = input[k],mean = DISH[[k]][[g]][,2],sd = sqrt(DISH[[k]][[g]][,3]^2+unknowns$var$U[[m]]*unknowns$var$V[k])))+ALPHA.sup[[k]][g]*p.L(L = input[k],var.m = unknowns$var$U[[m]]*unknowns$var$V[k],p.00 = P.00[[k]][g],taus =Taus[[k]],p.taus = P.Taus[[k]],p.0 = mean(P.0[[k]]),w = W[[k]],p.w = P.W[[k]]))/(sum(DISH.assign[[k]][[g]])+ALPHA.sup[[k]][g])})))) }))+log(alpha.meta)+log(length(subpathways.j))-log(length(subpathways))
      }
    }
    post.sub <- c(post.sub,post.new)
    sub.append <- c(subpathways,"new subpathway")

    post.sub <- post.sub-max(post.sub)
    post.sub <- exp(post.sub)
    post.sub <- post.sub/sum(post.sub)

    x <- c(post.sub[1:length(subpathways)],sum(post.sub[(length(subpathways)+1):length(post.sub)]))

  }else{
    sub.append <- subpathways
    post.sub <- post.sub-max(post.sub)
    post.sub <- exp(post.sub)
    x <- post.sub/sum(post.sub)
  }

  #sub.append[which.max(x)]
  if(plt==T){
    if(which.max(x)==(length(subpathways)+1)){
      labels <- c(rep("",length(subpathways)),"new subpathway")
      plot(x,pch=19,xlab = "Subpathway index",ylab = "Bayes posterior",main = unknowns$metabolites$BIOCHEMICAL[m],col=c(rep("black",length(subpathways)),"red"),ylim=c(min(x)*0.9,max(x)*1.1))
      text(1:(length(subpathways)+1),x*c(rep(0.95,length(subpathways)),1.25),labels)
      p.j <- post.sub[(length(subpathways)+1):length(post.sub)]
      plot(p.j,pch=19,xlab = "Superpathway index",ylab = "Bayes posterior",main = unknowns$metabolites$BIOCHEMICAL[m],col=ifelse(p.j==p.j[which.max(p.j)],"red","black"),ylim=c(min(p.j)*0.9,max(p.j)*1.1))
      if(isTRUE(plot.true)){
        plot(rep(factors,each=sum(data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m])),sapply(1:K, function(k){data$L[data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m],k]}),xlab = "Factor index",ylab = "Metabolite loading",main = paste0("True pathway:",data$metabolites$SUB_PATHWAY[m]),pch=6)
        abline(h=0,lty="dashed")
        lines(factors,input,pch=25,col="red", lwd=2,type = "o")
        lines(x = factors,y = sapply(factors, function(xxx){ mean( MU[[xxx]][,which(subpathways==data$metabolites$SUB_PATHWAY[m])] ) }), lwd=2,col="blue")
      }
    }else{
      labels <- ifelse(sub.append==sub.append[which.max(x)],sub.append[which.max(x)],"");#labels[length(labels)] <- "new subpathway"
      #pdf(file = paste0(unknowns$metabolites$BIOCHEMICAL[m]," classification result.pdf"),width = 8,height = 6)
      #par(mfrow=c(2,1))
      #plot(x,pch=19,xlab = "Subpathway index",ylab = "Posterior probability",main = unknowns$metabolites$BIOCHEMICAL[m],col=ifelse(subpathways==subpathways[which.max(post.sub)],"red","black"),xlim=c(1:length(x)),ylim=c(min(x)*0.9,max(x)*1.1))

      #add vetical lines indicating superpathways
      subs.supers <- sapply(subpathways, function(xxx){find_sup_for_sub(IDs = data$metabolites,sub = xxx)})
      plot(x,pch=19,xlab = "Subpathway index",ylab = "Posterior probability",main = unknowns$metabolites$BIOCHEMICAL[m],col=ifelse(subpathways==subpathways[which.max(post.sub)],"red","black"),ylim=c(0,1),xlim=c(1,length(x)))
      abline(v = which(subs.supers[-length(subs.supers)]!=subs.supers[-1])+0.5,lty="dashed" )
      text(1:length(x),x*c(rep(0.95,length(x))),labels)
      #dev.off()
      if(isTRUE(plot.true)){
        b <- which.max(post.sub)
        p1 <- sapply(factors, function(k){ log10(mean( P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+(1-P.out[[k]])*dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) )) })
        b <- which(subpathways ==data$metabolites$SUB_PATHWAY[m] )
        p2 <- sapply(factors, function(k){ log10(mean( P.out[[k]]*( pnorm((B.u[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) - pnorm((B.l[k]-input[k])/sqrt(unknowns$var$U[[m]][k,k]*unknowns$var$V[k])))/(B.u[k]-B.l[k])+(1-P.out[[k]])*dnorm(x = input[k],mean = MU[[k]][,b],sd = sqrt(SIGMA[[k]][,b]^2+unknowns$var$U[[m]][k,k]*unknowns$var$V[k])) )) })

        pdf(file = paste0("figs/COPSAC_old_",unknowns$metabolites$COMP_ID[m],"_loading_plot.pdf"),width = 8,height = 10)
        par(mfrow=c(3,1))
        plot(rep(c(factors),each=sum(data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)])),sapply(1:K, function(k){data$L[data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)],k]}),xlab = "Factor index",ylab = "Metabolite loading",main = paste0("Pred pathway:",subpathways[which.max(post.sub)]),pch=1,ylim=c( 1.25*min(sapply(1:K,function(k){ min(data$L[data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)],k]) })),1.25*max(sapply(1:K,function(k){ max(data$L[data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)],k]) })) ))
        abline(h=0,lty="dashed")
        lines(factors,input,pch=25,col="red", lwd=2,type = "o")
        lines(x = factors,y = sapply(factors, function(xxx){ mean( MU[[xxx]][,which(subpathways==subpathways[which.max(post.sub)])] ) }), lwd=2,col="black")
        lim.up <- 0.1+sapply(factors,function(k){ max(data$L[data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)],k]) })
        lim.up[lim.up<0] <- 0.1
        text(factors,lim.up,round(p1-p2,2),col=ifelse(p1-p2>0,"blue","red"),cex = 0.8)

        plot(rep(factors,each=sum(data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m])),sapply(1:K, function(k){data$L[data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m],k]}),xlab = "Factor index",ylab = "Metabolite loading",main = paste0("True pathway:",data$metabolites$SUB_PATHWAY[m]),pch=1,ylim=c( 1.25*min(sapply(1:K,function(k){ min(data$L[data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m],k]) })),1.25*max(sapply(1:K,function(k){ max(data$L[data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m],k]) })) ))
        abline(h=0,lty="dashed")
        lines(factors,input,pch=25,col="red", lwd=2,type = "o")
        lines(x = factors,y = sapply(factors, function(xxx){ mean( MU[[xxx]][,which(subpathways==data$metabolites$SUB_PATHWAY[m])] ) }), lwd=2,col="black")
        lim.up <- 0.1+sapply(factors,function(k){ max(data$L[data$metabolites$SUB_PATHWAY==data$metabolites$SUB_PATHWAY[m],k]) })
        lim.up[lim.up<0] <- 0.1
        text(factors,lim.up,round(p1-p2,2),col=ifelse(p1-p2>0,"blue","red"),cex = 0.8)

        plot(factors,p1-p2,col=ifelse(p1-p2>0,"blue","red"),type="h",main=paste0(round(sum((p1-p2)[p1-p2 >0]),3),"",round(sum((p1-p2)[p1-p2 <0]),3)));abline(h = 0)
        dev.off()
      }else{
        plot(rep(factors,each=sum(data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)])),sapply(factors, function(k){data$L[data$metabolites$SUB_PATHWAY==subpathways[which.max(post.sub)],k]}),xlab = "Factor index",ylab = "Metabolite loading",main = paste0("Pred pathway:",subpathways[which.max(post.sub)]),pch=1)
        abline(h=0,lty="dashed")
        lines(factors,input,pch=25,col="red", lwd=2,type = "o")
        lines(x = factors,y = sapply(factors, function(xxx){ mean( MU[[xxx]][,which(subpathways==subpathways[which.max(post.sub)])] ) }), lwd=2,col="black")
      }
    }
  }
  result <- list(sub.append[which.max(x)],x,unknowns$metabolites$BIOCHEMICAL[m],unknowns$metabolites$SUB_PATHWAY[m])
  names(result) <- c("pred.path","Bayes.post","BIOCHEMICAL","true.label")
  return(result)
}
