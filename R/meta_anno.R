#' Title
#'
#' @param mtb_data Output of 'preProcess'
#' @param size.clbr The number of metabolites selected as calibration set
#' @param Iter Total number of Gibbs samplers to draw
#' @param burnin The number of burn-ins to discard
#' @param new.pathway If a metabolite does not fit well with any existing pathway, determine whether to assign it to an unobserved pathway. Defaults to False. Setting this to True may significantly increase computation time
#' @param seed The seed for selecting calibration metabolites
#'
#' @return
#' * `annt.result` A data frame of three columns: metabolites biochemical ids, annotated pathway labels and confidence scores after calibration.
#'
#'* `calibration.para` The fitted parameters for two-parameter Beta calibration, a for slop and c for intercept.
#'* `pred.details` The prediction result for unknown metabolites before calibration.
#' @export
#'
#' @examples
metabAnnotate = function(mtb_data,size.clbr=30,Iter=100,burnin=50,new.pathway = F,seed=NULL){
  data = mtb_data$data
  data.unknown = mtb_data$data.unknown
  #estimate concentration parameter
  alpha.meta <- estAlpha(IDs = data$metabolites)
  #randomly select calibration set
  set.seed(seed = seed)
  clbr <- sample(x = 1:nrow(data$metabolites),size = size.clbr,replace = F)
  train <- setdiff(1:nrow(data$metabolites),clbr)

  Var <- list(data$var$U[train],data$var$V);names(Var) <- c("U","V")
  data.train <- list(data$L[train,],Var,data$metabolites[train,]);names(data.train) <- c("L","var","metabolites")
  Var <- list(data$var$U[clbr],data$var$V);names(Var) <- c("U","V")
  data.clbr <- list(data$L[clbr,],Var,data$metabolites[clbr,]);names(data.clbr) <- c("L","var","metabolites")

  Gibbs_samples = HDP_Gibbs(data = data.train,Iter=Iter,record=burnin+1)

  input=Input_Gibbs(Gibbs_Samples = Gibbs_samples)

  MU = input$MU;SIGMA = input$SIGMA;ALPHA.sub = input$ALPHA.sub;ALPHA.sup = input$ALPHA.sup;DISH = input$DISH;DISH.assign = input$DISH.assign;P.00 = input$P.00;P.out = input$P.out;P.0 = input$P.0;Taus = input$Taus;P.Taus = input$P.Taus;W = input$W;P.W = input$P.W;alpha.out = input$p.outlier.alpha;beta.out = input$p.outlier.beta;PCH = input$PCH;Iter = input$Iter;record = input$record;B.u = input$B.u;B.l = input$B.l
  factors = c(1:length(MU))
  Result.clbr = list()
  for (m in 1:nrow(data.clbr$metabolites)) {
    Result.clbr[[m]] = pred.sub(unknowns = data.clbr,data = data.train,m = m,factors = factors,Iter = Iter,record = record,alpha.meta = alpha.meta,MU = MU,SIGMA = SIGMA,ALPHA.sub = ALPHA.sub,ALPHA.sup = ALPHA.sup,DISH = DISH,DISH.assign = DISH.assign,P.00 = P.00,P.out = P.out,P.0 = P.0,Taus = Taus,P.Taus = P.Taus,W = W,P.W = P.W,B.u = B.u,B.l = B.l,alpha.out = p.outlier.alpha,beta.out = p.outlier.beta,PCH = PCH,plt = F,new.pathway = new)
  }
  s.clbr = sapply(Result.clbr, function(xxx){max(xxx$Bayes.post)})
  y.clbr = sapply(Result.clbr, function(xxx){xxx$pred.path==xxx$true.label})

  s = log(s.clbr)-log(1-s.clbr)
  y.clbr = y.clbr[!is.infinite(s)];s = s[!is.infinite(s)]
  fit.logit = glm(y.clbr~s, family = "binomial")
  a = fit.logit$coefficients[2];c = fit.logit$coefficients[1]
  Result.unknown = list()
  for (m in 1:nrow(data.unknown$metabolites)) {
    Result.unknown[[m]] = pred.sub(unknowns = data.unknown,data = data.train,m = m,factors = factors,Iter = Iter,record = record,alpha.meta = alpha.meta,MU = MU,SIGMA = SIGMA,ALPHA.sub = ALPHA.sub,ALPHA.sup = ALPHA.sup,DISH = DISH,DISH.assign = DISH.assign,P.00 = P.00,P.out = P.out,P.0 = P.0,Taus = Taus,P.Taus = P.Taus,W = W,P.W = P.W,B.u = B.u,B.l = B.l,alpha.out = p.outlier.alpha,beta.out = p.outlier.beta,PCH = PCH,plt = F,new.pathway = new)
  }
  scores.raw = sapply(Result.unknown, function(xxx){max(xxx$Bayes.post)})
  scores = 1/(1+1/(exp(c)*scores.raw^a/(1-scores.raw)^a))
  pred.pathway = sapply(Result.unknown, function(xxx){xxx$pred.path})
  calibration.para = list(unname(a),unname(c));names(calibration.para) = c("a","c")
  annt.result = cbind(data.unknown$metabolites$BIOCHEMICAL,pred.pathway,scores)
  colnames(annt.result) = c("BIOCHEMICAL","pred.pathway","confidence")
  annt.result = as.data.frame(annt.result)
  results = list(annt.result,calibration.para,Result.unknown)
  names(results) = c("annt.result","calibration.para","pred.details")
  return(results)
}
