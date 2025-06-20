# MetabAnnotate
A factor analysis framework for annotating pathway labels of unknown metabolites

## Load example data.
```
library("MetabAnnotate")
data("exampleData")
```

## Step 1. Estimate the number of latent factors

We recommend use BEMA for raw data, and dBEMA for summary-level data.

The BEMA function is part of the package and was downloaded from the GitHub repository maintained by ZhengTracyKe, available at: https://github.com/ZhengTracyKe/BEMA.
```
Y.copy = na.omit(Y)# remove rows including NAs in Y
K = BEMA(data = t(Y.copy),n = ncol(Y.copy),p = nrow(Y.copy),alpha = 0.1)
```
The dBEMA function is part of the package and was obtained from the GitHub repository maintained by Weiqiong Huang, available at: https://github.com/WeiqiongHuang/HiGSS.
```{r}
est.K = dBEMA(Stand.B = B,N = N,alpha = 0.1,n.avg = 5)# 5 simulations to estimate the bulk eigenvalue distribution for time efficiency, could use more for real tasks.
K = est.K$nFactors
```

## Step 2. Transform the input data into the format required by MetabAnnotate
### Raw data.
```
mtb_data = preProcess(IDs = IDs,Y = Y,K = K,impute = T,imputation.methods = "Mean",COMD.column = "CHEMICAL_ID")
```
### Summary-level data.
```
mtb_data = preProcess(IDs = IDs,B = B,N = N,K = K,COMD.column = "CHEMICAL_ID")
```

## Step 3. Implement the annotation function
Implement the annotation function by randomly selecting 20 metabolites for calibration. Generate 20 Gibbs samples, discarding the first 10 iterations as burn-in. It is recommended to run 100 iterations and 50 burn-in for real annotation tasks (Iter = 100, burnin = 50).
```{r}
test = metabAnnotate(mtb_data = mtb_data,size.clbr = 20,Iter = 20,burnin = 10,new.pathway = F,seed = 1)
```

Check the annotation results.
```{r}
head(test$annt.result)
```
