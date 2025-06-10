# MetabAnnotate
A factor analysis framework for annotating pathway labels of unknown metabolites

Load necessary packages
```
library(MSNIMBLE)
```

Randomly select a set of metabolites for calibration
```
clbr <- sample(x = 1:nrow(data$metabolites),size = size.clbr,replace = F)
train <- setdiff(1:nrow(data$metabolites),clbr)
```

Estimate concentration parameter
```
alpha.meta <- estAlpha(IDs = data.train$metabolites)
```
Draw Gibbs samples from the HDP process
```
Result = fitHDP(data = data.train,Iter=10,record=5)
```
Transfer the Gibbs samples to the data format that metabAnnotate requires
```
classifierInput=summarizeResult(fitOutput = Result,K = K)
```
Run the annotation algorithm
```
result=metabAnnotate(data.test = data.test1,data.unknown = data.test2,input = classifierInput,factors = 1:K,alpha.meta = alpha.meta,test = T,new = T)
```
  
