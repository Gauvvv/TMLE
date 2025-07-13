# Install Packages
# install.packages("survey")
# install.packages("SuperLearner")
# install.packages("tmle")
# install.packages("randomForest")

# Load Packages
library(survey)
library(SuperLearner)
library(tmle)
library(randomForest)

# Set random seed to replicate results
set.seed(123)

N=100
num.sim=100
results.tmle.SL1 <- results.tmle.GLM1 <- results.tmle.GLM2 <- results.tmle.GLM3 <-
  NULL
results.mle.SL1 <- results.mle.GLM1 <- results.mle.GLM2 <- NULL
results.psw.SL1 <- results.psw.GLM2 <- NULL
SL.library <- c("SL.glm", "SL.randomForest")

for(i in 1:num.sim){
  
  ######################################################################################
  ### GENERATE DATA
  ######################################################################################
  ### X1=Gender; X2=Therapy; X3=Antidepressant use
  X1 <- rbinom(N, 1, prob=.55)
  X2 <- rbinom(N, 1, prob=.30)
  X3 <- rbinom(N, 1, prob=.25)
  ### Exposure=regular physical exercise
  A <- rbinom(N, 1, plogis(-0.5 + 0.75*X1 + 1*X2 + 1.5*X3))
  ### Outcome=CES-D score
  Y <- 24 - 3*A + 3*X1 -4*X2 - 6*X3 - 1.5*A*X3 + rnorm(N,mean=0,sd=4.5)
  data <- data.frame(cbind(A,X1,X2,X3,Y))
  
  
  ######################################################################################
  ### TMLE approach: Super Learning
  W=cbind(X1,X2,X3)
  tmleSL1 <- tmle(Y, A, W, Q.SL.library = SL.library, g.SL.library = SL.library)
  results.tmle.SL1 <- c(results.tmle.SL1, tmleSL1$estimates$ATE$psi)
  
  ######################################################################################
  ### TMLE approach: GLM, MT misspecification of outcome
  # Misspecified outcome regression: Y ~ A + X1 + X2 + X3
  W=cbind(X1,X2,X3)
  tmleGLM1 <- tmle(Y, A, W, Qform="Y~A+X1+X2+X3", gform="A~X1+X2+X3")
  results.tmle.GLM1 <- c(results.tmle.GLM1, tmleGLM1$estimates$ATE$psi)
  
  ######################################################################################
  ### TMLE approach: GLM, OV misspecification of outcome (X3)
  # Misspecified outcome regression: Y ~ A + X1 + X2
  W=cbind(X1,X2,X3)
  tmleGLM2 <- tmle(Y, A, W, Qform="Y~A+X1+X2", gform="A~X1+X2+X3")
  results.tmle.GLM2 <- c(results.tmle.GLM2, tmleGLM2$estimates$ATE$psi)
  
  ######################################################################################
  ### TMLE approach: GLM, OV misspecification of exposure (X3)
  # Misspecified exposure regression: A ~ X1 + X2
  W=cbind(X1,X2,X3)
  tmleGLM3 <- tmle(Y, A, W, Qform="Y~A+X1+X2+X3+A:X3", gform="A~X1+X2")
  results.tmle.GLM3 <- c(results.tmle.GLM3, tmleGLM3$estimates$ATE$psi)
  
  
  ######################################################################################
  ### G-comp: Super Learning
  newData <- rbind(cbind(A=1, data[,2:4]),
                   cbind(A=0, data[,2:4]))
  SL.fit1 <- SuperLearner(Y=data[,5], X=data[,1:4], SL.library=SL.library,
                          family="gaussian",method="method.NNLS", newX=newData, verbose=TRUE)
  data$Y1.pred <- SL.fit1$SL.predict[1:100]
  data$Y0.pred <- SL.fit1$SL.predict[101:200]
  mle.ATE.SL1 <- mean(data$Y1.pred - data$Y0.pred)
  results.mle.SL1 <- c(results.mle.SL1, mle.ATE.SL1)
  
  ######################################################################################
  ### G-comp: GLM, MT misspecification of outcome
  ### Misspecified outcome regression: Y ~ A + X1 + X2 + X3
  out.reg.M1 <- glm(Y ~ A + X1 + X2 + X3, data = data)
  data$Y_0.M1 <- predict(out.reg.M1, newdata = data.frame(A = 0, X1, X2, X3))
  data$Y_1.M1 <- predict(out.reg.M1, newdata = data.frame(A = 1, X1, X2, X3))
  mle.ATE.GLM1 <- mean(data$Y_1.M1 - data$Y_0.M1)
  results.mle.GLM1 <- c(results.mle.GLM1, mle.ATE.GLM1)
  
  ######################################################################################
  ### G-COMP: GLM, OV misspecification of outcome (X3)
  ### Misspecified outcome regression: Y ~ A + X1 + X2
  out.reg.M2 <- glm(Y ~ A + X1 + X2, data = data)
  data$Y_0.M2 <- predict(out.reg.M2, newdata = data.frame(A = 0, X1, X2, X3))
  data$Y_1.M2 <- predict(out.reg.M2, newdata = data.frame(A = 1, X1, X2, X3))
  mle.ATE.GLM2 <- mean(data$Y_1.M2 - data$Y_0.M2)
  results.mle.GLM2 <- c(results.mle.GLM2, mle.ATE.GLM2)
  
  
  ######################################################################################
  ### IPW: Super Learning
  ### Exposure regression includes X1, X2, X3
  SL.fit3 <- SuperLearner(Y=data[,1], X=data[,2:4], SL.library=SL.library,
                          family="binomial",method="method.NNLS", verbose=TRUE)
  predictions <- cbind(SL.fit3$SL.predict,SL.fit3$library.predict)
  data$pi_1 <- SL.fit3$SL.predict
  data$ps[data$A==1] <- 1/data$pi_1[data$A==1]
  data$ps[data$A==0] <- 1/(1-data$pi_1[data$A==0])
  PS.ATE.SL1 <- weighted.mean(data$Y[data$A==1],w=data$ps[data$A==1])-
    weighted.mean(data$Y[data$A==0],w=data$ps[data$A==0])
  results.psw.SL1 <- c(results.psw.SL1, PS.ATE.SL1)
  
  ######################################################################################
  ### IPW: GLM, OV misspecification of exposure (X3)
  ### Misspecified exposure regression: A ~ X1 + X2
  exp.reg <- glm(A ~ X1 + X2 , family = binomial, data = data)
  data$pi_1 <- predict(exp.reg, type = "response")
  data$ps[data$A==1] <- 1/data$pi_1[data$A==1]
  data$ps[data$A==0] <- 1/(1-data$pi_1[data$A==0])
  PS.ATE.GLM2 <- weighted.mean(data$Y[data$A==1],w=data$ps[data$A==1])-
    weighted.mean(data$Y[data$A==0],w=data$ps[data$A==0])
  results.psw.GLM2 <- c(results.psw.GLM2, PS.ATE.GLM2)
}

######################################################################################
### Results
######################################################################################
### Calculate mean estimates
tmle.SL1 <- mean(results.tmle.SL1)
tmle.GLM1 <- mean(results.tmle.GLM1)
tmle.GLM2 <- mean(results.tmle.GLM2)
tmle.GLM3 <- mean(results.tmle.GLM3)
mle.SL1 <- mean(results.mle.SL1)
mle.GLM1 <- mean(results.mle.GLM1)
mle.GLM2 <- mean(results.mle.GLM2)
psw.SL1 <- mean(results.psw.SL1)
psw.GLM2 <- mean(results.psw.GLM2)
estimates <- rbind(tmle.SL1, tmle.GLM1, tmle.GLM2, tmle.GLM3,
                   mle.SL1, mle.GLM1, mle.GLM2, psw.SL1, psw.GLM2)
### Calculate SEs
tmle.SL1.se <- sd(results.tmle.SL1)
tmle.GLM1.se <- sd(results.tmle.GLM1)
tmle.GLM2.se <- sd(results.tmle.GLM2)
tmle.GLM3.se <- sd(results.tmle.GLM3)
mle.SL1.se <- sd(results.mle.SL1)
mle.GLM1.se <- sd(results.mle.GLM1)
mle.GLM2.se <- sd(results.mle.GLM2)
psw.SL1.se <- sd(results.psw.SL1)
psw.GLM2.se <- sd(results.psw.GLM2)
se <- rbind(tmle.SL1.se, tmle.GLM1.se, tmle.GLM2.se, tmle.GLM3.se,
            mle.SL1.se, mle.GLM1.se, mle.GLM2.se,psw.SL1.se, psw.GLM2.se)
### Calculate 95% CI
tmle.SL1.ci <- quantile(results.tmle.SL1, c(.025, .975))
tmle.GLM1.ci <- quantile(results.tmle.GLM1, c(.025, .975))
tmle.GLM2.ci <- quantile(results.tmle.GLM2, c(.025, .975))
tmle.GLM3.ci <- quantile(results.tmle.GLM3, c(.025, .975))
mle.SL1.ci <- quantile(results.mle.SL1, c(.025, .975))
mle.GLM1.ci <- quantile(results.mle.GLM1, c(.025, .975))
mle.GLM2.ci <- quantile(results.mle.GLM2, c(.025, .975))
psw.SL1.ci <- quantile(results.psw.SL1, c(.025, .975))
psw.GLM2.ci <- quantile(results.psw.GLM2, c(.025, .975))
ci <- rbind(tmle.SL1.ci, tmle.GLM1.ci, tmle.GLM2.ci, tmle.GLM3.ci,
            mle.SL1.ci, mle.GLM1.ci, mle.GLM2.ci,psw.SL1.ci, psw.GLM2.ci)
