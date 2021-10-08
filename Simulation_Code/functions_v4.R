library(survival)
library(rms)
library(mvtnorm)

# data generation function
data_gen_censor = function(n=2000, censor = "50%", overlap = "strong") {
  
  #n=50000; censor = "50%"; overlap = "weak"
  
  # generate continuous covariates X1 X2 X3 
  rho=0.5
  R <- (1-rho) * diag(3) + rho * matrix(1,3,3)
  Xc <- rmvnorm(n, rep(0,3), R)
  # generate binary covariates X4 X5 X6
  Xb = cbind(rbinom(n=n,size=1,p=0.5),rbinom(n=n,size=1,p=0.5),rbinom(n=n,size=1,p=0.5))
  X=cbind(Xc,Xb)
  # regression coefficients
  betac = c(0.2, 0.3, 0.3)
  betab = c(-0.2,-0.3,-0.2)
  # propensity score
  if (overlap=="strong") {
    e = c(plogis(0.4 + Xc %*% c(1*betac) + Xb %*% c(1*betab)))
  }
  
  if (overlap=="medium") {
    e = c(plogis(1.1 + Xc %*% c(3.0*betac) + Xb %*% c(3.0*betab)))
  }
  
  if (overlap=="weak") {
    e = c(plogis(1.8 + Xc %*% c(5*betac) + Xb %*% c(5*betab)))
  }
  # generate intervention status 
  z <- rbinom(n, 1, e)
  
  # generate counterfactural outcome T1 and T0
  U=runif(n)
  scale1 = 0.95
  h=0.4
  LP1 = -1.0-h + c(Xc %*% c(0.4-h,0.2-h,0.1-h)) + c(Xb %*% c(-0.1-h,-0.2-h,-0.3-h))
  T1= (-log(1-U)/(exp(LP1)) )^(1/scale1)
  
  scale0 = 0.95
  LP0 = -1 + c(Xc %*% c(0.4,0.2,0.1)) + c(Xb %*% c(-0.1,-0.2,-0.3))
  T0 = (-log(1-U)/(exp(LP0)))^(1/scale0)
  
  #independent censoring
  if (censor == "50%"){
    scalec = 1
    LPc = -2.1 + c(Xc %*% c(-0.1,0.1,0.2)) + c(Xb %*% c(0.2,-0.2,-0.1))
    C = (-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
  }
  
  if (censor == "25%"){
    scalec = 1
    LPc = -3.3 + c(Xc %*% c(-0.1,0.1,0.2)) + c(Xb %*% c(0.2,-0.2,-0.1))
    C = (-log(1-runif(n))/(exp(LPc)) )^(1/scalec)
  }
  # generate Time
  Time = T1*z+T0*(1-z)
  # indicator I(Time<C)
  Event = as.numeric(Time<C)
  Time = ifelse(Time<C,Time,C)
  # Survival probability for censoring process
  Kc=exp(-exp(c(LPc)) *Time)
  
  return(
    list(
      LP1 = LP1,
      LP0 = LP0,
      e=e,
      C=C,
      X = as.data.frame(cbind(x1=X[,1],x2=X[,2],x3=X[,3],x4=X[,4],x5=X[,5],x6=X[,6], 
                              z, Time = Time, Event =Event)),
      T1 = T1,
      T0 = T0,
      scale1=scale1,scale0=scale0,Kc=Kc+0.000001
    )
  )
  
}

# True WATE values calculated by the survival function
S.f1=function(u) {
  # survival function for the treatment group
  S1.f=function(u) {
    exp(-exp(dat$LP1)*u^dat$scale1)
  }
  # survival function for the control group
  S0.f=function(u) {
    exp(-exp(dat$LP0)*u^dat$scale0)
  }
  # IPW
  Delta.ipw=mean(S1.f(u)-S0.f(u))
  # OW
  Delta.ow = mean(S1.f(u)*dat$e*(1-dat$e)-S0.f(u)*dat$e*(1-dat$e))/mean(dat$e*(1-dat$e))
  # function to obtain WATE by IPW symmetric trimming
  ipwc.f=function(q) {
    ps = dat$e
    keep <- ((ps>= q) & (ps <= (1-q)))
    LP0=dat$LP0[keep]
    LP1=dat$LP1[keep]
    mean(exp(-exp(LP1)*u^dat$scale1)-exp(-exp(LP0)*u^dat$scale0))
  }
  # IPWC
  Delta.ipwc1=ipwc.f(q=0.05)
  Delta.ipwc2=ipwc.f(q=0.1)
  Delta.ipwc3=ipwc.f(q=0.15)
  # function to obtain WATE by IPW asymmetric trimming
  ipwa.f=function(q) {
    ps = dat$e
    ps0 <- ps[dat$X$z == 0]
    ps1 <- ps[dat$X$z == 1]
    lps <- max(min(ps0), min(ps1))
    ups <- min(max(ps0), max(ps1))
    keep <- rep(NA, length(dat$X$z))
    alpha0 <- as.numeric(quantile(ps0, 1-q))
    alpha1 <- as.numeric(quantile(ps1, q))
    keep[dat$X$z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
    keep[dat$X$z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
    keep <- ((ps>= q) & (ps <= (1-q)))
    LP0=dat$LP0[keep]
    LP1=dat$LP1[keep]
    mean(exp(-exp(LP1)*u^dat$scale1)-exp(-exp(LP0)*u^dat$scale0))
  }
  # IPWA
  Delta.ipwa1=ipwa.f(q=0)
  Delta.ipwa2=ipwa.f(q=0.01)
  Delta.ipwa3=ipwa.f(q=0.05)
  
  c(ow=Delta.ow,ipw=Delta.ipw,
    ipwc1=Delta.ipwc1,ipwc2=Delta.ipwc2,ipwc3=Delta.ipwc3,
    ipwa1=Delta.ipwa1,ipwa2=Delta.ipwa2,ipwa3=Delta.ipwa3)
}
# True WATE values calculated by the survival function
S.f2=function(u,dat) {
  X = as.matrix(cbind(1,dat$X[,1:6]))
  colnames(X) = c("int","x1","x2","x3","x4","x5","x6")
  Z = dat$X$z
  Time = dat$X$Time
  Event = dat$X$Event
  # IPW
  ps=dat$e
  w.trt=1/ps
  w.con=1/(1-ps)
  Delta.ipw=(sum(w.trt*Z*Event*as.numeric(Time>u)/dat$Kc)/sum(w.trt*Z*Event/dat$Kc)
             -sum(w.con*(1-Z)*Event*as.numeric(Time>u)/dat$Kc)/sum(w.con*(1-Z)*Event/dat$Kc))
  # OW
  w.trt=1-ps
  w.con=ps
  Delta.ow=(sum(w.trt*Z*Event*as.numeric(Time>u)/dat$Kc)/sum(w.trt*Z*Event/dat$Kc)
            -sum(w.con*(1-Z)*Event*as.numeric(Time>u)/dat$Kc)/sum(w.con*(1-Z)*Event/dat$Kc))
  # IPW with symmetric trimming
  ipwc.f=function(q) {
    ps = dat$e
    w.trt=1/ps
    w.con=1/(1-ps)
    keep <- ((ps>= q) & (ps <= (1-q)))
    # trimming
    LP0=dat$LP0[keep]
    LP1=dat$LP1[keep]
    Z=Z[keep]
    Time=Time[keep]
    Event=Event[keep]
    w.trt=w.trt[keep]
    w.con=w.con[keep]
    Kc=dat$Kc[keep]
    (sum(w.trt*Z*Event*as.numeric(Time>u)/Kc)/sum(w.trt*Z*Event/Kc)
               -sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc)/sum(w.con*(1-Z)*Event/Kc))
  }
  Delta.ipwc1=ipwc.f(q=0.05)
  Delta.ipwc2=ipwc.f(q=0.1)
  Delta.ipwc3=ipwc.f(q=0.15)
  # IPW with Asymmetric trimming
  ipwa.f=function(q) {
    ps = dat$e
    ps0 <- ps[dat$X$z == 0]
    ps1 <- ps[dat$X$z == 1]
    lps <- max(min(ps0), min(ps1))
    ups <- min(max(ps0), max(ps1))
    keep <- rep(NA, length(dat$X$z))
    alpha0 <- as.numeric(quantile(ps0, 1-q))
    alpha1 <- as.numeric(quantile(ps1, q))
    keep[dat$X$z == 0] <- ((ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
    keep[dat$X$z == 1] <- ((ps1 >= alpha1) & (ps1 >= lps) & (ps1 <= ups))
    keep <- ((ps>= q) & (ps <= (1-q)))
    # trimming
    LP0=dat$LP0[keep]
    LP1=dat$LP1[keep]
    Z=Z[keep]
    Time=Time[keep]
    Event=Event[keep]
    w.trt=w.trt[keep]
    w.con=w.con[keep]
    Kc=dat$Kc[keep]
    (sum(w.trt*Z*Event*as.numeric(Time>u)/Kc)/sum(w.trt*Z*Event/Kc)
    -sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc)/sum(w.con*(1-Z)*Event/Kc))
  }
  # IPWA
  Delta.ipwa1=ipwa.f(q=0)
  Delta.ipwa2=ipwa.f(q=0.01)
  Delta.ipwa3=ipwa.f(q=0.05)
  
  c(ow=Delta.ow,ipw=Delta.ipw,
    ipwc1=Delta.ipwc1,ipwc2=Delta.ipwc2,ipwc3=Delta.ipwc3,
    ipwa1=Delta.ipwa1,ipwa2=Delta.ipwa2,ipwa3=Delta.ipwa3)
}

# estimating WATE with an overlap weight 
OW.f=function(u=6,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1-ps
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = mean(ps*(1-ps))
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>u) - tau1))/Kc.trt) %*% (-1*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = ps
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Etau = mean(ps*(1-ps))
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>u) - tau0))/Kc.con) %*% (ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# estimating WATE with an inverse propensity score weight
IPW.f=function(u=6,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# estimating WATE with symmetric triming
IPWC.f=function(u=6,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  keep <- ((ps>= q) & (ps <= (1-q)))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

# estimating WATE with asymmetric triming
IPWA.f=function(u=6,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  ps0 <- ps[Z == 0]
  ps1 <- ps[Z == 1]
  lps <- max(min(ps0), min(ps1))
  ups <- min(max(ps0), max(ps1))
  
  # PS Asymmetric trimming
  keep <- rep(NA, length(Z))
  alpha0 <- as.numeric(quantile(ps0, 1-q))
  alpha1 <- as.numeric(quantile(ps1, q))
  keep[Z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
  keep[Z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}


# estimating WATE with an overlap weight 
OW.f2=function(u=6,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1-ps
  tau1 = sum(w.trt*Z*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z)
  Etau = mean(ps*(1-ps))
  Htheta1 = 1/n * t((w.trt * Z * (Time>u) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * (Time>u) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*(Time>u))/Kc.trt) %*% (-1*ps*(1-ps)*X) - 1/n* t(tau1*Z) %*% (-1*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = ps
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>u) ))/Kc.con) %*% (ps*(1-ps)*X) - 1/n* t(tau0*(1-Z)) %*% (ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>u))/Kc.trt - w.con*(1-Z)*((Time>u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# estimating WATE with an inverse propensity score weight
IPW.f2=function(u=6,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>u)))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>u))/Kc.trt - w.con*(1-Z)*((Time>u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# estimating WATE with symmetric triming
IPWC.f2=function(u=6,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  keep <- ((ps>= q) & (ps <= (1-q)))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>u) ))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>u))/Kc.trt - w.con*(1-Z)*((Time>u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

# estimating WATE with asymmetric triming
IPWA.f2=function(u=6,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  ps0 <- ps[Z == 0]
  ps1 <- ps[Z == 1]
  lps <- max(min(ps0), min(ps1))
  ups <- min(max(ps0), max(ps1))
  
  # PS Asymmetric trimming
  keep <- rep(NA, length(Z))
  alpha0 <- as.numeric(quantile(ps0, 1-q))
  alpha1 <- as.numeric(quantile(ps1, q))
  keep[Z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
  keep[Z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * X) / n
  #Sbeta = (Z-ps)*X
  # Taylor Expansion for censoring score in treatment group
  Etheta1 = crossprod(sqrt(Z*(Time^(gamma1.est))*exp(c(X %*% theta1.est))) * X)/n
  Egamma1 = 1/n * sum( Z* ((1-Event)/(gamma1.est)^2 + Time^gamma1.est *(log(Time))^2 * exp(c(X %*% theta1.est))  )  )
  #Stheta1 = Z * (X*( (1-Event) - Time^gamma1.est * exp(c(X %*% theta1.est))))
  #Sgamma1 = Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est *log(Time) *exp(c(X %*% theta1.est))   )
  # Taylor Expansion for censoring score in control group
  Etheta0 = crossprod(sqrt((1-Z)*(Time^(gamma0.est))*exp(c(X %*% theta0.est))) * X)/n
  Egamma0 = 1/n * sum( (1-Z) *( (1-Event)/(gamma0.est)^2 + Time^gamma0.est *(log(Time))^2 * exp(c(X %*% theta0.est)) )   )
  #Stheta0 = (1-Z) * (X*( (1-Event) - Time^gamma0.est * exp(c(X %*% theta0.est))))
  #Sgamma0 = (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est *log(Time) *exp(c(X %*% theta0.est))   )
  # Taylor Expansion for the Survival function in treatment group
  w.trt = 1/ps
  tau1 = sum(w.trt*Z*as.numeric(Time>u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*X) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*X)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>u) ))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*X) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*X)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>u))/Kc.trt - w.con*(1-Z)*((Time>u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*X) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

# main estimating function combining OW, IPW, IPWC, IPWA
estimation.f = function(X,Z,Time,Event,uvec) {
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  #Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  #Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  OW.est = rbind(OW.f(u=uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f(u=uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f(u=uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f(u=uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model))
  IPW.est = rbind(IPW.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model))
  IPWC1.est=rbind(IPWC.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05))
  IPWC2.est=rbind(IPWC.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1))
  IPWC3.est=rbind(IPWC.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15))
  IPWA1.est=rbind(IPWA.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0))
  IPWA2.est=rbind(IPWA.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01))
  IPWA3.est=rbind(IPWA.f(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05))
  
  
  
  OW.est2 = rbind(OW.f2(u=uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f2(u=uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f2(u=uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                 OW.f2(u=uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model))
  IPW.est2 = rbind(IPW.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model),
                  IPW.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model))
  IPWC1.est2=rbind(IPWC.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWC.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05))
  IPWC2.est2=rbind(IPWC.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1),
                  IPWC.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1))
  IPWC3.est2=rbind(IPWC.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15),
                  IPWC.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15))
  IPWA1.est2=rbind(IPWA.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0),
                  IPWA.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0))
  IPWA2.est2=rbind(IPWA.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01),
                  IPWA.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01))
  IPWA3.est2=rbind(IPWA.f2(uvec[1],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f2(uvec[2],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f2(uvec[3],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05),
                  IPWA.f2(uvec[4],X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05))
  list(ow=OW.est,ipw=IPW.est,
       ipwc1=IPWC1.est,ipwc2=IPWC2.est,ipwc3=IPWC3.est,
       ipwa1=IPWA1.est,ipwa2=IPWA2.est,ipwa3=IPWA3.est,
       sow=OW.est2,sipw=IPW.est2,
       sipwc1=IPWC1.est2,sipwc2=IPWC2.est2,sipwc3=IPWC3.est2,
       sipwa1=IPWA1.est2,sipwa2=IPWA2.est2,sipwa3=IPWA3.est2)
}


# estimating function for drawing the efficiency-by-uvec plot
estimation.f.plot = function(X,Z,Time,Event,u) {
  data = as.data.frame(cbind(X,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  ps.formula = as.formula(paste("Z~",paste(colnames(X),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T)
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  #Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T)
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  #Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  
  ## estimators
  OW.est = OW.f(u=u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
  IPW.est = IPW.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
  IPWC1.est=IPWC.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05)
  IPWC2.est=IPWC.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1)
  IPWC3.est=IPWC.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15)
  IPWA1.est=IPWA.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0)
  IPWA2.est=IPWA.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01)
  IPWA3.est=IPWA.f(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05)
  
  OW.est2 = OW.f2(u=u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
  IPW.est2 = IPW.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
  IPWC1.est2=IPWC.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05)
  IPWC2.est2=IPWC.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.1)
  IPWC3.est2=IPWC.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.15)
  IPWA1.est2=IPWA.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0)
  IPWA2.est2=IPWA.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01)
  IPWA3.est2=IPWA.f2(u,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.05)
  
  data.frame(ow=OW.est,ipw=IPW.est,
       ipwc1=IPWC1.est[1:2],ipwc2=IPWC2.est[1:2],ipwc3=IPWC3.est[1:2],
       ipwa1=IPWA1.est[1:2],ipwa2=IPWA2.est[1:2],ipwa3=IPWA3.est[1:2],
       sow=OW.est2,sipw=IPW.est2,
       sipwc1=IPWC1.est2[1:2],sipwc2=IPWC2.est2[1:2],sipwc3=IPWC3.est2[1:2],
       sipwa1=IPWA1.est2[1:2],sipwa2=IPWA2.est2[1:2],sipwa3=IPWA3.est2[1:2])
}



