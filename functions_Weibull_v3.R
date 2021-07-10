library(survival)

# Type 1 overlap weight with asymtotic variance estimator
OW.f=function(u=6,W,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
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
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = mean(ps*(1-ps))
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>=u) - tau1))/Kc.trt) %*% (-1*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = ps
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Etau = mean(ps*(1-ps))
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>=u) - tau0))/Kc.con) %*% (ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>=u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>=u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# Type 1 IPW with asymtotic variance estimator
IPW.f=function(u=6,W,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
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
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>=u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>=u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>=u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>=u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# Type 1 symmetric triming with asymtotic variance estimator
IPWC.f=function(u=6,W,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  keep <- ((ps>= q) & (ps <= (1-q)))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  W=W[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,W,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>=u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>=u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>=u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>=u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

# Type 1 asymmetric triming with asymtotic variance estimator
IPWA.f=function(u=6,W,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  ps0 <- ps[Z == 0]
  ps1 <- ps[Z == 1]
  lps <- max(min(ps0), min(ps1))
  ups <- min(max(ps0), max(ps1))
  
  # PS Asymmetric trimming
  keep <- rep(NA, length(Z))
  alpha0 <- as.numeric(quantile(ps0, 1-q/100))
  alpha1 <- as.numeric(quantile(ps1, q/100))
  keep[Z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
  keep[Z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  W=W[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,W,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*Event*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * Event * ((Time>=u) - tau1) )/Kc.trt * Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*Event*((Time>=u) - tau1))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
  Htheta0 = 1/n * t((w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * Event * ((Time>=u) - tau0) )/Kc.con * Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*Event*((Time>=u) - tau0))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*Event*((Time>=u)-tau1)/Kc.trt - w.con*(1-Z)*Event*((Time>=u)-tau0)/Kc.con
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}


# Type 2 overlap weight with asymtotic variance estimator
OW.f2=function(u=6,W,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
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
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z)
  Etau = mean(ps*(1-ps))
  Htheta1 = 1/n * t((w.trt * Z * (Time>=u) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * (Time>=u) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*(Time>=u))/Kc.trt) %*% (-1*ps*(1-ps)*W) - 1/n* t(tau1*Z) %*% (-1*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = ps
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>=u) ))/Kc.con) %*% (ps*(1-ps)*W) - 1/n* t(tau0*(1-Z)) %*% (ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>=u))/Kc.trt - w.con*(1-Z)*((Time>=u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# Type 2 IPW with asymtotic variance estimator
IPW.f2=function(u=6,W,X,Z,Time,Event,ps.model,cen.trt.model,cen.con.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
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
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>=u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>=u)))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>=u))/Kc.trt - w.con*(1-Z)*((Time>=u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est) 
}

# Type 2 symmetric triming with asymtotic variance estimator
IPWC.f2=function(u=6,W,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  keep <- ((ps>= q) & (ps <= (1-q)))
  ptrim <- 1 - mean(keep)
  # trim the dataset
  X=X[keep,]
  W=W[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,W,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>=u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>=u) ))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>=u))/Kc.trt - w.con*(1-Z)*((Time>=u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

# Type 2 asymmetric triming with asymtotic variance estimator
IPWA.f2=function(u=6,W,X,Z,Time,Event,q,ps.model,...) {
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  # propensity score
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
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
  W=W[keep,]
  Z=Z[keep]
  Time=Time[keep]
  Event=Event[keep]
  data = as.data.frame(cbind(X,W,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  # censoring function
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  Kc.trt = exp(-exp(c(X %*% theta1.est)) *u^gamma1.est)
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  Kc.con = exp(-exp(c(X %*% theta0.est)) *u^gamma0.est)
  # variance estimate
  # Tayler Expansion for PS model
  beta.est = ps.model$coefficients
  n=length(Z)
  Ebeta = crossprod(sqrt(ps*(1-ps)) * W) / n
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
  tau1 = sum(w.trt*Z*as.numeric(Time>=u)/Kc.trt)/sum(w.trt*Z)
  Etau = 1
  Htheta1 = 1/n * t((w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * exp(c(X %*% theta1.est))) %*% X 
  Hgamma1 = 1/n * sum( (w.trt * Z * ((Time>=u) ) )/Kc.trt * u^gamma1.est * log(u) * exp(c(X %*% theta1.est)))
  Hbeta1  = 1/n * t((Z*((Time>=u) ))/Kc.trt) %*% (-1/(ps^2)*ps*(1-ps)*W) - 1/n* t(tau1*Z) %*% (-1/(ps^2)*ps*(1-ps)*W)
  # Taylor Expansion for the Survival function in control group
  w.con = 1/(1-ps)
  tau0 = sum(w.con*(1-Z)*as.numeric(Time>=u)/Kc.con)/sum(w.con*(1-Z))
  Htheta0 = 1/n * t((w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * exp(c(X %*% theta0.est))) %*% X 
  Hgamma0 = 1/n * sum( (w.con * (1-Z) * ((Time>=u) ) )/Kc.con * u^gamma0.est * log(u) * exp(c(X %*% theta0.est)))
  Hbeta0  = 1/n * t(((1-Z)*((Time>=u) ))/Kc.con) %*% (1/((1-ps)^2) *ps*(1-ps)*W) - 1/n* t(tau0*(1-Z)) %*% (1/((1-ps)^2) *ps*(1-ps)*W)
  # obtain variance estimator
  Itau = w.trt*Z*((Time>=u))/Kc.trt - w.con*(1-Z)*((Time>=u))/Kc.con - w.trt*tau1*Z + w.con*tau0*(1-Z)
  Itheta1 = ( Z*( ((1-Event) - Time^gamma1.est * exp( c(X %*% theta1.est)))*X  ) ) %*% t(Htheta1 %*% solve(Etheta1))
  Itheta0 = ( (1-Z)*( ((1-Event) - Time^gamma0.est * exp( c(X %*% theta0.est)))*X  ) ) %*% t(Htheta0 %*% solve(Etheta0))
  Igamma1 = Hgamma1 * 1/Egamma1 * Z * ( (1-Event)*(1/gamma1.est + log(Time)) - Time^gamma1.est * log(Time) * exp(c(X %*% theta1.est))  )
  Igamma0 = Hgamma0 * 1/Egamma0 * (1-Z) * ( (1-Event)*(1/gamma0.est + log(Time)) - Time^gamma0.est * log(Time) * exp(c(X %*% theta0.est))  )
  Ibeta = ((Z-ps)*W) %*% t((Hbeta1 -Hbeta0) %*% solve(Ebeta))
  I = 1/Etau *(c(Itau)+c(Itheta1)-c(Itheta0)+c(Igamma1)-c(Igamma0)+c(Ibeta))
  var.est=sum(I^2)/(n^2)
  # point estimate
  p.est = tau1-tau0
  c(point=p.est,variance=var.est,ptrim=ptrim) 
}

SurvEffectWeibull=function(Data,
                           t=60,
                           Treatment,
                           SurvTime,
                           Status,
                           PS.formula,
                           Censor.formula,
                           Type=1,
                           Method="IPW",
                           alpha=0.1,
                           q=0.01) {
  # Data=rhc
  # t=60
  # Treatment="swang1"
  # SurvTime="survtime"
  # Status="death"
  # PS.formula=as.formula(paste(Treatment,"~",paste(Covariates,collapse="+"),sep=""))
  # Censor.formula=as.formula(paste("Surv(",SurvTime,",I(1-", Status,"))","~",paste(Covariates[-c(1:2)],collapse="+"),sep=""))
  # Type=1
  # Method="IPW"
  # alpha=0.1
  # q=0.01
  # identify the baseline covariates, treatment indicator, survival time, 
  # and censoring indicator
  W=model.matrix(PS.formula,data=Data)
  colnames(W)=paste("W",1:dim(W)[2],sep="")
  X=model.matrix(Censor.formula,data=Data)
  colnames(X)=paste("X",1:dim(X)[2],sep="")
  Z=Data[,Treatment]
  Time=Data[,SurvTime]
  Event=Data[,Status]
  # reconstruct the dataset
  data = as.data.frame(cbind(X,W,Z,Time,Event))
  data.trt = subset(data,Z==1)
  data.con = subset(data,Z==0)
  # propensity score
  ps.formula = as.formula(paste("Z~",paste(colnames(W),collapse = "+"),"-1",sep=""))
  ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
  ps = 1/(1+exp(-c(W %*% ps.model$coefficients))) 
  # censoring function
  surv.formula = as.formula(paste("Surv(Time, I(1-Event)) ~ -1 +",paste(colnames(X),collapse = "+"),sep=""))
  cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
  gamma1.est = 1/cen.trt.model$scale
  cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=350))
  theta0.est = -cen.con.model$coefficients/cen.con.model$scale
  gamma0.est = 1/cen.con.model$scale
  if (Method=="IPW") {
    if (Type==1) {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*Event*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
      # return(tau1-tau0)
      res=IPW.f(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
    } else {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *t^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *t^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z))
      # return(tau1-tau0)
      res=IPW.f2(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
    }
  } else if (Method=="OW") {
    if (Type==1) {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
      # w.trt = 1-ps
      # tau1 = sum(w.trt*Z*Event*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
      # w.con = ps
      # tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
      # return(tau1-tau0)
      res=OW.f(u=t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
    } else {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *t^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *t^gamma0.est)
      # w.trt = 1-ps
      # tau1 = sum(w.trt*Z*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z)
      # w.con = ps
      # tau0 = sum(w.con*(1-Z)*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z))
      # return(tau1-tau0)
      res=OW.f2(u=t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,cen.trt.model=cen.trt.model,cen.con.model=cen.con.model)
    }
  } else if (Method=="Symmetric") {
    # keep <- ((ps>= alpha) & (ps <= (1-alpha)))
    # X=X[keep,]
    # Z=Z[keep]
    # Time=Time[keep]
    # Event=Event[keep]
    # data = as.data.frame(cbind(X,Z,Time,Event))
    # data.trt = subset(data,Z==1)
    # data.con = subset(data,Z==0)
    # ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
    # ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
    # # censoring function
    # cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=300))
    # theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
    # gamma1.est = 1/cen.trt.model$scale
    # cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=300))
    # theta0.est = -cen.con.model$coefficients/cen.con.model$scale
    # gamma0.est = 1/cen.con.model$scale
    if (Type==1) {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*Event*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
      # return(tau1-tau0)
      res=IPWC.f(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=alpha)
    } else {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *t^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *t^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z))
      # return(tau1-tau0)
      res=IPWC.f2(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=alpha)
    }
  } else if (Method=="Asymmetric") {
    # ps0 <- ps[Z == 0]
    # ps1 <- ps[Z == 1]
    # lps <- max(min(ps0), min(ps1))
    # ups <- min(max(ps0), max(ps1))
    # keep <- rep(NA, length(Z))
    # alpha0 <- as.numeric(quantile(ps0, 1-q))
    # alpha1 <- as.numeric(quantile(ps1, q))
    # keep[Z == 0] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
    # keep[Z == 1] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))
    # # trim the dataset
    # X=X[keep,]
    # Z=Z[keep]
    # Time=Time[keep]
    # Event=Event[keep]
    # data = as.data.frame(cbind(X,Z,Time,Event))
    # data.trt = subset(data,Z==1)
    # data.con = subset(data,Z==0)
    # # propensity score
    # ps.model = glm(ps.formula,data=data,family=binomial(link="logit"))
    # ps = 1/(1+exp(-c(X %*% ps.model$coefficients))) 
    # # censoring function
    # cen.trt.model = survreg(surv.formula, data=data.trt,dist='weibull',score=T,control=survreg.control(maxiter=150))
    # theta1.est = -cen.trt.model$coefficients/cen.trt.model$scale
    # gamma1.est = 1/cen.trt.model$scale
    # cen.con.model = survreg(surv.formula, data=data.con,dist='weibull',score=T,control=survreg.control(maxiter=150))
    # theta0.est = -cen.con.model$coefficients/cen.con.model$scale
    # gamma0.est = 1/cen.con.model$scale
    if (Type==1) {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *Time^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *Time^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*Event*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z*Event/Kc.trt)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*Event*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z)*Event/Kc.con)
      # return(tau1-tau0)
      res=IPWA.f(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01)
    } else {
      # Kc.trt = exp(-exp(c(X %*% theta1.est)) *t^gamma1.est)
      # Kc.con = exp(-exp(c(X %*% theta0.est)) *t^gamma0.est)
      # w.trt = 1/ps
      # tau1 = sum(w.trt*Z*as.numeric(Time>=t)/Kc.trt)/sum(w.trt*Z)
      # w.con = 1/(1-ps)
      # tau0 = sum(w.con*(1-Z)*as.numeric(Time>=t)/Kc.con)/sum(w.con*(1-Z))
      # return(tau1-tau0)
      res=IPWA.f2(t,W=W,X=X,Z=Z,Time=Time,Event=Event,ps.model=ps.model,q=0.01)
    }
  }
  output=c(Estimate=res[1],SE=sqrt(res[2]),CI.lower=res[1]-1.96*sqrt(res[2]),CI.upper=res[1]+1.96*sqrt(res[2]))
  names(output)=c("Estimate","SE","CI.lower","CI.upper")
  output
}

















