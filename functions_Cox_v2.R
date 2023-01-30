SurvEffectWithCox=function(Data,
                           t=60,
                           Treatment,
                           SurvTime,
                           Status,
                           PS.formula,
                           Censor.formula,
                           Type=1,
                           Method="IPTW",
                           alpha=0.05,
                           q=0.01) {
  CensorScoreFun=function(CoxModel,TimeVec) {
    LinearPredictor=CoxModel$linear.predictors
    BaselineHazardForm=basehaz(CoxModel,centered=T) # baseline hazard data frame
    BaselineHazard=sapply(TimeVec, function(x) { # obtain the baseline hazard for Time Vec
      BaselineHazardForm[which.min(abs(BaselineHazardForm$time-x))[1],"hazard"]
    })
    CensorScore=  exp(-BaselineHazard*exp(LinearPredictor)) # obtain censoring probability
    return(CensorScore)
  }
  
  # estimate propensity scores
  PS.model=glm(PS.formula,data=Data,family=binomial(link="logit"))
  PS = 1/(1+exp(-c(PS.model$linear.predictors)))
  # Cox Model for the treatment group
  Data.trt=subset(Data,Data[,Treatment]==1)
  assign("Data.trt", Data.trt, envir = .GlobalEnv)
  Censor.trt.model = survival::coxph(Censor.formula,data=Data.trt) 
  # Cox Model for the control group
  Data.con=subset(Data,Data[,Treatment]==0)
  assign("Data.con", Data.con, envir = .GlobalEnv)
  Censor.con.model = survival::coxph(Censor.formula,data=Data.con) 
  if (Method=="IPTW") {
    w.trt = 1/PS[which(Data[,Treatment]==1)]
    w.con = 1/(1-PS[which(Data[,Treatment]==0)])
    if (Type==1) {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=Data.trt[,SurvTime])
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=Data.con[,SurvTime])
      delta1 = 1-sum(w.trt*Data.trt[,Status]*as.numeric(Data.trt[,SurvTime]<=t)/Kc.trt)/sum(w.trt)
      delta0 = 1-sum(w.con*Data.con[,Status]*as.numeric(Data.con[,SurvTime]<=t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    } else {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=rep(t,length(Censor.trt.model$y)))
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=rep(t,length(Censor.con.model$y)))
      delta1 = sum(w.trt*as.numeric(Data.trt[,SurvTime]>t)/Kc.trt)/sum(w.trt)
      delta0 = sum(w.con*as.numeric(Data.con[,SurvTime]>t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    }
  } else if (Method=="OW") {
    w.trt = 1-PS[which(Data[,Treatment]==1)]
    w.con = PS[which(Data[,Treatment]==0)]
    if (Type==1) {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=Data.trt[,SurvTime])
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=Data.con[,SurvTime])
      delta1 = 1-sum(w.trt*Data.trt[,Status]*as.numeric(Data.trt[,SurvTime]<=t)/Kc.trt)/sum(w.trt)
      delta0 = 1-sum(w.con*Data.con[,Status]*as.numeric(Data.con[,SurvTime]<=t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    } else {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=rep(t,length(Censor.trt.model$y)))
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=rep(t,length(Censor.con.model$y)))
      delta1 = sum(w.trt*as.numeric(Data.trt[,SurvTime]>t)/Kc.trt)/sum(w.trt)
      delta0 = sum(w.con*as.numeric(Data.con[,SurvTime]>t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    }
  } else if (Method=="Symmetric") {
    keep <- ((PS>= alpha) & (PS <= (1-alpha)))
    Data.trim=Data[keep,]
    ## refit propensity and censoring models
    PS.model=glm(PS.formula,data=Data.trim,family=binomial(link="logit"))
    PS = 1/(1+exp(-c(PS.model$linear.predictors)))
    Data.trim.trt=subset(Data.trim,Data.trim[,Treatment]==1)
    Censor.trt.model = coxph(Censor.formula,data=Data.trim.trt) 
    Data.trim.con=subset(Data.trim,Data.trim[,Treatment]==0)
    Censor.con.model = coxph(Censor.formula,data=Data.trim.con) 
    assign("Data.trim.trt", Data.trim.trt, envir = .GlobalEnv)
    assign("Data.trim.con", Data.trim.con, envir = .GlobalEnv)
    if (Type==1) {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=Data.trim.trt[,SurvTime])
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=Data.trim.con[,SurvTime])
      w.trt = 1/PS[which(Data.trim[,Treatment]==1)]
      delta1 = 1-sum(w.trt*Data.trim.trt[,Status]*as.numeric(Data.trim.trt[,SurvTime]<=t)/Kc.trt)/sum(w.trt)
      w.con = 1/(1-PS[which(Data.trim[,Treatment]==0)])
      delta0 = 1-sum(w.con*Data.trim.con[,Status]*as.numeric(Data.trim.con[,SurvTime]<=t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    } else {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=rep(t,length(Censor.trt.model$y)))
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=rep(t,length(Censor.con.model$y)))
      w.trt = 1/PS[which(Data.trim[,Treatment]==1)]
      delta1 = sum(w.trt*as.numeric(Data.trim.trt[,SurvTime]>t)/Kc.trt)/sum(w.trt)
      w.con = 1/(1-PS[which(Data.trim[,Treatment]==0)])
      delta0 = sum(w.con*as.numeric(Data.trim.con[,SurvTime]>t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    }
  } else if (Method=="Asymmetric") {
    PS0 <- PS[Data[,Treatment] == 0]
    PS1 <- PS[Data[,Treatment] == 1]
    LowerPS <- max(min(PS0), min(PS1))
    UpperPS <- min(max(PS0), max(PS1))
    keep <- rep(NA, length(PS))
    alpha0 <- as.numeric(quantile(PS0, 1-q))
    alpha1 <- as.numeric(quantile(PS1, q))
    keep[Data[,Treatment] == 0] <- ((PS0 >= alpha1) & (PS0 <= alpha0) & (PS0 >= LowerPS) & (PS0 <= UpperPS))
    keep[Data[,Treatment] == 1] <- ((PS1 >= alpha1) & (PS1 <= alpha0) & (PS1 >= LowerPS) & (PS1 <= UpperPS))
    Data.trim=Data[keep,]
    ## refit propensity and censoring models
    PS.model=glm(PS.formula,data=Data.trim,family=binomial(link="logit"))
    PS = 1/(1+exp(-c(PS.model$linear.predictors)))
    Data.trim.trt=subset(Data.trim,Data.trim[,Treatment]==1)
    Censor.trt.model = coxph(Censor.formula,data=Data.trim.trt) 
    Data.trim.con=subset(Data.trim,Data.trim[,Treatment]==0)
    Censor.con.model = coxph(Censor.formula,data=Data.trim.con) 
    assign("Data.trim.trt", Data.trim.trt, envir = .GlobalEnv)
    assign("Data.trim.con", Data.trim.con, envir = .GlobalEnv)
    if (Type==1) {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=Data.trim.trt[,SurvTime])
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=Data.trim.con[,SurvTime])
      w.trt = 1/PS[which(Data.trim[,Treatment]==1)]
      delta1 = 1-sum(w.trt*Data.trim.trt[,Status]*as.numeric(Data.trim.trt[,SurvTime]<=t)/Kc.trt)/sum(w.trt)
      w.con = 1/(1-PS[which(Data.trim[,Treatment]==0)])
      delta0 = 1-sum(w.con*Data.trim.con[,Status]*as.numeric(Data.trim.con[,SurvTime]<=t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    } else {
      Kc.trt=CensorScoreFun(CoxModel=Censor.trt.model,TimeVec=rep(t,length(Censor.trt.model$y)))
      Kc.con=CensorScoreFun(CoxModel=Censor.con.model,TimeVec=rep(t,length(Censor.con.model$y)))
      w.trt = 1/PS[which(Data.trim[,Treatment]==1)]
      delta1 = sum(w.trt*as.numeric(Data.trim.trt[,SurvTime]>t)/Kc.trt)/sum(w.trt)
      w.con = 1/(1-PS[which(Data.trim[,Treatment]==0)])
      delta0 = sum(w.con*as.numeric(Data.trim.con[,SurvTime]>t)/Kc.con)/sum(w.con)
      return(delta1-delta0)
    }
  } else {
    stop("The Method should be one of the following four choices: IPTW, OW, Symmetric, and Asymmetric.")
  }
}


AllSurvEffect=function(Data=rhc,t=60,Treatment="swang1",SurvTime="survtime",Status="death",
                       PS.formula,Censor.formula,alpha=0.01,q=0.01) {
  Delta1=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=1,Method="IPTW")
  Delta2=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=1,Method="OW")
  Delta3=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=1,Method="Symmetric",alpha=alpha)
  Delta4=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=1,Method="Asymmetric",q=q)
  Delta5=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=2,Method="IPTW")
  Delta6=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=2,Method="OW")
  Delta7=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=2,Method="Symmetric",alpha=alpha)
  Delta8=SurvEffect(Data=Data,t=t,Treatment=Treatment,SurvTime=SurvTime,Status=Status,
                    PS.formula=PS.formula,Censor.formula=Censor.formula,Type=2,Method="Asymmetric",q=q)
  output=c(Delta1,Delta2,Delta3,Delta4,Delta5,Delta6,Delta7,Delta8)
  names(output)=paste(rep(c("IPTW","OW","Symmetric","Asymmetric"),2),rep(c(1,2),each=4),sep="")
  output
}
