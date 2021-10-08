rm(list=ls())
set.seed(2021)
source("functions_v4.R")

# define uvec
uvec.mat=matrix(0,ncol=4,nrow=20)
for (i in (1:20)) {
  dat = data_gen_censor(n=500000,overlap="strong",censor = "50%")
  uvec.mat[i,]=quantile(dat$X$Time,prob=c(0.2,0.4,0.6,0.8)) 
}
uvec=apply(uvec.mat,2,mean)

##### calculate "true" WATEs
dat = data_gen_censor(n=500000,overlap="strong",censor = "50%")
Delta.true=rbind(S.f1(u=uvec[1]),S.f1(u=uvec[2]),S.f1(u=uvec[3]),S.f1(u=uvec[4]))
for (i in (1:19)) {
  dat = data_gen_censor(n=500000,overlap="strong",censor = "50%")
  Delta.true=Delta.true+rbind(S.f1(u=uvec[1]),S.f1(u=uvec[2]),S.f1(u=uvec[3]),S.f1(u=uvec[4]))
}
Delta.true=Delta.true/20
name.D=c(colnames(Delta.true),paste("s",colnames(Delta.true),sep=""))
Delta.true=cbind(Delta.true,Delta.true)
colnames(Delta.true)=name.D

#### simulations
replication=2000
p1=p2=p3=p4=v1=v2=v3=v4=matrix(NA,ncol=16,nrow=replication)
for (i in (1:replication)) {
  dat=data_gen_censor(n=2000,overlap="strong",censor = "50%")
  X = as.matrix(cbind(1,dat$X[,1:6]))
  colnames(X) = c("int","x1","x2","x3","x4","x5","x6")
  Z = dat$X$z
  Time = dat$X$Time
  Event = dat$X$Event
  est=estimation.f(X=X,Z=Z,Time=Time,Event=Event,uvec=uvec)
  # 1st time point
  p1[i,]=c(est$ow[1,1],est$ipw[1,1],est$ipwc1[1,1],est$ipwc2[1,1],est$ipwc3[1,1],est$ipwa1[1,1],est$ipwa2[1,1],est$ipwa3[1,1],
           est$sow[1,1],est$sipw[1,1],est$sipwc1[1,1],est$sipwc2[1,1],est$sipwc3[1,1],est$sipwa1[1,1],est$sipwa2[1,1],est$sipwa3[1,1])
  v1[i,]=c(est$ow[1,2],est$ipw[1,2],est$ipwc1[1,2],est$ipwc2[1,2],est$ipwc3[1,2],est$ipwa1[1,2],est$ipwa2[1,2],est$ipwa3[1,2],
           est$sow[1,2],est$sipw[1,2],est$sipwc1[1,2],est$sipwc2[1,2],est$sipwc3[1,2],est$sipwa1[1,2],est$sipwa2[1,2],est$sipwa3[1,2])
  # 2nd time point
  p2[i,]=c(est$ow[2,1],est$ipw[2,1],est$ipwc1[2,1],est$ipwc2[2,1],est$ipwc3[2,1],est$ipwa1[2,1],est$ipwa2[2,1],est$ipwa3[2,1],
           est$sow[2,1],est$sipw[2,1],est$sipwc1[2,1],est$sipwc2[2,1],est$sipwc3[2,1],est$sipwa1[2,1],est$sipwa2[2,1],est$sipwa3[2,1])
  v2[i,]=c(est$ow[2,2],est$ipw[2,2],est$ipwc1[2,2],est$ipwc2[2,2],est$ipwc3[2,2],est$ipwa1[2,2],est$ipwa2[2,2],est$ipwa3[2,2],
           est$sow[2,2],est$sipw[2,2],est$sipwc1[2,2],est$sipwc2[2,2],est$sipwc3[2,2],est$sipwa1[2,2],est$sipwa2[2,2],est$sipwa3[2,2])
  # 3rd time point
  p3[i,]=c(est$ow[3,1],est$ipw[3,1],est$ipwc1[3,1],est$ipwc2[3,1],est$ipwc3[3,1],est$ipwa1[3,1],est$ipwa2[3,1],est$ipwa3[3,1],
           est$sow[3,1],est$sipw[3,1],est$sipwc1[3,1],est$sipwc2[3,1],est$sipwc3[3,1],est$sipwa1[3,1],est$sipwa2[3,1],est$sipwa3[3,1])
  v3[i,]=c(est$ow[3,2],est$ipw[3,2],est$ipwc1[3,2],est$ipwc2[3,2],est$ipwc3[3,2],est$ipwa1[3,2],est$ipwa2[3,2],est$ipwa3[3,2],
           est$sow[3,2],est$sipw[3,2],est$sipwc1[3,2],est$sipwc2[3,2],est$sipwc3[3,2],est$sipwa1[3,2],est$sipwa2[3,2],est$sipwa3[3,2])
  # 4th time point
  p4[i,]=c(est$ow[4,1],est$ipw[4,1],est$ipwc1[4,1],est$ipwc2[4,1],est$ipwc3[4,1],est$ipwa1[4,1],est$ipwa2[4,1],est$ipwa3[4,1],
           est$sow[4,1],est$sipw[4,1],est$sipwc1[4,1],est$sipwc2[4,1],est$sipwc3[4,1],est$sipwa1[4,1],est$sipwa2[4,1],est$sipwa3[4,1])
  v4[i,]=c(est$ow[4,2],est$ipw[4,2],est$ipwc1[4,2],est$ipwc2[4,2],est$ipwc3[4,2],est$ipwa1[4,2],est$ipwa2[4,2],est$ipwa3[4,2],
           est$sow[4,2],est$sipw[4,2],est$sipwc1[4,2],est$sipwc2[4,2],est$sipwc3[4,2],est$sipwa1[4,2],est$sipwa2[4,2],est$sipwa3[4,2])
}
# time point 1
PerBias1  = apply((p1-t(replicate(replication, Delta.true[1,])))/t(replicate(replication, Delta.true[1,])),2,mean,na.rm=T)
RelEff1 = apply(p1,2,var,na.rm=T)[2]/apply(p1,2,var,na.rm=T)
twosd1=1.96*sqrt(v1)
CR1 = apply(((p1-twosd1)<t(replicate(replication, Delta.true[1,]))) & ((p1+twosd1)>t(replicate(replication, Delta.true[1,]))),2,mean,na.rm=T)

Tab1=as.data.frame(cbind(names(CR1),
                         round(Delta.true[1,],3),
                         round(PerBias1,3)*100,
                         round(RelEff1,3),
                         round(CR1,3)*100))
colnames(Tab1)=c("Method","WATE","PerBias","RelEff","CR")
rownames(Tab1)=names(CR1)
Tab1

# time point 2
PerBias2  = apply((p2-t(replicate(replication, Delta.true[2,])))/t(replicate(replication, Delta.true[2,])),2,mean,na.rm=T)
RelEff2 = apply(p2,2,var,na.rm=T)[2]/apply(p2,2,var,na.rm=T)
twosd2=1.96*sqrt(v2)
CR2 = apply(((p2-twosd2)<t(replicate(replication, Delta.true[2,]))) & ((p2+twosd2)>t(replicate(replication, Delta.true[2,]))),2,mean,na.rm=T)

Tab2=as.data.frame(cbind(names(CR2),round(Delta.true[2,],3),
                         round(PerBias2,3)*100,
                         round(RelEff2,3),
                         round(CR2,3)*100))
colnames(Tab2)=c("Method","WATE","PerBias","RelEff","CR")
rownames(Tab2)=names(CR2)
Tab2

# time point 3
PerBias3  = apply((p3-t(replicate(replication, Delta.true[3,])))/t(replicate(replication, Delta.true[3,])),2,mean,na.rm=T)
RelEff3 = apply(p3,2,var,na.rm=T)[2]/apply(p3,2,var,na.rm=T)
twosd3=1.96*sqrt(v3)
CR3 = apply(((p3-twosd3)<t(replicate(replication, Delta.true[3,]))) & ((p3+twosd3)>t(replicate(replication, Delta.true[3,]))),2,mean,na.rm=T)

Tab3=as.data.frame(cbind( names(CR3),round(Delta.true[3,],3),
                          round(PerBias3,3)*100,
                          round(RelEff3,3),
                          round(CR3,3)*100))
colnames(Tab3)=c("Method","WATE","PerBias","RelEff","CR")
rownames(Tab3)=names(CR3)
Tab3



# time point 4
PerBias4  = apply((p4-t(replicate(replication, Delta.true[4,])))/t(replicate(replication, Delta.true[4,])),2,mean,na.rm=T)
RelEff4 = apply(p4,2,var,na.rm=T)[2]/apply(p4,2,var,na.rm=T)
twosd4=1.96*sqrt(v4)
CR4 = apply(((p4-twosd4)<t(replicate(replication, Delta.true[4,]))) & ((p4+twosd4)>t(replicate(replication, Delta.true[4,]))),2,mean,na.rm=T)

Tab4=as.data.frame(cbind(names(CR4),round(Delta.true[4,],3),
                         round(PerBias4,3)*100,
                         round(RelEff4,3),
                         round(CR4,3)*100))
colnames(Tab4)=c("Method","WATE","PerBias","RelEff","CR")
rownames(Tab4)=names(CR4)
Tab4

library(rtf)
rtffile <- RTF("result_strong_50cen.doc")  # this can be an .rtf or a .doc
addParagraph(rtffile, paste("Results for time point 1 at t=", round(uvec[1],3),":\n",sep=""))
addTable(rtffile,Tab1)
addParagraph(rtffile, paste("\n\nResults for time point 2 at t=", round(uvec[2],3),":\n",sep=""))
addTable(rtffile,Tab2)
addParagraph(rtffile, paste("\n\nResults for time point 3 at t=", round(uvec[3],3),":\n",sep=""))
addTable(rtffile,Tab3)
addParagraph(rtffile, paste("\n\nResults for time point 4 at t=", round(uvec[4],3),":\n",sep=""))
addTable(rtffile,Tab4)
done(rtffile)

myTab=cbind(Tab1,Tab2,Tab3,Tab4)
write.csv(myTab,"result_strong_50cen.csv")



