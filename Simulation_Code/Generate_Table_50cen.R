rm(list=ls())
{
setwd("/Users/chaocheng/Dropbox/Mediation/Causal/simulation/HPC_Sim_210706")
dir()
library(rtf)

dat1=read.csv("result_strong_50cen.csv")[,-1]
dat2=read.csv("result_medium_50cen.csv")[,-1]
dat3=read.csv("result_weak_50cen.csv")[,-1]

###################################################
# bias table
###################################################
i1=c(1:8-1)*3+1
i2=c(1:8-1)*3+2
i3=c(1:8-1)*3+3
nameT1="PerBias";nameT2="PerBias.1";nameT3="PerBias.2";nameT4="PerBias.3"
dat_tab=data.frame(matrix(0,ncol=11,nrow=24))
colnames(dat_tab)=c("Method","Overlap","True","T1","T2","T3","T4","ST1","ST2","ST3","ST4")
dat_tab$Method=rep(c("OW","IPW","IPWC1","IPWC2","IPWC3","IPWA1","IPWA2","IPWA3"),each=3)
dat_tab$Overlap=rep(c("strong","medium","weak"),8)
# True
dat_tab$True[i1]=dat1$WATE[1:8]
dat_tab$True[i2]=dat2$WATE[1:8]
dat_tab$True[i3]=dat3$WATE[1:8]
# T1
dat_tab$T1[i1]=dat1[1:8,nameT1]
dat_tab$T1[i2]=dat2[1:8,nameT1]
dat_tab$T1[i3]=dat3[1:8,nameT1]
dat_tab$ST1[i1]=dat1[9:16,nameT1]
dat_tab$ST1[i2]=dat2[9:16,nameT1]
dat_tab$ST1[i3]=dat3[9:16,nameT1]
# T2
dat_tab$T2[i1]=dat1[1:8,nameT2]
dat_tab$T2[i2]=dat2[1:8,nameT2]
dat_tab$T2[i3]=dat3[1:8,nameT2]
dat_tab$ST2[i1]=dat1[9:16,nameT2]
dat_tab$ST2[i2]=dat2[9:16,nameT2]
dat_tab$ST2[i3]=dat3[9:16,nameT2]
# T3
dat_tab$T3[i1]=dat1[1:8,nameT3]
dat_tab$T3[i2]=dat2[1:8,nameT3]
dat_tab$T3[i3]=dat3[1:8,nameT3]
dat_tab$ST3[i1]=dat1[9:16,nameT3]
dat_tab$ST3[i2]=dat2[9:16,nameT3]
dat_tab$ST3[i3]=dat3[9:16,nameT3]
# T4
dat_tab$T4[i1]=dat1[1:8,nameT4]
dat_tab$T4[i2]=dat2[1:8,nameT4]
dat_tab$T4[i3]=dat3[1:8,nameT4]
dat_tab$ST4[i1]=dat1[9:16,nameT4]
dat_tab$ST4[i2]=dat2[9:16,nameT4]
dat_tab$ST4[i3]=dat3[9:16,nameT4]

dat_tab_bias = dat_tab

###################################################
# efficiency table
###################################################

i1=c(1:8-1)*3+1
i2=c(1:8-1)*3+2
i3=c(1:8-1)*3+3
nameT1="RelEff";nameT2="RelEff.1";nameT3="RelEff.2";nameT4="RelEff.3"
dat_tab=data.frame(matrix(0,ncol=11,nrow=24))
colnames(dat_tab)=c("Method","Overlap","True","T1","T2","T3","T4","ST1","ST2","ST3","ST4")
dat_tab$Method=rep(c("OW","IPW","IPWC1","IPWC2","IPWC3","IPWA1","IPWA2","IPWA3"),each=3)
dat_tab$Overlap=rep(c("strong","medium","weak"),8)
# True
dat_tab$True[i1]=dat1$WATE[1:8]
dat_tab$True[i2]=dat2$WATE[1:8]
dat_tab$True[i3]=dat3$WATE[1:8]
# T1
dat_tab$T1[i1]=dat1[1:8,nameT1]
dat_tab$T1[i2]=dat2[1:8,nameT1]
dat_tab$T1[i3]=dat3[1:8,nameT1]
dat_tab$ST1[i1]=dat1[9:16,nameT1]
dat_tab$ST1[i2]=dat2[9:16,nameT1]
dat_tab$ST1[i3]=dat3[9:16,nameT1]
# T2
dat_tab$T2[i1]=dat1[1:8,nameT2]
dat_tab$T2[i2]=dat2[1:8,nameT2]
dat_tab$T2[i3]=dat3[1:8,nameT2]
dat_tab$ST2[i1]=dat1[9:16,nameT2]
dat_tab$ST2[i2]=dat2[9:16,nameT2]
dat_tab$ST2[i3]=dat3[9:16,nameT2]
# T3
dat_tab$T3[i1]=dat1[1:8,nameT3]
dat_tab$T3[i2]=dat2[1:8,nameT3]
dat_tab$T3[i3]=dat3[1:8,nameT3]
dat_tab$ST3[i1]=dat1[9:16,nameT3]
dat_tab$ST3[i2]=dat2[9:16,nameT3]
dat_tab$ST3[i3]=dat3[9:16,nameT3]
# T4
dat_tab$T4[i1]=dat1[1:8,nameT4]
dat_tab$T4[i2]=dat2[1:8,nameT4]
dat_tab$T4[i3]=dat3[1:8,nameT4]
dat_tab$ST4[i1]=dat1[9:16,nameT4]
dat_tab$ST4[i2]=dat2[9:16,nameT4]
dat_tab$ST4[i3]=dat3[9:16,nameT4]

dat_tab[,4:11]=round(dat_tab[,4:11],2)

dat_tab_releff = dat_tab

###################################################
# 95% CI coverage rate table
###################################################

i1=c(1:8-1)*3+1
i2=c(1:8-1)*3+2
i3=c(1:8-1)*3+3
nameT1="CR";nameT2="CR.1";nameT3="CR.2";nameT4="CR.3"
dat_tab=data.frame(matrix(0,ncol=11,nrow=24))
colnames(dat_tab)=c("Method","Overlap","True","T1","T2","T3","T4","ST1","ST2","ST3","ST4")
dat_tab$Method=rep(c("OW","IPW","IPWC1","IPWC2","IPWC3","IPWA1","IPWA2","IPWA3"),each=3)
dat_tab$Overlap=rep(c("strong","medium","weak"),8)
# True
dat_tab$True[i1]=dat1$WATE[1:8]
dat_tab$True[i2]=dat2$WATE[1:8]
dat_tab$True[i3]=dat3$WATE[1:8]
# T1
dat_tab$T1[i1]=dat1[1:8,nameT1]
dat_tab$T1[i2]=dat2[1:8,nameT1]
dat_tab$T1[i3]=dat3[1:8,nameT1]
dat_tab$ST1[i1]=dat1[9:16,nameT1]
dat_tab$ST1[i2]=dat2[9:16,nameT1]
dat_tab$ST1[i3]=dat3[9:16,nameT1]
# T2
dat_tab$T2[i1]=dat1[1:8,nameT2]
dat_tab$T2[i2]=dat2[1:8,nameT2]
dat_tab$T2[i3]=dat3[1:8,nameT2]
dat_tab$ST2[i1]=dat1[9:16,nameT2]
dat_tab$ST2[i2]=dat2[9:16,nameT2]
dat_tab$ST2[i3]=dat3[9:16,nameT2]
# T3
dat_tab$T3[i1]=dat1[1:8,nameT3]
dat_tab$T3[i2]=dat2[1:8,nameT3]
dat_tab$T3[i3]=dat3[1:8,nameT3]
dat_tab$ST3[i1]=dat1[9:16,nameT3]
dat_tab$ST3[i2]=dat2[9:16,nameT3]
dat_tab$ST3[i3]=dat3[9:16,nameT3]
# T4
dat_tab$T4[i1]=dat1[1:8,nameT4]
dat_tab$T4[i2]=dat2[1:8,nameT4]
dat_tab$T4[i3]=dat3[1:8,nameT4]
dat_tab$ST4[i1]=dat1[9:16,nameT4]
dat_tab$ST4[i2]=dat2[9:16,nameT4]
dat_tab$ST4[i3]=dat3[9:16,nameT4]

dat_tab_cr = dat_tab



rtffile <- RTF("tab_50cen.doc")  # this can be an .rtf or a .doc
addParagraph(rtffile, "Percent Bias")
addTable(rtffile,dat_tab_bias)
addParagraph(rtffile, "Relative Efficiency")
addTable(rtffile,dat_tab_releff)
addParagraph(rtffile, "95% CI Coverage Rate")
addTable(rtffile,dat_tab_cr)
done(rtffile)
}
