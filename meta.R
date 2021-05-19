###########################################################################################
# @Title: The association between SYCP3 c.657T>C polymorphism and recurrent pregnancy loss
# @Author: Jianhua Chen 
# @Date: 2021-05-19   
###########################################################################################

### WARNING：THIS IS A WORKING CODE ###

### set working path ###
setwd("G:\\meta")
rm(list = ls())
data = "data extraction.csv"

### 1 load package ###
  #install.packages("meta")
  library(meta)
  #install.packages("genetics")
  library(genetics)

### 2 data preparation ###
  rt=read.csv(data,sep=",",header = T)

### 3 caculate OR value & forest plot ###
  
  # 3.1 Allele model: C vs. T
  metabin(case.C,case.T+case.C,control.C,control.T+control.C,data=rt[-2,],sm="OR")
  alle <- metabin(case.C,case.T+case.C,control.C,control.T+control.C,data=rt[-2,],sm="OR",
                  studlab=paste(first.author,year),comb.random=FALSE,
                  label.e="SYCP3 T657C",label.c = "Control")
  pdf("forest_Allele.pdf",width = 15,height = 10)
  forest(alle) 
  dev.off()
  
  # 3.2 Dominant model: CC+CT vs. TT
  metabin(case.CC+case.CT,case,control.CC+control.CT,control,data=rt[-2,],sm="OR")
  domi <- metabin(case.CC+case.CT,case,control.CC+control.CT,control,data=rt[-2,],sm="OR",
                  studlab=paste(first.author,year),comb.random=FALSE,
                  label.e="SYCP3 T657C",label.c = "Control")
  pdf("forest_Dominant.pdf",width = 15,height = 10)
  forest(domi)
  dev.off()
  
  # 3.3 Recessive model: CC vs. TT+CT
  metabin(case.CC,case,control.CC,control,data=rt[3,],sm="OR")
  rece <- metabin(case.CC,case,control.CC,control,data=rt[3,],sm="OR",
                  studlab=paste(first.author,year),comb.random=FALSE,
                  label.e="SYCP3 T657C",label.c = "Control")
  pdf("forest_recessive.pdf",width = 15,height = 10)
  forest(rece)
  dev.off()
  
  # 3.4 co-dominant model: CT vs. TT
  metabin(case.CT,case.CT+case.TT,control.CT,control.CT+control.TT,data=rt[-2,],sm="OR")
  codo <- metabin(case.CT,case.CT+case.TT,control.CT,control.CT+control.TT,data=rt[-2,],sm="OR",
                  studlab=paste(first.author,year),comb.random=FALSE,
                  label.e="SYCP3 T657C",label.c = "Control")
  pdf("forest_Co-Dominant.pdf",width = 15,height = 10)
  forest(codo)
  dev.off()
  
  # 3.5 homozygote model: CC vs. TT
  metabin(case.CC,case.CC+case.TT,control.CC,control.CC+control.TT,data=rt[3],sm="OR")
  homo <- metabin(case.CC,case.CC+case.TT,control.CC,control.CC+control.TT,data=rt[3,],sm="OR",
                  studlab=paste(first.author,year),comb.random=FALSE,
                  label.e="SYCP3 T657C",label.c = "Control")
  pdf("forest_homozygote.pdf",width = 15,height = 10)
  forest(homo)
  dev.off()
  
### 5 publish bias ###
  
  # as k=4 is too small, choose trim and filled to correct
  metabias(alle,method.bias = "Begg",k.min=5)
  metabias(alle,method.bias = "Egger",k.min=5)
  tf1_alle <- trimfill(alle,comb.fixed=TRUE)
  summary(tf1_alle)
  pdf("funnel_Allele.pdf")
  funnel(tf1_alle)
  dev.off()
  
  metabias(domi,method.bias = "Begg",k.min=5)
  metabias(domi,method.bias = "Egger",k.min=5)
  tf1_domi <- trimfill(domi,comb.fixed=TRUE)
  summary(tf1_domi)
  pdf("funnel_dominant.pdf")
  funnel(tf1_domi)
  dev.off()
  
  metabias(codo,method.bias = "Begg",k.min=5)
  metabias(codo,method.bias = "Egger",k.min=5)
  tf1_codo <- trimfill(codo,comb.fixed=TRUE)
  summary(tf1_codo)
  pdf("funnel_Co-Dominant.pdf")
  funnel(tf1_codo)
  dev.off()
  
  tf1_homo <- trimfill(homo,comb.fixed=TRUE)
  summary(tf1_homo)
  pdf("funnel_homozygote.pdf")
  funnel(tf1_homo)
  dev.off()
  
  tf1_rece <- trimfill(rece,comb.fixed=TRUE)
  summary(tf1_rece)
  pdf("funnel_recessive.pdf")
  funnel(tf1_rece)
  dev.off()
  
### 6 sensitivity  analysis ###
  metainf(alle, pooled="fixed")
  forest(metainf(alle), comb.fixed=TRUE)
  
  metainf(codo, pooled="fixed")
  forest(metainf(codo), comb.fixed=TRUE)
  
  metainf(rece, pooled="fixed")
  forest(metainf(rece), comb.fixed=TRUE)
  
  metainf(domi, pooled="fixed")
  forest(metainf(domi), comb.fixed=TRUE)
  
  metainf(homo, pooled="fixed")
  forest(metainf(homo), comb.fixed=TRUE)
  
 ### 7 HWE test in cotrol group ###
  #install.packages("genetics")
  library(genetics)
  
  k <- 1
  hwe1 <- genotype(c(rep("T/T",rt[k,9]),rep("C/T",rt[k,11]),rep("C/C",rt[k,13])))
  HWE.chisq(hwe1)
  
  k <- 3
  hwe3 <- genotype(c(rep("T/T",rt[k,9]),rep("C/T",rt[k,11]),rep("C/C",rt[k,13])))
  HWE.chisq(hwe3)
  
  k <- 4
  hwe4 <- genotype(c(rep("T/T",rt[k,9]),rep("C/T",rt[k,11]),rep("C/C",rt[k,13])))
  HWE.chisq(hwe4)
  
  k <- 5
  hwe5 <- genotype(c(rep("T/T",rt[k,9]),rep("C/T",rt[k,11]),rep("C/C",rt[k,13])))
  HWE.chisq(hwe5)
  
  k <- 6
  hwe6 <- genotype(c(rep("T/T",rt[k,9]),rep("C/T",rt[k,11]),rep("C/C",rt[k,13])))
  HWE.chisq(hwe6)
  
  
  
  
  
  
  
  
  
  
  