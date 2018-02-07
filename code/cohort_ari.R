## COHORT ARI MODEL
# use structure in tbaricon

### Libraries
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(reshape2)
theme_set(theme_bw())
library(plyr)

### home 
home <- "~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/"
code <- "~/Documents/Latent-MDR_git/code/"

##************************************************************************************###

ARIs <- c(matrix(0.01, 1, (1950-1800)), seq(0.01,0.001,length.out = (2015-1950)))
ARIr <- c(matrix(1, 1, (1950-1800)), seq(0.001,0.01,length.out = (2015-1950)))
steps <- 2015-1800
ages <- 101

birthrate <- as.data.frame(t(rbind(seq(1800,2014,1), matrix(0.01,1,steps))))
colnames(birthrate) <- c("years","value")
deathrate <- as.data.frame(t(rbind(seq(1800,2014,1), matrix(0.01,1,steps))))
colnames(deathrate) <- c("years","value")

dis_rec <- c(matrix(4.06,1,10),seq(4.06,13.8,length.out=11),matrix(13.80,1,80))
dis_react <- c(matrix(0,1,11),seq(0,0.03,length.out=9),matrix(0.03,1,81))
dis_reinf <- c(matrix(6.89,1,10),seq(6.89,8.25,length.out=11),matrix(8.25,1,80))

prop_inf <- c(matrix(0.1,1,10),seq(0.1,0.65,length.out=10),seq(0.653125,0.9,length.out=81))

Ds <- as.data.frame(matrix(0,ages*steps, 11))
colnames(Ds) <- c("year","age","psize","LTBI","new_inf","old_inf","reinf","dis_rec","dis_reac","dis_reinf","dis")
Dr <- as.data.frame(matrix(0,ages*steps, 11))
colnames(Dr) <- c("year","age","psize","LTBI","new_inf","old_inf","reinf",
                  "dis_rec","dis_reac","dis_reinf","dis")

Ds$age <- seq(1,ages,1); Dr$age <- seq(1,ages,1)

## Initial
Ds[1:(ages),c("psize")] <- 1000 # all 1000 to start

## Run
for(i in 2:steps){
  # what year?
  year <- seq(1800,2015)[i]
  # what ARI this year? 
  aris <- ARIs[i]
  arir <- ARIr[i]
  # births
  births = Ds[(a-ages):(b-ages),"psize"] * birthrate[which(birthrate$years == year),"value"]
  # deaths
  deaths = Ds[(a-ages):(b-ages),"psize"] * deathrate[which(birthrate$years == year),"value"]
  # fill in 
  a <- ((i-1)*(ages)+1); b <- (i*ages)
  Ds[a:b,"year"] = year; Dr[a:b,"year"] = year
  # new infections
  Ds[a:b,c("new_inf","old_inf","reinf")] = 
    c(aris*(1 - Ds[(a-ages):(b-ages),"LTBI"])*Ds[(a-ages):(b-ages),"psize"],
      (Ds[(a-ages):(b-ages),"LTBI"])*Ds[(a-ages):(b-ages),"psize"],
      aris*Ds[(a-ages):(b-ages),"LTBI"]*Ds[(a-ages):(b-ages),"psize"])
  # new infections
  Ds[a:b,c("new_inf","old_inf","reinf")] = 
    c(arir*(1 - Dr[(a-ages):(b-ages),"LTBI"])*Ds[(a-ages):(b-ages),"psize"],
      (Dr[(a-ages):(b-ages),"LTBI"])*Dr[(a-ages):(b-ages),"psize"],
      arir*Dr[(a-ages):(b-ages),"LTBI"]*Ds[(a-ages):(b-ages),"psize"])
  # new disease
  Ds[a:b,c("dis_rec","dis_reac","dis_reinf","dis")] = 
    c(Ds[a:b,"new_inf"]*dis_rec*prop_inf,
      Ds[a:b,"old_inf"]*dis_react*prop_inf,
      Ds[a:b,"reinf"]*dis_reinf*prop_inf,
      rowSums(Ds[a:b,c("dis_rec","dis_reac","dis_reinf")]))
  # LTBI - infected whether disease or not
  Ds[a:b,"LTBI"] = Ds[(a-ages):(b-ages),"LTBI"] + Ds[a:b,"new_inf"] 
  # psize 
  Ds[a:b,"psize"] = Ds[(a-ages):(b-ages),"psize"] - births - deaths
}




