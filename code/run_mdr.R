### RUN_MDR model

# Libraries
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(reshape2)
theme_set(theme_bw())
library(plyr)

# home 
home <- "~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/"
code <- "~/Documents/Latent-MDR_git/code/"

# get functions
setwd(code)
source("Nat_Hist_fn.R")
# time step determines treatment matrix - if change have to change parameters.R
dt <- 0.5
source("parameters.R")

# Variable parameters para_v
para_v <-         c(202,    0.3,  0.33)
names(para_v) <- c("beta", "f", "x")

year1 <- 1800
yearend <- 2015
steps <- 1 + (yearend-year1) / dt

X <- nat_hist(para_v, para_s, mort, birth, c(year1, yearend, dt), initial)

X$prev[431]
X$inc[431] ## India = 211

##*** ANALYSIS *** ###
# summary
summ_row <- as.data.frame(rbind(rowSums(X$U),rowSums(X$LS),rowSums(X$LR),rowSums(X$AS),rowSums(X$AR),rowSums(X$TR),rowSums(X$TS),
                                rowSums(X$LS_p),rowSums(X$LR_p),rowSums(X$AS_p),rowSums(X$AR_p),rowSums(X$TR_p),rowSums(X$TS_p)))
colnames(summ_row)<-seq(1:steps);
summ_row$pop<-c("U","LS","LR","AS","AR","TR","TS","LS_p","LR_p","AS_p","AR_p","TR_p","TS_p")
m_summ_row <- melt(summ_row[-1,],id.var = "pop")

ggplot(m_summ_row,aes(x=variable, y = value, colour = pop, group = pop)) + geom_line(size = 1.2) + 
  scale_color_manual(values= palette(rainbow(12))) # + scale_y_continuous(limit = c(0,1000))

ggplot(m_summ_row,aes(x = variable, y = value, colour = pop, group = pop)) + geom_line() + facet_wrap(~ pop, scales = "free")

par(mfrow=c(3,2))
# psize change
plot(seq(1:steps),X$psize, type = "l")

# incidence change
if(max(X$inc[2,]) > max(X$inc[1,])){
plot (seq(1:steps),X$inc[2,],type = "l", col ="red")
lines(seq(1:steps),X$inc[1,], col = "blue")
} else {plot (seq(1:steps),X$inc[1,],type = "l", col ="red")
  lines(seq(1:steps),X$inc[2,], col = "blue")}

# prev change
if(max(X$prev[2,]) > max(X$prev[1,])){
  plot (seq(1:steps),X$prev[2,],type = "l", col ="red")
  lines(seq(1:steps),X$prev[1,], col = "blue")
} else {plot (seq(1:steps),X$prev[1,],type = "l", col ="red")
  lines(seq(1:steps),X$prev[2,], col = "blue")}

# percentage of new TB cases that are MDR (should be ~ 20%)
plot (seq(1:steps),X$ratio_mdr,type = "l", ylab = "Perc. new TB = MDR")

# proportion of population latently infected - varies by country? 
w<-c(which(summ_row$pop == "LS"), which(summ_row$pop == "LR"), which(summ_row$pop == "LS_p"), which(summ_row$pop == "LR_p"))
plot( seq(1:steps), 100*colSums(summ_row[w,1:steps])/X$psize, type="l", ylab = "Perc. latently infected")


# proportion change
p_summ_row <- summ_row
for(ii in 1:length(X$psize)){p_summ_row[,ii]<-p_summ_row[,ii]/X$psize[ii]}
m_p_summ_row<-melt(p_summ_row,id.vars = "pop")
ggplot(m_p_summ_row,aes(x = variable, y = value,fill=pop)) + geom_bar(stat="identity")

# age distribution
mAS <- as.data.frame(X$AS)
colnames(mAS)<-seq(1:Mnage);
mAS$time <- seq(1:steps)
mAS$population <- "AS"
m_mAS <- melt(mAS,id.vars = c("time","population"))

mAR <- as.data.frame(X$AR)
colnames(mAR)<-seq(1:Mnage);
mAR$time <- seq(1:steps)
mAR$population <- "AR"
m_mAR <- melt(mAR,id.vars = c("time","population"))
m_mAR$value <- -m_mAR$value # flip for axis

m_pyr <- rbind(m_mAR,m_mAS)
m_pyr1<-m_pyr[which(m_pyr$time == 5),]

ggplot(m_pyr1, aes(x = variable, y = value, fill = population)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  scale_y_continuous("Population size",breaks = seq(-1,1,0.1)) + scale_x_discrete("Age",breaks = seq(0,100,10))

m_pyr2<-m_pyr[which(m_pyr$time == steps),]

ggplot(m_pyr2, aes(x = variable, y = value, fill = population)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  scale_y_continuous("Population size",breaks = seq(-1,2,0.1)) + scale_x_discrete("Age",breaks = seq(0,100,10))

mpsz_age <- as.data.frame(X$psize_age)
colnames(mpsz_age)<-seq(1:Mnage);
mpsz_age$time <- seq(1:steps)
m_psz_age <- melt(mpsz_age,id.vars = c("time"))

w<-c(which(m_psz_age$time == 2),which(m_psz_age$time == round(steps / 3,0)),which(m_psz_age$time == round(2*steps/3,0)),which(m_psz_age$time == steps))
ggplot(m_psz_age[w,], aes(x = variable, y = value )) + geom_bar(stat = "identity") + coord_flip() + facet_wrap(~time)

