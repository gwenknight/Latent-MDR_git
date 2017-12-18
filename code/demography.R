### Demography checks

###*** background to run
# Libraries
library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(reshape2)
theme_set(theme_bw())
library(plyr)

# home 
home <- "~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/"
code <- "~/Documents/Latent-MDR_git/code/"
plots <- "~/Documents/Latent-MDR_git/plots"

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
yearend <- 2014
steps <- 1 + (yearend-year1) / dt

###*** Run with no TB
initial[2:7,] <- 0 
X <- nat_hist(para_v, para_s, mort, birth, c(year1, yearend, dt), initial)

# Normalise by population size 
mpsz_age <- as.data.frame(matrix(0,dim(X$psize_age)[1],dim(X$psize_age)[2]))
for(i in 1:dim(X$psize_age)[1]){mpsz_age[i,] <- 100 * X$psize_age[i,] / X$psize[i] }
# Label columns by age and add in times, melt
colnames(mpsz_age)<-seq(1:Mnage);
mpsz_age$time <- seq(1:steps)
mpsz_age$year <- seq(year1, yearend, dt)
m_psz_age <- melt(mpsz_age,id.vars = c("year","time"))

setwd(plots)
# prior to 1950 (data starts)
times <- round(seq(1800,1950,length.out = 9),0)
w<-c()
for(i in 1:length(times)){w<-c(w,which(m_psz_age$year == times[i]))}
ggplot(m_psz_age[w,], aes(x = variable, y = value )) + geom_bar(stat = "identity") + coord_flip() + facet_wrap(~year) + scale_y_continuous("Percentage")
ggsave(paste(country,"_prior1950.pdf",sep=""))

# near to 1950 - stable? yes
times <- seq(1935,1950,dt)
w<-c()
for(i in 1:length(times)){w<-c(w,which(m_psz_age$year == times[i]))}
ggplot(m_psz_age[w,], aes(x = variable, y = value )) + geom_bar(stat = "identity") + coord_flip() + facet_wrap(~year)
ggsave(paste(country,"_prior1935_1950.pdf",sep=""))

# 1950 onwards?
times <- round(seq(1950,2014,length.out = 20),0)
w<-c()
for(i in 1:length(times)){w<-c(w,which(m_psz_age$year == times[i]))}
g <- ggplot(m_psz_age[w,], aes(x = variable, y = value )) + geom_bar(stat = "identity") + coord_flip() + 
   facet_wrap(~year) + scale_x_discrete(breaks = seq(0,100,5))
g
ggsave(paste(country,"_1950onwards.pdf",sep=""))
