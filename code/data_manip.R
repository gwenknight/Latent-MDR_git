#### Data manipulation
# Mortality - make 1 year ranges

library(reshape2)
library(ggplot2)

#### Where is the data?
data_home <- "~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/Data/"
setwd(data_home)

#### mortality rates
mort_orig <- read.csv("mort.csv")
# rep rows for each time period
mort <- mort_orig[rep(seq_len(nrow(mort_orig)), each=5),]
mort$year <- rep(seq(1950, 2014, 1),4)

# rep columns for each age
mort <- mort[,rep(seq_len(ncol(mort)), each=5)]
mort <- mort[,-c(2,3,4,5,7,8,9,10)]
mort <- mort[,-c(107,106,105,104)]
colnames(mort) <- c("country","time",seq(0,99,1),"year")

# melt 
m_mort <- melt(mort, id.vars = c("country","time","year"))
colnames(m_mort) <-c("country","time","year","age","value")

# add back to 1800s
w<-which(m_mort$year == 1950)
m_1950 <- m_mort[w,]
for(i in 1949:1800){
  m_1950$year <- i
  m_mort <- rbind(m_1950,m_mort)
}

# save
write.csv(m_mort, "m_mort.csv")


# plot
ggplot(m_mort,aes(x=year,y = value, group = age)) + geom_line(aes(colour=age)) + facet_wrap(~country)

w<-which(as.numeric(m_mort$age) > 80)
ggplot(m_mort[w,],aes(x=year,y = value, group = age)) + geom_line(aes(colour=age)) + facet_wrap(~country)
