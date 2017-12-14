#### Data manipulation
# Mortality - make 1 year ranges

library(reshape2)
library(ggplot2)
options(stringsAsFactors = FALSE)

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

# interpolate only age 5-80
m_mort$in_value <- m_mort$value
countries <- c("South Africa","India","China","United Kingdom")
for(cc in 1:length(countries)){
  c <- which(m_mort$country == countries[cc])
  w<-intersect(c,intersect(which(as.numeric(m_mort$age) > 6), which(as.numeric(m_mort$age) < 100))) # when as.numeric goes to level number
  for(j in 1800:2014){
    w1<-intersect(w, which(m_mort$year==j))
    aa<-approx(m_mort[w1,"value"][seq(1,dim(m_mort)[1],5)],n = length(w1)) # jumps every 5 yrs: interpolate between
    m_mort[w1,"in_value"] <- aa$y
  }
}

plot((as.numeric(m_mort[c[1:100],"age"])-1),m_mort[c[1:100],"value"])
points((as.numeric(m_mort[c[1:100],"age"])-1),m_mort[c[1:100],"in_value"], col="red")


plot((as.numeric(m_mort[c[20000:21500],"age"])-1),m_mort[c[20000:21500],"value"])
points((as.numeric(m_mort[c[20000:21500],"age"])-1),m_mort[c[20000:21500],"in_value"], col="red")




# save
write.csv(m_mort, "m_mort.csv")


# plot
ggplot(m_mort,aes(x=year,y = value, group = age)) + geom_line(aes(colour=age)) + facet_wrap(~country)

w<-which(as.numeric(m_mort$age) > 80)
ggplot(m_mort[w,],aes(x=year,y = value, group = age)) + geom_line(aes(colour=age)) + facet_wrap(~country)


###**** Births ***###
birth_orig <- read.csv("birth.csv")

# rep rows for each time period
birth <- birth_orig[rep(seq_len(nrow(birth_orig)), each=5),]
birth$year <- rep(seq(2019, 1950, -1),4)

# add back to 1800s
w<-which(birth$year == 1950)
m_1950 <- birth[w,]
for(i in 1949:1800){
  m_1950$year <- i
  birth <- rbind(birth,m_1950)
}

# save
write.csv(birth, "m_birth.csv")

# plot
ggplot(birth,aes(x=year,y = births)) + geom_line(aes(colour=Country)) 
