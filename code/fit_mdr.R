#### Attempt to fit

library("lhs")


# home 
home <- "~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/"
code <- "~/Documents/Latent-MDR_git/code/"

# get functions
setwd(code)
source("Nat_Hist_fn.R")
# time step determines treatment matrix - if change have to change parameters.R
dt <- 0.5 # half a year
source("parameters.R")


### select parameters
runs <- 100000 # 3 sec for 10, 
npara <- 9
Xp <- randomLHS(runs, npara)
Y <- matrix(0, nrow=runs, ncol=npara)

para <- as.data.frame(matrix(0,npara,3))
colnames(para) <- c("names","min","max")
para$names <- c("wr", "ws", "eps", "ma","sigma","p","beta", "f", "x")
para$min <- c(1/3,  1, 0.04, 0.1, 0.0005, 0.05, 10, 0, 0.1)  
para$max <- c(1, 4, 0.08, 0.3, 0.002, 0.4, 200, 1, 0.7)

for(i in 1:npara){
  Y[,i] <- qunif(Xp[,i], min=para[i,"min"], max=para[i,"max"])  
}

# SET
year1 <- 1800
yearend <- 2014
steps <- 1 + (yearend-year1) / dt

# Country
country <- "India"

# RUN
setwd("~/Dropbox/MRC SD Fellowship/Research/MDR/Latent MDR/fits/")
err_i <- c()
Ystore<-c()

for(i in 1:runs){
  # Empty store places
  store_para <- c() # STORE
  store_para1 <- c() # STORE
  store_para2 <- c() # STORE
  store_para3 <- c() # STORE
  store_X <- c() # STORE
  
  para_s <- Y[i,1:6]
  para_v <- Y[i,7:9]
  names(para_s) <- c("wr", "ws", "eps", "ma","sigma","p")
  names(para_v) <- c("beta", "f", "x")
  tryCatch({
    X <- nat_hist(para_v, para_s, mort, birth, c(year1, yearend, dt), initial)
    Ystore <- rbind(Ystore,Y[i,]) # only store if nat_hist works
  if(X$prev[1,steps] > 20 & X$prev[1,steps] < 1200){ ## AS Between 20-1200 (Kendall 1st)
    store_para1 <- rbind(store_para,c(para_s, para_v))
    if(X$inc[1,steps] > 100 & X$inc[1,steps] < 500){ ## AS Between 20-1400 per year (Kendall 1st). India = 211 (109–345)
      store_para2 <- rbind(store_para,c(para_s, para_v))
      if(X$ratio_mdr[1,steps] > 2 & X$ratio_mdr[1,steps] < 3.5){
        store_para3 <- rbind(store_para,c(para_s, para_v))
        if(X$ratio_mdr[2,steps] > 10 & X$ratio_mdr[2,steps] < 13){# new: 5%, prev.treat: 25% India: 2.8% (2–3.5) 12% (10–13)
          
          store_para <- rbind(store_para,c(para_s, para_v))
          
        }
      }  
    }
  }
  
  
  # every run
  store_X <- rbind(store_X, c(X$prev[1,steps],X$inc[1,steps],X$ratio_mdr[1,steps]))
  }, error=function(e,err_i){cat("ERROR :",conditionMessage(e), " ",i," \n")})
  
  
  
  # every 1,000 runs write out 
  if(i %% 1000 == 0){ 
    
    # read in
    astore_para <- read.csv("store_para.csv")[,-1]
    astore_para1 <- read.csv("store_para1.csv")[,-1]
    astore_para2 <- read.csv("store_para2.csv")[,-1]
    astore_para3 <- read.csv("store_para3.csv")[,-1]
    astore_X <- read.csv("store_x.csv")[,-1]
    
    write.csv(rbind(astore_para,store_para), "store_para.csv")
    write.csv(rbind(astore_para1,store_para1), "store_para1.csv")
    write.csv(rbind(astore_para2,store_para2), "store_para2.csv")
    write.csv(rbind(astore_para3,store_para3), "store_para3.csv")
    write.csv(rbind(astore_X,store_X), "store_x.csv")
    
  }
  
}

write.csv(Y, "store_all_para.csv") # store parameters

if(runs < 1000){
  write.csv(store_para, "store_para.csv")
  write.csv(store_para1, "store_para1.csv")
  write.csv(store_para2, "store_para2.csv")
  write.csv(store_para3, "store_para3.csv")
  write.csv(store_X, "store_x.csv")
}
