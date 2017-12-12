### Parameters
# ?? make births / deaths non constant?
# ?? yearly changes in treatment success

# Highest age group
Mnage <- 100
upp = Mnage - 1

# Standard parameters para_s
para_s        <- c(0.0005,        1/2,   1,  0.08, 0.187, 0.0015, 0.14)
names(para_s) <- c("birth_rate",  "wr", "ws", "eps", "ma","sigma","p")
# assume 12months to S detection, 2 years to MDR detection

# proportion infectious by age
p_i <- matrix(1,1,Mnage) 
p_i[1:15] <- 0.02; p_i[16:Mnage] <- 0.15;

# mortality rates
#m <- matrix(0.0125,1,Mnage) # LE of 80yrs
m[1:10] <- 400/100000; # 40/1000 high from USAID
m[10:40]<- 30/100000; # from US https://www.google.co.uk/imgres?imgurl=https://www.cdc.gov/nchs/images/databriefs/1-50/db26_Fig_2.png&imgrefurl=https://www.cdc.gov/nchs/products/databriefs/db26.htm&h=686&w=960&tbnid=0ZQ2wBdlGuXr1M:&tbnh=150&tbnw=211&usg=___lRLJd8LVdpoKcbwF_EMZRQe9Bo%3D&vet=10ahUKEwih2a7CtM_XAhUrCsAKHeQuBZ8Q9QEILDAA..i&docid=0AfpBMnU6MIZFM&client=firefox-b-ab&sa=X&ved=0ahUKEwih2a7CtM_XAhUrCsAKHeQuBZ8Q9QEILDAA 
m[40:Mnage] <- 450 / 100000
m[Mnage] <- 1;
plot(seq(1,99,1),m[,1:(Mnage-1)])
1/mean(m) # want to be 80 - currently 76

###**** treatment parameters
# Length of treatment
rx_s_length <- 0.5 / dt
rx_r_length <- 1.5 / dt

# fail
fails <- matrix(0.01,rx_s_length,upp)
failr <- matrix(0.1,rx_r_length,upp)
# cure
cures <- matrix(0.01,rx_s_length,upp)
curer <- matrix(0,rx_r_length,upp)
# mortality 
mts <- matrix(0.01,rx_r_length,upp)
mtr <- matrix(0.01,rx_r_length,upp)

###*** Initial conditions
initial <- matrix(0,7,Mnage)
# U[1] <- init[1,]; LS[1] <- init[2,]; LS[1] <- init[3,]
# AS[1] <- init[4,]; AR[1] <- init[5,]; 
# TS[1] <- init[6,]; TR[1] <- init[7,]; 
pop_prop <- c(rep(0.02,16), rep(0.01,54), rep(0.14/30,30))
#plot(pop_prop)
initial[1,] <- (100000 - 5) * pop_prop # population of 100,000
initial[4,18:25] <- 5 # 5 18yos with TB 



