
##########
ep = epsilon
para <- c(beta, birth_rate, pchild, v, x, u[i], wr,ws, mts, mtr, rxs, rxr)*dt
# proportion infectious by age
p_i <- matrix(1,1,Mnage) 
p_i[1:15] <- 0.1; p_i[16:Mnage] <- 0.5;

mts
rxs
# change by age (?) and time since start
failr #fail prop MDR

