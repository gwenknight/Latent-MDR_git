### Check treatment matrix works

# need
steps = 10
rx_s_length = 1
rx_r_length = 3
Mnage = 10
upp = Mnage - 1
TR <- matrix(0,steps,Mnage); TS <- matrix(0,steps,Mnage);
TS_m <- matrix(0,rx_s_length,Mnage); TS_m_new <- matrix(0,rx_s_length,Mnage)
TR_m <- matrix(0,rx_r_length,Mnage); TR_m_new <- matrix(0,rx_r_length,Mnage)

# para
fails <- matrix(0.1,rx_s_length,Mnage)
failr <- matrix(0.1,rx_r_length,Mnage)
cures <- matrix(0.2,rx_s_length,Mnage)
curer <- matrix(0.2,rx_r_length,Mnage)
mts <- matrix(0.01,rx_r_length,Mnage)
mtr <- matrix(0.01,rx_r_length,Mnage)
eps <- 0.01 / rx_r_length * (1/12)

save_TS_M <- c(); save_TR_M <- c()
for(i in 1:steps){
  
  TR[i,2:Mnage] = colSums(TR_m[,1:upp])
  TS[i,2:Mnage] = colSums(TR_m[,1:upp])
  
  # leave treatment now are:
  new_AS_from_rx <- colSums(fails*TS_m) # fail treatment 
  new_AR_from_rx <- colSums(failr*TR_m) + eps * TS_m[rx_s_length,] # fail treatment + acquisitions
  new_LS_from_rx <- colSums(cures*TS_m) # cured
  new_LR_from_rx <- colSums(curer*TR_m) # cured
  
  # new
  TS_m_new[1,] = 10
  TR_m_new[1,] = 10
  # update - here bottom row (rx_r/s_length) lost - to failr / cure / death / acquisition(s-r)
  if(rx_s_length > 1){
  for(i in 2:rx_s_length){ TS_m_new[i,] = TS_m[i-1,] - mts[i-1,]*TS_m[i-1,] - cures[i-1,]*TS_m[i-1,] - fails[i-1,]*TS_m[i-1,] }
  }
  for(i in 2:rx_r_length){ TR_m_new[i,] = TR_m[i-1,] - mtr[i-1,]*TR_m[i-1,] - curer[i-1,]*TR_m[i-1,] - failr[i-1,]*TR_m[i-1,] }
  
  print(TS_m)
  print(TR_m)
  print(new_AR_from_rx)
  # save 
  save_TS_M <- rbind(save_TS_M,rowSums(TS_m))
  save_TR_M <- rbind(save_TR_M,rowSums(TR_m))
  
  # new
  TR_m = TR_m_new
  TS_m = TS_m_new
  
}
