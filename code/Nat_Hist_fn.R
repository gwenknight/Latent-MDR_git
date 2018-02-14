#### All functions needed 
#(1) nat_hist:  main function for MDR tb dynamics
# ?? : issues / questions / thoughts

### ****** MAIN Function to generate model output ******************************************************** 
nat_hist <- function(para_v, para_s, mort, birth, times_v, init){
  # para_v <- parameters that vary
  # para_s <- standard nat hist parameters
  # times <- c(year1, yearend, dt)
  # initial <- initial conditions. (U, LS, LR, AS, AR, TS, TR) rows. Age in columns
  
  ####**** INITIALISE ***####
  ### Parameter inputs
  # assign parameters to values inputted
  for(i in 1:length(para_v)){assign(names(para_v)[i],para_v[i])}
  for(i in 1:length(para_s)){assign(names(para_s)[i],para_s[i])}
  # correct timestep for rates (all originally with 1 year denominator)
  mort$in_value <- mort$in_value*dt; ma <- ma*dt; # mortality
  sigma <- sigma*dt # reactivation rate
  beta<-beta*dt # transmission rate per year
  wr<-wr*dt; ws<-ws*dt # treatment rate per year
  
  ### Time inputs
  # Timesteps with inputted start end and timesteps
  year1<-times_v[1]; yearend <- times_v[2]
  dt <- times_v[3]
  times<-seq(year1, yearend, dt);
  steps<-length(times);
  if(dt > 0.5){print("ERROR: dt too big")} # need to be at least as long as DS treatment
  
  # Length of treatment
  rx_s_length <- 0.5 / dt
  rx_r_length <- 1.5 / dt
  
  ### Initialise populations (zeros, row = time, col = age)
  U  <- matrix(0,steps,Mnage); LS <- matrix(0,steps,Mnage); AS <- matrix(0,steps,Mnage)
  LR <- matrix(0,steps,Mnage); AR <- matrix(0,steps,Mnage)
  TR <- matrix(0,steps,Mnage); TS <- matrix(0,steps,Mnage);
  U_p  <- matrix(0,steps,Mnage); LS_p <- matrix(0,steps,Mnage); AS_p <- matrix(0,steps,Mnage)
  LR_p <- matrix(0,steps,Mnage); AR_p <- matrix(0,steps,Mnage)
  TR_p <- matrix(0,steps,Mnage); TS_p <- matrix(0,steps,Mnage);
  # From inputted initial conditions
  U[1,] <- init[1,]; LS[1,] <- init[2,]; LS[1,] <- init[3,]
  AS[1,] <- init[4,]; AR[1,] <- init[5,]; 
  TS[1,] <- init[6,]; TR[1,] <- init[7,];
  # none in previously treated!
  # transmission
  lambdaS <- matrix(0,1,steps)
  lambdaR <- matrix(0,1,steps)
  # new infections / reactivations
  new_AR_inf <- matrix(0,steps,Mnage);  new_AS_inf <- matrix(0,steps,Mnage);
  new_AR_rea <- matrix(0,steps,Mnage);  new_AS_rea <- matrix(0,steps,Mnage);
  new_AR_inf_p <- matrix(0,steps,Mnage);  new_AS_inf_p <- matrix(0,steps,Mnage);
  new_AR_rea_p <- matrix(0,steps,Mnage);  new_AS_rea_p <- matrix(0,steps,Mnage);
  # treatment matrix
  TS_m <- matrix(0,rx_s_length,Mnage); TS_m_new <- matrix(0,rx_s_length,Mnage)
  TR_m <- matrix(0,rx_r_length,Mnage); TR_m_new <- matrix(0,rx_r_length,Mnage)
  TS_m_p <- matrix(0,rx_s_length,Mnage); TS_m_new_p <- matrix(0,rx_s_length,Mnage)
  TR_m_p <- matrix(0,rx_r_length,Mnage); TR_m_new_p <- matrix(0,rx_r_length,Mnage)
  
  ## build output matrices
  psize <- matrix(0,1,steps); inc <- matrix(0,2,steps); prev <- matrix(0,4,steps); ratio_mdr <- matrix(0,2,steps);
  psize_age <- matrix(0,steps,Mnage)
  
  ### At start set up year tracker
  year = year1
  psize[1]<-sum(U[1,],LS[1,],LR[1,],AS[1,],AR[1,],TS[1,],TR[1,]) # no previously treated
  psize_age[1,] <- U[1,] + LS[1,] + LR[1,] + AS[1,] + AR[1,] + TS[1,] + TR[1,] + LS_p[1,] + LR_p[1,] + AS_p[1,] + AR_p[1,] + TS_p[1,] + TR_p[1,]
  #print(c("psize",psize[1]))
  
  ######################## ******** RUN ******** ######################################################################################
  ### If start of year
  for (i in 2:steps){
    print(c(i,year))
    #print(i)
    year <- year + dt
    
    ###**** Birth and deaths ***####
    birth_rate <- birth[which(birth$year == round(year,0)), "in_births"] # "in" gives interpolated
    births = dt * birth_rate * psize[i-1] # occur over the year not just at start
    m <- mort[which(mort$year == round(year,0)),"in_value"]
    #print(c("births",births, dt, birth_rate, psize[i-1], year))
    
    ####***** TB model ****####
    
    # Transmission
    lambdaS[i-1] =     (beta/psize[i-1]) * sum(p_i*(AS[i-1,] + AS_p[i-1,])); 
    lambdaR[i-1] = f * (beta/psize[i-1]) * sum(p_i*(AR[i-1,] + AR_p[i-1,])) # else 0 (in initial conditions)
    #print(c("beta",beta,psize[i-1]))
    
    ###** Treatment matrix bit special **####
    # In treatment now is sum from last time
    TR[i,2:Mnage] = colSums(TR_m[,2:Mnage])
    TS[i,2:Mnage] = colSums(TR_m[,2:Mnage])
    TR_p[i,2:Mnage] = colSums(TR_m_p[,2:Mnage])
    TS_p[i,2:Mnage] = colSums(TR_m_p[,2:Mnage])
    
    # leave treatment now are:
    new_AS_from_rx <- colSums(fails*TS_m[,2:Mnage]) # leave + not cured (if during treatment) = fail treatment 
    new_AR_from_rx <- colSums(failr*TR_m[,2:Mnage]) + eps * TS_m[rx_s_length,2:Mnage] # leave + not cured = fail treatment + acquisitions
    new_LS_from_rx <- colSums(cures*TS_m[,2:Mnage]) # leave + cured
    new_LR_from_rx <- colSums(curer*TR_m[,2:Mnage]) # leave + cured
    
    new_AS_from_rx_p <- colSums(fails*TS_m_p[,2:Mnage]) # fail treatment 
    new_AR_from_rx_p <- colSums(failr*TR_m_p[,2:Mnage]) + eps * TS_m_p[rx_s_length,2:Mnage] # fail treatment + acquisitions
    new_LS_from_rx_p <- colSums(cures*TS_m_p[,2:Mnage]) # cured
    new_LR_from_rx_p <- colSums(curer*TR_m_p[,2:Mnage]) # cured
    
    # new
    TS_m_new[1,]   = ws * AS[i-1,1:Mnage]
    TR_m_new[1,]   = wr * AR[i-1,1:Mnage]
    TS_m_new_p[1,] = ws * AS_p[i-1,1:Mnage]
    TR_m_new_p[1,] = wr * AR_p[i-1,1:Mnage]
    
    ## update TR/TS matrix
    # lose bottom row (rx_r/s_length) lost - to failr / cure / death / acquisition(s-r)
    # first column = age 0 = not changed but none? and not read out. 
    # end column = dead = drops off
    if(rx_s_length > 1){
      for(jj in 2:Mnage){
        for(ii in 2:rx_s_length){ TS_m_new[ii,jj] = TS_m[ii-1,jj-1]*(1 - mts[ii-1,jj-1] - cures[ii-1,jj-1] - fails[ii-1,jj-1]);
                                  TS_m_new_p[ii,jj] = TS_m_p[ii-1,jj-1]*(1 - mts[ii-1,jj-1] - cures[ii-1,jj-1] - fails[ii-1,jj-1])}
      }
    }
    for(jj in 2:Mnage){
      for(ii in 2:rx_r_length){ TR_m_new[ii,jj] = TR_m[ii-1,jj-1]*(1 - mtr[ii-1,jj-1] - curer[ii-1,jj-1] - failr[ii-1,jj-1]);
                                TR_m_new_p[ii,jj] = TR_m_p[ii-1,jj-1]*(1 - mtr[ii-1,jj-1] - curer[ii-1,jj-1] - failr[ii-1,jj-1])}
    }
    # save 
    #save_TS_M <- rbind(save_TS_M,rowSums(TS_m))
    #save_TR_M <- rbind(save_TR_M,rowSums(TR_m))
    
    # new
    TR_m = TR_m_new; TR_m_p = TR_m_new_p
    TS_m = TS_m_new; TS_m_p = TS_m_new_p
    # END TREATMENT MATRIX
    
    ####**** Standard dynamics ***######
    upp <- Mnage - 1 # top age - 1
    
    U [i,1] = births # spread out over year so new ones in each time step
    U [i,2:Mnage] = U [i-1,1:upp] - (m[1:upp]+lambdaS[i-1]+lambdaR[i-1])*U[i-1,1:upp] # start = 1 only at begin of year
    LS[i,2:Mnage] = LS[i-1,1:upp] + lambdaS[i-1]*(1 - p)*(U[i-1,1:upp] + x*LR[i-1,1:upp]) - (sigma + m[1:upp] + x*(lambdaS[i-1]*p + lambdaR[i-1]) )*LS[i-1,1:upp]
    LR[i,2:Mnage] = LR[i-1,1:upp] + lambdaR[i-1]*(1 - p)*(U[i-1,1:upp] + x*LS[i-1,1:upp]) - (sigma + m[1:upp] + x*(lambdaR[i-1]*p + lambdaS[i-1]) )*LR[i-1,1:upp]
    #print(c(m[1:upp],lambdaS[i-1],lambdaR[i-1]))
    
    new_AR_inf[i,2:Mnage] = lambdaR[i-1]*p*( U[i-1,1:upp] + x*(LS[i-1,1:upp] + LR[i-1,1:upp]) )
    new_AR_rea[i,2:Mnage] = sigma * LR[i-1,1:upp] 
    
    new_AS_inf[i,2:Mnage] = lambdaS[i-1]*p*( U[i-1,1:upp] + x*(LS[i-1,1:upp] + LR[i-1,1:upp]) )
    new_AS_rea[i,2:Mnage] = sigma * LS[i-1,1:upp] 
    
    AR[i,2:Mnage] = AR[i-1,1:upp] + new_AR_inf[i,2:Mnage] + new_AR_rea[i,2:Mnage] - (wr + m[1:upp] + ma)*AR[i-1,1:upp] 
    AS[i,2:Mnage] = AS[i-1,1:upp] + new_AS_inf[i,2:Mnage] + new_AS_rea[i,2:Mnage] - (ws + m[1:upp] + ma)*AS[i-1,1:upp] 
    
    # previous_treated
    LS_p[i,2:Mnage] = LS_p[i-1,1:upp] + lambdaS[i-1]*(1 - p)*(x*LR_p[i-1,1:upp]) + new_LS_from_rx + new_LS_from_rx_p - (sigma + m[1:upp] + x*(lambdaS[i-1]*p + lambdaR[i-1]) )*LS_p[i-1,1:upp]
    LR_p[i,2:Mnage] = LR_p[i-1,1:upp] + lambdaR[i-1]*(1 - p)*(x*LS_p[i-1,1:upp]) + new_LR_from_rx + new_LR_from_rx_p - (sigma + m[1:upp] + x*(lambdaR[i-1]*p + lambdaS[i-1]) )*LR_p[i-1,1:upp]
    
    new_AR_inf_p[i,2:Mnage] = lambdaR[i-1]*p*(x*(LS_p[i-1,1:upp] + LR_p[i-1,1:upp]))
    new_AR_rea_p[i,2:Mnage] = sigma * LR_p[i-1,1:upp] 
    
    new_AS_inf_p[i,2:Mnage] = lambdaS[i-1]*p*(x*(LS_p[i-1,1:upp] + LR_p[i-1,1:upp]))
    new_AS_rea_p[i,2:Mnage] = sigma * LS_p[i-1,1:upp] 
    
    AR_p[i,2:Mnage] = AR_p[i-1,1:upp] + new_AR_inf_p[i,2:Mnage] + new_AR_rea_p[i,2:Mnage] - (wr + m[1:upp] + ma)*AR_p[i-1,1:upp] + new_AR_from_rx + new_AR_from_rx_p
    AS_p[i,2:Mnage] = AS_p[i-1,1:upp] + new_AS_inf_p[i,2:Mnage] + new_AS_rea_p[i,2:Mnage] - (ws + m[1:upp] + ma)*AS_p[i-1,1:upp] + new_AS_from_rx + new_AS_from_rx_p
    
    ### look at output
    # print(c("U","LS","LR","AS","AR","TR","TS","LS_p","LR_p","AS_p","AR_p","TR_p","TS_p"))
    # print(c(sum(U[i,1:Mnage]), sum(LS[i,1:Mnage]), sum(LR[i,1:Mnage]), sum(AS[i,1:Mnage]),
    #       sum(AR[i,1:Mnage]), sum(TR[i,1:Mnage]), sum(TS[i,1:Mnage]), sum(LS_p[i,1:Mnage]),
    #       sum(LR_p[i,1:Mnage]), sum(AS_p[i,1:Mnage]), sum(AR_p[i,1:Mnage]),
    #       sum(TR_p[i,1:Mnage]), sum(TS_p[i,1:Mnage])))
    
        
    ####**** Summary Output ****#####
    ## POPULATION SIZE
    psize[i]<-sum(U[i,],LS[i,],LR[i,],AS[i,],AR[i,],TS[i,],TR[i,],
                  LS_p[i,],LR_p[i,],AS_p[i,],AR_p[i,],TS_p[i,],TR_p[i,])
    psize_age[i,] <- U[i,] + LS[i,] + LR[i,] + AS[i,] + AR[i,] + TS[i,] + TR[i,] + 
                         LS_p[i,] + LR_p[i,] + AS_p[i,] + AR_p[i,] + TS_p[i,] + TR_p[i,]
    #
    ## Incidence per year
    inc[1,i] <-  1/dt*100000*sum(new_AS_inf[i,1:Mnage] + new_AS_rea[i,1:Mnage] + new_AS_inf_p[i,1:Mnage] + new_AS_rea_p[i,1:Mnage])/psize[i]
    inc[2,i] <-  1/dt*100000*sum(new_AR_inf[i,1:Mnage] + new_AR_rea[i,1:Mnage] + new_AR_inf_p[i,1:Mnage] + new_AR_rea_p[i,1:Mnage])/psize[i]
    
    ## Prevalence of active (1: S, 2: R) and latent [3: S, 4: R]
    prev[1,i] <- 100000*sum(AS[i,1:Mnage] + AS_p[i,1:Mnage])/psize[i]
    prev[2,i] <- 100000*sum(AR[i,1:Mnage] + AR_p[i,1:Mnage])/psize[i]
    prev[3,i] <- 100000*sum(LS[i,1:Mnage] + LS_p[i,1:Mnage])/psize[i]
    prev[4,i] <- 100000*sum(LR[i,1:Mnage] + LR_p[i,1:Mnage])/psize[i]
    
    ## Ratio of MDR: 
    # 1) how many new are MDR? 
    # 2) how many prev treat are MDR? 
    if(sum(new_AR_inf_p[i,1:Mnage] + new_AR_rea_p[i,1:Mnage] + 
           new_AS_inf_p[i,1:Mnage] + new_AS_rea_p[i,1:Mnage]) > 0){
    ratio_mdr[1,i] <- 100*sum(new_AR_inf[i,1:Mnage] + new_AR_rea[i,1:Mnage]) / 
      sum(new_AR_inf[i,1:Mnage] + new_AR_rea[i,1:Mnage] + 
            new_AS_inf[i,1:Mnage] + new_AS_rea[i,1:Mnage])
    ratio_mdr[2,i] <- 100*sum(new_AR_inf_p[i,1:Mnage] + new_AR_rea_p[i,1:Mnage]) / 
      sum(new_AR_inf_p[i,1:Mnage] + new_AR_rea_p[i,1:Mnage] + 
            new_AS_inf_p[i,1:Mnage] + new_AS_rea_p[i,1:Mnage])
    }
    
    # # Number of TB deaths in HIV-, in HIV+, all form HIV deaths
    # TBDeaths[i,]=ma * (AR[i-1,1:upp] + AS[i-1,]) + ma[1:upp] * (U[i-1,] + LR[i-1,] + LS[i-1,]);
    
  }
  
  ###*** What outputting ***####
  return(list(U = U, LS = LS, LR = LR, AS = AS, AR = AR, TS = TS, TR = TR, 
              LS_p = LS_p, LR_p = LR_p, AS_p = AS_p, AR_p = AR_p, TS_p = TS_p, TR_p = TR_p, 
              lambdaR = lambdaR, lambdaS = lambdaS, TR_m = TR_m, TS_m = TS_m,
              psize = psize, inc = inc, prev = prev, ratio_mdr = ratio_mdr, psize_age = psize_age))
}