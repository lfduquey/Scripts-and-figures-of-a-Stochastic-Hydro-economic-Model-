################################################################################################################
# Script for exploring the economic effectiveness of river flood warning systems
# Author: Luis F. Duque
# Date: 12 January 2026

################################################################################################################
################################################################################################################
################################################################################################################
# Analytical framework parameters

# Statistical moments of the River Wansbeck at Mitford (river system parameters)
loguy<--1.06
logsdy<-0.8
uy<-exp(loguy+0.5*logsdy^2)
sdy<-(exp(2*loguy+logsdy^2)*(exp(logsdy^2)-1))^0.5
CVy<-sdy/uy
yo<-2
m<-1.60
L<-6 # (assumed)

# Flood forecasting system parameters
IPrho<-0.85
uyhat<-uy
sdyhat<-sdy
yohat<-yo
LT=as.list(c(6))

# At-risk community and Monte Carlo integration parameters
n<-50000
Tf<-5
nThouses<-1000
Vtheta<-1
alpha<-0.5
gama<-0.1
################################################################################################################
################################################################################################################
# Calculations



################################################################################################################
# Calculate the flood threshold
Pyf<-(1/(Tf*m))
yT<-qlnorm(Pyf,loguy,logsdy,lower.tail = FALSE)+yo
ybankfull<-qlnorm((1/(2.5*m)),loguy,logsdy,lower.tail = FALSE)+yo


################################################################################################################
# Impact curve estimation: affected properties as a percentage of the 200-year flood impact

Function_AHouses<-function(loguy,logsdy,yo,m,nThouses){
  
  AhousesPer<-c(0,5,10,25,80,93,100)
  Tr<-c(2.5,5,10,25,50,100,200)
  Ahouses=AhousesPer*nThouses/100
  #Affected_Houses<-data.frame(AhousesPer,Ahouses=AhousesPer*nThouses/100,Tr) 
  
  Py2_5<-(1/(2.5*m))
  Py5<-(1/(5*m))
  Py10<-(1/(10*m))
  Py25<-(1/(25*m))
  Py50<-(1/(50*m))
  Py100<-(1/(100*m))
  Py200<-(1/(200*m))
  
  y2_5<-qlnorm(Py2_5,loguy,logsdy,lower.tail = FALSE)+yo
  y5<-qlnorm(Py5,loguy,logsdy,lower.tail = FALSE)+yo
  y10<-qlnorm(Py10,loguy,logsdy,lower.tail = FALSE)+yo
  y25<-qlnorm(Py25,loguy,logsdy,lower.tail = FALSE)+yo
  y50<-qlnorm(Py50,loguy,logsdy,lower.tail = FALSE)+yo
  y100<-qlnorm(Py100,loguy,logsdy,lower.tail = FALSE)+yo
  y200<-qlnorm(Py200,loguy,logsdy,lower.tail = FALSE)+yo
  
  yTr<-c(y2_5,y5,y10,y25,y50,y100,y200)
  
  TrInt<-seq(2.5,200,by=.5)
  AhousesPerInt<-approx(Tr,AhousesPer,TrInt)$y
  AhousesInt<-approx(Tr,Ahouses,TrInt)$y
  yTrInt<-approx(Tr,yTr,TrInt)$y
  
  Affected_Houses<-data.frame(AhousesPer=AhousesPerInt,Ahouses=AhousesInt,yTr=yTrInt,Tr=TrInt)
  
  return(Affected_Houses)
}   
Affected_Houses<-Function_AHouses(loguy,logsdy,yo,m,nThouses)
Affected_Houses$AhousesPer[Affected_Houses$Tr<Tf]<-0    # Set all values below the flood threshold to zero
Affected_Houses$Ahouses[Affected_Houses$Tr<Tf]<-0       ## Set all values below the flood threshold to zero
#---------------------------------------------------------------
#---------------------------------------------------------------
# Function to perform SA of the forecasting system performance

LT_PWS<-function(uy,sdy,yo,uyhat,sdyhat,yohat,LT,n,IPrho,L,yT,ybankfull,Affected_Houses,Vtheta,alpha,gama){
  
  
  ################################################################################################################
  # Monte Carlo flood and forecast generator (BLND)
  
  MCFG<-function(uy,sdy,yo,uyhat,sdyhat,yohat,LT,n,IPrho,L){
    
    
    #--------------------------------------------
    #--------------------------------------------
    # # Calculate summary statistics for log(y) and log(yhat)
    
    # observed values y
    uz <- 2*log(uy)-0.5*log(sdy^2 + uy^2)
    vz <- -2*log(uy)+log(sdy^2 + uy^2)
    sdz <- sqrt(vz)
    
    # predicted values yhat
    uzhat <- 2*log(uyhat)-0.5*log(sdyhat^2 + uyhat^2)
    vzhat <- -2*log(uyhat)+log(sdyhat^2 + uyhat^2)
    sdzhat <- sqrt(vzhat)
    
    
    #--------------------------------------------
    #--------------------------------------------
    #  Calculate the correlation coefficient in real space
    
    Lead_Time<- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
    rhof<-ifelse(Lead_Time<=L,-0.015*(Lead_Time-L)+IPrho,-0.03*(Lead_Time-L)+IPrho)
    rho<-approx(Lead_Time, rhof, LT, rule = 2)$y
    #--------------------------------------------
    #--------------------------------------------
    # Calculate the correlation coefficient in normal space
    a<-((exp(vz)-1)*(exp(vzhat)-1))^0.5
    rhoN<-log(a*rho+1)/(sdz*sdzhat)
    
    #--------------------------------------------
    #--------------------------------------------
    # Calculate standardized values of zhat
    
    etahat<- rnorm(n,0,1)
    
    #--------------------------------------------
    #--------------------------------------------
    # Observed values conditional on etahat and associated with rhoN
    Cmeaneta<-sapply(etahat,function(x){rhoN*x})
    Csdeta<-sqrt(1-rhoN^2)
    eta<-sapply(Cmeaneta, function(x){rnorm(1,mean=x,sd=Csdeta)})
    
    #--------------------------------------------
    #--------------------------------------------
    # Calculate the logarithms of the observed values and their forecasts
    z<-uz+eta*sdz
    zhat<-uzhat+etahat*sdzhat
    
    #--------------------------------------------
    #--------------------------------------------
    # Obtain the observed values and their forecasts
    
    y<-exp(z)+yo
    yhat<-exp(zhat)+yohat
    
    
    #--------------------------------------------
    #--------------------------------------------
    # Obtain EP for the observed values and their forecasts
    
    ep_y<-plnorm(y-yo,meanlog =uz, sdlog = sdz, lower.tail = FALSE)
    ep_yhat<-plnorm(yhat-yohat,meanlog =uzhat, sdlog = sdzhat, lower.tail = FALSE)
    
    
    #--------------------------------------------
    #--------------------------------------------
    # Joint pairs of values
    Pairsyyhat<-data.frame(yhat,y,ep_y,ep_yhat)
    
    
    #--------------------------------------------
    #--------------------------------------------
    # Compute predictive uncertainty PU
    
    p<-seq(0.001,0.999,1/1000) # Probabilities of conditional quantiles
    MPU_eta_N<-sapply(Cmeaneta,function(x){qnorm(p,mean=x,sd=Csdeta)})
    MPU_z_N<-uz+MPU_eta_N*sdz
    MPU_y<-exp(MPU_z_N)+yo
    
    MPU_PDF_y<-apply(MPU_y,2,function (x){density(x)$y})
    MPU_y_grid<-apply(MPU_y,2,function (x){density(x)$x})
    
    cdfy<-function(x){
      CDFy_function<-ecdf(x)
      y<-CDFy_function(x)
    }
    MPU_CDF_y<-apply(MPU_y_grid,2,cdfy)
    
    PU<-list()
    for (i in 1:n){
      PU[[paste0("Cond_",i)]]<- data.frame(y=MPU_y_grid[,i],PDF=MPU_PDF_y[,i],CDF=MPU_CDF_y[,i])
    }
    
    #--------------------------------------------
    #--------------------------------------------
    # Export results
    
    Results<-list()
    Results[["Pairsyyhat"]]<-Pairsyyhat
    Results[["PU"]]<-PU
    Results[["rho"]]<-rho
    return(Results)
    
  }
  Obs_For_Peaks<-MCFG(uy=uy,sdy=sdy,yo=yo,uyhat=uyhat,sdyhat=sdyhat,yohat=yohat,LT=LT,n=n,IPrho=IPrho,L=L)
  
  ################################################################################################################
  # Deterministic decision making
  
  df<-Obs_For_Peaks$Pairsyyhat
  DDM<-function(df,yT,ybankfull,Affected_Houses,Vtheta,alpha,gama,LT){
    
    
    #---------------------------
    #---------------------------
    # Function of losses
    
    
    
    Losses<-function(df_w,Vtheta,alpha,gama,LT){
      
      
      #------------------------------------------
      #------------------------------------------
      #------------------------------------------
      # Percentage contents damage for a single house under different 
      # mitigation times and flood depths (Carsell curves)
      
      depthMT<-c(-0.91,-0.61,-0.3,0,0.3,0.61,0.91,1.22,1.52,1.83,2.44,3.05,4.57,6.1)  
      #MT0
      D0<- c(0,0,0,10,17,23,29,35,40,45,55,60,60,60) 
      #MT1
      D1<- c(0,0,0,8,15,21,26,32,37,42,51,55,55,55)
      #MT6
      D6<-c(0 ,0 ,0 ,6 ,12 ,17 ,22 ,27,33,36 ,44 ,48 ,48,48) 
      #MT12
      D12<-c(0 ,0 ,0 ,5 ,11 ,16 ,21 ,27 ,32 ,36 ,43 ,46,47,47) 
      #MT24
      D24<-c(0,0,0,4,9,13,16,23,27,31,38,40,40,40) 
      #MT36
      D36<-c(0,0,0,4,9,11,15,20,25,27,33,35,35,35) 
      #MT48
      D48<-c(0,0,0,4,8,11,14,19,23,25,30,31,31,31) 
      
      #------------------------------------------
      #------------------------------------------
      #------------------------------------------
      # Calculate unmitigated damage (Lu) for a single house and all houses
      
      DMT0<- approx(depthMT,D0,df_w$depth,rule = 2)$y 
      df_w$Lu<- round(DMT0*Vtheta/100*df_w$housesF,2)
      
      
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #----------------------------------------------------------------------------------------------------------- 
      # Calculate mitigated damage with lead time (LT) for a single house
      
      if(unlist(LT)==1) {DMT<-approx(depthMT,D1,df_w$depth,rule = 2)$y
      
      } else if (unlist(LT)==6) {
        
        DMT<-approx(depthMT,D6,df_w$depth,rule = 2)$y} else if (unlist(LT)==12){
          
          DMT<-approx(depthMT,D12,df_w$depth,rule = 2)$y} else if (unlist(LT)==24){
            
            DMT<-approx(depthMT,D24,df_w$depth,rule = 2)$y} else if(unlist(LT)==36){
              
              DMT<-approx(depthMT,D36,df_w$depth,rule = 2)$y} else {
                
                DMT<-approx(depthMT,D48,df_w$depth,rule = 2)$y
              }
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #----------------------------------------------------------------------------------------------------------- 
      # Obtain the forecast umitigated damage / single house
      DMT0hat<- approx(depthMT,D0,df_w$depth_hat,rule = 2)$y 
      
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #----------------------------------------------------------------------------------------------------------- 
      # Calculate forecast-based mitigated damage with LT for a single house
    
      if(unlist(LT)==1) {DMThat<-approx(depthMT,D1,df_w$depth_hat,rule = 2)$y
      
      } else if (unlist(LT)==6) {DMThat<-approx(depthMT,D6,df_w$depth_hat,rule = 2)$y} else if (unlist(LT)==12){
        
        DMThat<-approx(depthMT,D12,df_w$depth_hat,rule = 2)$y} else if (unlist(LT)==24){
          
          DMThat<-approx(depthMT,D24,df_w$depth_hat,rule = 2)$y} else if(unlist(LT)==36){
            
            DMThat<-approx(depthMT,D36,df_w$depth_hat,rule = 2)$y} else {
              
              DMThat<-approx(depthMT,D48,df_w$depth_hat,rule = 2)$y
            }
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      # Compute losses for y
      
      # Misses
      df_w$D_misses<-ifelse(df_w$TypeE=="M",round(DMT0/100*Vtheta*df_w$housesF,2),0)
      
      # Damage Hits 1
      df_w$Dm_H1<-ifelse(df_w$TypeE=="H1",round(Vtheta*df_w$housesF*((DMT0/100)-alpha*((DMT0-DMT)/100)),2),0)
      
      # Damage Hits 2
      df_w$Du_H2<-ifelse(df_w$TypeE=="H2",round(Vtheta*DMT0/100*df_w$housesNWF,2),0)
      df_w$Dm_H2<-ifelse(df_w$TypeE=="H2",round(Vtheta*df_w$housesWF*(DMT0/100-alpha*((DMT0-DMT)/100)),2),0)
      df_w$Dum_H2<-df_w$Du_H2+df_w$Dm_H2               
      
      # Perfect forecast
      df_w$Lmc_pf<-Vtheta*df_w$housesF*(DMT0/100-alpha*(DMT0-DMT)/100)+ gama*alpha*Vtheta*df_w$housesF*((DMT0-DMT)/100)
      
      
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      # Compute economic consequences for forecasts of y
      
      # Cw false alarms
      
      df_w$Cw_FA<-ifelse(df_w$TypeE=="FA",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
      
      # Cw hits H1
      
      df_w$Cw_H1<-ifelse(df_w$TypeE=="H1",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
      
      
      # Cw hits H2
      df_w$Cw_H2<-ifelse(df_w$TypeE=="H2",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
      
      # Net damage hits 1
      df_w$NetD_H1<-df_w$Cw_H1+df_w$Dm_H1
      
      
      # Net damage hits 2
      df_w$NetD_H2<-df_w$Cw_H2+df_w$Dum_H2
      
      
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      # Economic consequences of imperfect forecasts
      df_w$EC_Imf<-ifelse(df_w$TypeE=="H1",df_w$NetD_H1,ifelse(df_w$TypeE=="H2",df_w$NetD_H2,
                                                               ifelse(df_w$TypeE=="M",df_w$D_misses,
                                                                      ifelse(df_w$TypeE=="FA",df_w$Cw_FA,0))))
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      # Computation of EAD
      
      
      
      #-------------------------------------------
      #-------------------------------------------
      # Non-warning
      depth_NWPS<-df_w[df_w$depth>=0,]
      histogram_fd<-function(x,minx,max,delta){
        
        B1x<-seq(minx,max,by=delta)
        B2x<-c(B1x[-1],max)
        
        Counts<-matrix(0, nrow =length(x),ncol = length(B1x))
        
        
        for (i in 1:length(B1x)) {
          
          if (i==length(B1x)) { 
            Counts[,i]<-ifelse(x>=B1x[i],1,0)} else {
              Counts[,i]<-ifelse(x>=B1x[i] & x < B2x[i],1,0)
            }
          
          
        }
        
        Count<-colSums(Counts)
        rFrequency<-Count/sum(Count)
        hist_Counts<-data.frame(B1x,B2x,Count,rFrequency)
        hist_Counts$Midpoints<-(B1x+B2x)/2
        return(hist_Counts)
        
      }
      hist_depth<-histogram_fd(depth_NWPS$depth,minx=0,max=15,delta=0.25)
      depth_mI<-hist_depth$Midpoints
      damage_NW<-approx(depth_NWPS$depth,depth_NWPS$Lu,depth_mI)$y
      depth_freq<-hist_depth$Count/n
      EAD_nw<-sum(depth_freq*damage_NW, na.rm = TRUE)
      n<-length(df_w$Lu)
      
      
      #-------------------------------------------
      #-------------------------------------------
      # Perfect forecast
      
      # Perfect forecast
      
      damage_PS<-approx(depth_NWPS$depth,depth_NWPS$Lmc_pf,depth_mI)$y
      EAD_pf<-sum(depth_freq*damage_PS, na.rm = TRUE)
      
      #-------------------------------------------
      #-------------------------------------------
      # Deterministic forecast
      # Misses
      depth_misses<-df_w[df_w$m==1,]
      hist_depth_misses<-histogram_fd(depth_misses$depth,minx=0,max=15,delta=0.25)
      depth_mI_misses<-hist_depth_misses$Midpoints
      damage_misses<-approx(depth_misses$depth,depth_misses$D_misses,depth_mI_misses)$y
      depth_freq_misses<-hist_depth_misses$Count/n
      D_misses<-sum(depth_freq_misses*damage_misses, na.rm = TRUE)
      
      #-----------------
      # Hits 
      
      depth_hits1<-df_w[df_w$TypeE=="H1",]
      
      # Damage hits 1
      
      hist_depth_hits1<-histogram_fd(depth_hits1$depth,minx=0,max=15,delta=0.25)
      depth_mI_hits1<-hist_depth_hits1$Midpoints
      damage_hits1<-approx(depth_hits1$depth,depth_hits1$Dm_H1,depth_mI_hits1)$y
      depth_freq_hits1<-hist_depth_hits1$Count/n
      D_hits1<-sum(depth_freq_hits1*damage_hits1, na.rm = TRUE)
      
      
      # Damage hits 2
      
      depth_hits2<-df_w[df_w$TypeE=="H2",]
      hist_depth_hits2<-histogram_fd(depth_hits2$depth,minx=0,max=15,delta=0.25)
      depth_mI_hits2<-hist_depth_hits2$Midpoints
      damage_hits2<-approx(depth_hits2$depth,depth_hits2$Dum_H2,depth_mI_hits2)$y
      depth_freq_hits2<-hist_depth_hits2$Count/n
      D_hits2<-sum(depth_freq_hits2*damage_hits2, na.rm = TRUE)
      
      
      
      # Warning cost
      
      # Cw of Hits H1
      
      hist_depthat_hits1<-histogram_fd(depth_hits1$depth_hat,minx=0,max=15,delta=0.25)
      depthat_mI_hits1<-hist_depthat_hits1$Midpoints
      Cw_H1<-approx(depth_hits1$depth_hat,depth_hits1$Cw_H1,depthat_mI_hits1)$y
      depthat_freq_hits1<-hist_depthat_hits1$Count/n
      Cw_hits1<-sum(depthat_freq_hits1*Cw_H1, na.rm = TRUE)
      
      
      # Cw of Hits H2
      
      hist_depthat_hits2<-histogram_fd(depth_hits2$depth_hat,minx=0,max=15,delta=0.25)
      depthat_mI_hits2<-hist_depthat_hits2$Midpoints
      Cw_H2<-approx(depth_hits2$depth_hat,depth_hits2$Cw_H2,depthat_mI_hits2)$y
      depthat_freq_hits2<-hist_depthat_hits2$Count/n
      Cw_hits2<-sum(depthat_freq_hits2*Cw_H2, na.rm = TRUE)
      
      
      #-----------------
      
      # False alarms
      depthat_fa<-df_w[df_w$f==1,]
      hist_depthat_fa<-histogram_fd(depthat_fa$depth_hat,minx=0,max=15,delta=0.25)
      depthat_mI_fa<-hist_depthat_fa$Midpoints
      Cw_f<-approx(depthat_fa$depth_hat,depthat_fa$Cw_FA,depthat_mI_fa)$y
      depthat_freq_fa<-hist_depthat_fa$Count/n
      Cw_fa<-sum(depthat_freq_fa*Cw_f, na.rm = TRUE)
      #-----------------
      
      EAD_wdet<-D_misses+D_hits1+D_hits2+Cw_hits1+Cw_hits2+Cw_fa
      
    
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      #-----------------------------------------------------------------------------------------------------------
      # Export results
      
      Results<-list()
      Results[['df_EC']]<-df_w
      Results[['EAD_nw']]<-EAD_nw
      Results[['EAD_pf']]<-EAD_pf
      Results[['EAD_wdet']]<-EAD_wdet
      return(Results)
      
    }
    
    
    # Analysis of the deterministic decision making
    #---------------------------
    #---------------------------
    # Characterization of the flood events
    
    # Binary indicator of warning
    df$w<-ifelse(df$yhat>yT,1,0)
    
    # Binary indicator of flood
    df$flood<-ifelse(df$y>yT,1,0)
    
    # State of the flood event
    df$h<-ifelse(df$w==1 & df$flood==1,1,0)
    df$m<-ifelse(df$w==0 & df$flood==1,1,0)
    df$f<-ifelse(df$w==1 & df$flood==0,1,0)
    df$cn<-ifelse(df$w==0 & df$flood==0,1,0)
    
    # Types of flood events 
    df$TypeE<- rep("H",length(df$y))
    df$TypeE<-ifelse(df$m==1,'M',df$TypeE)
    df$TypeE<-ifelse(df$f==1,'FA',df$TypeE)
    df$TypeE<-ifelse(df$cn==1,'CN',df$TypeE)
    df$TypeE<-ifelse(df$h==1 & df$y<df$yhat,"H1",df$TypeE)
    df$TypeE<-ifelse(df$h==1 & df$y>df$yhat,"H2",df$TypeE)
    
    #---------------------------
    #---------------------------
    #  Analysis of affected houses 
    
    # Houses warned
    df$housesW<-ifelse(df$w==0,0,round(approx(Affected_Houses$yTr,Affected_Houses$Ahouses,df$yhat,rule = 2)$y,0))
    
    # Houses flooded
    df$housesF<-ifelse(df$y<yT,0,round(approx(Affected_Houses$yTr,Affected_Houses$Ahouses,df$y,rule = 2)$y,0))
    
    # Houses warned and flooded
    
    df$housesWF<-ifelse(df$TypeE=='H1',df$housesF,0)
    df$housesWF<-ifelse(df$TypeE=='H2',df$housesW,df$housesWF)
    
    # Houses not warned and flooded
    
    df$housesNWF<-ifelse(df$TypeE=='M',df$housesF,0)
    df$housesNWF<-ifelse(df$TypeE=='H2',df$housesF-df$housesW,df$housesNWF)
    
    # Houses warned and not flooded
    
    df$housesWNF<-ifelse(df$TypeE=='FA',df$housesW,0)
    df$housesWNF<-ifelse(df$TypeE=='H1',df$housesW-df$housesF,df$housesWNF)
    
    # Houses not warned and not flooded
    
    df$housesNWNF<-max(Affected_Houses$Ahouses)-(df$housesWF+df$housesNWF+df$housesWNF)
    
    #---------------------------
    #---------------------------
    #  Analysis of flood depth
    
    # Observed flood depth
    
    df$depth<-ifelse(df$flood==1,df$y-ybankfull-0.25,0)
    df$depth_hat<-ifelse(df$w==1,df$yhat-ybankfull-0.25,0)
    
    #---------------------------
    #---------------------------
    #  Economic consequences Flood warning scenarios
    
    Losses_Scenarios<-Losses(df_w=df,Vtheta,alpha,gama,LT)
    
    #---------------------------
    #---------------------------
    #  Effectiveness
    
    # No warning ED
    EAD_nw<-Losses_Scenarios$EAD_nw        
    # Perfect forecast ED
    EAD_PF<-Losses_Scenarios$EAD_pf        
    
    # Imperfect forecast ED
    EAD_DF<-Losses_Scenarios$EAD_wdet     
    
    # Perfect forecast
    Ew_PF<-(EAD_nw-EAD_PF)/EAD_nw*100     
    # Imperfect forecast
    Ew_DF<-(EAD_nw-EAD_DF)/EAD_nw*100     
    
    
    #---------------------------
    #---------------------------
    #  Results of warning issues
    

    
    # Affected property-based approach
    Hits_houses<-sum(df$housesWF)
    Hits1_houses<-sum(df$housesWF[df$TypeE=="H1"])
    Hits2_houses<-sum(df$housesWF[df$TypeE=="H2"])
    Hits1_fa_houses<-sum(df$housesWNF[df$TypeE=="H1"])
    Hits2_m_houses<-sum(df$housesNWF[df$TypeE=="H2"])
    Miss_houses<-sum(df$housesNWF)
    FA_houses<-sum(df$housesWNF)
    CN_houses<-sum(df$housesNWNF)
    #---------------------------
    #---------------------------
    # Compute skill scores FAR, POD and CSI
    

    
    # Affected property-based approach
    POD_houses<-Hits_houses/(Hits_houses+Miss_houses)
    FAR_houses<-FA_houses/(FA_houses+Hits_houses)
    
    
    #---------------------------
    #---------------------------
    # Export Results
    Results<-list()
    RDDM<- data.frame(Total=length(df$TypeE),Hits_houses,Miss_houses,FA_houses,
                      CN_houses,Totalhouses=max(Affected_Houses$Ahouses)*length(df$TypeE),POD_houses,
                      FAR_houses,Hits1_houses,Hits2_houses,Hits1_fa_houses,Hits2_m_houses,
                      EAD_nw,EAD_PF,EAD_DF,Ew_PF,Ew_DF,stringsAsFactors=FALSE)
    Results[['Events']]<-Losses_Scenarios$df_EC
    Results[['RDDM']]<-RDDM
    return(Results)
  }
  Results_DDM<-DDM(df,yT=yT,ybankfull=ybankfull,Affected_Houses=Affected_Houses,Vtheta=Vtheta,alpha=alpha,gama=gama,LT=LT)
  
  ################################################################################################################
  # Probabilistic decision making 
  
  #---------------------------
  #---------------------------
  # Function to determine the optimal PT
  
  PU<-Obs_For_Peaks$PU
  
  Opt_PFC<-function(df,PU,yT,ybankfull,Affected_Houses,Vtheta,alpha,gama,LT){
    
    
    df$yhat<-NULL
    
    #---------------------------
    #---------------------------
    # Define the the value of PT
    # to be analysed
    
    PT<-seq(0.05,0.95,by=0.01)
    
    
    #---------------------------
    #---------------------------
    # Function to represent Pro. warning decisions
    
    PFC<-function(df,yT,ybankfull,Affected_Houses,Vtheta,alpha,gama,LT){
      
    #---------------------------
    #---------------------------
    # Function of losses
      
      Losses<-function(df_w,Vtheta,alpha,gama,LT){
        
        
        #------------------------------------------
        #------------------------------------------
        #------------------------------------------
        # Damage to content for a single house as a percentage for a different mitigation time and flood depths
        
        depthMT<-c(-0.91,-0.61,-0.3,0,0.3,0.61,0.91,1.22,1.52,1.83,2.44,3.05,4.57,6.1)  # Original Carsell Curves
        #MT0
        D0<- c(0,0,0,10,17,23,29,35,40,45,55,60,60,60) # Original Carsell Curves
        #MT1
        D1<- c(0,0,0,8,15,21,26,32,37,42,51,55,55,55)
        #MT6
        D6<-c(0 ,0 ,0 ,6 ,12 ,17 ,22 ,27,33,36 ,44 ,48 ,48,48) # Original Carsell Curves
        #MT12
        D12<-c(0 ,0 ,0 ,5 ,11 ,16 ,21 ,27 ,32 ,36 ,43 ,46,47,47) # Original Carsell Curves
        #MT24
        D24<-c(0,0,0,4,9,13,16,23,27,31,38,40,40,40) # Original Carsell Curves
        #MT36
        D36<-c(0,0,0,4,9,11,15,20,25,27,33,35,35,35) # Original Carsell Curves
        #MT48
        D48<-c(0,0,0,4,8,11,14,19,23,25,30,31,31,31) # Original Carsell Curves
       
        
        
        #------------------------------------------
        #------------------------------------------
        #------------------------------------------
        # # Obtain Lu/ single and all houses
        DMT0<- approx(depthMT,D0,df_w$depth,rule = 2)$y # Interpolation with original Carsell curves
        df_w$Lu<- round(DMT0*Vtheta/100*df_w$housesF,2)
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------- 
        # Obtain umitigated damage with LT/ single house
        
        # Interpolation with original Carsell curves
        if(unlist(LT)==1) {DMT<-approx(depthMT,D1,df_w$depth,rule = 2)$y
        
        } else if (unlist(LT)==6) {
          
          DMT<-approx(depthMT,D6,df_w$depth,rule = 2)$y} else if (unlist(LT)==12){
            
            DMT<-approx(depthMT,D12,df_w$depth,rule = 2)$y} else if (unlist(LT)==24){
              
              DMT<-approx(depthMT,D24,df_w$depth,rule = 2)$y} else if(unlist(LT)==36){
                
                DMT<-approx(depthMT,D36,df_w$depth,rule = 2)$y} else {
                  
                  DMT<-approx(depthMT,D48,df_w$depth,rule = 2)$y
                }
        
        
        
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------- 
        # Obtain the forecast umitigated damage / single house
        DMT0hat<- approx(depthMT,D0,df_w$depth_hat,rule = 2)$y # Interpolation with original Carsell curves
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #----------------------------------------------------------------------------------------------------------- 
        # Obtain the forecast umitigated damage with LT/ single house
        
        # Interpolation with original Carsell curves
        if(unlist(LT)==1) {DMThat<-approx(depthMT,D1,df_w$depth_hat,rule = 2)$y
        
        } else if (unlist(LT)==6) {DMThat<-approx(depthMT,D6,df_w$depth_hat,rule = 2)$y} else if (unlist(LT)==12){
          
          DMThat<-approx(depthMT,D12,df_w$depth_hat,rule = 2)$y} else if (unlist(LT)==24){
            
            DMThat<-approx(depthMT,D24,df_w$depth_hat,rule = 2)$y} else if(unlist(LT)==36){
              
              DMThat<-approx(depthMT,D36,df_w$depth_hat,rule = 2)$y} else {
                
                DMThat<-approx(depthMT,D48,df_w$depth_hat,rule = 2)$y
              }
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        # Computation of losses for y
        
        # Misses
        df_w$D_misses<-ifelse(df_w$TypeE=="M",round(DMT0/100*Vtheta*df_w$housesF,2),0)
        # Damage Hits 1
        df_w$Dm_H1<-ifelse(df_w$TypeE=="H1",round(Vtheta*df_w$housesF*((DMT0/100)-alpha*((DMT0-DMT)/100)),2),0)
        
        # Damage Hits 2
        df_w$Du_H2<-ifelse(df_w$TypeE=="H2",round(Vtheta*DMT0/100*df_w$housesNWF,2),0)
        df_w$Dm_H2<-ifelse(df_w$TypeE=="H2",round(Vtheta*df_w$housesWF*(DMT0/100-alpha*((DMT0-DMT)/100)),2),0)
        df_w$Dum_H2<-df_w$Du_H2+df_w$Dm_H2               
        
        # Perfect forecast
        df_w$Lmc_pf<-Vtheta*df_w$housesF*(DMT0/100-alpha*(DMT0-DMT)/100)+ gama*alpha*Vtheta*df_w$housesF*((DMT0-DMT)/100)
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        # Computation of economic consequences for the forecasts of y
        
        # Cw false alarms
        
        df_w$Cw_FA<-ifelse(df_w$TypeE=="FA",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
        
        # Cw hits H1
        
        df_w$Cw_H1<-ifelse(df_w$TypeE=="H1",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
        
        
        # Cw hits H2
        df_w$Cw_H2<-ifelse(df_w$TypeE=="H2",round(gama*alpha*Vtheta*df_w$housesW*((DMT0hat-DMThat)/100),2),0)
        
        # Net damage hits 1
        df_w$NetD_H1<-df_w$Cw_H1+df_w$Dm_H1
        
        
        # Net damage hits 2
        df_w$NetD_H2<-df_w$Cw_H2+df_w$Dum_H2
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        # Economic consequences Imperfect forecast
        df_w$EC_Imf<-ifelse(df_w$TypeE=="H1",df_w$NetD_H1,ifelse(df_w$TypeE=="H2",df_w$NetD_H2,
                                                                 ifelse(df_w$TypeE=="M",df_w$D_misses,
                                                                        ifelse(df_w$TypeE=="FA",df_w$Cw_FA,0))))
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        # Computation of EAD
        
        
        
        #-------------------------------------------
        #-------------------------------------------
        # Non warning
        
        
        EAD_nw<-round(mean(df_w$Lu),2)
        n<-length(df_w$Lu)
        
        #-------------------------------------------
        #-------------------------------------------
        # Perfect forecast
        
        # Perfect forecast
        EAD_pf<-round(mean(df_w$Lmc_pf),2)
        
        #-------------------------------------------
        #-------------------------------------------
        # Deterministic forecast
        
        EAD_wdet<-mean(df_w$EC_Imf)
        
        
        
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        #-----------------------------------------------------------------------------------------------------------
        # Export results
        
        Results<-list()
        Results[['df_EC']]<-df_w
        Results[['EAD_nw']]<-EAD_nw
        Results[['EAD_pf']]<-EAD_pf
        Results[['EAD_wdet']]<-EAD_wdet

        
        
        return(Results)
        
      }
      
      #---------------------------
      #---------------------------
      # Characterization of the flood events
      
      # binary indicator of warning
      df$w<-ifelse(df$yw>yT,1,0)
      
      # binary indicator of flood
      df$flood<-ifelse(df$y>yT,1,0)
      
      #  The state of the flood event
      df$h<-ifelse(df$w==1 & df$flood==1,1,0)
      df$m<-ifelse(df$w==0 & df$flood==1,1,0)
      df$f<-ifelse(df$w==1 & df$flood==0,1,0)
      df$cn<-ifelse(df$w==0 & df$flood==0,1,0)
      
      
      # Characterization of the flood event  
      df$TypeE<- rep("H",length(df$y))
      df$TypeE<-ifelse(df$m==1,'M',df$TypeE)
      df$TypeE<-ifelse(df$f==1,'FA',df$TypeE)
      df$TypeE<-ifelse(df$cn==1,'CN',df$TypeE)
      df$TypeE<-ifelse(df$h==1 & df$y<df$yw,"H1",df$TypeE)
      df$TypeE<-ifelse(df$h==1 & df$y>df$yw,"H2",df$TypeE)
      
      #---------------------------
      #---------------------------
      #  Analysis of affected houses
      
      # Houses warned
      df$housesW<-ifelse(df$yw==0,0,round(approx(Affected_Houses$yTr,Affected_Houses$Ahouses,df$yw,rule = 2)$y,0))
      
      # Houses flooded
      df$housesF<-ifelse(df$y<yT,0,round(approx(Affected_Houses$yTr,Affected_Houses$Ahouses,df$y,rule = 2)$y,0))
      
      # Houses warned and flooded
      
      df$housesWF<-ifelse(df$TypeE=='H1',df$housesF,0)
      df$housesWF<-ifelse(df$TypeE=='H2',df$housesW,df$housesWF)
      
      # Houses not warned and flooded
      
      df$housesNWF<-ifelse(df$TypeE=='M',df$housesF,0)
      df$housesNWF<-ifelse(df$TypeE=='H2',df$housesF-df$housesW,df$housesNWF)
      
      # Houses warned and not flooded
      
      df$housesWNF<-ifelse(df$TypeE=='FA',df$housesW,0)
      df$housesWNF<-ifelse(df$TypeE=='H1',df$housesW-df$housesF,df$housesWNF)
      
      # Houses not warned and not flooded
      df$housesNWNF<-max(Affected_Houses$Ahouses)-(df$housesWF+df$housesNWF+df$housesWNF)
      
      
      #---------------------------
      #---------------------------
      #  Analysis of flood depth
      
      # Observed flood depth
      
      df$depth<-ifelse(df$flood==1,df$y-ybankfull-0.25,0)
    
      # forecast flood depth
      
      df$depth_hat<-ifelse(df$w==1,df$yw-ybankfull-0.25,0)
      #---------------------------
      #---------------------------
      #  Economic consequences Flood warning scenarios
      
      
      Losses_Scenarios<-Losses(df_w=df,Vtheta,alpha,gama,LT)
      
      #---------------------------
      #---------------------------
      #  Effectiveness
      
      # No warning ED
      EAD_nw<-Losses_Scenarios$EAD_nw        
      # Perfect forecast ED
      EAD_PF<-Losses_Scenarios$EAD_pf        
      # Imperfect forecast ED
      EAD_PrF<-Losses_Scenarios$EAD_wdet     
      # EFECTIVNESS
      # Perfect forecast
      Ew_PF<-(EAD_nw-EAD_PF)/EAD_nw*100             
      # Imperfect forecast
      Ew_PrF<-(EAD_nw-EAD_PrF)/EAD_nw*100                
      
      #---------------------------
      #---------------------------
      #  Results of warning issues
      
      
      
      # Affected property-based approach
      Hits_houses<-sum(df$housesWF)
      Hits1_houses<-sum(df$housesWF[df$TypeE=="H1"])
      Hits2_houses<-sum(df$housesWF[df$TypeE=="H2"])
      Hits1_fa_houses<-sum(df$housesWNF[df$TypeE=="H1"])
      Hits2_m_houses<-sum(df$housesNWF[df$TypeE=="H2"])
      Miss_houses<-sum(df$housesNWF)
      FA_houses<-sum(df$housesWNF)
      CN_houses<-sum(df$housesNWNF)
      #---------------------------
      #---------------------------
      # Compute skill score FAR and POD, and FR
      
      # Affected property-based approach
      POD_houses<-Hits_houses/(Hits_houses+Miss_houses)
      FAR_houses<-FA_houses/(FA_houses+Hits_houses)
      
      #---------------------------
      #---------------------------
      # Export Results
      Results<-list()
      CT_FPC<- data.frame(Total=length(df$TypeE),Hits_houses,Miss_houses,FA_houses,
                          CN_houses,Totalhouses=max(Affected_Houses$Ahouses)*length(df$TypeE),POD_houses,
                          FAR_houses,Hits1_houses,Hits2_houses,Hits1_fa_houses,Hits2_m_houses,
                          EAD_nw,EAD_PF,EAD_PrF,Ew_PF,Ew_PrF,stringsAsFactors=FALSE)
      
      Results[['Events']]<-Losses_Scenarios$df_EC
      Results[['CT_FPC']]<-CT_FPC
      return(Results)
    }
    
    #---------------------------
    #---------------------------
    # Performance analysis
    
    Results_PT<-list()
  
    for (i in 1:length(PT)){
      
      df$yw<-unlist(lapply(PU,function(x){approx(1-x[,3],x[,1],PT[i],rule = 2)$y}),use.names = FALSE)
      Results_PT[[paste0("PT_",PT[i])]]<-PFC(df,yT,ybankfull,Affected_Houses,Vtheta,alpha,gama,LT)
    }
    
    POD_houses<-unlist(lapply(Results_PT,function(x){x$CT_FPC$POD_houses}), use.names = FALSE)
    FAR_houses<-unlist(lapply(Results_PT,function(x){x$CT_FPC$FAR_houses}), use.names = FALSE)
    Ew_PDR<-round(unlist(lapply(Results_PT,function(x){x$CT_FPC$Ew_PrF}), use.names = FALSE),2)
    Performance_PT<-data.frame(PT,POD_houses,FAR_houses,Ew_PDR)
    
    #---------------------------
    #---------------------------
    # Define the optimal PT
    Index_OP<-which.max(Ew_PDR)
    Opt_PT<-Results_PT[[Index_OP]]$CT_FPC
    Opt_PT$PT<-PT[Index_OP]
    
    #---------------------------
    #---------------------------
    # Export results
    Results<-list()
    Results[["Performance_PT"]]<-Performance_PT
    Results[["Opt_PT"]]<-Opt_PT
    return(Results)
  }
  
  #---------------------------
  #---------------------------
  # Obtain the optimal PT and its performance
  
  Results_PFC<-Opt_PFC(df,PU=PU,yT=yT,ybankfull=ybankfull,Affected_Houses=Affected_Houses,Vtheta=Vtheta,
                       alpha=alpha,gama=gama,LT=LT)
  
  
  ###############################################################################################################
  # Export results
  Results<-list()
  Results[["Results_DDM"]]<-Results_DDM
  Results[["Results_PDR"]]<-Results_PFC
  Results[["Obs_For_Peaks"]]<-Obs_For_Peaks
  return(Results)
  
}



Results_LT<-LT_PWS(uy=uy,sdy=sdy,yo=yo,uyhat=uyhat,sdyhat=sdyhat,yohat=yohat,
                   LT=LT,n=n,IPrho=IPrho,L=L,yT=yT,ybankfull=ybankfull,
                   Affected_Houses=Affected_Houses,Vtheta=Vtheta,alpha=alpha,gama=gama)

#---------------------------------------------------------------
#---------------------------------------------------------------
# Get values of contingency tables

Performance_DFDR<-Results_LT$Results_DDM$RDDM
Op_Performance_PDR<-Results_LT$Results_PDR$Opt_PT
Performance_PDR<-Results_LT$Results_PDR$Performance_PT[,c(1,4)]

################################################################################################################
# Plots

#--------------------------------------------
#-------------------------------------------- 

library(ggplot2)
library(reshape2)

Plot_EW_BS<-ggplot()+
  theme_classic()+
  geom_line(data = Performance_PDR,aes(x=PT,y=Ew_PDR,color='PrFS'),size=1.1)+
  geom_hline(aes(yintercept=Op_Performance_PDR$Ew_PF,color='PFS'),size=1.1)+
  geom_hline(aes(yintercept=Performance_DFDR$Ew_DF,color='DFS'),size=1.1)+
  geom_point(data = Op_Performance_PDR,aes(x=PT,y=Ew_PrF,shape='PT*'),size=3,color='red')+
  scale_color_manual("Scenarios",values = c("DFS"="dodgerblue4","PFS"="black","PrFS"="gray"))+
  scale_shape_manual("Optimal PT",values = c("PT*"=3))+
  labs(y=expression(E[w]~"[%]"))+
  theme(text=element_text(size=8, family="Times New Roman"))


CT_Re<-data.frame(POD=c(Performance_DFDR$POD_houses,Op_Performance_PDR$POD_houses),
                  FAR=c(Performance_DFDR$FAR_houses,Op_Performance_PDR$FAR_houses),
                  Scenario=c('DFS','PrFS'))

CT_Re<-melt(data =CT_Re,
            id.vars = c("Scenario"),
            variable.name = "Skill_score",
            value.name = "Value")

Plot_reliability_BS<-ggplot(data = CT_Re, aes(x=Skill_score, y=Value, fill=Scenario))+
  theme_classic()+
  geom_bar(stat="identity",position=position_dodge())+
  scale_fill_manual('',values=c("dodgerblue4","gray"), labels=c("DFS"="DFDR","PrFS"="PDR"))+
  theme(text=element_text(size=8, family="Times New Roman"))

# Export figures
library(ggpubr)
png(filename="Plot_EW_BS.png", 
    units="cm", 
    width=12, 
    height=7.5, 
    pointsize=8, 
    res=400)
ggpubr::ggarrange(Plot_EW_BS,ncol = 1,nrow = 1,
                  font.label = list(size = 10, color = "black", face = "bold", family = "Times New Roman"))
dev.off()


png(filename="Plot_Re_BS.png", 
    units="cm", 
    width=12, 
    height=7.5, 
    pointsize=8, 
    res=400)
ggpubr::ggarrange(Plot_reliability_BS,ncol = 1,nrow = 1,
                  font.label = list(size = 10, color = "black", face = "bold", family = "Times New Roman"),
                  common.legend = TRUE,legend="bottom")
dev.off()



#--------------------------------------------
#-------------------------------------------- 
# Export



save(Plot_EW_BS,Plot_reliability_BS,Performance_PDR,Op_Performance_PDR,Performance_DFDR,CT_Re,file = "Results_BS.Rdata")
