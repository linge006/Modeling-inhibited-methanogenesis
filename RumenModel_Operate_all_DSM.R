
library(FME)

##########################################################################################################################

# Set working directory 

setwd("~/UCDavis/RumenModel_H2/CodePaperDSM")

# The 'NAD_FME' function needs to be uploaded

### Set simulation time...

## Timer will set a run time sequence with non-equal step sizes (saves calculation time)
# Times_nf will contain a periodic time sequence with period of 12 h;
# p is the number of 12-h-periods, should be an integer;
# delt is stepsize from 0 to t_switch h of every period;
# mult_delt*delt is stepsize from t_switch to 12 h of every period;
# times_nf lasts for 'final' h, with 'final' a multiple of 12.

timer <- function(t_switch_1=0.5, t_switch_2=12, t_switch_3=12.5, final=72, delt_t=2.5e-3, mult_delt=2){
  
  times_nf <- NULL
  
  for (p in 1:(final/24)){ 
    times_nf <- c(times_nf,((p-1)*24+
                              c(seq(0,t_switch_1,delt_t),
                                seq(t_switch_1,t_switch_2-(mult_delt*delt_t),mult_delt*delt_t)[-1],
                                seq(t_switch_2,t_switch_3,delt_t),
                                seq(t_switch_3,24-(mult_delt*delt_t),mult_delt*delt_t)[-1])
    ))
  } # End of for loop
  return(c(times_nf,final))
} # End of 'timer' function

final <- 120
delt <- 2e-3 # time step to be used
times_nf <- timer(final=final,delt_t=delt)

### Initial values for parameters to be estimated 
Parms_7 <- c(k_HyEm=7.490259,v_HeAc=0.05277139,v_HeAP=0.007528171,v_HeBu=0.002435758,v_HyMe=0.2935437,M_Hy_HyMe=6.289578e-07,v_VfAb=0.8384857)

### DSM plot function

plot_Dyn_DSM <- function(data,var_y,ylim,lgrthm=F,ylab,xlim=c(final-24,final)){
  
  cols <- c("magenta","limegreen",'dodgerblue','orange','grey','red','brown','navy')
  
  if (lgrthm==F){
    plot(as.formula(paste(var_y,'~time')),data=data[[1]][data[[1]]$time >= xlim[1],],type='l',lwd=2,xlab='',ylab=ylab,
         cex.lab=1.5,cex.axis=1.5,cex.lab=1.5,ylim=ylim,xlim=xlim,axes=F)
    for (i in 2:length(data)){
      lines(as.formula(paste(var_y,'~time')),data=data[[i]][data[[i]]$time >= xlim[1],],type='l',lwd=2,col=cols[i-1],lty=1)
    } # End of FOR loop
    box()
  } # End of IF statement; plot with non-logarithmic y-axis
  
  if (lgrthm==T){
    plot(as.formula(paste(var_y,'~time')),data=data[[1]][data[[1]]$time >= xlim[1],],type='l',lwd=2,xlab='',ylab=ylab,
         cex.lab=1.5,cex.axis=1.5,cex.lab=1.5,ylim=ylim,xlim=xlim,log="y",axes=F)
    for (i in 2:length(data)){
      lines(as.formula(paste(var_y,'~time')),data=data[[i]][data[[i]]$time >= xlim[1],],type='l',lwd=2,col=cols[i-1])
    } # End of FOR loop
    box()
  } # End of IF statement; plot with logarithmic y-axis
  
} # End of plot_Dyn_DSM

###

source("RumenModel_Operate_HNO3.R")
source("RumenModel_Operate_HNO32.R")
source("RumenModel_Operate_3NOP.R")
source("RumenModel_Operate_3NOP_HNO32.R")
