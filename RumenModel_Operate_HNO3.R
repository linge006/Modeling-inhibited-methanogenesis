##########################################################################################################################

# Set run time sequence
times_nf <- timer(final=final,delt=2e-3, mult_delt=1)

### Initialize Q(t=0) values and assign names

initial <- as.numeric(read.table("initial_rumen.txt",header=T))
names(initial) <- names(read.table("initial_rumen.txt",header=T))

initial <- c(initial[c("Q_Fd","Q_Sd","Q_Wr")],1e-5,initial[c("Q_He","Q_Mi","Q_Ac","Q_Pr","Q_Bu","Q_Hy","Q_Me","f_NADH")],rep(0,5))

more_names <- c("Q_Na","Q_Mt","f_NaIn","f_NaAb","f_NaEx","f_NaAm")
names(initial)[names(initial) %in% ""] <- more_names


### Load feed composition and feed intake rate input data 
feedcomp	<- read.table("Feed_comp_VL2017.txt",header=T)
feedcomp[c('St','Sr','WSC')] <- c(feedcomp['St']+0.5*feedcomp['Sr'],Sr=0,feedcomp['WSC']+0.5*feedcomp['Sr']) # Equally divide soluble starch (Sr) over Sd and Wr

FI_rate_2m <- as.data.frame(read.table("FI_rate_2m.txt",header=T)) # Read feed intake rate throughout 0.5 day
FI_rate_2m[,'FI_rate'] <- feedcomp$c_DM*FI_rate_2m[,'FI_rate'] # Convert to fresh product to DM
FI_rate_2m <- rbind(FI_rate_2m,cbind(12+FI_rate_2m[,c('T_begin','T_eind')],FI_rate=FI_rate_2m[,'FI_rate'])) # Set feed intake rate for a whole day

feedcomp <- feedcomp[!names(feedcomp) %in% c('c_DM','Sr')] # Remove dry matter and soluble starch contents


### Initialize constants and parameters and FME model functions/code...

source("RumenModel_HNO3_Code.R")
source("RumenModel_ConstPar_DSM.R")

Na_par <- c(Parms_7,k_NaAm=6.987)

# Run model with three different doses of HNO3
out_HNO3 <- list(out_000 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_HNO3=0.00),Parms=Na_par,Constants=Constants,
                                        stpsz=delt,Yields=Yields,feedintake=FI_rate_2m),
                 out_010 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_HNO3=0.16),Parms=Na_par,Constants=Constants,
                                        stpsz=delt,Yields=Yields,feedintake=FI_rate_2m), # 10 g HNO3/kg DM
                 out_020 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_HNO3=0.32),Parms=Na_par,Constants=Constants,
                                        stpsz=delt,Yields=Yields,feedintake=FI_rate_2m)) # 20 g HNO3/kg DM

cat('Methane emissions is',
    (max(out_HNO3[[3]]$Q_Mt) - out_HNO3[[3]]$Q_Mt[out_HNO3[[3]]$time==final-24])*16,',',
    (max(out_HNO3[[2]]$Q_Mt) - out_HNO3[[2]]$Q_Mt[out_HNO3[[2]]$time==final-24])*16,'and',
    (max(out_HNO3[[1]]$Q_Mt) - out_HNO3[[1]]$Q_Mt[out_HNO3[[1]]$time==final-24])*16,
    'g/d after 0, 0.16 and 0.32 mol/kg DMI of Nitrate supplementation, respectively.')

pdf("HNO3_vs_Time.pdf",height=7.25,width=7.25)
#tiff("HNO3_vs_Time.tiff",height=7.25,width=7.25,units='in',res=600,compression='lzw')
par(mfrow=c(3,3),mar=c(2,5,0,0),oma=c(2,0.1,0.5,0.5))
plot_Dyn_DSM(data=out_HNO3,var_y='C_Na',ylim=c(0,15e-3),ylab=expression(italic(C)[NO[3]^{'-'}]~' (mM)'),lgrthm=F)
axis(side=2,at=c(0,5e-3,10e-3,15e-3,20e-3,25e-3),labels=c('0',5,'10',15,'20','25'),las=2,cex.axis=1.1) 
legend("topleft",c(expression(paste('0.00',~mol~'(kg DMI)'^{-1})),expression(paste('0.16',~mol~'(kg DMI)'^{-1})),expression(paste('0.32',~mol~'(kg DMI)'^{-1}))),
       cex=1.05,col=c("black","magenta","limegreen"),lty=1,lwd=2)
plot_Dyn_DSM(data=out_HNO3,var_y='C_VFA',ylim=c(0.05,0.130),ylab=expression(italic(C)[VFA]~ '(mM)'),lgrthm=F)
axis(side=2,at=c(0.05,0.07,0.09,0.11,0.13),labels=c(50,70,90,110,130),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_HNO3,var_y='P_Hy',ylim=c(1e-4,0.3),ylab=expression(italic(p)[H[2]]~' (atm)'),lgrthm=T)
axis(side=2,at=c(1e-4,1e-3,1e-2,1e-1),labels=c('1e-4','1e-3','1e-2','1e-1'),las=2,cex.axis=1.1) 
plot_Dyn_DSM(data=out_HNO3,var_y='P_Mt_HyMe',ylim=c(0,2),ylab=expression(CH[4]~'emission rate'~paste('(',mol~h^{-1},')')),lgrthm=F)
axis(side=2,at=c(0,0.5,1,1.5,2),labels=c('0.0',0.5,'1.0',1.5,'2.0'),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_HNO3,var_y='F_T_Fd',ylim=c(-3,1),ylab=expression(italic(F)['T']),lgrthm=F)
axis(side=2,at=c(-3,-2,-1,0,1),labels=c(-3,-2,-1,0,1),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_HNO3,var_y='r_NAD',ylim=c(0,4),ylab=expression(italic(r)[NAD]),lgrthm=F)
axis(side=2,at=c(0,1,2,3,4),labels=c(0,1,2,3,4),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_HNO3,var_y='p_Ac',ylim=c(0.4,0.8),ylab=expression(Ac^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2) 
axis(side=2,at=c(0.4,0.5,0.6,0.7,0.8),labels=c(0.4,0.5,0.6,0.7,0.8),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_HNO3,var_y='p_Pr',ylim=c(0.1,0.4),ylab=expression(Pr^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2)
axis(side=2,at=c(0.1,0.2,0.3,0.4),labels=c(0.1,0.2,0.3,0.4),las=2,cex.axis=1.2)
plot_Dyn_DSM(data=out_HNO3,var_y='p_Bu',ylim=c(0.11,0.17),ylab=expression(Bu^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2)
axis(side=2,at=c(0.11,0.13,0.15,0.17),labels=c(0.11,0.13,0.15,0.17),las=2,cex.axis=1.15) 
mtext(outer=T,side=1, line=0.75, 'Time (h)',cex=1.2)
dev.off()


#########################################
### Simulations using literature data ###
#########################################

times_nf <- timer(final=240,delt=2e-3, mult_delt=1)

# initialize Q_i(t=0)
initial <- read.table("initial_HNO3.txt",header=T)
initial[,more_names] <- matrix(0,nrow=nrow(initial),ncol=length(more_names)) # Set Q_zeros at zero !!! Should more_names be used here, or fewer names?


# initialize feed composition
feedcomp <- read.table("Feed_comp_HNO32.csv",header=T,sep=',')
feedcomp[,c('St','Sr','WSC')] <- cbind(feedcomp[,'St']+0.5*feedcomp[,'Sr'],Sr=0,feedcomp[,'WSC']+0.5*feedcomp[,'Sr']) # Equally divide soluble starch (Sr) over degradable starch (St) and soluble sugars (WSC)
feedcomp <- feedcomp[,-c(match(c('c_DM','Sr'),colnames(feedcomp)))]


# initialize feed intake rate profile
FI_rate_2m <- read.table("FI_rate_Olijhoek.csv",header=T,sep=',') # Load feed intake rate
FI_rate_2m <- list(FI_rate_2m[,c("T_begin","T_eind","FI_rate")],
                   FI_rate_2m[,c("T_begin","T_eind","FI_rate_med")],
                   FI_rate_2m[,c("T_begin","T_eind","FI_rate_high")],
                   as.data.frame(read.table("FI_rate_SvZ.csv",header=T,sep=',')),
                   as.data.frame(read.table("FI_rate_Veneman.csv",header=T,sep=',')))


# Rename FI_rate_xxx to FI_rate
for (j in 2:3) colnames(FI_rate_2m[[j]])[3] <- 'FI_rate'


### Model initial runs after setting time for five runs (3 times Olijhoek, 1 VZ, 1 Veneman) !!! Why do I need this?

times_nf <- list(seq(0,final,5e-3),seq(0,final,5e-3),seq(0,final,5e-3),
                 seq(0,final,2e-3),seq(0,final,2e-3))

##################################################################################################################################################


### initialize Control feed composition
feedcomp_Ctrl <- read.table("Feed_comp_HNO32_cntrls.csv",header=T,sep=',')
feedcomp_Ctrl[,c('St','Sr','WSC')] <- cbind(feedcomp_Ctrl[,'St']+0.5*feedcomp_Ctrl[,'Sr'],Sr=0,feedcomp_Ctrl[,'WSC']+0.5*feedcomp_Ctrl[,'Sr']) # Equally divide soluble starch (Sr) over degradable starch (St) and soluble sugars (WSC)
feedcomp_Ctrl <- feedcomp_Ctrl[,!names(feedcomp_Ctrl) %in% c('c_DM','Sr')]

# initialize feed intake rate profile
FI_rate_2m_Ctrl <- list(as.data.frame(read.table("FI_rate_Olijhoek_cntrl.csv",header=T,sep=',')), # Load feed intake rate,
                        as.data.frame(read.table("FI_rate_SvZ.csv",header=T,sep=',')),
                        as.data.frame(read.table("FI_rate_Veneman_cntrl.csv",header=T,sep=',')))

colnames(FI_rate_2m_Ctrl[[1]])[3] <- 'FI_rate' # 

# Running Ctrl treatments
out_HNO3_Ctrl <- NULL
for (i in 1:nrow(feedcomp_Ctrl)){
  
  Qt_0 <- as.numeric(initial[i,]); names(Qt_0) <- colnames(initial)
  
  out_HNO3_Ctrl[[length(out_HNO3_Ctrl)+1]] <- NAD_confurc(times=times_nf[[i]],initial=Qt_0,feed=feedcomp_Ctrl[i,],
                                                            Parms=c(Parms_7,Na_par),Constants=Constants,
                                                            Yields=Yields,feedintake=FI_rate_2m_Ctrl[[i]],stpsz=delt)
} # End of for loop

# Compute methane production (g/d) for final day that was simulated
CH4_Ctrl_P <- 16*c((max(out_HNO3_Ctrl[[1]]$Q_Mt) - out_HNO3_Ctrl[[1]]$Q_Mt[out_HNO3_Ctrl[[1]]$time==final-24]),
                   (max(out_HNO3_Ctrl[[2]]$Q_Mt) - out_HNO3_Ctrl[[2]]$Q_Mt[out_HNO3_Ctrl[[2]]$time==final-24]),
                   (max(out_HNO3_Ctrl[[3]]$Q_Mt) - out_HNO3_Ctrl[[3]]$Q_Mt[out_HNO3_Ctrl[[3]]$time==final-24]))

CH4_Ctrl_O <- c(341,368,394) # Observed values per 24 h (From Olijhoek 2016, Van Zijderveld 2011 en Veneman 2011)

cat('Observed values are: ',CH4_Ctrl_O[1],',', CH4_Ctrl_O[2],'and',CH4_Ctrl_O[3],'g/d, respectively')

#######################################################################################################################


# Nitrate observed means... Required for HNO3 model
obsdata <- c(0.71*c(468,424,396),mean(c(283,313,318,326)),317)* # These Ctrl treatment methane values are from the three papers
  c((CH4_Ctrl_P/CH4_Ctrl_O)[1],
    (CH4_Ctrl_P/CH4_Ctrl_O)[1],
    (CH4_Ctrl_P/CH4_Ctrl_O)[1],
    (CH4_Ctrl_P/CH4_Ctrl_O)[2],
    (CH4_Ctrl_P/CH4_Ctrl_O)[3])


### Costfunc daily average CH_4

Costfunc <- function(parms, Cnstnts, feedcomp, obsdata, times_nf, delt=delt, verbose=TRUE) { 
  
  print("Parameter value is:"); print(parms)
  
  initial <- read.table("initial_HNO3.txt",header=T)
  initial[,more_names] <- matrix(0,nrow=nrow(initial),ncol=length(more_names)) # Set Q_zeros at zero
  
  out <- NULL; tailrows <- NULL
  for (i in 1:nrow(feedcomp)){
    Qt_0 <- as.numeric(initial[i,])
    names(Qt_0) <- colnames(initial); #print(Qt_0)
    out[[length(out)+1]] <- NAD_confurc(times=times_nf[[i]],initial=Qt_0,feed=feedcomp[i,],Parms=parms,Constants=Cnstnts,
                                        Yields=Yields,
                                        feedintake=FI_rate_2m[[i]],stpsz=delt,output='FINAL_24h')
    
    plot(P_Mt_HyMe~time,data=out[[i]],type='l',lwd=2,col='dodgerblue') # print(tail(out[[i]]$Q_Wr))
    tailrows <- rbind(tailrows,out[[i]][nrow(out[[i]]),2:14]) # Extract final row per treatment and state variable columns
  } # End of for loop
  
  print(tailrows); #write.table(tailrows,"initial_HNO3.txt",sep="\t",row.names=F)
  
  if (length(obsdata) != length(out)) {stop('Lengths of obsdata and out are different')} # Check is obsdata and out have same length
  
  for (q in 1:length(obsdata)){
    out[[q]] <- 16*diff(out[[q]][which(match(out[[q]]$time, c(final,final-24)) >= 0),"Q_Mt"]) # CH_4 production g/d
  } # End of for loop q; CH_4 production in g/d
  
  cost <- modCost(model = cbind(time=1:length(out),CH4=unlist(out)), 
                  obs = cbind(time=1:length(obsdata),CH4=obsdata), weight="mean")
  cat('SSR is:',cost$model,'\n')
  return(cost) 
  
} # End function Costfunc


# Estimate par
parFit_HNO3q <- modFit(p=c(k_NaAm=7.5),lower=0.1*c(k_NaAm=7),upper=10*c(k_NaAm=7),
                      f=Costfunc,Cnstnts=c(Constants,Parms_7),feedcomp=feedcomp,times_nf=times_nf,
                      obsdata=obsdata,method="bobyqa",delt=delt) # no reltol for bobyqa
print(parFit_HNO3q$par); print(parFit_HNO3q$par); print(parFit_HNO3q$convergence) # '0' indicates successful convergence, '1' count 'max-eval' was reached
write.table(summary(parFit_HNO3q)$par,"pars_HNO3.csv",sep=",")
