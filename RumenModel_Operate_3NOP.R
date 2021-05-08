##########################################################################################################################

# Set run time sequence
final <- 120
times_nf <- timer(final=final)


### Initialize Q(t=0) values and assign names

initial <- as.numeric(read.table("initial_rumen.txt",header=T))
names(initial) <- names(read.table("initial_rumen.txt",header=T))

initial <- c(initial[c("Q_Fd","Q_Sd","Q_Wr")],1e-5,initial[c("Q_He","Q_Mi","Q_Ac","Q_Pr","Q_Bu","Q_Hy","Q_Me","f_NADH")],0)
names(initial)[names(initial) %in% ""] <- c("Q_NP","Q_Mt")


### load feed composition and feed intake rate input data 

feedcomp	<- read.table("Feed_comp_VL2017.txt",header=T)
feedcomp[c('St','Sr','WSC')] <- c(feedcomp['St']+0.5*feedcomp['Sr'],Sr=0,feedcomp['WSC']+0.5*feedcomp['Sr']) # Equally divide soluble starch (Sr) over Sd and Wr

FI_rate_2m <- as.data.frame(read.table("FI_rate_2m.txt",header=T)) # Read feed intake rate throughout 0.5 day
FI_rate_2m[,'FI_rate'] <- feedcomp$c_DM*FI_rate_2m[,'FI_rate'] # Convert to fresh product to DM
FI_rate_2m <- rbind(FI_rate_2m,cbind(12+FI_rate_2m[,c('T_begin','T_eind')],FI_rate=FI_rate_2m[,'FI_rate'])) # Set feed intake rate for a whole day

feedcomp	<- feedcomp[!names(feedcomp) %in% c('c_DM','Sr')] # Remove dry matter and soluble starch contents

### Initialize constants and parameters and FME model code
### Initial values for parameters to be estimated 

#NP_par <- c(Parms_7,J_NP_HyMe=2.05e-5)
NP_par <- c(Parms_7,J_NP_HyMe=1.93e-5)
source("RumenModel_3NOP_Code.R")


# Run model with three different doses of 3NOP

out_3NOP <- list(out_00 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_3NOP=0),Parms=unlist(NP_par),Constants=Constants,
                                       stpsz=delt,Yields=Yields,feedintake=FI_rate_2m),
                 out_05 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_3NOP=0.495e-3),Parms=unlist(NP_par),Constants=Constants,
                                       stpsz=delt,Yields=Yields,feedintake=FI_rate_2m), # 60 mg 3NOP/kg DM
                 out_10 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_3NOP=0.991e-3),Parms=unlist(NP_par),Constants=Constants,
                                       stpsz=delt,Yields=Yields,feedintake=FI_rate_2m)) # 120 mg 3NOP/kg DM

cat('Methane emissions is',
    (max(out_3NOP[[3]]$Q_Mt) - out_3NOP[[3]]$Q_Mt[out_3NOP[[3]]$time==final-24])*16,',',
    (max(out_3NOP[[2]]$Q_Mt) - out_3NOP[[2]]$Q_Mt[out_3NOP[[2]]$time==final-24])*16,'and',
    (max(out_3NOP[[1]]$Q_Mt) - out_3NOP[[1]]$Q_Mt[out_3NOP[[1]]$time==final-24])*16,
    'g/d after 0, 0.5 and 1.0 mmol/kg DMI of 3NOP supplementation, respectively.')

pdf("3NOP_vs_Time.pdf",height=7.25,width=7.25)
#tiff("3NOP_vs_Time.tiff",height=7.25,width=7.25,units='in',res=600,compression='lzw')
par(mfrow=c(3,3),mar=c(2,5,0,0),oma=c(2,0.1,0.5,0.5))
plot_Dyn_DSM(data=out_3NOP,var_y='C_NP',ylim=c(0,0.10e-3),ylab=expression(italic(C)['3NOP']~' (mM)'),lgrthm=F)
axis(side=2,at=c(0,0.5e-4,1e-4,1.5e-4,2e-4,2.5e-4),labels=c('0.00',0.05,'0.10',0.15,'0.20','0.25'),las=2,cex.axis=1.1)
legend("topleft",c(expression(paste('0.0',~mmol~'(kg DMI)'^{-1})),expression(paste('0.5',~mmol~'(kg DMI)'^{-1})),expression(paste('1.0',~mmol~'(kg DMI)'^{-1}))),
       cex=1.05,col=c("black","magenta","limegreen"),lty=1,lwd=2)
plot_Dyn_DSM(data=out_3NOP,var_y='C_VFA',ylim=c(0.05,0.130),ylab=expression(italic(C)[VFA]~ '(mM)'),lgrthm=F)
axis(side=2,at=c(0.05,0.07,0.09,0.11,0.13),labels=c(50,70,90,110,130),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_3NOP,var_y='P_Hy',ylim=c(1e-4,0.3),ylab=expression(italic(p)[H[2]]~' (atm)'),lgrthm=T)
axis(side=2,at=c(1e-4,1e-3,1e-2,1e-1),labels=c('1e-4','1e-3','1e-2','1e-1'),las=2,cex.axis=1.1) 
plot_Dyn_DSM(data=out_3NOP,var_y='P_Mt_HyMe',ylim=c(0,2),ylab=expression(CH[4]~'emission rate'~paste('(',mol~h^{-1},')')),lgrthm=F)
axis(side=2,at=c(0,0.5,1,1.5,2),labels=c('0.0',0.5,'1.0',1.5,'2.0'),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_3NOP,var_y='F_T_Fd',ylim=c(-12,1),ylab=expression(italic(F)['T']),lgrthm=F)
axis(side=2,at=c(-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1),labels=c(-11,'','',-8,'','',-5,'','',-2,'','',1),las=2,cex.axis=1.2)
plot_Dyn_DSM(data=out_3NOP,var_y='r_NAD',ylim=c(0,4),ylab=expression(italic(r)[NAD]),lgrthm=F)
axis(side=2,at=c(0,1,2,3,4),labels=c(0,1,2,3,4),las=2,cex.axis=1.2)
plot_Dyn_DSM(data=out_3NOP,var_y='p_Ac',ylim=c(0.4,0.8),ylab=expression(Ac^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2)
axis(side=2,at=c(0.4,0.5,0.6,0.7,0.8),labels=c(0.4,0.5,0.6,0.7,0.8),las=2,cex.axis=1.2) 
plot_Dyn_DSM(data=out_3NOP,var_y='p_Pr',ylim=c(0.1,0.4),ylab=expression(Pr^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2)
axis(side=2,at=c(0.1,0.2,0.3,0.4),labels=c(0.1,0.2,0.3,0.4),las=2,cex.axis=1.2)
plot_Dyn_DSM(data=out_3NOP,var_y='p_Bu',ylim=c(0.11,0.17),ylab=expression(Bu^{'-'}~'proportion (-)'),lgrthm=F)
axis(side=1,at=final-c(24,18,12,6,0),labels=c(0,6,12,18,24),las=1,cex.axis=1.2)
axis(side=2,at=c(0.11,0.13,0.15,0.17),labels=c(0.11,0.13,0.15,0.17),las=2,cex.axis=1.15) 
mtext(outer=T, side=1, line=0.75, 'Time (h)', cex=1.2)
dev.off()


#########################################
### Simulations using literature data ###
#########################################

final <- 240
times_nf <- timer(final=final)

# initialize Q_i(t=0)
initial <- read.table("initial_3NOP.txt",header=T)
initial$Q_Mt <- rep(0,nrow(initial)) # Set Q_Mt at zero

# initialize feed composition
feedcomp <- read.table("Feed_comp_3NOP.csv",header=T,sep=',')
feedcomp[,c('St','Sr','WSC')] <- cbind(feedcomp[,'St']+0.5*feedcomp[,'Sr'],Sr=0,feedcomp[,'WSC']+0.5*feedcomp[,'Sr']) # Equally divide soluble starch (Sr) over degradable starch (St) and soluble sugars (WSC)
feedcomp <- feedcomp[,-c(match(c('c_DM','Sr'),colnames(feedcomp)))]


# initialize feed intake rate profile
FI_rate_2m <- lapply(1:nrow(feedcomp), function(x) mean(x)) # Make list of length nrow(feedcomp)

FI_rate_2m[[1]] <- as.data.frame(read.table("FI_rate_3NOP.csv",header=T,sep=',')) # Assign feed intake rate profile to list[[1]]

for (i in nrow(feedcomp):1){FI_rate_2m[[i]] <- cbind(FI_rate_2m[[1]][,c('T_begin','T_eind')], # Compute FI rate based on DMI and assign to list[[i]]
                                                     FI_rate=feedcomp$DMI[i]*FI_rate_2m[[1]][,'FI_rate'])} 

FI_rate_2m_Ctrl <- FI_rate_2m[feedcomp$c_3NOP == 0]
FI_rate_2m <- FI_rate_2m[feedcomp$c_3NOP != 0]

feedcomp_Ctrl <- subset(feedcomp, c_3NOP==0)
feedcomp <- subset(feedcomp, c_3NOP != 0)

##################################################################################################################################################

out_3NOP_Ctrl <- NULL
for (i in 1:nrow(feedcomp_Ctrl)){
  Qt_0 <- as.numeric(initial[i,])
  names(Qt_0) <- c("Q_Fd","Q_Sd","Q_Wr","Q_NP","Q_He","Q_Mi","Q_Ac","Q_Pr","Q_Bu","Q_Hy","Q_Me","f_NADH",'Q_Mt')
  out_3NOP_Ctrl[[length(out_3NOP_Ctrl)+1]] <- NAD_confurc(times=times_nf,initial=Qt_0,feed=feedcomp_Ctrl[i,],
                                                          Parms=unlist(c(Parms_7,NP_par)),
                                                          Constants=Constants,Yields=Yields,feedintake=FI_rate_2m_Ctrl[[i]],stpsz=delt)
} # End of for loop


CH4_Ctrl_P <- NULL
for (j in 1:length(out_3NOP_Ctrl)) {
  CH4_Ctrl_P <- c(CH4_Ctrl_P ,16*c((max(out_3NOP_Ctrl[[j]]$Q_Mt) - out_3NOP_Ctrl[[j]]$Q_Mt[out_3NOP_Ctrl[[j]]$time==final-24])))}


cat('Predicted control methane emissions are',
    CH4_Ctrl_P[1:4],'and',
    CH4_Ctrl_P[5],
    'g/d, respectively.') # 490, 430, 551, 638 and 483 g/d

# Van Wesemael 2019, Haisan 2017, Lopes 2016, Hristov 2015, Haisan 2014
CH4_Ctrl_O <- c(525,378,487,481,372) # Observed values per 24 h


cat('Observed control methane emissions are: ',CH4_Ctrl_O,'g/d, respectively')

# Van Wesemael 2019 (2x), Haisan 2017 (2x), Lopes 2016 (1x), Hristov 2015 (3x), Haisan 2014 (1x)
obsdata <- c(380,403,275,231,335,363,333,319,132) # Observed values per 24 h


obsdata <- obsdata+
  c((CH4_Ctrl_P-CH4_Ctrl_O)[1],
    (CH4_Ctrl_P-CH4_Ctrl_O)[1],
    (CH4_Ctrl_P-CH4_Ctrl_O)[2],
    (CH4_Ctrl_P-CH4_Ctrl_O)[2],
    (CH4_Ctrl_P-CH4_Ctrl_O)[3],
    (CH4_Ctrl_P-CH4_Ctrl_O)[4],
    (CH4_Ctrl_P-CH4_Ctrl_O)[4],
    (CH4_Ctrl_P-CH4_Ctrl_O)[4],
    (CH4_Ctrl_P-CH4_Ctrl_O)[5]) # All

### Costfunc daily average CH_4


Costfunc <- function(parms, Cnstnts, feedcomp, obsdata, times_nf, delt=delt, verbose=TRUE) { 
  
  #parms <- unlist(c(Parms_7,J_NP_HyMe=1e-4)); Cnstnts <- Constants
  print("Parameter value is:"); print(parms)
  
  initial <- read.table("initial_3NOP.txt",header=T)
  initial <- cbind(initial[,c(1:12)],Q_Mt=rep(0,nrow(initial))) # Set Q_Mt at zero
  
  out <- NULL; tailrows <- NULL
  for (i in 1:nrow(feedcomp)){
    Qt_0 <- as.numeric(initial[i,])
    names(Qt_0) <- colnames(initial); #print(Qt_0)
    out[[length(out)+1]] <- NAD_confurc(times=times_nf,initial=Qt_0,feed=feedcomp[i,],Parms=parms,Constants=Cnstnts,
                                        Yields=Yields,
                                        feedintake=FI_rate_2m[[i]],stpsz=delt,output='FINAL_24h')
    
    plot(P_Mt_HyMe~time,data=out[[i]],type='l',lwd=2,col='dodgerblue',ylim=c(0,2.5)) # print(tail(out[[i]]$Q_Wr))
    tailrows <- rbind(tailrows,out[[i]][nrow(out[[i]]),2:14]) # Extract final row per treatment and state variable columns
  } # End of for loop
  
  print(tailrows); # write.table(tailrows,"initial_3NOP.txt",sep="\t",row.names=F)
  
  if (length(obsdata) != length(out)) {stop('Lengths of obsdata and out are different')} # Check is obsdata and out have same length
  
  for (q in 1:length(obsdata)){
    out[[q]] <- 16*diff(out[[q]][which(match(out[[q]]$time, c(final,final-24)) >= 0),"Q_Mt"]) # CH_4 production g/d
  } # End of for loop q; CH_4 production in g/d
  
  cost <- modCost(model = cbind(time=1:length(out),CH4=unlist(out)), 
                  obs = cbind(time=1:length(obsdata),CH4=obsdata), weight="mean")
  write.table(c(parms,cost$model),'par_SSR.csv',sep=',',append=T, row.names=F, col.names=F)
  cat('SSR is:',cost$model,'\n')
  return(cost) 
  
} # End function Costfunc

########################################################################################################


parFit_3NOP <- modFit(p=c(J_NP_HyMe=2e-5),lower=0.5*c(J_NP_HyMe=2e-5),upper=2*c(J_NP_HyMe=2e-5),
                      f=Costfunc,Cnstnts=c(Constants,Parms_7),feedcomp=feedcomp,times_nf=times_nf,
                      obsdata=obsdata,method="bobyqa",delt=delt) # Faster than BFGS, needs about 2/3 of time

print(summary(parFit_3NOP))
print(parFit_3NOP$convergence) # '0' indicates successful convergence. '1' indicates function evaluation count 'max-eval' was reached
print(parFit_3NOP$par)
write.table(summary(parFit_3NOP)$par,"pars_3NOP.csv",sep=",")
