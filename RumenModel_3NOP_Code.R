#############################################################################################
### This code represents a rumen model in which thermodynamic driving force of hydrogen  	###
### on the microbial metabolism, more specifically the VFA formation pathways, is 		###
### incorporated; the model is used for prediction of methane emission from cows and 	###
### require feed intake and composition as input. 							###
#############################################################################################

NAD_confurc <- function(times, initial, feed, Parms, Constants, Yields, feedintake, stpsz, output="CLEAN") {

# times is a 'periodic' time sequence with periods of 12 h;
# intial is a sequence containing all Q(t=0)
# feed is a sequence containing feed composition and hydrolysis rates of carbohydrates
# Parms is a sequence containing all model parameters that are potentially optimized
# Constants is a sequence containing all model constants
# Yield is a sequence containing all model yield factors
# feedintake is a dataframe containing the feed intake rate per period of 24 h
  
   Constants <- c(Constants, V_fluid=47.86+1.759*feed$DMI, k_FlEx=(-3.40+1.224*DMI-0.030*DMI^2+5.93*feed$c_For)/100, 
                 k_SoEx=(feed$c_For*(1.74+0.15*DMI)+(1-feed$c_For)*(10.1-0.96*DMI+0.037*DMI^2))/100,
                 k_LiEx=1/(-0.2*feed$DMI+13))

######################################################################################################################################

NAD_confurc_FME <- function(times,initial,feed,Parms,Constants,Yields,feedintake,delt) {

with(as.list(c(initial,Parms,feed, feedintake, Constants, Yields, delt)), 
	{ 

### Calculation of the concentrations

C_NP <- Q_NP/V_fluid									# 3-NOP oncentration	  
C_He <- Q_He/V_fluid									# Hexose concentration
C_Mi <- Q_Mi/V_fluid									# Concentration of microbes
C_Ps <- 3e-3 										# Soluble protein concentration
C_Am <- 5e-3 										# Ammonia concentration
C_Ac <- Q_Ac/V_fluid									# Acetate concentration
C_Pr <- Q_Pr/V_fluid									# Propionate concentration
C_Bu <- Q_Bu/V_fluid									# Butyrate concentration
C_VFA <- sum(C_Ac,C_Pr,C_Bu)								# total VFA concentration
p_VFA <- c(C_Ac,C_Pr,C_Bu)/sum(C_Ac,C_Pr,C_Bu)					# proportions of Ac, Pr and Bu of total VFA
P_Hy <- Q_Hy * (1e3*8.3145*312/1.01325e5) / V_headspace			# Hydrogen partial pressure in rumen headspace
C_Hy <- P_Hy / k_H #(Q_Hy * 8.3145 * 312) / (k_H * 1.01325e5 * V_fluid) # Dissolved H2 concentration based on ideal gas law and Henry's law 
pH = 7.73 - (0.014*1e3*sum(C_Ac,C_Pr,C_Bu))					# rumen fluid pH
pH_cell <- 6.43+3.62e-8*exp(2.4*pH)								# intracellular pH calculated based on extracellulair/fluid pH
r_NAD <- (1-f_NADH)/f_NADH									# NAD+ to NADH ratio

### Calculation of flows and states

# when first line and second line below have a TRUE match, 
# then select FI_rate of the 14 choices in FI_rate_dynmodel.txt,
# calculate DMI intake rate multiplying by dry matter content, c_DM.

D_DMI <- feedintake[(feedintake[,"T_begin"] <= (times %% 24)) & 
				((times %% 24) < feedintake[,"T_eind"]),
				"FI_rate"] # dry matter intake rate [kg/h]


############
# Q_Fd (g) #
############

P_Fd_InFd <- D_DMI*NDF # Inflow from feed [kg/h] * [g/kg] = [g/h]

U_Fd_FdHe <- k_FdHe*C_Mi/C_starMi*Q_Fd # Hydrolysis to hexose
U_Fd_FdEx <- k_SoEx*Q_Fd # Fractional outflow

############
# Q_Sd (g) #
############

P_Sd_InSd <- D_DMI*St # Inflow from feed

U_Sd_SdHe <- k_SdHe*C_Mi/C_starMi*Q_Sd # Hydrolysis to hexose
U_Sd_SdEx <- k_SoEx*Q_Sd # Fractional outflow

############
# Q_Wr (g) #
############

P_Wr_InWr <- D_DMI*WSC # Inflow from feed

U_Wr_WrHe <- k_WrHe*C_Mi/C_starMi*Q_Wr # Oliggosaccharide hydrolysis
U_Wr_WrEx <- k_FlEx*Q_Wr # Fractional outflow

##############
# Q_NP (mol) #
##############

P_NP_InNP <- D_DMI*c_3NOP # Inflow from feed

U_NP_NPAb <- k_NPAb*Q_NP # Absorption of 3-NOP
U_NP_NPEx <- k_FlEx*Q_NP # Fractional outflow


##############
# Q_He (mol) #
##############

P_He_LaHe <- Y_He_LaHe * D_DMI*Lain/89 # Lactate to be fermented in He equivalents

P_He_XdHe <- (U_Fd_FdHe+U_Sd_SdHe+U_Wr_WrHe)/162 # hexose yield from fiber and starch hydrolysis in the rumen [mol/h]

###

U_AcHe_HeMi <- (v_HeAc/q_GM) * Q_Mi / ((1+M_He_HeMi/C_He)*(1+M_NAD_Ac/r_NAD))	# He uptake for microbial growth yielding acetate 
U_BuHe_HeMi <- (v_HeBu/q_GM) * Q_Mi / (1+M_He_HeMi/C_He)				# He uptake for microbial growth yielding butyrate 
U_APHe_HeMi <- (v_HeAP/q_GM) * Q_Mi / ((1+M_He_HeMi/C_He)*(1+r_NAD/J_NAD_AP))	# He uptake for microbial growth yielding acetate and propionate 

U_He_HeAc <- v_HeAc * Q_Mi / ((1+M_He_HeVa/C_He)*(1+C_Am/J_Am_HeVf+C_Ps/J_Ps_HeVf)*(1+M_NAD_Ac/r_NAD))	# He uptake for microbial maintenance yielding acetate 
U_He_HeBu <- v_HeBu * Q_Mi / ((1+M_He_HeVa/C_He)*(1+C_Am/J_Am_HeVf+C_Ps/J_Ps_HeVf)) 				# He uptake for microbial maintenance yielding butyrate
U_He_HeAP <- v_HeAP * Q_Mi / ((1+M_He_HeVa/C_He)*(1+C_Am/J_Am_HeVf+C_Ps/J_Ps_HeVf)*(1+r_NAD/J_NAD_AP))	# He uptake for microbial maintenance yielding acetate and propionate

###

U_He_HeEx <- k_FlEx * Q_He # Fractional outflow of hexose with fluid fraction

############
# Q_Mi (g) #
############

P_Mi_HeMi <- Y_AcMi_HeMi*U_AcHe_HeMi + Y_BuMi_HeMi*U_BuHe_HeMi + Y_APMi_HeMi*U_APHe_HeMi # Microbial growth

U_Mi_MiEx <- (0.15*k_FlEx+0.65*k_SoEx) * Q_Mi # Fractional outflow of microbes

##############
# Q_Ac (mol) #
##############

P_Ac_in 	<- D_DMI*Acin/59 # Acetate intake from feed  
P_Ac_HeAc 	<- Y_Ac_HeAc * (f_He_HeAc*U_AcHe_HeMi+U_He_HeAc) + Y_Ac_HeAP * (f_He_HeAP*U_APHe_HeMi+U_He_HeAP) # Acetate yield from hexose fermentation

U_Ac_AcAb <- (v_VfAb*v_AcAb*Q_Ac*(V_fluid^0.75)) / ((1+(M_Ac_AcAb/C_Ac)^1.17)*(1+(pH/6.02)^3.91))	# Van Lingen 2017

U_Ac_AcEx <- k_FlEx * Q_Ac # Fractional outflow of acetate

k_AcAb <- U_Ac_AcAb/Q_Ac

##############
# Q_Pr (mol) #
##############

P_Pr_in 	<- D_DMI*Prin/73 # Propionate intake from feed  
P_Pr_HePr 	<- Y_Pr_HeAP * (f_He_HeAP*(U_APHe_HeMi) + U_He_HeAP) # Propionate yield from hexose fermentation

U_Pr_PrAb <- (v_VfAb*v_PrAb*Q_Pr*(V_fluid^0.75)) / ((1+(M_Pr_PrAb/C_Pr)^0.95)*(1+(pH/6.02)^4.61))	# Van Lingen 2017

U_Pr_PrEx 	<- k_FlEx * Q_Pr # Fractional outflow of propionate 

k_PrAb <- U_Pr_PrAb/Q_Pr

##############
# Q_Bu (mol) #
##############

P_Bu_in <- D_DMI*Buin/87 # Butyrate intake from feed  
P_Bu_HeBu <- Y_Bu_HeBu * (f_He_HeBu*U_BuHe_HeMi+U_He_HeBu) # Butyrate yield from hexose fermentation

U_Bu_BuAb <- (v_VfAb*v_BuAb*Q_Bu*(V_fluid^0.75)) / ((1+(M_Bu_BuAb/C_Bu)^0.99)*(1+(pH/6.02)^5.13))	# Van Lingen 2017

U_Bu_BuEx <- k_FlEx * Q_Bu # Fractional outflow of butyrate

k_BuAb <- U_Bu_BuAb/Q_Bu

##############
# Q_Hy (mol) #
##############

P_Hy_HeAc <- Y_Hy_HeAc * (f_He_HeAc*U_AcHe_HeMi+U_He_HeAc) # Hydrogen yield from microbial acetate production 
P_Hy_HeBu <- Y_Hy_HeBu * (f_He_HeBu*U_BuHe_HeMi+U_He_HeBu) # Hydrogen yield from microbial butyrate production 

U_Hy_HyMe <- (v_HyMe*Q_Me) / (1+M_Hy_HyMe/C_Hy+C_NP/J_NP_HyMe) # Methanogenesis 

U_Hy_HyEm <- (k_HyEm + f_b * V_mol / (V_headspace * k_H)) * Q_Hy 		# Eructation/Exhalation of hydrogen + release of hydrogen via blood/lungs
U_Hy_HyEx <- V_fluid * k_FlEx * Q_Hy * V_mol / (V_headspace * k_H)	# Passage to hindgut

############
# Q_Me (g) #
############

P_Me_HyMt <- Y_Me_HyMt * Y_Mt_HyMe * U_Hy_HyMe # Methanogenic growth on hydrogen

U_Me_MeEx <- (0.4*k_SoEx+0.4*k_FlEx) * Q_Me # Methanogen outflow

############
# Q_LI (g) #
############

# Fermentable OM that reaches the hindgut is calculated from fiber leaving the rumen, starch that is still not degraded in the small intestine,
# protein fraction of bacteria and archaea; the latter is taken twice to get a rough estimate of the intestinally degraded rumen undegradable protein.

P_FmLi <- U_Fd_FdEx * k_FdHe/(k_FdHe + k_LiEx) +
  U_Sd_SdEx * f_SdSi * k_SdHe/(k_SdHe + k_LiEx) +
  2*(U_Mi_MiEx+U_Me_MeEx)*c_PdMi * f_PdSi * k_PdPs/(k_PdPs + k_LiEx)

##############
# A_Mt (mol) #
##############

P_Mt_HyMe <- Y_Mt_HyMe * U_Hy_HyMe # Ruminal methane production rate
P_Mt_MtLi <-  P_FmLi * Y_Mt_FmLi # Hindgut methane production rate

###############
# Q_NAD (mol) #
###############

P_NADH_Ac <- 2*(f_He_HeAc*U_AcHe_HeMi+U_He_HeAc) # NAD^+ reduction to NADH associated with hexose fermentation to 2 acetate

U_NADH_AP <- 0.67*(f_He_HeAP*U_APHe_HeMi+U_He_HeAP) # NADH oxidation to NAD^+ associated with hexose fermentation to 4/3 propionate + 2/3 acetate

F_T_Fd <- 1 - (((9*P_Hy^2)/(9*(10^(-1*pH_cell))^3)))^(1/2) * exp( (-102e3)/(2*8.3145*312) ) # thermodynamic potential factor for hydrogenase catalyzed NADH oxidation via confurcation

if (F_T_Fd > 0) {U_NADH_Fe <- k_NADH_Fe * f_NADH * c_NAD * Q_Mi * F_T_Fd}
	else {U_NADH_Fe <- 0} # hydrogenase catalyzed NADH oxidation via confurcation; U_NADH_Fe is set zero for F_T_Fd < 0

###

# Calculation of all dQ/dt

dFd <- P_Fd_InFd - U_Fd_FdHe - U_Fd_FdEx
dSd <- P_Sd_InSd - U_Sd_SdHe - U_Sd_SdEx
dWr <- P_Wr_InWr - U_Wr_WrHe - U_Wr_WrEx
dNP <- P_NP_InNP - U_NP_NPAb - U_NP_NPEx

dHe <- P_He_LaHe + P_He_XdHe - U_AcHe_HeMi - U_APHe_HeMi - U_BuHe_HeMi - U_He_HeAc - U_He_HeAP - U_He_HeBu - U_He_HeEx
dMi <- P_Mi_HeMi - U_Mi_MiEx

dAc <- P_Ac_HeAc + P_Ac_in - U_Ac_AcAb - U_Ac_AcEx
dPr <- P_Pr_HePr + P_Pr_in - U_Pr_PrAb - U_Pr_PrEx
dBu <- P_Bu_HeBu + P_Bu_in - U_Bu_BuAb - U_Bu_BuEx

dNAD <- P_NADH_Ac - U_NADH_AP - U_NADH_Fe

dHy <- P_Hy_HeAc + P_Hy_HeBu - U_Hy_HyMe - U_Hy_HyEm - U_Hy_HyEx
dMe <- P_Me_HyMt - U_Me_MeEx
dMt <- P_Mt_HyMe + P_Mt_MtLi

f_NADH <- (f_NADH*c_NAD*Q_Mi + dNAD*delt)/(c_NAD*Q_Mi) # NADH as a fraction of NAD^+ + NADH

### list all variables to be returned

list(c(dFd,dSd,dWr,dNP,dHe,dMi,dAc,dPr,dBu,dHy,dMe,dNAD,dMt),c(D_DMI=D_DMI,P_Fd_InFd=P_Fd_InFd,U_Fd_FdHe=U_Fd_FdHe,U_Fd_FdEx=U_Fd_FdEx,
	P_Sd_InSd=P_Sd_InSd,U_Sd_SdHe=U_Sd_SdHe,U_Sd_SdEx=U_Sd_SdEx,P_Wr_InWr=P_Wr_InWr,U_Wr_WrHe=U_Wr_WrHe,U_Wr_WrEx=U_Wr_WrEx,
	P_NP_InNP=P_NP_InNP,U_NP_NPEx=U_NP_NPEx,C_NP=C_NP,
	C_He=C_He,p_Ac=p_VFA[1],p_Pr=p_VFA[2],p_Bu=p_VFA[3],P_Ac_in=P_Ac_in,U_Ac_AcAb=U_Ac_AcAb,U_Ac_AcEx=U_Ac_AcEx,
	P_Pr_in=P_Pr_in,U_Pr_PrAb=U_Pr_PrAb,U_Pr_PrEx=U_Pr_PrEx,
	P_Bu_in=P_Bu_in,U_Bu_BuAb=U_Bu_BuAb,U_Bu_BuEx=U_Bu_BuEx,
	P_He_LaHe=P_He_LaHe,P_He_XdHe=P_He_XdHe,U_AcHe_HeMi=U_AcHe_HeMi,U_APHe_HeMi=U_APHe_HeMi,U_BuHe_HeMi=U_BuHe_HeMi,
	U_He_HeAc=U_He_HeAc,U_He_HeAP=U_He_HeAP,U_He_HeBu=U_He_HeBu,U_He_HeEx=U_He_HeEx,
	P_Hy_HeAc=P_Hy_HeAc,P_Hy_HeBu=P_Hy_HeBu,U_Hy_HyMe=U_Hy_HyMe,U_Hy_HyEm=U_Hy_HyEm,U_Hy_HyEx=U_Hy_HyEx,
	P_Me_HyMt=P_Me_HyMt,U_Me_MeEx=U_Me_MeEx,P_Mi_HeMi=P_Mi_HeMi,U_Mi_MiEx=U_Mi_MiEx,
	P_NADH_Ac=P_NADH_Ac,U_NADH_AP=U_NADH_AP,U_NADH_Fe=U_NADH_Fe,
	C_VFA=C_VFA,C_Ac=C_Ac,C_Pr=C_Pr,C_Bu=C_Bu,
	P_Ac_HeAc=P_Ac_HeAc,P_Pr_HePr=P_Pr_HePr,P_Bu_HeBu=P_Bu_HeBu,k_AcAb=k_AcAb,k_PrAb=k_PrAb,k_BuAb=k_BuAb,pH=pH,
	pH_cell=pH_cell,F_T_Fd=F_T_Fd,C_Hy=C_Hy,P_Hy=P_Hy,r_NAD=r_NAD,
	P_NADH_Ac=P_NADH_Ac,U_NADH_AP=U_NADH_AP,U_NADH_Fe=U_NADH_Fe,U_Hy_HyMe=U_Hy_HyMe,U_Hy_HyEm=U_Hy_HyEm,P_Mt_HyMe=P_Mt_HyMe+P_Mt_MtLi,P_Mt_MtHg=P_Mt_MtLi,
	v_HeAP=v_HeAP,M_NAD_Ac=M_NAD_Ac)) 
		} ) } 

################################################################################################

out = as.data.frame(lsoda(y=initial,times=times,func=NAD_confurc_FME,feed=feed,feedintake=feedintake,parms=Parms,
                          Constants=Constants,Yields=Yields,delt=stpsz))
cat("Run at:",as.character.Date(Sys.time()),"\n") # Print date and time
if (output=="FINAL_24h") {
  out <- out[out$time %in% seq(0,final,10*stpsz),]
  return(out[out$time >= max(out$time) - 24,])} # return out
if (output=="CLEAN") {return(out[out$time %in% seq(0,final,10*stpsz),])} # return out while cleaning up the out file
}
