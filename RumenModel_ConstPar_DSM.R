
Constants <- c(V_headspace=40,V_mol=25,C_starMi=12.5,k_H=1382,f_b=703,M_He_PsMi=2.24e-2,
               M_He_AmMi=1.35e-3,J_Am_HeVf=8.61e-3,J_Ps_HeVf=1.465e-2,
               v_AcAb=0.07271622,v_PrAb=0.22636199,v_BuAb=1.01266850,k_NPAb=0.3,
               M_Ac_AcAb=0.0791,M_Pr_PrAb=0.112,M_Bu_BuAb=0.4934,
               M_He_HeMi=0.02,M_He_HeVa=0.055,q_GM=3,M_NAD_Ac=9,J_NAD_AP=1,
               k_NADH_Fe=202,k_NaAb=0.3)

DMI <- sum((FI_rate_2m[,2] - FI_rate_2m[,1]) * FI_rate_2m[,3]) # DMI for DM 24 h DMI (no total product!) input 

Yields <- c(Y_Ac_HeAc=2,Y_Bu_HeBu=1,Y_Ac_HeAP=2/3,Y_Pr_HeAP=4/3,Y_AcMi_HeMi=84.25,
            Y_BuMi_HeMi=68.95,Y_APMi_HeMi=77.86,Y_He_LaHe=2.78e-3,Y_Hy_PsMi=0.58e-3,X_Hy_AmMi=0.41e-3,
            Y_Hy_HeAc=4,Y_Hy_HeBu=2,Y_Hy_NaAm=4,
            Y_Hy_NaNi=1,Y_Hy_NiAm=3,
            c_PdMi=0.55,c_NAD=7e-6,Y_Mt_HyMe=0.25,Y_Me_HyMt=2,Y_Mt_FmLi=2e-3,
            f_He_HeAc=0.645,f_He_HeBu=0.706,f_He_HeAP=0.67,
            f_SdSi=0.1,f_PdSi=0.25)
