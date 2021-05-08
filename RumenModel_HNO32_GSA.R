
GSA_constants <- unlist(c(Parms_7,Constants[which(names(Constants) != 'k_NaAb')]))
GSA_pars <- c(Na_par,Constants['k_NaAb'])


LH_par <- cbind(min=0.75*GSA_pars,max=1.25*GSA_pars)
n_smpl <- 1000

# Points used for GSA
t_GSA <- final - c(12,11.5,11,10,8,6,2)

Na_pars_LHS <- list(as.data.frame(Latinhyper(LH_par,num=n_smpl)),as.data.frame(matrix(ncol=length(t_GSA),nrow=n_smpl)))

# k_NiAm = 1.067163e-03 

for (i in 1:n_smpl){
  out_HNO32_GSA <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_HNO3=0.32),Parms=Na_pars_LHS[[1]][i,],
                         Constants=GSA_constants,
                         stpsz=delt,Yields=Yields,feedintake=FI_rate_2m,output="FINAL_24h")
  Na_pars_LHS[[2]][i,] <- out_HNO32_GSA$P_Mt_HyMe[out_HNO32_GSA$time %in% t_GSA]
} # End of for loop i


methane_times <- paste0('CH4_',c('00','05','10','20','40','60','100'))

colnames(Na_pars_LHS[[2]]) <- methane_times
Na_pars_LHS <- cbind(Na_pars_LHS[[1]],Na_pars_LHS[[2]])


write.table(Na_pars_LHS,paste0('GSA_HNO32_',n_smpl,'.csv'),sep=',',row.names=F)
colnames(Na_pars_LHS)

library(corrplot)

HNO32_cor <- cor(Na_pars_LHS)
write.table(HNO32_cor,paste0('GSA_HNO32_cor_',n_smpl,'.csv'),sep=',',row.names=F)

pdf(paste0('GSA_HNO32_',n_smpl,'.pdf'),width=8.5,height=11)
par(mar=c(0,1,4,1),mfrow=c(2,1),oma=c(1,1,1,3))
corrplot(HNO32_cor[1:4,methane_times], method = 'shade',
         addgrid.col = 'lightgrey',tl.col = "black",main='Nitrate model')
dev.off()

# 0, 0.5, 1, 2, 4, 6, 10