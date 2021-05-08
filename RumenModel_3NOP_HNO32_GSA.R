
# Assign constants to Constants object
GSA_constants <- unlist(c(NP_pars[c("k_HyEm","v_HeAc","v_HeAP","v_HeBu","v_HyMe","M_Hy_HyMe","v_VfAb")],
                          Constants[which(names(Constants) != 'k_NPAb')])) 

# Assign pars to pars object
GSA_pars <- unlist(c(NP_pars[c("J_NP_HyMe","k_NPNi")],k_NPNa=as.numeric(3.5*NP_pars['k_NPNi']),Constants['k_NPAb'])) 


# Prepare parameter matrix
LH_par <- cbind(min=0.75*GSA_pars,max=1.25*GSA_pars)
n_smpl <- 1000

# Points used for GSA
t_GSA <- final - c(12,11.5,11,10,8,6,2)


NP_pars_LHS <- list(as.data.frame(Latinhyper(LH_par,num=n_smpl)),as.data.frame(matrix(ncol=length(t_GSA),nrow=n_smpl)))

methane_times <- paste0('CH4_',c('00','05','10','20','40','60','100'))
colnames(NP_pars_LHS[[2]]) <- methane_times

# Run model over parameter space

for (i in 1:n_smpl){
  cat('i ran to',i,' ')
  out_099 <- NAD_confurc(times=times_nf,initial=initial,feed=c(feedcomp,c_3NOP=0.991e-3),Parms=NP_pars_LHS[[1]][i,],
                         Constants=GSA_constants,
                         stpsz=delt,Yields=Yields,feedintake=FI_rate_2m,output="FINAL_24h")
    NP_pars_LHS[[2]][i,] <- out_099$P_Mt_HyMe[out_099$time %in% t_GSA]
} # End of for loop i


NP_pars_LHS <- cbind(NP_pars_LHS[[1]],NP_pars_LHS[[2]])


# Write and plot

write.table(NP_pars_LHS,paste0('GSA_3NOP_HNO32_',n_smpl,'.csv'),sep=',',row.names=F)
#NP_pars_LHS <- read.table(paste0('GSA_3NOP_HNO32_',n_smpl,'.csv'),sep=',',header=T)

library(corrplot)

NP_cor <- cor(NP_pars_LHS)
NP_cor # shade circle square

write.table(NP_cor,paste0('GSA_3NOP_HNO32_cor_',n_smpl,'.csv'),sep=',',row.names=F)


pdf('GSA_3NOP_HNO2.pdf',width=7,height=5)
par(mar=c(0,0,0,0))
corrplot(NP_cor[1:4,methane_times], method = 'shade',
         addgrid.col ='lightgrey', tl.col = "black")
dev.off()


GSA_nitrate <- read.table('GSA_HNO32_1000.csv',sep=',',header=T)
NaNi_cor <- cor(GSA_nitrate)
# 0, 0.5, 1, 2, 4, 6, 10

#NP_cor <- read.table(paste0('GSA_3NOP_HNO32_cor_',n_smpl,'.csv'),sep=',',header=T)

colnames(NP_cor)[5:11] <- c('0.0','0.5','1.0','2.0','4.0','6.0','10.0')
colnames(NaNi_cor)[5:11] <- c('0.0','0.5','1.0','2.0','4.0','6.0','10.0')

NaNi_cor <- NaNi_cor[c(3,1,2,4),]

pdf('GSA_3NOP_Na.pdf',width=8,height=5)
#tiff('GSA_3NOP_Na.tiff',width=8,height=5,units='in',res=600,compression='lzw')
par(mfrow=c(1,2),mar=c(2,2,2,0),oma=c(2,2,0,0.25))
barplot(as.matrix(NP_cor[1:4,c('0.0','0.5','1.0','2.0','4.0','6.0','10.0')]),beside=T,ylim=c(-0.2,1.4),cex.lab=1.25,cex.axis=1.1,cex.names = 1.1,
        legend=F, col=c('dodgerblue','magenta','orange','limegreen'),xlab='',ylab='',main='3-NOP+nitrite model',yaxt='n')
axis(2,at=seq(-0.2,1,0.2),labels = c('-0.2','0.0','0.2','0.4','0.6','0.8','1.0'),cex.lab=1.25,cex.axis=1.1)
legend(23.5,1.325,legend=c(expression(italic(J)[paste('MCR;',H[2],',Me')]),
                        expression(italic(k)[paste('3NOP,',NO[2]^{'-'})]),
                        expression(italic(k)[paste('3NOP,',NO[3]^{'-'})]),
                        expression(italic(k)['3NOP,Ab'])),
       cex=1.1,col = c('dodgerblue','magenta','orange','limegreen'),pch=15, bty='n')
abline(h=0)
barplot(NaNi_cor[1:4,c('0.0','0.5','1.0','2.0','4.0','6.0','10.0')],beside=T,ylim=c(-0.2,1.4),cex.lab=1.25,cex.axis=1.1,cex.names = 1.1,
        legend=F, col=c('dodgerblue','magenta','orange','limegreen'),xlab='',ylab='',main='Nitrate+nitrite model',yaxt='n')
legend(23.5,1.325,legend=c(expression(italic(J)[paste(NO[2]^{'-'},';',H[2],',Me')]),expression(italic(k)[paste(NO[3]^{'-'},',',NO[2]^{'-'})]),expression(italic(k)[paste(NO[2]^{'-'},',',NH[3])]),expression(italic(k)[paste(NO[italic(x)]^{'-'},',Ab')])),
       cex=1.1,col = c('dodgerblue','magenta','orange','limegreen'),pch=15, bty='n')
abline(h=0)
mtext(side=1,'Time (h)',outer=T,cex=1.15,line=0.5)
mtext(side=2,'Correlation coefficient',outer=T,cex=1.15,line=0.75)
dev.off()
