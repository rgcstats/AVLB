simresults <- NULL
for(k in c(1:184)){
  load(paste0("simres",k,".rdata"))
  simresults <- rbind(simresults,simres)
}

array.ratio.mse <- with(simresults,
                        tapply(error^2,list(n=n,method=method,bw=bw,power.selnprob=power.selnprob,maxratio.z=maxratio.z,work.var.index=work.var.index),mean)
                        /tapply(UB,list(n=n,method=method,bw=bw,power.selnprob=power.selnprob,maxratio.z=maxratio.z,work.var.index=work.var.index),mean))

save(simresults,array.ratio.mse,file="simresults_250119.rdata")


# Table 1: MSE-ratio, Best-BLUP, n=500, var correct
tab1 <- rbind(t(array.ratio.mse["500","best",,"0.5",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["500","best",,"1",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["500","best",,"2",c("1","1.5","2.5"),"1"]))
tab1 <- format(round(tab1,digits=3),nsmall=3)
tab1 <- cbind(power.selnprob=rep(c("0.5","1","2"),each=3),
              maxratio.z=as.numeric(rownames(tab1))-1,tab1)
tab1

# Table 1b: MSE-ratio, Best-BLUP, n=25000, var correct
tab1b <- rbind(t(array.ratio.mse["25000","best",,"0.5",c("1","1.5","2.5"),"1"]),
               t(array.ratio.mse["25000","best",,"1",c("1","1.5","2.5"),"1"]),
               t(array.ratio.mse["25000","best",,"2",c("1","1.5","2.5"),"1"]))
tab1b <- format(round(tab1b,digits=3),nsmall=3)
tab1b <- cbind(power.selnprob=rep(c("0.5","1","2"),each=3),
               maxratio.z=as.numeric(rownames(tab1b))-1,tab1b)
tab1b

# Table 2: MSE-ratio, Best-BLUP, n=500, var incorrect (power=2 should be 1)
tab2 <- rbind(t(array.ratio.mse["500","best",,"0.5",c("1","1.5","2.5"),"2"]),
              t(array.ratio.mse["500","best",,"1",c("1","1.5","2.5"),"2"]),
              t(array.ratio.mse["500","best",,"2",c("1","1.5","2.5"),"2"]))
tab2 <- format(round(tab2,digits=3),nsmall=3)
tab2 <- cbind(power.selnprob=rep(c("0.5","1","2"),each=3),
              maxratio.z=as.numeric(rownames(tab2))-1,tab2)
tab2

# Table 3: MSE-ratio, Best-BLUP, n=1500, var correct (power=1 should be 1)
tab3 <- rbind(t(array.ratio.mse["1500","best",,"0.5",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["1500","best",,"1",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["1500","best",,"2",c("1","1.5","2.5"),"1"]))
tab3 <- format(round(tab3,digits=3),nsmall=3)
tab3 <- cbind(power.selnprob=rep(c("0.5","1","2"),each=3),
              maxratio.z=as.numeric(rownames(tab3))-1,tab3)
tab3

# Table 4: MSE-ratio, simple BLUP, n=500, var correct
tab4 <- rbind(t(array.ratio.mse["500","lm",,"0.5",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["500","lm",,"1",c("1","1.5","2.5"),"1"]),
              t(array.ratio.mse["500","lm",,"2",c("1","1.5","2.5"),"1"]))
tab4 <- format(round(tab4,digits=3),nsmall=3)
tab4 <- cbind(power.selnprob=rep(c("0.5","1","2"),each=3),
              maxratio.z=as.numeric(rownames(tab4))-1,tab4)
tab4

paperdir <- "/Users/rclark/Google Drive/Model-Based AV Bound/Survey Methodology second submission January 2019/"
write.table(tab1,file=paste0(paperdir,"table1.txt"),sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(tab2,file=paste0(paperdir,"table2.txt"),sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(tab3,file=paste0(paperdir,"table3.txt"),sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(tab4,file=paste0(paperdir,"table4.txt"),sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(tab1b[c(1,4,7),],file=paste0(paperdir,"table1b.txt"),sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

#################################
# Plot examples of mu(x)
#################################

xvals <- seq(from=-1,to=3,by=0.01)
setEPS()
postscript(paste0(paperdir,"mu_x.eps"))
plot(xvals,4*xvals,type="n",main="",xlab="x1",ylab="E[Y|x1]",
     ylim=c(0,7.5),xlim=c(0,1.75))
all.bw <- c(0.5,1,2,5)
for(k in c(1:4)) lines(xvals,4*xvals + 1*sin(xvals*2*pi/all.bw[k]),lty=k,col=k)
legend("topleft",lty=c(1:4),col=c(1:4),legend=paste0("period=",all.bw))
dev.off()

# Plot MSE ratio vs bw and power.selnprob for best estimator
#    n=1500, var correct, working power=1, seln prob prop to x^0.5
fig2dat0.5 <- array.ratio.mse["1500","best","0.5","0.5",,"1"]
fig2dat1 <- array.ratio.mse["1500","best","1","0.5",,"1"]
fig2dat2 <- array.ratio.mse["1500","best","2","0.5",,"1"]
fig2dat5 <- array.ratio.mse["1500","best","5","0.5",,"1"]
setEPS()
postscript(paste0(paperdir,"MSEplot.eps"))
plot(as.numeric(names(fig2dat0.5))-1,fig2dat0.5,ylim=c(0.6,1.02),pch=1,lty=1,
     type="b",col=1,xlab="c",ylab="Ratio of MSE to GJLB")
points(as.numeric(names(fig2dat1))-1,fig2dat1,pch=2,lty=2,type="b",col=2)
points(as.numeric(names(fig2dat2))-1,fig2dat2,pch=3,lty=3,type="b",col=3)
points(as.numeric(names(fig2dat5))-1,fig2dat5,pch=4,lty=4,type="b",col=4)
legend("topright",lty=1:4,pch=1:4,col=1:4,
       legend=c("0.5","1","2","5"),title="period (h)")
dev.off()

########################################
# How often was correct model chosen?
########################################

models.with.z <- c(1,3,5,7,9,11,13,15,17,19,21,23,25)
models.sans.x1 <- c(1,2)
table(simresults$model.chosen[simresults$method=="best"],exclude=NULL)
table(simresults$model.chosen[simresults$method=="best"] %in% models.with.z,exclude=NULL)
mean(simresults$model.chosen[simresults$method=="best"] %in% models.with.z,exclude=NULL)
100*(1-with(simresults[simresults$method=="best",],
            tapply(model.chosen %in% models.with.z,k,FUN=mean)))
range(with(simresults[simresults$method=="best",],
           tapply(model.chosen %in% models.with.z,k,FUN=mean)))

