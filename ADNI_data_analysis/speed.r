cd output/;module load r/4.4.0;R
T0=10^(seq(4,7,length.out=10))
Tcost=matrix(NA, 10,5)
for(ii in 1:10)
{
print(10-ii)
temp=system.time(permu_pvalue_fast(Y, X, delta,B=T0[ii]))
Tcost[ii,1]=temp[3]
print(10-ii)
temp=system.time(permu_pvalue_fast(Y, X, delta,B=T0[ii],Fw=function(x){punif(x,min=-1,max=1)}))
Tcost[ii,2]=temp[3]
print(10-ii)
temp=system.time(ipcw.dcov.test(cbind(Y,delta), X,B=T0[ii]))
Tcost[ii,3]=temp[3]
print(10-ii)
temp=system.time(wild_bootstrap_test_logrank_covariates(array(as.matrix(X),c(length(Y),1)),array(Y),array(delta),'gau','gau', 
num_bootstrap_statistics=as.integer(T0[ii])))
Tcost[ii,4]=temp[3]
print(10-ii)
temp=system.time(cph_test(X=as.matrix(X),z=array(Y),d=array(delta)))
Tcost[ii,5]=temp[3]
}
write.csv(Tcost,file='timeT.csv',quote=F)

### make a time plot ###
Tcost=read.csv('timeT.csv');
pdf('timecost.pdf')
# Plot with custom colors/shapes and x-axis specification
matplot(Tcost[,1], Tcost[,-c(1,6)],  # Specify x-values explicitly
        col = 1:5,              # 5 distinct colors (for methods 1-5)
        pch = 1:5,              # 5 distinct point shapes
        type = "b",             # Both points and lines
        xaxt = "n",             # Suppress default x-axis
        xlab = "B (permutation/bootstrap times)", 
        ylab = "Time cost (in seconds)",
        main = "Time cost Comparison",lwd=2,cex=1.5,cex.lab=1.5,cex.axis=1.5)

# Custom x-axis with 10^ notation
T0=10^(seq(4,7,length.out=10))
axis(1, at = c(1,4,7,10), 
     labels = expression(10^4, 10^5, 10^6, 10^7),
     las = 1,cex.axis=1.5,cex.lab=1.5)  # Horizontal labels

legend("topleft", 
       legend = parse(text = c("CRC[N(0,1)]", "CRC[U(-1,1)]", "IPCW", "KLR")),
	   col = 1:5, 
       pch = 1:5,
	   lwd=2,title = "Methods",cex=1.5,title.cex=1.5)

dev.off()
	   
