cd output/;module load r/4.4.0;R
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
	   