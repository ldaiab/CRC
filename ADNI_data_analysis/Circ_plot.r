#cd /overflow/KernelLogrankTest-master/data/protein_result;module load r; R
library(Matrix)
library(qqman)
library(ggplot2) # Data visualization
library(data.table)
library(xgboost)
library(caret)
#library(qvalue)
library(circlize)
library(car)
require(sparseLDA)
require(sda)
require(gridExtra)
library(ComplexHeatmap)
A=read.csv('C:\\Users\\RealData\\output\\LOG_U1.csv',check.names = FALSE)
#A=read.csv('N01_ADNI1.csv',check.names = FALSE)
A=A[,-1]
M1=read.csv('C:\\coefficients\\coefficients\\RealData\\PseudoRealData\\protein_dic.csv');
temp=data.frame(list(gene=M1$Target,count=rep(0,length(M1$Target))))
UG=unique(M1$Target);L=length(UG)
temp3=duplicated(M1$Target)
for(ii in 1:L)
{
temp2=which(M1$Target==UG[ii])
temp[temp2,2]=cumsum(temp3[temp2])
}
temp1=paste0(M1$Target,'_',temp[,2])
temp1=gsub('_0','',temp1)
p=dim(A)[1]
modality=temp1
out='C:\\Users\\coefficients\\coefficients\\RealData\\Script'
Q0=p.adjust(A[,2],'fdr')
stat0=data.frame(list(COR=A[,4],pval=A[,2],qval=Q0))
rownames(stat0)=modality
PATH0='C:\\Users\\coefficients\\coefficients\\RealData\\Script\\circ.pdf'
temp0=sort(stat0$pval,ind=T)$ix
stat0=stat0[temp0,]
stat0=stat0[1:100,]
stat0[stat0[,1]>0.999,1]=0.999
rownames(stat0)=gsub('_1','',rownames(stat0))
##################################################
circplot<-function(stat0,PATH0='C:\\Users\\coefficients\\coefficients\\RealData\\Script\\circ.pdf')
{
	factorstemp=substr(as.character(stat0$COR+rnorm(length(stat0$COR),0,0.001)),1,7)
	factors2=factor(factorstemp,levels =factorstemp)
	gsub1<-function(i,vec,vec1){return(gsub(vec[i],'',vec1[i]))}
	modality=unlist(lapply(strsplit(rownames(stat0),'_'),'[',1))
	factors3=unlist(lapply(1:length(modality),gsub1,vec=paste0(modality,'_'),vec1=rownames(stat0)))
	grDevices::cairo_pdf(PATH0,width=15,height=12)
	par(mar=c(6,0,5,0), xpd=T)
	par(bg = "white")
	rbPal <- colorRampPalette(c("#00BFFF", "yellow", "#FF4500"))
	allh <- stat0$COR
	col_hc<- rbPal(25)[as.numeric(cut(allh,breaks = 25))]
	col_hf<- rbPal(25)[as.numeric(cut(-log(stat0$pval)/log(10),breaks = 25))]
	col_hm<- rbPal(25)[as.numeric(cut(-log(stat0$qval)/log(10),breaks = 25))]
	par(bg = "white")
	circos.par("track.height" = 0.2)
	circos.initialize(factors2,xlim=c(0,3))
	print('succeed!!!')
	circos.track(factors = factors2,ylim = c(0, 1), panel.fun = function(x, y) {
	  chr = get.cell.meta.data("sector.index")
	  xlim = get.cell.meta.data("xlim")
	  ylim = get.cell.meta.data("ylim")
	  circos.rect(0, xlim[1], xlim[2], 1,border = NA
	  )
	  circos.text(mean(xlim),mean(ylim)-0.3,substr(allh,1,5)[which(factorstemp==chr)], 
	  cex = 1.3, adj = c(0, degree(0)),facing = "clockwise", niceFacing = T)
	  circos.text(mean(xlim),mean(ylim)+0.7,factors3[which(factorstemp==chr)], cex = 1.3, 
	  adj = c(0, degree(0)),facing = "clockwise", niceFacing = T)
	}, bg.border = NA,bg.col=col_hc)
	print('succeed!!!')
	circos.track(factors = factors2, ylim = c(0, 1), panel.fun = function(x, y) {
	  chr = get.cell.meta.data("sector.index")
	  xlim = get.cell.meta.data("xlim")
	  ylim = get.cell.meta.data("ylim")
	  circos.rect(0, xlim[1], xlim[2], 1, border = NA)
	  circos.text(mean(xlim), mean(ylim) - 0.3, 
				  labels = formatC(stat0$pval[which(factorstemp == chr)], format = "f", digits = 4), 
				  cex = 1.3, adj = c(0, degree(0)), facing = "clockwise", niceFacing = T)
	}, bg.border = NA, bg.col = col_hf)
	circos.track(factors = factors2, ylim = c(0, 1), panel.fun = function(x, y) {
	  chr = get.cell.meta.data("sector.index")
	  xlim = get.cell.meta.data("xlim")
	  ylim = get.cell.meta.data("ylim")
	  circos.rect(0, xlim[1], xlim[2], 1, border = NA)
	  circos.text(mean(xlim), mean(ylim) - 0.3, 
				  labels = formatC(stat0$qval[which(factorstemp == chr)], format = "f", digits = 4), 
				  cex = 1.3, adj = c(0, degree(0)), facing = "clockwise", niceFacing = T)
	}, bg.border = NA, bg.col = col_hm)
	circos.info()
	tt1=-log(stat0$pval)/log(10)
	tt2=-log(stat0$qval)/log(10)
	A0=c(round(min(stat0$COR),2), round((min(stat0$COR)+max(stat0$COR))/2,2),round(max(stat0$COR),2))
	col_fun0 = colorRamp2(A0, c("#00BFFF", "yellow", "#FF4500"))
	#lgd_links = Legend(at = A0,labels=A0, col_fun=colorRamp2(A0, c("green", "yellow", "red")),grid_width=NULL,legend_width=1,title_position = "topleft", title = "Outer: correlation", direction = "horizontal")
	A1=round(10^(-c(round(min(tt1),1), round((min(tt1)+max(tt1))/2,1),round(max(tt1),1))),5)
	col_fun1 = colorRamp2(-log10(A1),c("#00BFFF", "yellow", "#FF4500"))
	#lgd_links2 = Legend(at = A0,labels=A0, col_fun=colorRamp2(A0, c("green", "yellow", "red")),grid_width=NULL,legend_width=1,title_position = "topleft", title = expression("Middle: -log"["10"]*"p"), direction = "horizontal")
	A2=round(10^(-c(round(min(tt2),1), round((min(tt2)+max(tt2))/2,1),round(max(tt2),1))),2)
	col_fun2 = colorRamp2(-log10(A2), c("#00BFFF", "yellow", "#FF4500"))
	#lgd_links3 = Legend(at = A0,labels=A0, col_fun=colorRamp2(A0, c("green", "yellow", "red")),grid_width=NULL,legend_width=1,title_position = "topleft", title = expression("Inner: -log"["10"]*"q"), direction = "horizontal")
	lgd_links = Legend(at = A0,labels=A0, col_fun=col_fun0,title_position = "topleft", title = "Outer: correlation", 
	title_gp = gpar(fontsize = 18, fontface = "bold"),labels_gp = gpar(fontsize = 16, fontface = "bold"),direction = "horizontal")
	lgd_links2 = Legend(at = -log10(A1),labels=A1, col_fun=col_fun1,title_position = "topleft", title = "Middle: p-value", 
	title_gp = gpar(fontsize = 18, fontface = "bold"),labels_gp = gpar(fontsize = 16, fontface = "bold"),direction = "horizontal")
	lgd_links3 = Legend(at = -log10(A2),labels=A2, col_fun=col_fun2,title_position = "topleft", title = "Inner: FDR", 
	title_gp = gpar(fontsize = 18, fontface = "bold"),labels_gp = gpar(fontsize = 16, fontface = "bold"),direction = "horizontal")
	lgd_list_vertical<-packLegend(lgd_links,lgd_links2,lgd_links3)
	lgd_list_vertical
	lgd_grob <- grid.grabExpr(grid.draw(lgd_list_vertical))

	# Position and draw
	pushViewport(viewport(
	  x = unit(50, "mm"), 
	  y = unit(40, "mm"),
	  width = grobWidth(lgd_grob)*1.5,
	  height = grobHeight(lgd_grob)*1.5,
	  just = c("left", "bottom")
	))
	grid.draw(lgd_grob)
	
	#pushViewport(viewport(x = unit(20, "mm"), y = unit(30, "mm"),   width = grobWidth(lgd_list_vertical),  
	#height = grobHeight(lgd_list_vertical), just = c("left", "bottom")))
	#grid.draw(lgd_list_vertical)
	upViewport()
	circos.clear()
	dev.off()
}
circplot(stat0)
