u12DensityPlotIntron<-function(x, type=c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), fcType="edgeR", sampleAnnotation=c(),  sampleAnnoCol=c(), group=c(), intExCol="int_ex", intTypeCol="int_type", strandCol="strand", geneIdCol="collapsed_transcripts", naUnstrand=FALSE, col=1, lty=1, lwd=1, plotLegend=TRUE, cexLegend=1, xLegend="topright", yLegend=NULL, legend=c(), randomSeed=NULL, xlab="", ...){

	u12DensityPlot (x=x, type=type, fcType=fcType, 
	sampleAnnotation=sampleAnnotation,  sampleAnnoCol=sampleAnnoCol, 
	group=group, intExCol=intExCol, intTypeCol=intTypeCol, 
	intronExon="intron", strandCol=strandCol, geneIdCol=geneIdCol, 
	naUnstrand=naUnstrand, col=col, lty=lty, lwd=lwd, 
	plotLegend=plotLegend, cexLegend=cexLegend, xLegend=xLegend, 
	yLegend=yLegend, legend=legend, randomSeed=randomSeed, 
	xlab=xlab, ...)

}


