u12BoxplotNb<-function(
	x, sampleAnnoCol=2, intExCol="int_ex", intTypeCol="int_type", 
	intronExon, strandCol="strand", geneIdCol, col=c(), names=c(), lasNames=1, 
	outline=FALSE, plotLegend=TRUE, cexLegend=1, xLegend="topright", 
	yLegend=NULL, legend=c(), ...)
{
	object=x
	group=as.vector(SummarizedExperiment::colData(object)[,sampleAnnoCol])
	uniGroup=unique(as.vector(group))
	if(length(col)!=1 & length(col)!=length(uniGroup) ){
		stop(paste("Length of col should be either 1 or ",length(uniGroup),
			", i.e. the same length as c('", 
			paste(uniGroup, collapse="','", sep=""), "')", collapse='', 
				sep=""))
	}
	u12NbInd= u12NbIndex(object, intExCol=intExCol, intTypeCol=intTypeCol, 
		strandCol=strandCol, geneIdCol=geneIdCol, naUnstrand=FALSE)
	u12Ind=u12Index(object, intExCol=intExCol, intTypeCol=intTypeCol)
	datFpkm=scaledRetention(object)
	if(length(col)==1){
		col=rep(col,length(uniGroup))
	}
	if(length(col)>1)
		names(col)=uniGroup	
	upPlot=c()
	u12Plot=c()
	dnPlot=c()
	updnPlot=c()
	upIntType=SummarizedExperiment::rowData(object)[u12NbInd$upIntron, 
		intTypeCol]
	dnIntType=SummarizedExperiment::rowData(object)[u12NbInd$downIntron, 
		intTypeCol]
	for(cnt in 1:ncol(datFpkm)){
		if(intronExon=="intron"){
			upTmp=datFpkm[u12NbInd$upIntron,cnt]
			dnTmp=datFpkm[u12NbInd$downIntron,cnt]
			upTmp[upIntType=="U12"]=NA
			dnTmp[dnIntType=="U12"]=NA
			upPlot=c(upPlot, list(upTmp))
			u12Plot=c(u12Plot, list(datFpkm[u12Ind,cnt]))
			dnPlot=c(dnPlot, list(dnTmp))
			updnPlot=c(updnPlot, list(c(upTmp,dnTmp)))
		} else if(intronExon=="exon"){
			#Fix for strand
			upIndU2Tmp1=u12NbInd$upExon-1
			upIndU2Tmp2=u12NbInd$upExon-2
			upIndU2Tmp1[u12NbInd$upExon>u12NbInd$downExon]=
				u12NbInd$upExon[u12NbInd$upExon>u12NbInd$downExon]+1
			upIndU2Tmp2[u12NbInd$upExon>u12NbInd$downExon]=
				u12NbInd$upExon[u12NbInd$upExon>u12NbInd$downExon]+2
			dnIndU2Tmp1=u12NbInd$downExon+1
			dnIndU2Tmp2=u12NbInd$downExon+2
			dnIndU2Tmp1[u12NbInd$upExon>u12NbInd$downExon]=
				u12NbInd$downExon[u12NbInd$upExon>u12NbInd$downExon]-1
			dnIndU2Tmp2[u12NbInd$upExon>u12NbInd$downExon]=
				u12NbInd$downExon[u12NbInd$upExon>u12NbInd$downExon]-2

			upIndsU2=unique(c(upIndU2Tmp1, upIndU2Tmp2))
			dnIndsU2=unique(c(dnIndU2Tmp1, dnIndU2Tmp1))
			upIndsU12=u12NbInd$upExon
			dnIndsU12=u12NbInd$downExon
			upIndsU2=unique(upIndsU2[which(!is.na(upIndsU2))])
			dnIndsU2=unique(dnIndsU2[which(!is.na(dnIndsU2))])
			upIndsU12=unique(upIndsU12[which(!is.na(upIndsU12))])
			dnIndsU12=unique(dnIndsU12[which(!is.na(dnIndsU12))])
			upTmp= datFpkm[upIndsU2[which(is.na(match(upIndsU2, 
				c(upIndsU12,dnIndsU12))))], cnt]
			dnTmp= datFpkm[dnIndsU2[which(is.na(match(dnIndsU2, 
				c(upIndsU12,dnIndsU12))))], cnt]
			upPlot= c(upPlot, list(upTmp))
			u12Plot= c(u12Plot, 
				list(datFpkm[unique(c(upIndsU12, dnIndsU12)),cnt]))
			dnPlot=c(dnPlot, list(dnTmp))
			updnPlot=c(updnPlot, list(c(upTmp,dnTmp)))			
		}else {
			msg<- 
'The intronExon parameter should be defined and set to either "intron" for 
plotting intron retention levels or "exon" for plotting exon-exon junction 
levels.'
			stop(msg)
		}
	}
	
	color=c(as.vector(col[group]), NA)
	
	strand=SummarizedExperiment::rowData(object)[u12Ind, strandCol]

	if(unique(strand=="*") & length(unique(strand))==1){
		typePlot=2
		plotList=c(updnPlot,list(NA),u12Plot)
	} else {
		typePlot=1
		plotList=c(upPlot,list(NA),u12Plot,list(NA),dnPlot)
	}
	if(typePlot==1){
		graphics::boxplot(plotList, names=c(), xaxt = "n", outline=outline, 
			col=color, ...)
		axisAt=seq(from= trunc((length(plotList)-2)/6), 
			by=trunc((length(plotList)-2)/3)+1, length.out=3)
		graphics::axis(1 ,at=axisAt, labels=c("Upstream U2 intron", 
			"U12 intron", "Downstream U2 intron"), las=lasNames, ...)
		if(plotLegend)
			graphics::legend(x=xLegend,y=yLegend, legend=uniGroup, fill=col, 
				cex=cexLegend)
	} else if(typePlot==2){
		graphics::boxplot(plotList, names=c(), xaxt = "n", outline=outline, 
			col=color, ...)
		axisAt=seq(from= trunc((length(plotList)-1)/4), 
			by=trunc((length(plotList)-1)/2)+1, length.out=2)
		graphics::axis(1 ,at=axisAt, 
			labels=c("(Up/Down)stream U2 intron", "U12 intron"), las=lasNames, 
			...)
		if(length(legend)==0 & plotLegend)
			legend=uniGroup
		if(plotLegend)
			graphics::legend(x=xLegend,y=yLegend, legend=legend, 
				fill=unique(color), cex=cexLegend)
	}

}

