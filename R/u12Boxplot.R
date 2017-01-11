u12Boxplot<-function(x, sampleAnnoCol=NA, intExCol="int_ex", 
	intTypeCol="int_type", intronExon, col="white", boxplotNames=c(), 
	lasNames=3, outline=FALSE, ...)
{
	object=x
	if(!is.na(sampleAnnoCol)){
		groups=object@sampleAnnotation[,sampleAnnoCol]
	} else {
		groups=x@sampleNames
	}
	uniGroup=unique(as.vector(groups))
	color=c()	
	plotList=c()
	axisAt=c()
	for(cnt in 1:length(uniGroup)){
		axisAt=c(axisAt,TRUE,TRUE,FALSE)
		if(length(col)>1)
			color=c(color, col[(cnt*2)-1],col[(cnt*2)], NA)
		if(length(col)==1)
			color=c(color, col, col, NA)
		if(intronExon=="intron"){
			plotList=c(plotList, 
				list( as.vector(
					unlist(object@scaledRetention[
						object@interestDf[,intExCol]=="intron" & 
						object@interestDf[,intTypeCol]=="U12", 
						which(groups==uniGroup[cnt])])), 
					as.vector(unlist(object@scaledRetention[
						object@interestDf[,intExCol]=="intron" & 
						object@interestDf[,intTypeCol]=="U2",
					which(groups==uniGroup[cnt])])), NA))
		} else if (intronExon=="exon") {
			indChooseIntU12=unique(which(object@interestDf[,intExCol]=="intron"
				& object@interestDf[,intTypeCol]=="U12"))
			indChooseIntU2=unique(which(object@interestDf[,intExCol]=="intron"
				& (object@interestDf[,intTypeCol]!="U12") ))

			indChooseExU12=unique(unlist(lapply(indChooseIntU12, function(tmp)
				return(c(tmp-1, tmp+1)))))
			indChooseExU2=unique(unlist(lapply(indChooseIntU2, function(tmp)
				return(c(tmp-1, tmp+1)))))
			indChooseExU2=indChooseExU2[which(is.na(match(indChooseExU2,
				indChooseExU12)))]
			plotList=c(plotList, 
				list( as.vector(unlist(object@scaledRetention[indChooseExU12,
					which(groups==uniGroup[cnt])])), 
					as.vector(unlist(object@interestDf[indChooseExU2,
					which(groups==uniGroup[cnt])])), NA))				
		}
		if(length(boxplotNames)==0){
			names(plotList)[(length(plotList)-2):length(plotList)]=
				c(paste("U12",uniGroup[cnt], sep=" "),
					paste("U2",uniGroup[cnt], sep=" "),"")
		} else {
			names(plotList)[(length(plotList)-2):length(plotList)]=
				c(boxplotNames[(cnt*2)-1],boxplotNames[(cnt*2)],"")
		}
	}
	plotList=plotList[-length(plotList)]
	axisAt=axisAt[-length(axisAt)]
	color=color[-length(color)]

	graphics::boxplot(plotList, names=c(), xaxt = "n", outline=outline, 
		col=color, ...)
	if(length(boxplotNames)==0)
		boxplotNames=names(plotList)
	graphics::axis(1,at=which(axisAt), labels=boxplotNames[boxplotNames!=""],
		las=lasNames, ...)

}

