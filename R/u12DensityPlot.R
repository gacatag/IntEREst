u12DensityPlot<-function(x, type=c("U12", "U2Up", "U2Dn", "U2UpDn", "U2Rand"), 
	fcType="edgeR", sampleAnnotation=c(),  sampleAnnoCol=c(), group=c(), 
	intExCol="int_ex", intTypeCol="int_type", intronExon, strandCol="strand", 
	geneIdCol="collapsed_transcripts", naUnstrand=FALSE, col=1, lty=1, lwd=1, 
	plotLegend=TRUE, cexLegend=1, xLegend="topright", yLegend=NULL, legend=c(), 
	randomSeed=NULL, xlab="", ...){
	#Check parameters
	object=x
	if(length(col)==1 ){
		col=rep(col, length(type))
	} else if (length(col)!=length(type)){
		msg<-
"The length of the col parameter should be either 1 or the same size as type."
		stop(msg)
	}
	if(length(lty)==1 ){
		lty=rep(lty, length(type))
	} else if (length(lty)!=length(type)){
		msg<-
"The length of the lty parameter should be either 1 or the same size as type."
		stop(msg)
	}
	if(length(lwd)==1 ){
		lwd=rep(lwd, length(type))
	} else if (length(lwd)!=length(type)){
		msg<-
"The length of the lwd parameter should be either 1 or the same size as type."
		stop(msg)
	}

	fcRes=lfc(object, fcType=fcType, sampleAnnoCol=sampleAnnoCol, 
		sampleAnnotation=sampleAnnotation, group=group)
	u12Ind=u12Index(object,intExCol=intExCol, intTypeCol=intTypeCol)
	u12NbInd=u12NbIndex(object,intExCol=intExCol, intTypeCol=intTypeCol, 
		strandCol=strandCol, geneIdCol=geneIdCol, naUnstrand=naUnstrand)
	uniStrand=unique(SummarizedExperiment::rowData(object)[,strandCol])
	definedStrand=uniStrand[uniStrand!="*"]
	listPlot=c()	
	if(intronExon=="intron"){
		if("U12" %in% type)
			listPlot=c(listPlot, list(U12_introns=fcRes[u12Ind]))
		if("U2Up" %in% type)
			listPlot=c(listPlot, 
				list(Upstream_U2_introns=fcRes[u12NbInd$upIntron]))
		if("U2Dn" %in% type)
			listPlot=c(listPlot, 
				list(Downstream_U2_introns=fcRes[u12NbInd$downIntron]))
		if("U2UpDn" %in% type)
			listPlot=c(listPlot, list(UpDown_stream_U2_introns=
				fcRes[c(u12NbInd$downIntron,u12NbInd$upIntron)]))
		if("U2Rand"%in% type){
			set.seed(randomSeed)
			allIntronsInd=which(SummarizedExperiment::rowData(object)[,
				intExCol]=="intron")
			randInd=sample(
				allIntronsInd[which(is.na(match(allIntronsInd,u12Ind)))], 
				length(u12Ind))
			listPlot=c(listPlot, list(random_U2_introns=fcRes[randInd]))
		}
	} else if(intronExon=="exon"){
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
		upDnIndsU2=unique(c(upIndsU2,dnIndsU2))[
			which(is.na(match(unique(c(upIndsU2,dnIndsU2)), 
				c(upIndsU12,dnIndsU12))))]

		if("U12" %in% type)
			listPlot=c(listPlot, list(U12_introns=
				fcRes[unique(c(upIndsU12, dnIndsU12))]))
		if("U2Up" %in% type)
			listPlot=c(listPlot, list(Upstream_U2_introns=
				fcRes[ upIndsU2[which(is.na(match(upIndsU2, 
					c(upIndsU12,dnIndsU12))))] ]))
		if("U2Dn" %in% type)
			listPlot=c(listPlot, list(Downstream_U2_introns=
				fcRes[ dnIndsU2[which(is.na(match(dnIndsU2, 
					c(upIndsU12,dnIndsU12))))] ]))
		if("U2UpDn" %in% type)
			listPlot=c(listPlot, list(UpDown_stream_U2_introns=
				fcRes[ upDnIndsU2 ]))
		if("U2Rand"%in% type){
			set.seed(randomSeed)
			allExonsInd=which(SummarizedExperiment::rowData(object)[,
				intExCol]=="exon")
			randInd=sample(allExonsInd[ which(is.na(match(allExonsInd, 
				unique(c(upIndsU12, dnIndsU12)) )))], length(u12Ind))
			listPlot=c(listPlot, list(random_U2_introns=fcRes[randInd]))
		}
	}
	# Temporary index variables for correctly assigning lty col, and lwd
	indOrdFull=match(c("U12","U2Up","U2Dn","U2UpDn","U2Rand"), type)
	indOrd=indOrdFull[which(!is.na(indOrdFull))]

	#Index to use for correct order of listPlot 
	#(based on order of type definition)
	indOrd2=1:length(indOrd)
	indOrd2[indOrd]=1:length(indOrd)

	for(i in 1:length(listPlot)){
		if(i==1){
			graphics::plot(stats::density(unlist(listPlot[[i]]), na.rm=TRUE), 
				type='l', lty=lty[indOrd[i]], lwd=lwd[indOrd[i]], 
				col=col[indOrd[i]], main="", xlab=xlab, ...)
		}else{
			graphics::points(stats::density(unlist(listPlot[[i]]), na.rm=TRUE),
				type='l', lty=lty[indOrd[i]], lwd=lwd[indOrd[i]], 
				col=col[indOrd[i]])
		}
	}

	if(plotLegend | length(legend)!=0){
		if(length(legend)==0){
			legendTmp=c("U12 introns", "U2 upstream introns", 
				"U2 downstream introns", "Up/Down stream U2 introns", 
				"Random U2 introns" )
			if(intronExon=="exon")
				legendTmp=c("Exons flanking U12 introns", 
					"Exons flanking upstream U2 introns", 
					"Exons flanking downstream U2 introns", 
					"Exons flanking up/down stream U2 introns", 
					"Exons flanking random U2 introns" )
			names(legendTmp)=c("U12","U2Up","U2Dn","U2UpDn","U2Rand")
			#Index to order the correct legend names
			indOrdLegend=match(type,c("U12","U2Up","U2Dn","U2UpDn","U2Rand"))
			legend=paste(paste(legendTmp[indOrdLegend], unlist(lapply(indOrd2,
				function(tmp)  
					return(length(unlist(listPlot[[tmp]])[
						which(!is.na(unlist(listPlot[[tmp]])))])) )), 
				sep=" (n="), ")", sep="")
		}
		graphics::legend(x=xLegend,y=yLegend, legend=legend, col=col, lty=lty,
			lwd=lwd, cex=cexLegend)
	}


}


