#Calling the constructor
InterestResult <- function(resultFiles=c(), counts, scaledRetention, 
	scaleLength, scaleFragment, sampleAnnotation, rowData){
	if(!missing(sampleAnnotation)){
		if(length(resultFiles)>0)
			sampleAnnotation<- data.frame(resultFiles=resultFiles, 
				sampleAnnotation, row.names= rownames(sampleAnnotation))
		res<- SummarizedExperiment::SummarizedExperiment(
			assays=list(counts=counts, 
			scaledRetention=scaledRetention), 
			rowData=rowData,
			colData=sampleAnnotation, 
			metadata=list( scaleFragment= scaleFragment,
			scaleLength= scaleLength)		
		)
	} else {
		if(length(resultFiles)>0)
			sampleAnnotation<- data.frame(resultFiles=resultFiles,
				row.names=colnames(counts))
		res<- SummarizedExperiment::SummarizedExperiment(
			assays=list(counts=counts, 
			scaledRetention=scaledRetention), 
			rowData=rowData, 
			metadata=list( scaleFragment= scaleFragment,
			scaleLength= scaleLength)
		)
	}
	return(res)
}

# plot method
plot.InterestResult<-function(x, summary="none", subsetRows=NULL, 
	what="scaled", intronExon="intron", logScaleBase=NULL, logPseudoCnt=1, 
	plotLoess=TRUE, loessCol="red",	loessLwd=1, loessLty=1, cexText=1, 
	marPlot=c(2,2,2,2), mgpPlot=c(1, 1, 0), cexAxis=1, writeCor=TRUE, corCex=1,
	corMethod="pearson", corCol="grey63", upperCorXY=c("topleft", NULL),
	lowerCorXY= c("topleft", NULL), na.rm=TRUE, cex=1, sampleAnnoCol=c(),
	lowerPlot=FALSE, upperPlot=TRUE, ...){
	object=x
	if(summary=="mean")
		summaryFun=mean
	if(summary=="median")
		summaryFun=stats::median
	if(is.null(subsetRows))
		subsetRows=1:nrow(SummarizedExperiment::rowData(object))
	subDat=SummarizedExperiment::rowData(object)[subsetRows, ]

	indRow=!is.na(match(subDat$int_ex, intronExon))

	if(ncol(SummarizedExperiment::colData(object))==2){
		sampleGroups=SummarizedExperiment::colData(object)[,2]
	} else if( ncol(SummarizedExperiment::colData(object))>2 & 
		length(sampleAnnoCol)>0 & summary!="none"){
		sampleGroups=SummarizedExperiment::colData(object)[,sampleAnnoCol]
	} else if (ncol(SummarizedExperiment::colData(object))>2 & 
		length(sampleAnnoCol)==0 & summary!="none"){
		stop('Please define the sampleAnnoCol.')
	}

	if(summary=="none"& what=="scaled"){
		plotDat=scaledRetention(object)[subsetRows,]
	}
	if(summary=="none"& what=="counts"){
		plotDat=counts(object)[subsetRows,]
	}
	if(summary!="none"&what=="scaled"){
		plotList= tapply(1:ncol(scaledRetention(object)), sampleGroups, 
			function(x) apply(scaledRetention(object)[subsetRows,x], 1, 
				summaryFun) )
		plotList= plotList[match(unique(sampleGroups),names(plotList))]
		plotDat= matrix(unlist(plotList), ncol=length(unique(sampleGroups)), 
			byrow=FALSE)
		colnames(plotDat)= names(plotList)
		plotDat= as.data.frame(plotDat)
	}
	if(summary!="none"&what=="counts"){
		plotList=tapply(1:ncol(counts(object)), sampleGroups, function(x) 
			apply(counts(object)[subsetRows,x], 1, summaryFun) )
		plotList=plotList[match(unique(sampleGroups),names(plotList))]
		plotDat=matrix(unlist(plotList), ncol=length(unique(sampleGroups)), 
			byrow=FALSE)
		colnames(plotDat)=names(plotList)
		plotDat=as.data.frame(plotDat)
	}
	if(!is.null(logScaleBase))
		plotDat=log(plotDat+logPseudoCnt, base=logScaleBase)
	if(na.rm)
		indRow=indRow& (!is.na(rowSums(plotDat)))

	
	if(dim(plotDat)[2]==2){	
		loessInd=indRow& (!is.infinite(rowSums(plotDat)))
		loessDat=plotDat[loessInd, ]
		mainText=""
		if(writeCor)
			mainText=paste("cor (",corMethod,") = ", 
				format(stats::cor(loessDat[, 1], loessDat[,2], method=corMethod)
				, digits=2), sep="")
		graphics::par(mfrow=c(1,1))	
		graphics::par(cex.main=corCex*(3/ncol(plotDat)))
		graphics::par(mar=marPlot*(3/ncol(plotDat)))
		graphics::plot(plotDat[indRow, ], cex=cex*(3/ncol(plotDat)), 
			main=mainText, ...)
		if(dim(plotDat)[2]==2)
			if(plotLoess)
				graphics::points(stats::loess.smooth(plotDat[indRow, 1], 
					plotDat[indRow, 2]), col=loessCol, type="l", 
					lwd= loessLwd,lty=loessLty)
	}else{
		graphics::par(mfrow=c(ncol(plotDat),ncol(plotDat)))
		graphics::par(cex.main=corCex*(3/ncol(plotDat)))
		for(j in 1: ncol(plotDat)){
			for(i in 1:ncol(plotDat)){
				graphics::par(mar=marPlot*(3/ncol(plotDat)))
				if(i==j){
					graphics::plot(c(0, 1), c(0, 1), ann = FALSE, type = 'n', 
						xaxt = 'n', yaxt = 'n')
					graphics::text(x = 0.5, y = 0.5, colnames(plotDat)[i], 
						cex=cexText*(3/ncol(plotDat)))
				# If upperPlot is requested
				} else if(j<i & upperPlot) {
					xAxt='n'
					yAxt='n'
					mainText=""
					if(i==ncol(plotDat))
						yAxt='s'
					if(j==1)
						xAxt='s'

					loessInd=indRow& (!is.infinite(rowSums(plotDat)))
					loessDat=plotDat[loessInd, ]
					if(writeCor)
						mainText=paste("cor = ", format(stats::cor(loessDat[, i]
							, loessDat[,j], method=corMethod), digits=2), 
							sep="")

					graphics::plot(plotDat[indRow, i], plotDat[indRow, j], 
						ylab="", xlab="", xaxt='n', yaxt='n', 
						cex.axis=cexAxis*(3/ncol(plotDat)), 
						main="", cex=cex*(3/ncol(plotDat)), 
						xlim=c(min(plotDat[indRow,], na.rm=TRUE), 
							max(plotDat[indRow,], na.rm=TRUE)), 
						ylim=c(min(plotDat[indRow,], na.rm=TRUE), 
							max(plotDat[indRow,], na.rm=TRUE)), ...)
					if(xAxt=='s')
						graphics::axis(3, cex.axis=cexAxis*(3/ncol(plotDat)), 
							mgp=mgpPlot*(3/ncol(plotDat)))
					if(yAxt=='s')
						graphics::axis(4, cex.axis=cexAxis*(3/ncol(plotDat)), 
							mgp=mgpPlot*(3/ncol(plotDat)))

					graphics::legend(upperCorXY[1], upperCorXY[2], mainText, 
						bty="n", text.col=corCol, text.font=1, 
						cex=corCex*(3/ncol(plotDat)))
					if(plotLoess)
						graphics::points(stats::loess.smooth(loessDat[,i], 
							loessDat[,j]), col=loessCol, type="l", 
							lwd= loessLwd,lty=loessLty)
				} else if(j>i & lowerPlot) {
					xAxt='n'
					yAxt='n'
					mainText=""
					if(j==ncol(plotDat))
						xAxt='s'

					if(i==1)
						yAxt='s'						

					loessInd=indRow& (!is.infinite(rowSums(plotDat)))
					loessDat=plotDat[loessInd, ]
					if(writeCor)
						mainText=paste("cor = ", format(stats::cor(loessDat[, i]
							, loessDat[,j], method=corMethod) , digits=2), 
							sep="")

					graphics::plot(plotDat[indRow, i], plotDat[indRow, j], 
						ylab="", xlab="", xaxt='n', yaxt='n', 
						cex.axis=cexAxis*(3/ncol(plotDat)), 
						main="", cex=cex*(3/ncol(plotDat)), 
						xlim=c(min(plotDat[indRow,], na.rm=TRUE), 
							max(plotDat[indRow,], na.rm=TRUE)),
						ylim=c(min(plotDat[indRow,], na.rm=TRUE), 
							max(plotDat[indRow,], na.rm=TRUE)), ...)
					if(xAxt=='s')
						graphics::axis(1, cex.axis=cexAxis*(3/ncol(plotDat)), 
							mgp=mgpPlot*(3/ncol(plotDat)))
					if(yAxt=='s')
						graphics::axis(2, cex.axis=cexAxis*(3/ncol(plotDat)), 
							mgp=mgpPlot*(3/ncol(plotDat)))

					graphics::legend(lowerCorXY[1], lowerCorXY[2], mainText, 
						bty="n", text.col=corCol, text.font=1, 
						cex=corCex*(3/ncol(plotDat)))
					if(plotLoess)
						graphics::points(stats::loess.smooth(loessDat[,i], 
							loessDat[,j]), col=loessCol, type="l", 
							lwd= loessLwd,lty=loessLty)
				} else {
					graphics::par(mar=c(0,0,0,0))
					graphics::plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', 
						type = 'n', xaxt = 'n', yaxt = 'n')
					graphics::rect(par("usr")[1],par("usr")[3],par("usr")[2],
						par("usr")[4],col = "grey",border="grey")	
				}
			}
		}
	}
	
}

methods::setMethod("plot", "SummarizedExperiment", plot.InterestResult)

# Functions getting methods  
getAnnotation<-function(x){
	return(SummarizedExperiment::colData(x))
}

getRowData<-function(x){
	return(SummarizedExperiment::rowData(x))
}

addAnnotation<-function(x, sampleAnnotationType, sampleAnnotation){
	tmp=cbind(as.data.frame(SummarizedExperiment::colData(x)), 
		sampleAnnotation)
	colnames(tmp)[ncol(tmp)]= sampleAnnotationType
	SummarizedExperiment::colData(x)<-
		S4Vectors::DataFrame(tmp)
	return(x)
}


scaledRetention<-function(x){
	return(SummarizedExperiment::assays(x)$scaledRetention)
}

counts.InterestResults<-function(object){
	return(SummarizedExperiment::assays(object)$counts)
}
methods::setMethod("counts", "SummarizedExperiment", counts.InterestResults)
