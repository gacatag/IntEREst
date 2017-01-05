# Class definition
#scaledRetentionColIndex -> scaledRetention
#readFreqColIndex -> readFreq
#scaledRetentionColIndex -> scaledRetention
methods::setClass("interestResult", 
	methods::representation(resultFiles="character", readFreq="matrix", scaledRetention="matrix", sampleNames="character",
		scaleLength="logical", scaleFragment="logical", sampleAnnotation="data.frame", interestDf="data.frame") 
)

#Constructor
interestResult <- function(resultFiles=c(), readFreq, scaledRetention, sampleNames, scaleLength, scaleFragment, 
	sampleAnnotation, interestDf){
	methods::new("interestResult", resultFiles=resultFiles, readFreq=readFreq, scaledRetention=scaledRetention, 
		sampleNames=sampleNames, scaleLength= scaleLength, scaleFragment= scaleFragment, sampleAnnotation=sampleAnnotation, interestDf=interestDf)
}

#view largedf
viewLargeDf<-function(datf, nrow=4, nchar=40, col.names=TRUE, row.names=FALSE, sep=' '){
if(col.names & !row.names)
	cat(substr(paste(unlist(lapply(colnames(datf), as.character)), collapse=sep), 1, nchar), "...", sep="")
if(col.names & row.names)
	cat(sep, substr(paste(unlist(lapply(colnames(datf), as.character)), collapse=sep), 1, nchar), "...", sep="")
if(!row.names)
	cat("\n", substr(paste(unlist(lapply(datf[1,], as.character)), collapse=sep), 1, nchar), "...", sep="")
if(row.names)
	cat("\n", substr(paste(unlist(lapply(rownames(datf)[1], as.character)), datf[1,], collapse=sep), 1, nchar), "...", sep="")

if(!row.names & nrow>1)
	lapply(2:min(nrow,nrow(datf)), function(n) cat("\n", substr(paste(unlist(lapply(datf[n,], as.character)), collapse=sep), 1, nchar), sep=""))
if(row.names & nrow>1)
	lapply(2:min(nrow,nrow(datf)), function(n) cat("\n", substr(paste(unlist(lapply(rownames(datf)[n], as.character)), unlist(lapply(datf[n,], as.character)), collapse=sep), 1, nchar), sep=""))
cat("\n...\n")
}

#show method
methods::setMethod("show", "interestResult",
	function(object){
		cat("  @resultFiles: ", strtrim(paste(object@resultFiles, collapse=", "), 40), "... ",
			length(object@resultFiles), " files", sep="" )
		cat("\n", "  @sampleNames: ", strtrim(paste(object@sampleNames, collapse=", "), 40), "... ",
			length(object@sampleNames), " samples", sep="" )
		cat("\n", "  @scaleLength: ", object@scaleLength, sep="" )
		cat("\n", "  @scaleFragment: ", object@scaleFragment, sep="" )
		cat("\n\n", "  @sampleAnnotation: \n", sep="")
		viewLargeDf(object@sampleAnnotation)
		cat(nrow(object@sampleAnnotation), " rows and ", ncol(object@sampleAnnotation), " columns.\n", sep="")
		cat("\n\n", "  @interestDf: \n", sep="")
		viewLargeDf(object@interestDf)
		cat(nrow(object@interestDf), " rows and ", ncol(object@interestDf), " columns.\n\n", sep="")
		cat("Use nread() to get the raw reteniton levels (Number fo mapped fragments) and \nscaledRetention() to get the scaled retention levels.\n")}
)

# plot method
plot.interestResult<-function(x, summary="none", subsetRows=NULL, what="scaled", intronExon="intron", logScaleBase=NULL, logPseudoCnt=1, plotLoess=TRUE, 
	loessCol="red",	loessLwd=1, loessLty=1, cexText=1, marPlot=c(2,2,2,2), mgpPlot=c(1, 1, 0), cexAxis=1, writeCor=TRUE, corCex=1, corMethod="pearson", 	corCol="grey63", upperCorXY=c("topleft", NULL), lowerCorXY= c("topleft", NULL), na.rm=TRUE, cex=1, sampleAnnoCol=c(), lowerPlot=FALSE, 
	upperPlot=TRUE, ...){
	object=x
	if(summary=="mean")
		summaryFun=mean
	if(summary=="median")
		summaryFun=stats::median
	if(is.null(subsetRows))
		subsetRows=1:nrow(object@interestDf)
	subDat=object@interestDf[subsetRows, ]

	indRow=!is.na(match(subDat$int_ex, intronExon))

	if(ncol(object@sampleAnnotation)==2){
		sampleGroups=object@sampleAnnotation[,2]
	} else if( ncol(object@sampleAnnotation)>2 & length(sampleAnnoCol)>0 & summary!="none"){
		sampleGroups=object@sampleAnnotation[,sampleAnnoCol]
	} else if (ncol(object@sampleAnnotation)>2 & length(sampleAnnoCol)==0 & summary!="none"){
		stop('Please define the sampleAnnoCol.')
	}

	if(summary=="none"&what=="scaled"){
		plotDat=object@scaledRetention[subsetRows,]
	}
	if(summary=="none"&what=="readFreq"){
		plotDat=object@readFreq[subsetRows,]
	}
	if(summary!="none"&what=="scaled"){
		plotList=tapply(1:ncol(object@scaledRetention), sampleGroups, function(x) apply(object@scaledRetention[subsetRows,x], 1, summaryFun) )
		plotList=plotList[match(unique(sampleGroups),names(plotList))]
		plotDat=matrix(unlist(plotList), ncol=length(unique(sampleGroups)), byrow=FALSE)
		colnames(plotDat)=names(plotList)
		plotDat=as.data.frame(plotDat)
	}
	if(summary!="none"&what=="readFreq"){
		plotList=tapply(1:ncol(object@readFreq), sampleGroups, function(x) apply(object@readFreq[subsetRows,x], 1, summaryFun) )
		plotList=plotList[match(unique(sampleGroups),names(plotList))]
		plotDat=matrix(unlist(plotList), ncol=length(unique(sampleGroups)), byrow=FALSE)
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
			mainText=paste("cor (",corMethod,") = ", format(stats::cor(loessDat[, 1], loessDat[,2], method=corMethod)
				, digits=2), sep="")
		graphics::par(mfrow=c(1,1))	
		graphics::par(cex.main=corCex*(3/ncol(plotDat)))
		graphics::par(mar=marPlot*(3/ncol(plotDat)))
		graphics::plot(plotDat[indRow, ], cex=cex*(3/ncol(plotDat)), main=mainText, ...)
		if(dim(plotDat)[2]==2)
			if(plotLoess)
				graphics::points(stats::loess.smooth(plotDat[indRow, 1], plotDat[indRow, 2]), col=loessCol, type="l", 
					lwd= loessLwd,lty=loessLty)
	}else{
		graphics::par(mfrow=c(ncol(plotDat),ncol(plotDat)))
		graphics::par(cex.main=corCex*(3/ncol(plotDat)))
		for(j in 1: ncol(plotDat)){
			for(i in 1:ncol(plotDat)){
				graphics::par(mar=marPlot*(3/ncol(plotDat)))
				if(i==j){
					graphics::plot(c(0, 1), c(0, 1), ann = FALSE, type = 'n', xaxt = 'n', yaxt = 'n')
					graphics::text(x = 0.5, y = 0.5, colnames(plotDat)[i], cex=cexText*(3/ncol(plotDat)))
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
						mainText=paste("cor = ", format(stats::cor(loessDat[, i], loessDat[,j], method=corMethod)
							, digits=2), sep="")

					graphics::plot(plotDat[indRow, i], plotDat[indRow, j], ylab="", xlab="", xaxt='n', yaxt='n', cex.axis=cexAxis*(3/ncol(plotDat)), 
					main="", cex=cex*(3/ncol(plotDat)), xlim=c(min(plotDat[indRow,], na.rm=TRUE), max(plotDat[indRow,], na.rm=TRUE)), 
					ylim=c(min(plotDat[indRow,], na.rm=TRUE), max(plotDat[indRow,], na.rm=TRUE)), ...)
					if(xAxt=='s')
						graphics::axis(3, cex.axis=cexAxis*(3/ncol(plotDat)), mgp=mgpPlot*(3/ncol(plotDat)))
					if(yAxt=='s')
						graphics::axis(4, cex.axis=cexAxis*(3/ncol(plotDat)), mgp=mgpPlot*(3/ncol(plotDat)))

					graphics::legend(upperCorXY[1], upperCorXY[2], mainText, bty="n", text.col=corCol, text.font=1, cex=corCex*(3/ncol(plotDat)))
					if(plotLoess)
						graphics::points(stats::loess.smooth(loessDat[,i], loessDat[,j]), col=loessCol, type="l", 
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
						mainText=paste("cor = ", format(stats::cor(loessDat[, i], loessDat[,j], method=corMethod)
							, digits=2), sep="")

					graphics::plot(plotDat[indRow, i], plotDat[indRow, j], ylab="", xlab="", xaxt='n', yaxt='n', cex.axis=cexAxis*(3/ncol(plotDat)), 
					main="", cex=cex*(3/ncol(plotDat)), xlim=c(min(plotDat[indRow,], na.rm=TRUE), max(plotDat[indRow,], na.rm=TRUE)),
					ylim=c(min(plotDat[indRow,], na.rm=TRUE), max(plotDat[indRow,], na.rm=TRUE)), ...)
					if(xAxt=='s')
						graphics::axis(1, cex.axis=cexAxis*(3/ncol(plotDat)), mgp=mgpPlot*(3/ncol(plotDat)))
					if(yAxt=='s')
						graphics::axis(2, cex.axis=cexAxis*(3/ncol(plotDat)), mgp=mgpPlot*(3/ncol(plotDat)))

					graphics::legend(lowerCorXY[1], lowerCorXY[2], mainText, bty="n", text.col=corCol, text.font=1, cex=corCex*(3/ncol(plotDat)))
					if(plotLoess)
						graphics::points(stats::loess.smooth(loessDat[,i], loessDat[,j]), col=loessCol, type="l", 
							lwd= loessLwd,lty=loessLty)
				} else {
					graphics::par(mar=c(0,0,0,0))
					graphics::plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
					graphics::rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey",border="grey")	
				}
			}
		}
	}
	
}

methods::setMethod("plot", "interestResult", plot.interestResult)

# Functions getting methods  
getAnnotation<-function(x){
	return(x@sampleAnnotation)
}

addAnnotation<-function(x, sampleAnnotationType, sampleAnnotation){
	x@sampleAnnotation=cbind(x@sampleAnnotation, sampleAnnotation)
	colnames(x@sampleAnnotation)[ncol(x@sampleAnnotation)]=sampleAnnotationType
	return(x)
}


scaledRetention<-function(x){
	return(x@scaledRetention)
}

nread<-function(x){
	return(x@readFreq)
}

interestDf<-function(x){
	return(x@interestDf)
}


