eBayesInterest<- function(x, sampleAnnoCol=c(), sampleAnnotation=c(), 
	group=c(), design=c(), logBase=2, ...)
{

	if(length(group)==0 & length(sampleAnnoCol)>0 & length(sampleAnnotation)>0 &
		length(design)==0){
		colInd<- !is.na(match(x@sampleAnnotation[,sampleAnnoCol],
			sampleAnnotation))
		group<- x@sampleAnnotation[colInd,sampleAnnoCol]
	} else if (length(group)>0 & length(design)==0){
		#Check if sampleAnnotation param is set correctly
		if(length(sampleAnnotation)>0)
			if(length(which(is.na(match(sampleAnnotation,group))))>0 & 
				length(sampleAnnotation)>0)
				stop(
'The sampleAnnotation parameter should be a vector of size 2 which cotains
values from group; e.g. if group=c("test", "test", "ctrl","ctrl", ...), and the
goal is to compare "test" and "ctrl" samples, sampleAnnotation should either be
c("test","ctrl") or c("ctrl","test").')

		if(length(sampleAnnotation)==0)
			sampleAnnotation<- unique(group)
		colInd<- !is.na(match(group,sampleAnnotation))
	} else if(length(design)==0) {
		stop(
'Either group, desin or the sampleAnnotation and sampleAnnoCol parameters need 
to be set.')
	}

	
	if(length(design)==0)
		design<- stats::model.matrix(~factor(group))

	lmMat<- scaledRetention(subInterestResult(x, interestDfRow= 
		which(interestDf(x)[,"int_ex"]=="intron") ))
	if(!is.na(logBase))
		lmMat<- log(lmMat+1, base=2)

	y <- limma::lmFit(lmMat, design)
	ebayRes<- limma::eBayes(y, trend=TRUE, ...)

	return(ebayRes)
}
