lfc<-function(x, fcType="edgeR", sampleAnnoCol=c(), sampleAnnotation=c(), 
	silent=TRUE, group=c(), rejection.region="doubletail", pseudoCnt=1, 
	log2=TRUE, ... )
{
	object=x
	if(as.character(class(object))=="SummarizedExperiment"){
		if(fcType=="edgeR"){
			testRes=exactTestInterest(x=object, sampleAnnoCol=sampleAnnoCol, 
				sampleAnnotation=sampleAnnotation, silent=TRUE, group=group,
				rejection.region="doubletail", ...)
			fcRes=testRes$table[,"logFC"]
			if(!log2)
				fcRes=2^fcRes
		} else if(fcType=="scaledRetention"){
			if(length(sampleAnnoCol)>0){
				colInd=!is.na(match(
					as.character(SummarizedExperiment::colData(x)[,
						sampleAnnoCol]),
					sampleAnnotation))
				group=as.character(SummarizedExperiment::colData(x)[colInd,
					sampleAnnoCol])
			}
			if(length(sampleAnnotation)==0& length(unique(group))==2){
				sampleAnnotation=unique(group)
			} else if(length(sampleAnnotation)!=2 | length(unique(group))!=2) {
				msg<-
'The sampleAnnotation parameter should be a vector of size 2 which cotains 
values from colData(x); e.g. if colData(x)[,sampleAnnoCol] =
c("test", "test", "ctrl","ctrl", ...), and the goal is to compare "test" and 
"ctrl" samples, sampleAnnotation should either be c("test","ctrl") or 
c("ctrl","test").'
				stop(msg)
			}
			if(length(group)==0){
				ind1= (1:ncol(scaledRetention(object)))[getAnnotation(object)[
					,sampleAnnoCol]==sampleAnnotation[1]]
				ind2= (1:ncol(scaledRetention(object)))[getAnnotation(object)[
					,sampleAnnoCol]==sampleAnnotation[2]]
			} else {
				ind1= (1:ncol(scaledRetention(object)))[
					group==sampleAnnotation[1]]
				ind2= (1:ncol(scaledRetention(object)))[
					group==sampleAnnotation[2]]
			}

			fcRes= (apply(scaledRetention(object)[,ind2], 1,mean, na.rm=TRUE)+
				pseudoCnt)/ (apply(scaledRetention(object)[,ind1], 1,mean, 
					na.rm=TRUE)+ pseudoCnt)
			if(log2)
				fcRes=log(fcRes, base=2)
		}
	return(fcRes)	
	}	
}

