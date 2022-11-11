exactTestInterest<- function(x, sampleAnnoCol=c(), sampleAnnotation=c(),
	geneIdCol, silent=TRUE, group=c(), rejection.region="doubletail", 
	big.count=900, prior.count=0.125, disp="common", ...)
{

	if(length(group)==0 & length(sampleAnnoCol)>0 & length(sampleAnnotation)>0){
		colInd=!is.na(match(as.character(SummarizedExperiment::colData(x)[,
			sampleAnnoCol]), sampleAnnotation))
		group=as.character(SummarizedExperiment::colData(x)[colInd,
			sampleAnnoCol])
	} else if (length(group)>0){
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
			sampleAnnotation=unique(group)
		colInd=!is.na(match(group,sampleAnnotation))
	} else {
		stop(
'Either group or the sampleAnnotation and sampleAnnoCol parameters need to be
set.')
	}
	y <- edgeR::DGEList(counts=counts(x)[,colInd], group=group )
	if(is.character(disp)){
		if((is.na(match(disp, c("tagwise", "trended", "common", "genewise", 
			"auto")))) & length(disp)==1){
				stop(
'The disp parameter should be set "tagwise", "trended", "common", or 
"genewise"'
				)
		}
		if(disp!="genewise")
			dispTmp=edgeR::estimateDisp(y, ...)
	}
	if(disp=="tagwise" & length(disp)==1){
		dispersionType="tagwise"
		dispersion=dispTmp$tagwise.dispersion
	} else if(disp=="trended"& length(disp)==1){
		dispersionType="trended"
		dispersion=dispTmp$trended.dispersion
	} else if(disp=="common"& length(disp)==1){
		dispersionType="common"
		dispersion=dispTmp$common.dispersion
	} else if(disp=="genewise"& length(disp)==1){
		dispersionType="genewise"
		dispTmp=edgeR::estimateExonGenewiseDisp(counts(x)[,colInd], 
			geneID=SummarizedExperiment::rowData(x)[,geneIdCol], 
			group=group)
		dispersion=as.numeric(dispTmp[as.character(
			SummarizedExperiment::rowData(x)[,geneIdCol])])
	} else if(is.numeric(disp)){
		dispersionType="manualSet"
	}
	res <- edgeR::exactTest(y, dispersion=dispersion, pair=sampleAnnotation, 
		rejection.region=rejection.region, big.count=big.count, 
		prior.count=prior.count)
	if(!silent)
		print(edgeR::topTags(res))
	res$dispersionType=dispersionType
	res$dispersion=dispersion
	return(res)
}
