glmInterest<-function(x, design=c(), silent=TRUE, disp="common", coef=c(), contrast=NULL, ...){
	y=nread(x)
	if(length(design)==0){
		stop('Please provide a suitable design.')
	}
	
	if(disp=="tagwiseInitCommon" & length(disp)==1){
		yDisp <- edgeR::estimateGLMTagwiseDisp(y,design, 
			dispersion=edgeR::estimateGLMTrendedDisp(y,design))
		dispersionType=disp
	} else if(disp=="tagwiseInitTrended" & length(disp)==1){
		yDisp <- edgeR::estimateGLMTagwiseDisp(y,design, 
			dispersion=edgeR::estimateGLMCommonDisp(y,design))
		dispersionType=disp
	} else if(disp=="trended"& length(disp)==1){
		yDisp <- edgeR::estimateGLMTrendedDisp(y,design)
		dispersionType=disp
	} else if(disp=="common"& length(disp)==1){
		yDisp <- edgeR::estimateGLMCommonDisp(y,design)
		dispersionType=disp
	} else if(is.numeric(disp)){
		dispersionType="manualSet"
	}  else{
		stop('Argument disp should be either "common", "trended", "tagwiseInitCommon", "tagwiseInitTrended" or a number.')
	}

	fit <- edgeR::glmFit(y, design, dispersion=yDisp, ...)
	if(length(coef)==0)
		coef=ncol(fit$design)
	res <- edgeR::glmLRT(fit, coef=coef, contrast)
	if(!silent)
		print(edgeR::topTags(res))
	res$dispersionType=dispersionType
	res$dispersion=yDisp
	return(res)
}


