qlfInterest<-function(x, design=c(), silent=TRUE, disp="common", coef=c(),
	contrast=NULL, poisson.bound=TRUE, ...){
	y=counts(x)
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
		msg<-
'Argument disp should be either "common", "trended", "tagwiseInitCommon", 
"tagwiseInitTrended" or a number.'
		stop(msg)
	}

	fit <- edgeR::glmQLFit(y, design, dispersion=yDisp, ...)
	if(length(coef)==0)
		coef=ncol(fit$design)
	res <- edgeR::glmQLFTest(fit, coef=coef, contrast, 
		poisson.bound=poisson.bound)
	if(!silent)
		print(edgeR::topTags(res))
	res$dispersionType=dispersionType
	res$dispersion=yDisp
	return(res)
}

