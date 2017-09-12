DEXSeqIntEREst <- function(x, design, reducedModel= ~ sample + intex, 
fitExpToVar, intExCol, geneIdCol, bpparam, silent=TRUE, ...){
	if(missing(bpparam)){
		bpparam<- BiocParallel::MulticoreParam(workers=1)
	}

#Check which of the three 'intron' 'exon' or 'intex' variables are used
	parIntExVar<-c("intron", "exon", "intex")
	matchParIntExVar<-match(all.vars(design),parIntExVar)
	matchParIntExVar<-matchParIntExVar[which(!is.na(matchParIntExVar))]
	intEx<- parIntExVar[matchParIntExVar]

	if(length(intEx)==1 & (intEx=="intron"|intEx=="exon")){
		intronExon<-intEx
	} else if(length(intEx)==1 & intEx=="intex"){
		intronExon<-c("intron","exon")
	} else {
		stop(
"One of the variables 'intron', 'exon' or 'intex' must be used in the design.")
	}


# Replace 'intron' in the formula with 'exon' as expected by DEXSeq
	design<- eval(parse(text=paste("substitute(", deparse(design), 
		", list(intron=quote(exon)))", sep="")))
	design<- eval(parse(text=paste("substitute(", deparse(design), 
		", list(intex=quote(exon)))", sep=""))) 	
	reducedModel<- eval(parse(text=paste("substitute(", deparse(reducedModel),
		", list(intron=quote(exon)))", sep="")))
	reducedModel<- eval(parse(text=paste("substitute(", deparse(reducedModel),
		", list(intex=quote(exon)))", sep="")))
	design<-eval(parse(text=design))
	reducedModel<-eval(parse(text=reducedModel))

# Extract the introns or exon (or both) corresponding rows from the object
	x<- x[which(rowData(x)[,intExCol] %in% intronExon),]

	sampleTable<- as.data.frame(colData(x))
	dexObj<- DEXSeq::DEXSeqDataSet( countData=counts(x), 
		sampleData=sampleTable,
		design= design,
		featureID=as.character(1:nrow(counts(x))), 
		groupID=as.character(rowData(x)[,geneIdCol]), ...)

	dexRes<- DEXSeq::DEXSeq(dexObj, BPPARAM=bpparam, 
		fitExpToVar=fitExpToVar, reducedModel=reducedModel, quiet=silent)
	return(dexRes)
}

