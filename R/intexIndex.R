intexIndex<-function(x, intExCol="int_ex", what="intron"){
	return(which(as.character(SummarizedExperiment::rowData(x)[,intExCol])==
		what))
}

