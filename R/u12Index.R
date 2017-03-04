u12Index<-function(x, intExCol="int_ex", intTypeCol="int_type"){
	return(which(SummarizedExperiment::rowData(x)[,intExCol]=="intron" & 
		SummarizedExperiment::rowData(x)[,intTypeCol]=="U12"))
}

