u12Index<-
	function(x, intExCol="int_ex", intTypeCol="int_type", intronExon="intron"){
	return(which(SummarizedExperiment::rowData(x)[,intExCol]==intronExon & 
		SummarizedExperiment::rowData(x)[,intTypeCol]=="U12"))
}

