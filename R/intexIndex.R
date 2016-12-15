intexIndex<-function(x, intExCol="int_ex", what="intron"){
	return(which(x@interestDf[,intExCol]==what))
}

