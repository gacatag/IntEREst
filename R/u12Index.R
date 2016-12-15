u12Index<-function(x, intExCol="int_ex", intTypeCol="int_type"){
	return(which(x@interestDf[,intExCol]=="intron" & x@interestDf[,intTypeCol]=="U12"))
}

