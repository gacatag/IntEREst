u12NbIndex<- function(
	x, intExCol="int_ex", intTypeCol="int_type", strandCol="strand", 
	geneIdCol="collapsed_transcripts", naUnstrand=FALSE) {
	object=x
	dat=interestDf(object)
	minInd=tapply(1:nrow(dat), dat[,geneIdCol], min)
	maxInd=tapply(1:nrow(dat), dat[,geneIdCol], max)
	u12Ind=u12Index(object, intExCol=intExCol, intTypeCol=intTypeCol)
	u12Gen=dat[u12Ind,geneIdCol]
	u12Str=dat[u12Ind,strandCol]
	# Analyzing positive strand first
	upEx=u12Ind-1
	upInt=u12Ind-2
	upInt[upInt<minInd[u12Gen]]=NA
	upEx[upEx<minInd[u12Gen]]=NA
	dnEx=u12Ind+1
	dnInt=u12Ind+2
	dnInt[dnInt>maxInd[u12Gen]]=NA
	dnEx[dnEx>maxInd[u12Gen]]=NA
	#Correcting for minus strand
	upInt[u12Str=="-"]=dnInt[u12Str=="-"]
	upEx[u12Str=="-"]=dnEx[u12Str=="-"]
	# if naUnstrand is FALSE, NA is returned for Unstranded (strand=='*'); 
	#it would be analyzed like '+' strand otherwise
	if(naUnstrand){
			upInt[u12Str=="*"]=NA
			upEx[u12Str=="*"]=NA
			dnInt[u12Str=="*"]=NA
			dnEx[u12Str=="*"]=NA
	}

	# Correct in case intron and exon don't match!
	
	upInt[dat[upInt,intExCol]=="exon" & !is.na(upInt) ]=NA
	upEx[dat[upEx,intExCol]=="intron" & !is.na(upInt)]=NA
	dnInt[dat[dnInt,intExCol]=="exon" & !is.na(dnInt)]=NA
	dnEx[dat[dnEx,intExCol]=="intron" & !is.na(dnInt)]=NA
	
	res=c()
	res$upIntron=upInt
	res$downIntron=dnInt
	res$upExon=upEx
	res$downExon=dnEx
	return(res)
}
