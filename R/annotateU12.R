annotateU12<-function(
	pwmU12U2=c(), 
	pwmSsIndex=c(), 
	referenceChr, referenceBegin, referenceEnd, referenceIntronExon, 
	intronExon="intron",  matchWindowRelativeUpstreamPos=c() ,
	matchWindowRelativeDownstreamPos=c(), minMatchScore="80%", refGenome="", 	
	setNaAs="U2", annotateU12Subtype=TRUE, includeMatchScores=FALSE, 
	ignoreHybrid=TRUE, filterReference)
{

	if(!missing(filterReference)){
		if(as.character(class(filterReference))=="SummarizedExperiment"){
			filterReference<- SummarizedExperiment::rowData(query)[,1:3]
			filterReference<- 
				GenomicRanges::GRanges (as.character(filterReference[,1]), 
					IRanges::IRanges(as.numeric(filterReference[,2]), 
						as.numeric(filterReference[,3])))
		} else if (as.character(class(filterReference))!="GRanges"){
			stop(paste("The filterReference parameter should be either of",
				"'SummarizedExperiment' or 'GRanges' type. For the former,", 
				"i.e. SummarizedExperiment, the first 3 columns of the", 
				"rowData must be: chr name, start and end of the coordinates.",
				sep="\n"))
		}
	}

	if(length(pwmU12U2)==0){
		pwmU12db<-0
		data("pwmU12db")
		pwmU12U2= 
			list(pwmU12db[[1]][,11:17],
			pwmU12db[[2]],pwmU12db[[3]][,38:40],pwmU12db[[4]][,11:17],
			pwmU12db[[5]][,38:40])
	}
	if(length(pwmSsIndex)==0){
		pwmSsIndex=list(indexDonU12=1, indexBpU12=1, indexAccU12=3, 
			indexDonU2=1, indexAccU2=3)
	}
	if(length(matchWindowRelativeUpstreamPos)==0)
		matchWindowRelativeUpstreamPos=c(NA,-29,NA,NA,NA)
	if(length(matchWindowRelativeDownstreamPos)==0)
		matchWindowRelativeDownstreamPos=c(NA,-9,NA,NA,NA)
	if(length(minMatchScore)==1)
		minMatchScore=rep(minMatchScore, 5)
	u12DonMatch=suppressWarnings( matchPwmToReference(pwm=pwmU12U2[[1]], 
		pwmSsIndex=unlist(pwmSsIndex[1]), referenceChr=referenceChr, 
		referencePos=referenceBegin, intronExon="intron", 
		referenceIntronExon=referenceIntronExon, 
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[1]), 
			unlist(matchWindowRelativeDownstreamPos[1])), 
		minMatchScore=minMatchScore[1], 
		refGenome=refGenome) )
	u12DonCompMatch=suppressWarnings(
		matchPwmToReference(pwm=Biostrings::reverseComplement(pwmU12U2[[1]]), 
			pwmSsIndex=ncol(pwmU12U2[[1]])-unlist(pwmSsIndex[1])+1, 
			referenceChr=referenceChr, referencePos=referenceEnd, 
			intronExon="intron", referenceIntronExon=referenceIntronExon,
			matchWindowRelativePos=
				(-1)*c(unlist(matchWindowRelativeDownstreamPos[1]), 
					unlist(matchWindowRelativeUpstreamPos[1])), 
			minMatchScore=minMatchScore[1], 
			refGenome=refGenome))

	u12BpMatch= suppressWarnings(matchPwmToReference(pwm=pwmU12U2[[2]], 
		pwmSsIndex=unlist(pwmSsIndex[2]), referenceChr=referenceChr, 
		referencePos=referenceEnd, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[2]), 
			unlist(matchWindowRelativeDownstreamPos[2])), 
		minMatchScore=minMatchScore[2], 
		refGenome=refGenome))
	u12BpCompMatch= suppressWarnings(matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[2]]), 
		pwmSsIndex=ncol(pwmU12U2[[2]])-unlist(pwmSsIndex[2])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[2]), 
				unlist(matchWindowRelativeUpstreamPos[2])), 
		minMatchScore=minMatchScore[2], 
		refGenome=refGenome))

	u12AccMatch= suppressWarnings(matchPwmToReference(pwm=pwmU12U2[[3]], 
		pwmSsIndex=unlist(pwmSsIndex[3]), 
		referenceChr=referenceChr, 
		referencePos=referenceEnd, 
		intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[3]), 
			unlist(matchWindowRelativeDownstreamPos[3])), 
		minMatchScore=minMatchScore[3], 
		refGenome=refGenome))
	u12AccCompMatch= suppressWarnings(matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[3]]), 
		pwmSsIndex=ncol(pwmU12U2[[3]])-unlist(pwmSsIndex[3])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[3]), 
				unlist(matchWindowRelativeUpstreamPos[3])), 
		minMatchScore=minMatchScore[3], 
		refGenome=refGenome))

	u2DonMatch= suppressWarnings(matchPwmToReference(pwm=pwmU12U2[[4]], 
		pwmSsIndex=unlist(pwmSsIndex[4]), referenceChr=referenceChr, 
		referencePos=referenceBegin, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[4]), 
			unlist(matchWindowRelativeDownstreamPos[4])), 
		minMatchScore=minMatchScore[4], 
		refGenome=refGenome))
	u2DonCompMatch= suppressWarnings(matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[4]]), 
		pwmSsIndex=ncol(pwmU12U2[[4]])-unlist(pwmSsIndex[4])+1, 
		referenceChr=referenceChr, referencePos=referenceEnd, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[4]), 
				unlist(matchWindowRelativeUpstreamPos[4])), 
		minMatchScore=minMatchScore[4], 
		refGenome=refGenome))

	u2AccMatch= suppressWarnings(matchPwmToReference(pwm=pwmU12U2[[5]], 
		pwmSsIndex=unlist(pwmSsIndex[5]), referenceChr=referenceChr, 
		referencePos=referenceEnd, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[5]), 
			unlist(matchWindowRelativeDownstreamPos[5])), 
		minMatchScore=minMatchScore[5], refGenome=refGenome))
	u2AccCompMatch= suppressWarnings(matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[5]]), 
		pwmSsIndex=ncol(pwmU12U2[[5]])-unlist(pwmSsIndex[5])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[5]), 
				unlist(matchWindowRelativeUpstreamPos[5])), 
		minMatchScore=minMatchScore[5], refGenome=refGenome))
	resTmp= rep(NA,nrow(u12DonMatch))
	u12Bool= u12DonMatch[,1]& u12BpMatch[,1]& u12AccMatch[,1]
	u2Bool= u2DonMatch[,1]& u2AccMatch[,1]
	u12CompBool= u12DonCompMatch[,1]& u12BpCompMatch[,1]& u12AccCompMatch[,1]
	u2CompBool= u2DonCompMatch[,1]& u2AccCompMatch[,1]

	resTmp[u12Bool]="U12"
	resTmp[u2Bool]="U2"
	resTmp[u12Bool & u2Bool]="U12/U2"
	resTmp[u12CompBool]="U12"
	resTmp[u2CompBool]="U2"
	resTmp[u12CompBool & u2CompBool]="U12/U2"

	resTmp[is.na(resTmp)]=setNaAs

	strandMatchTmp=rep(NA, nrow(u12DonMatch))
	strandMatchTmp[u12Bool|u2Bool]="+"
	strandMatchTmp[u12CompBool|u2CompBool]="-"

	res=data.frame( int_type=rep(NA,length(referenceIntronExon)), 
		strandMatch=rep(NA,length(referenceIntronExon)))
	res[!is.na(match(referenceIntronExon,intronExon)),1]=resTmp
	res[!is.na(match(referenceIntronExon,intronExon)),2]=strandMatchTmp

	if(annotateU12Subtype & (length(which(res[,1]=="U12"))>0)){
		donPosSeq<- rep(NA, nrow(res))
		accPosSeq<- rep(NA, nrow(res))
		donNegSeq<- rep(NA, nrow(res))
		accNegSeq<- rep(NA, nrow(res))
		donNegSeqTmp<- rep(NA, nrow(res))
		accNegSeqTmp<- rep(NA, nrow(res))
		if (as.character(class(refGenome)) == "BSgenome"){
			if(length(which(res[,2]=="+"))>0){
				donPosSeq= Biostrings::getSeq(refGenome, 
					names=referenceChr[which(!is.na(res[,1])&
						res[,1]=="U12"&res[,2]=="+")],
					start=
						referenceBegin[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+")],
					end=
						referenceBegin[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+")]+
						1,
					as.character=TRUE)

				accPosSeq= Biostrings::getSeq(refGenome, 
					names=referenceChr[which(!is.na(res[,1])&
						res[,1]=="U12"&res[,2]=="+")],
					start=
						referenceEnd[which(!is.na(res[,1])&
							res[,1]=="U12"&res[,2]=="+")]-1,
					end=
						referenceEnd[which(!is.na(res[,1])&
							res[,1]=="U12"&res[,2]=="+")],
					as.character=TRUE)
			} 
			if(length(which(res[,2]=="-"))>0) {

				donNegSeqTmp=Biostrings::getSeq(refGenome, 
					names=referenceChr[which(!is.na(res[,1])&res[,1]=="U12"&
						res[,2]=="-")],
					start=
						referenceEnd[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-")]-1,
					end=
						referenceEnd[which(!is.na(res[,1])&
							res[,1]=="U12"&res[,2]=="-")],
					as.character=TRUE)

				accNegSeqTmp=Biostrings::getSeq(refGenome, 
					names=referenceChr[which(!is.na(res[,1])&res[,1]=="U12"&
						res[,2]=="-")],
					start=
						referenceBegin[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-")],
					end=
							referenceBegin[which(!is.na(res[,1])&
							res[,1]=="U12"&	res[,2]=="-")]+ 1,
					as.character=TRUE)
			}
		} else if (as.character(class(refGenome)) == "DNAStringSet"){

			if(length(which(res[,2]=="+"))>0){
				tmpGr<- GenomicRanges::GRanges (
						referenceChr[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+")], 
						IRanges(referenceBegin[which(!is.na(res[,1])& 
							res[,1]=="U12"& res[,2]=="+")],
							referenceBegin[which(!is.na(res[,1])&
								res[,1]=="U12"& res[,2]=="+")]+1))

				donPosSeq= as.character(Biostrings::getSeq(refGenome, 
					tmpGr))

				tmpGr<- GenomicRanges::GRanges (referenceChr[which(
					!is.na(res[,1])& res[,1]=="U12"&res[,2]=="+")], 
						IRanges(referenceEnd[which(!is.na(res[,1])&
								res[,1]=="U12"& res[,2]=="+")]-1,
							referenceEnd[which(!is.na(res[,1])&res[,1]=="U12"&
								res[,2]=="+")]))

				accPosSeq= as.character(Biostrings::getSeq(refGenome, 
					tmpGr))
			} 
			if (length(which(res[,2]=="-"))>0){
				tmpGr<- GenomicRanges::GRanges (
						referenceChr[which(!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-")], 
						IRanges(referenceEnd[which(!is.na(res[,1])&
								res[,1]=="U12"& res[,2]=="-")]-1,
							referenceEnd[which(!is.na(res[,1])&res[,1]=="U12"&
								res[,2]=="-")]
							))

				donNegSeqTmp=as.character(Biostrings::getSeq(refGenome, 
					tmpGr))

				tmpGr<- GenomicRanges::GRanges (referenceChr[which(
					!is.na(res[,1])& res[,1]=="U12"&res[,2]=="-")], 
						IRanges(referenceBegin[which(!is.na(res[,1])&
								res[,1]=="U12"& res[,2]=="-")],
							referenceBegin[which(!is.na(res[,1])&
								res[,1]=="U12"& res[,2]=="-")]+ 1
							))

				accNegSeqTmp=as.character(Biostrings::getSeq(refGenome, 
					tmpGr))
			}
		} else if (file.exists(refGenome)) {

			fa=seqinr::read.fasta(file = refGenome)
			if(length(which(res[,2]=="+"))>0){
				donPosSeq= unlist(lapply( 
					which(!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"), 
						function(temp) 
							paste(unlist(seqinr::getSequence(seqinr::getFrag(
								object=fa[as.character(referenceChr)[temp]], 
								begin=referenceBegin[temp], 
								end=referenceBegin[temp]+1, as.string=TRUE))),
								collapse="")))


				accPosSeq= unlist(lapply( 
					which(!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"), 
						function(temp) 
							paste(unlist(seqinr::getSequence(seqinr::getFrag(
								object=fa[as.character(referenceChr)[temp]], 
								begin=referenceEnd[temp]-1, 
								end=referenceEnd[temp], as.string=TRUE))),
								collapse="")))
			}
			if(length(which(res[,2]=="-"))>0){
				donNegSeqTmp=unlist(lapply( 
					which(!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"), 
						function(temp) 
							paste(unlist(seqinr::getSequence(seqinr::getFrag(
								object=fa[as.character(referenceChr)[temp]], 
								begin=referenceEnd[temp]-1, 
								end=referenceEnd[temp], as.string=TRUE))),
								collapse="")))

				accNegSeqTmp=unlist(lapply( 
					which(!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"), 
					function(temp) 
						paste(unlist(seqinr::getSequence(seqinr::getFrag(
							object=fa[as.character(referenceChr)[temp]], 
							begin=referenceBegin[temp], 
							end=referenceBegin[temp]+1, as.string=TRUE))),
							collapse="")))
			}
		}
	if(length(which(!is.na(donNegSeqTmp)))>0)
		donNegSeq= unlist(lapply(donNegSeqTmp, function(tmp) 
			as.character(Biostrings::reverseComplement(DNAString(tmp))) ))
	if(length(which(!is.na(accNegSeqTmp)))>0)
		accNegSeq= unlist(lapply(accNegSeqTmp, function(tmp) 
			as.character(Biostrings::reverseComplement(DNAString(tmp))) ))
		

	res$u12_subtype=rep(NA,nrow(res))
	res$u12_subtype[which((referenceIntronExon %in% intronExon) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="+")[which(donPosSeq=="GT" & 
			accPosSeq=="AG")]]="GT-AG"
	res$u12_subtype[which((referenceIntronExon %in% intronExon) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="-")[which(donNegSeq=="GT" & 
			accNegSeq=="AG")]]="GT-AG"
	res$u12_subtype[which((referenceIntronExon %in% intronExon) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="+")[which(donPosSeq=="AT" & 
			accPosSeq=="AC")]]="AT-AC"
	res$u12_subtype[which((referenceIntronExon %in% intronExon) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="-")[which(donNegSeq=="AT" & 
			accNegSeq=="AC")]]="AT-AC"
	if(ignoreHybrid)
		res$int_type[which(is.na(res$u12_subtype) & res$int_type=="U12")]<- 
			setNaAs

}

if (includeMatchScores){
	u12_donScoreTmp<- u12DonMatch$pwmMatchScore
	u12_donScoreTmp[u12CompBool]<- u12DonCompMatch$pwmMatchScore[u12CompBool]
	u12_bpScoreTmp<- u12BpMatch$pwmMatchScore
	u12_bpScoreTmp[u12CompBool]<- u12BpCompMatch$pwmMatchScore[u12CompBool]
	u12_accScoreTmp<- u12AccMatch$pwmMatchScore
	u12_accScoreTmp[u12CompBool]<- u12AccCompMatch$pwmMatchScore[u12CompBool]

	u2_donScoreTmp<- u2DonMatch$pwmMatchScore
	u2_donScoreTmp[u2CompBool]<- u2DonCompMatch$pwmMatchScore[u2CompBool]
	u2_accScoreTmp<- u2AccMatch$pwmMatchScore
	u2_accScoreTmp[u2CompBool]<- u2AccCompMatch$pwmMatchScore[u2CompBool]
	
	u12_donScore<- rep(NA, nrow(res))
	u12_donScore[which(!is.na(match(referenceIntronExon,intronExon)))]<- 
		u12_donScoreTmp
	u12_bpScore<- rep(NA, nrow(res))
	u12_bpScore[which(!is.na(match(referenceIntronExon,intronExon)))]<- 
		u12_bpScoreTmp
	u12_accScore<- rep(NA, nrow(res))
	u12_accScore[which(!is.na(match(referenceIntronExon,intronExon)))]<- 
		u12_accScoreTmp

	u2_donScore<- rep(NA, nrow(res))
	u2_donScore[which(!is.na(match(referenceIntronExon,intronExon)))]<- 
		u2_donScoreTmp
	u2_accScore<- rep(NA, nrow(res))
	u2_accScore[which(!is.na(match(referenceIntronExon,intronExon)))]<- 
		u2_accScoreTmp

	res<- 
		cbind(res, u12_donScore, u12_bpScore, u12_accScore, 
			u2_donScore, u2_accScore)


}

if(!missing(filterReference)){
	refGr<- 
		GenomicRanges::GRanges (as.character(referenceChr), 
			IRanges::IRanges(as.numeric(referenceBegin), 
				as.numeric(referenceEnd)))

	tmpMap<- GenomicRanges::findOverlaps(refGr, 
			filterReference, type= "equal")
	res$int_type[-c(unique(S4Vectors::queryHits(tmpMap)))]<-NA
	
}

print(table(res$int_type))

return(res)

}
