annotateU12<-function(
	pwmU12U2=c(), 
	pwmSsIndex=c(), 
	referenceChr, referenceBegin, referenceEnd, referenceIntronExon, 
	intronExon="intron",  matchWindowRelativeUpstreamPos=c() ,
	matchWindowRelativeDownstreamPos=c(), minMatchScore="80%", refGenome="", 	
	setNaAs="U2", annotateU12Subtype=TRUE)
{
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
	u12DonMatch=matchPwmToReference(pwm=pwmU12U2[[1]], 
		pwmSsIndex=unlist(pwmSsIndex[1]), referenceChr=referenceChr, 
		referencePos=referenceBegin, intronExon="intron", 
		referenceIntronExon=referenceIntronExon, 
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[1]), 
			unlist(matchWindowRelativeDownstreamPos[1])), 
		minMatchScore=minMatchScore[1], 
		refGenome=refGenome)
	u12DonCompMatch=
		matchPwmToReference(pwm=Biostrings::reverseComplement(pwmU12U2[[1]]), 
			pwmSsIndex=ncol(pwmU12U2[[1]])-unlist(pwmSsIndex[1])+1, 
			referenceChr=referenceChr, referencePos=referenceEnd, 
			intronExon="intron", referenceIntronExon=referenceIntronExon,
			matchWindowRelativePos=
				(-1)*c(unlist(matchWindowRelativeDownstreamPos[1]), 
					unlist(matchWindowRelativeUpstreamPos[1])), 
			minMatchScore=minMatchScore[1], 
			refGenome=refGenome)

	u12BpMatch= matchPwmToReference(pwm=pwmU12U2[[2]], 
		pwmSsIndex=unlist(pwmSsIndex[2]), referenceChr=referenceChr, 
		referencePos=referenceEnd, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[2]), 
			unlist(matchWindowRelativeDownstreamPos[2])), 
		minMatchScore=minMatchScore[2], 
		refGenome=refGenome)
	u12BpCompMatch= matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[2]]), 
		pwmSsIndex=ncol(pwmU12U2[[2]])-unlist(pwmSsIndex[2])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[2]), 
				unlist(matchWindowRelativeUpstreamPos[2])), 
		minMatchScore=minMatchScore[2], 
		refGenome=refGenome)

	u12AccMatch= matchPwmToReference(pwm=pwmU12U2[[3]], 
		pwmSsIndex=unlist(pwmSsIndex[3]), 
		referenceChr=referenceChr, 
		referencePos=referenceEnd, 
		intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[3]), 
			unlist(matchWindowRelativeDownstreamPos[3])), 
		minMatchScore=minMatchScore[3], 
		refGenome=refGenome)
	u12AccCompMatch= matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[3]]), 
		pwmSsIndex=ncol(pwmU12U2[[3]])-unlist(pwmSsIndex[3])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[3]), 
				unlist(matchWindowRelativeUpstreamPos[3])), 
		minMatchScore=minMatchScore[3], 
		refGenome=refGenome)

	u2DonMatch= matchPwmToReference(pwm=pwmU12U2[[4]], 
		pwmSsIndex=unlist(pwmSsIndex[4]), referenceChr=referenceChr, 
		referencePos=referenceBegin, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[4]), 
			unlist(matchWindowRelativeDownstreamPos[4])), 
		minMatchScore=minMatchScore[4], 
		refGenome=refGenome)
	u2DonCompMatch= matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[4]]), 
		pwmSsIndex=ncol(pwmU12U2[[4]])-unlist(pwmSsIndex[4])+1, 
		referenceChr=referenceChr, referencePos=referenceEnd, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[4]), 
				unlist(matchWindowRelativeUpstreamPos[4])), 
		minMatchScore=minMatchScore[4], 
		refGenome=refGenome)

	u2AccMatch= matchPwmToReference(pwm=pwmU12U2[[5]], 
		pwmSsIndex=unlist(pwmSsIndex[5]), referenceChr=referenceChr, 
		referencePos=referenceEnd, intronExon="intron", 
		referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=c(unlist(matchWindowRelativeUpstreamPos[5]), 
			unlist(matchWindowRelativeDownstreamPos[5])), 
		minMatchScore=minMatchScore[5], refGenome=refGenome)
	u2AccCompMatch= matchPwmToReference(
		pwm=Biostrings::reverseComplement(pwmU12U2[[5]]), 
		pwmSsIndex=ncol(pwmU12U2[[5]])-unlist(pwmSsIndex[5])+1, 
		referenceChr=referenceChr, referencePos=referenceBegin, 
		intronExon="intron", referenceIntronExon=referenceIntronExon,
		matchWindowRelativePos=
			(-1)*c(unlist(matchWindowRelativeDownstreamPos[5]), 
				unlist(matchWindowRelativeUpstreamPos[5])), 
		minMatchScore=minMatchScore[5], refGenome=refGenome)
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

		if (as.character(class(refGenome)) == "BSgenome"){
			donPosSeq= Biostrings::getSeq(refGenome, 
				names=referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"],
				start=
					referenceBegin[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"],
				end=referenceBegin[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"]+
					1,
				as.character=TRUE)

			accPosSeq= Biostrings::getSeq(refGenome, 
				names=referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"],
				start=
					referenceEnd[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"]-1,
				end=referenceEnd[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"],
				as.character=TRUE)

			donNegSeqTmp=Biostrings::getSeq(refGenome, 
				names=referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"],
				start=
					referenceEnd[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"]-1,
				end=referenceEnd[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"],
				as.character=TRUE)

			accNegSeqTmp=Biostrings::getSeq(refGenome, 
				names=referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"],
				start=
					referenceBegin[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"],
				end=referenceBegin[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"]+
					1,
				as.character=TRUE)
		} else if (as.character(class(refGenome)) == "DNAStringSet"){

			tmpGr<- GenomicRanges::GRanges (
					, 
					IRanges(,
						))

			tmpGr<- GenomicRanges::GRanges (
					referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="+"], 
					IRanges(referenceBegin[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+"],
						referenceBegin[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+"]+1))

			donPosSeq= as.character(Biostrings::getSeq(refGenome, 
				tmpGr))

			tmpGr<- GenomicRanges::GRanges (referenceChr[!is.na(res[,1])&
				res[,1]=="U12"&res[,2]=="+"], 
					IRanges(referenceEnd[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+"]-1,
						referenceEnd[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="+"]))

			accPosSeq= as.character(Biostrings::getSeq(refGenome, 
				tmpGr))

			tmpGr<- GenomicRanges::GRanges (
					referenceChr[!is.na(res[,1])&res[,1]=="U12"&res[,2]=="-"], 
					IRanges(referenceEnd[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-"]-1,
						referenceEnd[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-"]
						))

			donNegSeqTmp=as.character(Biostrings::getSeq(refGenome, 
				tmpGr))

			tmpGr<- GenomicRanges::GRanges (referenceChr[!is.na(res[,1])&
				res[,1]=="U12"&res[,2]=="-"], 
					IRanges(referenceBegin[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-"],
						referenceBegin[!is.na(res[,1])&res[,1]=="U12"&
							res[,2]=="-"]+ 1
						))

			accNegSeqTmp=as.character(Biostrings::getSeq(refGenome, 
				tmpGr))
		} else if (file.exists(refGenome)) {

			fa=seqinr::read.fasta(file = refGenome)

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
	donNegSeq= unlist(lapply(donNegSeqTmp, function(tmp) 
		as.character(Biostrings::reverseComplement(DNAString(tmp))) ))
	accNegSeq= unlist(lapply(accNegSeqTmp, function(tmp) 
		as.character(Biostrings::reverseComplement(DNAString(tmp))) ))
		

	res$u12_subtype=rep(NA,nrow(res))
	res$u12_subtype[which(!is.na(match(referenceIntronExon,intronExon)) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="+")[donPosSeq=="GT" & 
			accPosSeq=="AG"]]="GT-AG"
	res$u12_subtype[which(!is.na(match(referenceIntronExon,intronExon)) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="-")[donNegSeq=="GT" & 
			accNegSeq=="AG"]]="GT-AG"
	res$u12_subtype[which(!is.na(match(referenceIntronExon,intronExon)) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="+")[donPosSeq=="AT" & 
			accPosSeq=="AC"]]="AT-AC"
	res$u12_subtype[which(!is.na(match(referenceIntronExon,intronExon)) & 
		!is.na(res[,1])& res[,1]=="U12"& res[,2]=="-")[donNegSeq=="AT" & 
			accNegSeq=="AC"]]="AT-AC"

}

print(table(res))

return(res)

}
