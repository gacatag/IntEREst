matchPwmToReference<-function(pwm=pwmU12db$pwmDonU12, pwmSsIndex=11, referenceChr, referencePos, referenceIntronExon, intronExon="intron", matchWindowRelativePos=c(NA,NA), minMatchScore="80%", refGenome=""){
	print(refGenome)
	# Check if reference genome sequence is of BSgenome type
	if (as.character(class(refGenome)) == "BSgenome"){

		# Check if only intron sequences should be checked for U12 annotation 
		if(length(intronExon)==1){
			# Check whether an upstream or downstream window is NOT defined for the PWM matching
			if((!is.na(matchWindowRelativePos[1])) & (!is.na(matchWindowRelativePos[2]))){
				seqIntEx=Biostrings::getSeq(refGenome, names=referenceChr[referenceIntronExon==intronExon], 
					start=(referencePos+(matchWindowRelativePos[1]))[referenceIntronExon==intronExon], 
					end=(referencePos+matchWindowRelativePos[2])[referenceIntronExon==intronExon],
					as.character=TRUE)
			# Check whether an upstream or downstream window is defined for the PWM matching
			} else {
				seqIntEx=Biostrings::getSeq(refGenome, names=referenceChr[referenceIntronExon==intronExon], 
					start=(referencePos-(pwmSsIndex-1))[referenceIntronExon==intronExon], 
					end=(referencePos+(ncol(pwm)-pwmSsIndex))[referenceIntronExon==intronExon],
					as.character=TRUE)
			}
		# Check if all (intron and exon) sequences should be checked for U12 annotation  
		} else if ( length(which(intronExon=="intron"))==1 & length(which(intronExon=="exon"))==1 ) {
			
			if((!is.na(matchWindowRelativePos[1])) & (!is.na(matchWindowRelativePos[2]))){
				# run if upstream or downstream windows are NOT defined for the PWM matching
				seqIntEx=Biostrings::getSeq(refGenome, names=referenceChr[referenceIntronExon==intronExon[1] | referenceIntronExon==intronExon[2]], 
					start=(referencePos+(matchWindowRelativePos[1]))[referenceIntronExon==intronExon[1] |
						referenceIntronExon==intronExon[2]], 
					end=(referencePos+matchWindowRelativePos[2])[referenceIntronExon==intronExon[1] | 
						referenceIntronExon==intronExon[2]],
					as.character=TRUE)
			} else {
				# run if upstream or downstream windows are defined for the PWM matching
				seqIntEx=Biostrings::getSeq(refGenome, names=referenceChr[referenceIntronExon==intronExon[1] | referenceIntronExon==intronExon[2]], 
					start=(referencePos-(pwmSsIndex-1))[referenceIntronExon==intronExon[1] | referenceIntronExon==intronExon[2]], 
					end=(referencePos+(ncol(pwm)-pwmSsIndex))[referenceIntronExon==intronExon[1] | 
						referenceIntronExon==intronExon[2]],
					as.character=TRUE)
			}


		} else {
			stop("Invalid intronExon parameter. Check ?matchPwmToReference for more information.")
		}

		matchInd= lapply(seqIntEx, function(tmp) { 
			test=Biostrings::matchPWM(pwm, tmp, min.score=minMatchScore, with.score=TRUE)
			return( unlist(GenomicRanges::ranges(test)))
			}) 
		lenMatchInd=sapply(matchInd,length)

	} else if (file.exists(refGenome)) {
		#Run if reference genome sequence is a fasta file
		fa=seqinr::read.fasta(file = refGenome)
		# Check if only intron sequences should be checked for U12 annotation
		if(length(intronExon)==1){
			# Check whether an upstream or downstream window is NOT defined for the PWM matching
			if((!is.na(matchWindowRelativePos[1])) & (!is.na(matchWindowRelativePos[2]))){

				seqIntEx=unlist(lapply( which(referenceIntronExon==intronExon), 
					function(temp) seqinr::getSequence(seqinr::getFrag(object=fa[as.character(referenceChr)[temp]], 
						begin=(referencePos+(matchWindowRelativePos[1]))[temp], 
						end=(referencePos+matchWindowRelativePos[2])[temp]), as.string=TRUE)))
			
			# Check whether an upstream or downstream window is defined for the PWM matching
			} else {

				seqIntEx=unlist(lapply( which(referenceIntronExon==intronExon), 
					function(temp) seqinr::getSequence(seqinr::getFrag(object=fa[as.character(referenceChr)[temp]], 
						begin=(referencePos-(pwmSsIndex-1))[temp], 
						end=(referencePos+(ncol(pwm)-pwmSsIndex))[temp]), as.string=TRUE)))

			}
		# Check if all (intron and exon) sequences should be checked for U12 annotation
		} else if ( length(which(intronExon=="intron"))==1 & length(which(intronExon=="exon"))==1 ) {
			# Check Whether an upstream or downstream window is NOT defined for the PWM matching
			if((!is.na(matchWindowRelativePos[1])) & (!is.na(matchWindowRelativePos[2]))){

				seqIntEx=unlist(lapply( which(referenceIntronExon==intronExon[1]| referenceIntronExon==intronExon[2]), 
					function(temp) seqinr::getSequence(seqinr::getFrag(object=fa[as.character(referenceChr)[temp]], 
						begin=(referencePos+(matchWindowRelativePos[1]))[temp], 
						end=(referencePos+matchWindowRelativePos[2])[temp]), as.string=TRUE)))

			# Check whether an upstream or downstream window is defined for the PWM matching
			} else {

				seqIntEx=unlist(lapply( which(referenceIntronExon==intronExon[1]| referenceIntronExon==intronExon[2]), 
					function(temp) seqinr::getSequence(seqinr::getFrag(object=fa[as.character(referenceChr)[temp]], 
						begin=(referencePos-(pwmSsIndex-1))[temp], 
						end=(referencePos+(ncol(pwm)-pwmSsIndex))[temp]), as.string=TRUE)))

			}
		} else {
			stop("Invalid intronExon parameter. Check ?matchPwmToReference for more information.")
		}

		matchInd= lapply(seqIntEx, function(tmp) { 
			test=Biostrings::matchPWM(pwm, tmp, min.score=minMatchScore, with.score=TRUE)
			return( unlist(GenomicRanges::ranges(test)))
			}) 
		lenMatchInd=sapply(matchInd,length)


	}


	pwmMatchBool=lenMatchInd>0

	pwmMatchScore=unlist(sapply(matchInd, function(tmp) 
		if(length(tmp)>0){
			return(tmp[1])
		} else {
			return(0)
		}))
	res=data.frame(pwmMatchBool=pwmMatchBool, pwmMatchScore=pwmMatchScore)
	return(res)

}