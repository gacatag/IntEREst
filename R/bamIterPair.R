bamIterPair <- function(bf, isPairedDuplicate) {
	done <- FALSE
	if (!Rsamtools::isOpen( bf))
		open(bf)
	function() {
		if (done)
			return(NULL)
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(isPaired=TRUE,
				isDuplicate=isPairedDuplicate))
		yld=GenomicAlignments::readGAlignmentPairs(bf, 
			param=scParam)
		if (length(yld) == 0L) {
			close(bf)
			done <<- TRUE
			NULL
		} else yld
	}
}

