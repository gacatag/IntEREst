
bamIterSingle <- function(bf, isSingleReadDuplicate) {
	done <- FALSE
	if (!Rsamtools::isOpen( bf))
		open(bf)
	function() {
		if (done)
			return(NULL)
		scParam=Rsamtools::ScanBamParam(
			what=Rsamtools::scanBamWhat()[c(1,
				3,5,8,13,9, 10, 6, 4, 14, 15)], 
			flag=Rsamtools::scanBamFlag(hasUnmappedMate=TRUE,
				isPaired=TRUE, 
				isDuplicate=isSingleReadDuplicate))
		yld=GenomicAlignments::readGAlignmentsList(bf, 
			param=scParam)
		if (length(yld) == 0L) {
			close(bf)
			done <<- TRUE
			NULL
		} else yld
	}
}
