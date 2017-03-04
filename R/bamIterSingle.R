bamIterSingle <- function(bf, scParam) {
	done <- FALSE
	if (!Rsamtools::isOpen( bf))
		open(bf)
	function() {
		if (done)
			return(NULL)

		yld=GenomicAlignments::readGAlignmentsList(bf, 
			param=scParam)
		if (length(yld) == 0L) {
			close(bf)
			done <<- TRUE
			NULL
		} else yld
	}
}
