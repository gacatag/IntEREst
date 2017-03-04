\name{updateRowDataCol}
\alias{updateRowDataCol}

\title{
 Updating contents of \code{rowData} of \code{SummarizedExperiment} objects
}
\description{
Updates the values in a single column of the \code{rowData} of 
\code{SummarizedExperiment} objects.
}
\usage{
updateRowDataCol(x,  updateCol, value)
}
\arguments{
\item{x}{ 
	Object of type \code{SummarizedExperiment}.
}
\item{updateCol}{ 
Name or the number of the column in the \code{rowData} of \code{x} to be 
updated with the new values. if the \code{updateCol} does not match to any 
column names it will be added as a new column.
}
\item{value}{ 
	The new Replacing values.
}
}
\value{
	Returns an object of type \code{SummarizedExperiment}.
}
\seealso{
	\code{\link{annotateU12}}
}

\author{Ali Oghabian}

\examples{

test<- mdsChr22Obj
# See the the frequency of each intron type annotation
table(rowData(test)$int_type)

#Change the ambiguous (U12/U2) cases to U2
newIntType<- rowData(test)$int_type
newIntType[newIntType=="U12/U2" & 
	!is.na(newIntType=="U12/U2")]<- "U2"
#Updating values
test<- updateRowDataCol(test, updateCol="int_type", 
	value=newIntType)
#See the frequency of the updated intron type annotations
table(rowData(test)$int_type)

#Adding a new column
test<- updateRowDataCol(test, updateCol="new_column", 
	value=rep(NA, nrow(rowData(test))) )
rowData(test)[1,]

}

