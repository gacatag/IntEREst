\name{updateInterestDfCol}
\alias{updateInterestDfCol}

\title{
 Updating contents of \code{@interestDf} of \code{InterestResult} objects
}
\description{
Updates the values in a single column of the \code{@interestDf} dataframe of 
\code{InterestResult} objects.
}
\usage{
updateInterestDfCol(x,  updateCol, value)
}
\arguments{
\item{x}{ 
	Object of type \code{InterestResult}.
}
\item{updateCol}{ 
Name or the number of the column in the \code{@interestDf} dataframe attribute 
of \code{x} to be updated with the new values. if the column name does not 
match to any column names of the \code{@interestDf} data frame the new values 
will be added as a new column before the first column with the read counts (the
column name would also be the defined \code{updateCol}).
}
\item{value}{ 
	The new Replacing values.
}
}
\value{
	Returns an object of type \code{InterestResult}.
}
\seealso{
	\code{\link{annotateU12}}
}

\author{Ali Oghabian}

\examples{

test<- mdsChr22Obj
# See the the frequency of each intron type annotation
table(interestDf(test)$int_type)

#Change the ambiguous (U12/U2) cases to U2
newIntType<- interestDf(test)$int_type
newIntType[newIntType=="U12/U2" & 
	!is.na(newIntType=="U12/U2")]<- "U2"

test<- updateInterestDfCol(test, updateCol="int_type", 
	value=newIntType)
#See the frequency of the updated intron type annotations
table(interestDf(test)$int_type)

#Adding a new column
test<- updateInterestDfCol(test, updateCol="new_column", 
	value=rep(NA, nrow(interestDf(test))) )
interestDf(test)[1,]

}

