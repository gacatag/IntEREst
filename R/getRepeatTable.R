getRepeatTable<-function(
	dbUser="genome",
	dbHost="genome-mysql.cse.ucsc.edu",
	ucscGenome="hg19",
	ucscTable="rmsk",
	minLength=0,
	repFamilyFil="Alu",
	repFamilyCol="repFamily",
	repChrCol="genoName",
	repBegCol="genoStart",
	repEndCol="genoEnd",
	repStrandCol="strand",
	repNameCol="repName",
	repClassCol="repClass")
{
# Calling required functions
dbConnect<- DBI::dbConnect
dbGetQuery<- DBI::dbGetQuery
dbDisconnect<- DBI::dbDisconnect
MySQL<- RMySQL::MySQL

# Connecting to the online shared database and getting the data
con <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db="hg19")
tblRmsk <- dbGetQuery(con, paste("select * from ", ucscTable, ";", sep=""))
tblRmskFilType=tblRmsk[!is.na(match(tblRmsk[,repFamilyCol],repFamilyFil)),]
tblRmskFilTypeFilLen=tblRmskFilType[(tblRmskFilType[,repEndCol]-tblRmskFilType[,repBegCol]+1)>minLength,]
dbDisconnect(con)

# OR FOLLOWING CODES ALSO RETREIVES THE RMSK TABLE
#library(rtracklayer)
#session <- rtracklayer::browserSession(browserName, url=browserUrl)
#GenomeInfoDb::genome(session) <- ucscGenome
#tblRmsk <- rtracklayer::getTable(rtracklayer::ucscTableQuery(session, track=ucscTrackName, table=ucscTable))

res=tblRmskFilTypeFilLen[,c(repChrCol, repBegCol, repEndCol, repStrandCol, repNameCol, repClassCol, repFamilyCol)]
colnames(res)[1:3]=c("chr","begin","end")
return(res)

}
