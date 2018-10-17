buildSsTypePwms<-function(
	cexSeqLogo=1,
	pdfWidth=35,
	pdfHeight=10,
	tmpDir="./",
	u12dbSpecies="Homo_sapiens",
	pwmSource="U12DB",
	splicerackSsLinks=list(
		U12_AT_AC_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.25", 
		U12_AT_AC_branchpoint=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.26", 
		U12_AT_AC_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.29", 
		U12_GT_AG_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.22", 
		U12_GT_AG_branchpoint=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.27", 
		U12_GT_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.21", 
		U2_GC_AG_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.24", 
		U2_GC_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.30", 
		U2_GT_AG_donor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.23", 
		U2_GT_AG_acceptor=
			"http://katahdin.mssm.edu/splice/out/9606_logo_file.28"),
	u12dbLink="ftp://genome.imim.es/pub/software/u12/u12db_v1_0.sql.gz",
	u12dbDbName="u12db", u12dbDropDb=TRUE, pdfFileSeqLogos="",
	removeTempFiles=TRUE, ...)
{

	if (pwmSource=="SpliceRack"){
#from http://katahdin.mssm.edu/splice/splice_matrix.cgi?database=spliceNew
		utils::download.file(unlist(splicerackSsLinks[1]), 
			paste(tmpDir, "AT_AC_U12_Don.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[2]), 
			paste(tmpDir, "AT_AC_U12_BP.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[3]), 
			paste(tmpDir, "AT_AC_U12_Acc.txt", sep="/"))
		
		utils::download.file(unlist(splicerackSsLinks[8]), 
			paste(tmpDir, "GC_AG_U2_Acc.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[7]), 
			paste(tmpDir, "GC_AG_U2_Don.txt", sep="/"))

		utils::download.file(unlist(splicerackSsLinks[4]), 
			paste(tmpDir, "GT_AG_U12_Don.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[5]), 
			paste(tmpDir, "GT_AG_U12_BP.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[6]), 
			paste(tmpDir, "GT_AG_U12_Acc.txt", sep="/"))

		utils::download.file(unlist(splicerackSsLinks[9]), 
			paste(tmpDir, "GT_AG_U2_Don.txt", sep="/"))
		utils::download.file(unlist(splicerackSsLinks[10]), 
			paste(tmpDir, "GT_AG_U2_Acc.txt", sep="/"))


		donAtacU12=read.table(paste(tmpDir, "AT_AC_U12_Don.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		bpAtacU12=read.table(paste(tmpDir, "AT_AC_U12_BP.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		accAtacU12=read.table(paste(tmpDir, "AT_AC_U12_Acc.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')

		accGcagU2=read.table(paste(tmpDir, "GC_AG_U2_Acc.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		donGcagU2=read.table(paste(tmpDir, "GC_AG_U2_Don.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')

		donGtagU12=read.table(paste(tmpDir, "GT_AG_U12_Don.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		bpGtagU12=read.table(paste(tmpDir, "GT_AG_U12_BP.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		accGtagU12=read.table(paste(tmpDir, "GT_AG_U12_Acc.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')

		donGtagU2=read.table(paste(tmpDir, "GT_AG_U2_Don.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')
		accGtagU2=read.table(paste(tmpDir, "GT_AG_U2_Acc.txt", sep="/"), 
			skip=7, header=TRUE, stringsAsFactors=FALSE, sep='\t')



		rowDifDonU12<- nrow(donGtagU12)-nrow(donAtacU12)
		rowDifAccU12<- nrow(accGtagU12)-nrow(accAtacU12)
		rowDifBpU12<- nrow(bpGtagU12)-nrow(bpAtacU12)

		rowDifDonU2<- nrow(donGtagU2)-nrow(donGcagU2)
		rowDifAccU2<- nrow(accGtagU2)-nrow(accGcagU2)

		if(rowDifDonU12>0){
			donAtacU12<- rbind(donAtacU12, matrix(rep(.25, 4*rowDifDonU12), 
				nrow=rowDifDonU12) )
		} else if(rowDifDonU12<0){
			donGtagU12<- rbind(donGtagU12, matrix(rep(.25, 4*rowDifDonU12), 
				nrow=abs(rowDifDonU12)) )
		}

		if(rowDifAccU12>0){
			accAtacU12<- rbind(accAtacU12, matrix(rep(.25, 4*rowDifAccU12), 
				nrow=rowDifAccU12) )
		} else if(rowDifAccU12<0){
			accGtagU12<- rbind(accGtagU12, matrix(rep(.25, 4*rowDifAccU12), 
				nrow=abs(rowDifAccU12)) )
		}

		if(rowDifBpU12>0){
			bpAtacU12<- rbind(bpAtacU12, matrix(rep(.25, 4*rowDifBpU12), 
				nrow=rowDifBpU12) )
		} else if(rowDifBpU12<0){
			bpGtagU12<- rbind(bpGtagU12, matrix(rep(.25, 4*rowDifBpU12), 
				nrow=abs(rowDifBpU12)) )
		}
		if(rowDifDonU2>0){
			donGcagU2<- rbind(donGcagU2, matrix(rep(.25, 4*rowDifDonU2), 
				nrow=rowDifDonU2) )
		} else if(rowDifDonU2<0){
			donGtagU2<- rbind(donGtagU2, matrix(rep(.25, 4*rowDifDonU2), 
				nrow=abs(rowDifDonU2)) )
		}
		if(rowDifAccU2>0){
			accGcagU2<- rbind(accGcagU2, matrix(rep(.25, 4*rowDifAccU2), 
				nrow=rowDifAccU2) )
		} else if(rowDifAccU2<0){
			accGtagU2<- rbind(accGtagU2, matrix(rep(.25, 4*rowDifAccU2), 
				nrow=abs(rowDifAccU2)) )
		}


		donU12=t((donAtacU12+donGtagU12)/2)
		accU12=t((accAtacU12+accGtagU12)/2)
		bpU12=t((bpAtacU12+bpGtagU12)/2)

		donU2=t((donGtagU2+donGcagU2)/2)
		accU2=t((accGtagU2+accGcagU2)/2)
	
		# Plotting sequence logos

		if(pdfFileSeqLogos!=""){
			grDevices::pdf(pdfFileSeqLogos, width=pdfWidth, height=pdfHeight)
# Codes borrowed from 
#http://grokbase.com/t/r/bioconductor/1096ds3bnk/bioc-drawing-multiple-sequence
#-logos-with-seqlogo#2010090767dkq0ngtwg3fbx8cjx5mxwsn0 
			mySeqLogo = seqLogo::seqLogo

			bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | 
				sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
			body(mySeqLogo)[bad] = NULL


			#Defining dimensions
			ncol <- 3
			nrow <- 2

			#Creating plot

			grid::grid.newpage()

			listPwm=list(donU12, bpU12, accU12, donU2, c(), accU2)

			for(i in 1:nrow){
				for(j in 1:ncol){
					vp <- grid::viewport(x = (j-1)/ncol, y = 1-(i-1)/nrow, 
						w = 1/ncol, h = 1/nrow, just = c("left", "top"), 
						gp=grid::gpar(cex=cexSeqLogo))
					grid::pushViewport(vp)

					if(i!=2|j!=2){
						ind=(i-1)*3+j
						mySeqLogo(listPwm[[ind]])
						grid::grid.text(sprintf(c("U12 donor", 
							"U12 branch point", "U12 acceptor", "U2 donor", 
							"", "U2 acceptor" )[ind], row, col), x=0.5, y=0.95,
							just="top")
					}
					grid::upViewport()
				}
			}
			dev.off()
		}

		res=list(pwmDonorU12=Biostrings::unitScale(donU12), 
			pwmBpU12=Biostrings::unitScale(bpU12), 
			pwmAccU12=Biostrings::unitScale(accU12), 
			pwmDonU2=Biostrings::unitScale(donU2), 
			pwmAccU2=Biostrings::unitScale(accU2))

		if(removeTempFiles){
			file.remove(c(paste(tmpDir, "AT_AC_U12_Don.txt", sep="/"), 
				paste(tmpDir, "AT_AC_U12_BP.txt", sep="/"), 
				paste(tmpDir, "AT_AC_U12_Acc.txt", sep="/"), 
				paste(tmpDir, "GC_AG_U2_Acc.txt", sep="/"), 
				paste(tmpDir, "GC_AG_U2_Don.txt", sep="/"), 
				paste(tmpDir, "GT_AG_U12_Don.txt", sep="/"), 
				paste(tmpDir, "GT_AG_U12_BP.txt", sep="/"), 
				paste(tmpDir, "GT_AG_U12_Acc.txt", sep="/"),
				paste(tmpDir, "GT_AG_U2_Don.txt", sep="/"), 
				paste(tmpDir, "GT_AG_U2_Acc.txt", sep="/")))
		}

	} else if (pwmSource=="U12DB"){
		# Calling required functions
		dbConnect<- DBI::dbConnect
		dbSendQuery<- DBI::dbSendQuery
		dbGetQuery<- DBI::dbGetQuery
		dbDisconnect<- DBI::dbDisconnect
		MySQL<- RMySQL::MySQL

		# Downloading U12DB
		utils::download.file(u12dbLink, paste(tmpDir,"u12db.sql.gz",sep="/"))
		scriptsTmp=scan(gzfile(paste(tmpDir,"u12db.sql.gz",sep="/")), 
			what="character", sep='\n')

		dbc=dbConnect(MySQL(), ...)
		dbSendQuery(dbc, paste("CREATE DATABASE IF NOT EXISTS ", u12dbDbName, 
			";", sep=""))
		dbDisconnect(dbc)

		dbc=dbConnect( MySQL(), dbname = u12dbDbName, ...)

		scriptsTmp=scriptsTmp[-grep("^-",scriptsTmp)]
		scriptsTmp=scriptsTmp[-grep("\\/\\*",scriptsTmp)]

		indBeg= grep("(CREATE|INSERT|LOCK TABLES|UNLOCK TABLES|DROP TABLE)",
			scriptsTmp)
		indEnd=indBeg-1
		indBeg=c(1,indBeg)
		indBeg=indBeg[-length(indBeg)]
		if(length(scriptsTmp)!=indEnd[length(indEnd)]){
			indBeg=c(indBeg,length(scriptsTmp))
			indEnd=c(indEnd,length(scriptsTmp))
		}


		for(cnt in 1:length(indEnd)){
			test= paste(scriptsTmp[unique(indBeg[cnt]:indEnd[cnt])], 
				collapse="\n")
			dbSendQuery(dbc, test)
		}

		selCol= c("intron.gene_id", "intron.type", "splice_site.site_type", 
			"splice_site.seq", "splice_site.position_bp", 
			"splice_site.site_type")

		donor=dbGetQuery(dbc, paste('select ', paste(selCol, collapse=","), 
			' from intron,splice_site where intron.species_id="', u12dbSpecies,
			'" AND intron.donor_ss=splice_site.id', sep=""))
		selCol= c("intron.gene_id", "intron.type", "splice_site.site_type", 
			"splice_site.seq", "splice_site.position_bp", 
			"splice_site.position_ppt", "splice_site.site_type")
		acceptor=dbGetQuery(dbc, paste('select ', paste(selCol, collapse=","),
			' from intron,splice_site where intron.species_id="', u12dbSpecies,
			'" AND intron.acceptor_ss=splice_site.id', sep=""))

		bpStr= unlist(lapply(1:nrow(acceptor), 
			function(i) 
				substr(x=acceptor[i,"seq"], 
					start=nchar(acceptor[i,"seq"])+
						acceptor[i,"position_bp"]-3-9, 
					stop=nchar(acceptor[i,"seq"])+acceptor[i,"position_bp"]-4)))

		dbDisconnect(dbc)

		bpStrU12=bpStr[acceptor[,"type"]=="U12"]
		accStrU12=acceptor[acceptor[,"type"]=="U12","seq"]
		donStrU12=donor[donor[,"type"]=="U12","seq"]

		accStrU2=acceptor[acceptor[,"type"]=="U2","seq"]
		donStrU2=donor[donor[,"type"]=="U2","seq"]
#Filter Sequences includng characters other than A, C, G and T
		donStrU12<- donStrU12[grep("^[A|C|G|T]+$",donStrU12)]
		bpStrU12<- bpStrU12[grep("^[A|C|G|T]+$",bpStrU12)]
		accStrU12<- accStrU12[grep("^[A|C|G|T]+$",accStrU12)]

		donStrU2<- donStrU2[grep("^[A|C|G|T]+$",donStrU2)]
		accStrU2<- accStrU2[grep("^[A|C|G|T]+$",accStrU2)]
		pwmDonU12=Biostrings::PWM(donStrU12)
		pwmBpU12=Biostrings::PWM(bpStrU12)
		pwmAccU12=Biostrings::PWM(accStrU12)

		pwmDonU2=Biostrings::PWM(donStrU2)
		pwmAccU2=Biostrings::PWM(accStrU2)

		if(pdfFileSeqLogos!=""){
			grDevices::pdf(pdfFileSeqLogos, width=pdfWidth, height=pdfHeight)
# Codes borrowed from http://grokbase.com/t/r/bioconductor/1096ds3bnk/bioc-dra
#wing-multiple-sequence-logos-with-seqlogo#2010090767dkq0ngtwg3fbx8cjx5mxwsn0
			mySeqLogo = seqLogo::seqLogo

			bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | 
				sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
			body(mySeqLogo)[bad] = NULL


#Defining dimensions
			ncol <- 3
			nrow <- 2

#Creating plot

			grid::grid.newpage()

			listStr= list(donStrU12, bpStrU12, accStrU12, donStrU2, c(), 
				accStrU2)

			for(i in 1:nrow){
				for(j in 1:ncol){
					vp <- grid::viewport(x = (j-1)/ncol, y = 1-(i-1)/nrow, 
						w = 1/ncol, h =	1/nrow, just = c("left", "top"), 
						gp=grid::gpar(cex=cexSeqLogo))
					grid::pushViewport(vp)
			
					if(i!=2|j!=2){
						ind=(i-1)*3+j
						strVecTmp=listStr[[ind]]
						prTmp=lapply(strVecTmp, 
							function(tmp)unlist(strsplit(tmp, split="")))
						prTmpStrMat=matrix(unlist(prTmp),byrow=TRUE, 
							ncol=length(prTmp[[1]]))
		
						prTmpMat= matrix(unlist(lapply(1:ncol(prTmpStrMat), 
							function(tmp) 
								c(length(which(prTmpStrMat[,tmp]=="A")),
									length(which(prTmpStrMat[,tmp]=="C")),
									length(which(prTmpStrMat[,tmp]=="G")),
									length(which(prTmpStrMat[,tmp]=="T")))/
										nrow(prTmpStrMat))), 
							nrow=4, byrow=FALSE)
						rownames(prTmpMat)=c("A","C","G","T")
						mySeqLogo(prTmpMat)
						grid::grid.text(sprintf(c("U12 donor", 
							"U12 branch point", "U12 acceptor", "U2 donor", "",
							"U2 acceptor" )[ind], row, col), x=0.5, y=0.95, 
							just="top")
					}
					grid::upViewport()
				}
			}
			dev.off()
		}




		colnames(pwmDonU12)=as.character(1:ncol(pwmDonU12))
		colnames(pwmBpU12)=as.character(1:ncol(pwmBpU12))
		colnames(pwmAccU12)=as.character(1:ncol(pwmAccU12))
		colnames(pwmDonU2)=as.character(1:ncol(pwmDonU2))
		colnames(pwmAccU2)=as.character(1:ncol(pwmAccU2))

		res= list(pwmDonorU12=pwmDonU12, pwmBpU12=pwmBpU12, pwmAccU12=pwmAccU12,
			pwmDonU2=pwmDonU2, pwmAccU2=pwmAccU2)

		if(u12dbDropDb){
			dbc=dbConnect(MySQL())
			dbSendQuery(dbc, paste("DROP DATABASE ", u12dbDbName, ";", sep=""))
			dbDisconnect(dbc)
		}


		if(removeTempFiles){
			file.remove(paste(tmpDir,"u12db.sql.gz",sep="/"))
		}

	}

return (res)

}
