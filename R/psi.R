psi<- function(x, y, intCol, exCol, pseudoCnt=0){
	if(missing(y))
		y<-x

	xCnt<-
		counts(x)[,
			intCol]
	yCnt<-
		counts(y)[,
			exCol]

	psiRes<- (xCnt)/(xCnt+(yCnt*2)+pseudoCnt)
	return (psiRes)

}
