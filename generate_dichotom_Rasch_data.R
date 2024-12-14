#######################################
### Generate Dichotom Rasch Data    ###
#######################################
#theta~Norm(0,1) & beta~seq(-3.5,3.5,0.75), 10 item, 100 person
genresp<-function(b,theta){
	n=length(theta)   
	item=length(b)
	theta=matrix(nrow=n,ncol=item)
	for (i in 1:n){
		theta[i,]=rep(rnorm(1,mean=0,sd=1),item)
	}

	theta1=as.data.frame(theta[,1])
	
	beta=matrix(nrow=n,ncol=item)
	for (i in 1:item){
		beta[,i]=rep(b[i],n)
	}

	logit=theta-beta
	prob=1/(1+(exp(logit)^-1))

	dat_=matrix(prob,nrow=n,ncol=item)
	u=matrix(runif(n*item),nrow=n,ncol=item)
	dat=u<dat_
	dat=dat*1
	return(dat)

}
