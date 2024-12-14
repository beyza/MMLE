#################################
#################################
######## Anchored MLE ###########
#################################

#http://www.rasch.org/rmt/rmt102t.htm
#beta: numeric(betapar)
#dat:numeric(response data)

mle<-function(beta,resp,eps=0.01){
	L=length(beta)
	ifelse(sum(resp==1)==0,R<-0.5, R<-sum(resp==1))
	W=L-R
	if(sum(resp==0)==0) {
		W<-0.5; R<-L-W
	}

	M=mean(beta)+log(R/W)

	score=sum(exp(M-beta)/(1+exp(M-beta)))
	variance=sum(exp(M-beta)/(1+exp(M-beta))^2)
	M_new=M+((R-score)/variance)

	while(abs(M-M_new)>eps){
		M=max(c(min(c(M+1,M_new)),M-1))
		score=sum(exp(M-beta)/(1+exp(M-beta)))
		variance=sum(exp(M-beta)/(1+exp(M-beta))^2)
		M_new=M+((R-score)/variance)
	}
	M=max(c(min(c(M+1,M_new)),M-1))
	se_M=sqrt(1/variance)
return(c(M,se_M))

}
