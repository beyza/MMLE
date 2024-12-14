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


################
################
###### EAP #####
################
################

#beta: numeric(betapar)
#dat:numeric(response data)

eap<-function(beta,dat){
prior.mean=0
prior.sd=1
min.theta=-10
max.theta=10
inc=1
theta=seq(min.theta,max.theta,inc)
nq=length(theta)

# prior

prior=dnorm((theta-prior.mean)/prior.sd)

#likelihood
like=numeric(nq)
for(i in 1:nq){
 logit=theta[i]-beta
 prob=1/(1+(exp(logit)^-1))
         for (k in 1:length(dat)){
         if(dat[k]==0) prob[k]=1-prob[k]
        }
  like[i]=exp(sum(log(prob)))
}


# posterior

posterior=numeric(nq)
posterior=prior*like
posterior=posterior/(sum(posterior))


# EAP

eap=numeric(nq)
eap=theta*posterior
eap=sum(eap)
eap=round(eap,3)
     

# SE 

se=numeric(nq)
se=posterior*(theta-eap)^2
se=sqrt(sum(se))
se=round(se,3)       

    
return(cbind(eap,se))

}

########################
########################
######## MLE ###########
########################

#http://www.rasch.org/rmt/rmt102t.htm
#beta: numeric(betapar)
#dat:numeric(response data)

mle<-function(beta,resp,eps=0.01){
	L=length(beta)
	ifelse(sum(resp==1)==L,R<-0.5, R<-sum(resp==1))
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

########################
########################
##### MODIFIED MLE #####
########################
########################

mmle=function(beta,dat,eps=0.01){
n=length(beta) #kisiye sorulan toplam soru sayisi
d=seq(1:n) #d=1,...,n
data=cbind(dat,beta)

# 1) compute t_i
t=qlogis(d/(n+1))  # expected value of z_i

#2) compute alpha1_i, beta1_i, alpha2_i, beta2_i 
f=expression(exp(-u)/(1+exp(-u))^2) 
df=D(f,"u")  #derivative of standard logistic distribution function 
u=t
F=plogis(t)

g1=dlogis(t)/plogis(t) 
g2=dlogis(t)/(1-plogis(t))  

b1=(dlogis(t)^2-eval(df)*F)/(F^2)
a1= g1+b1*t

b2=(dlogis(t)^2+(1-F)*eval(df))/(1-F)^2
a2=g2-b2*t

#3&4)order statistics: beta(i) and u[i] & compute m_i and delta_i 

if(length(d)==1) {
m=b1*data[1]+b2*(1-data[1])
delta=a1*data[1]+a2*(1-data[1])
} else{
data=data[order(beta),] 
m=b1*data[,1]+b2*(1-data[,1])
delta=a1*data[,1]-a2*(1-data[,1])
}

#5) estimate theta using weighted average of m, delta and beta;

if(length(d)==1){
theta=(sum(delta)/sum(m) + sum(m*data[2])/sum(m)) 
} else {
theta=(sum(delta)/sum(m) + sum(m*data[,2])/sum(m))
}

#6)SE of theta
data=cbind(dat,beta)
z=theta+data[,2]
Q=dlogis(z)^2/(plogis(z)*(1-plogis(z)))
sd_theta=sqrt(1/sum(Q))

if(length(dat)>1&&sum(dat==1)!=length(dat)&&sum(dat==0)!=length(dat)){
#######################
## revised estimates ##
#######################

# 1) compute t_i

t=theta - sort(beta)  # expected value of z_i

#2)compute alpha1_i, beta1_i, alpha2_i, beta2_i
f=expression(exp(-u)/(1+exp(-u))^2) 
df=D(f,"u")  #standart lojistik dagilim fonksiyonu t revi
u=t
F=plogis(t)

g1=dlogis(t)/plogis(t) 
g2=dlogis(t)/(1-plogis(t))  

b1=(dlogis(t)^2-eval(df)*F)/(F^2)
a1= g1+b1*t

b2=(dlogis(t)^2+(1-F)*eval(df))/(1-F)^2
a2=g2-b2*t

#3&4)order statistics beta(i) and u[i] & compute m_i and delta_i 
if(length(d)==1) {
m=b1*data[1]+b2*(1-data[1])
delta=a1*data[1]+a2*(1-data[1])
} else{
data=data[order(beta),] 
m=b1*data[,1]+b2*(1-data[,1])
delta=a1*data[,1]-a2*(1-data[,1])
}

#5) revise theta using weighted average of m, delta and beta;

if(length(d)==1){
theta_new=(sum(delta)/sum(m) + sum(m*data[2])/sum(m)) 
} else {
theta_new=(sum(delta)/sum(m) + sum(m*data[,2])/sum(m))
}

if(abs(theta-theta_new)<eps){
	theta=theta_new
}

data=data=cbind(dat,beta)

while(abs(theta-theta_new)>eps){

	theta=theta_new

	# 1) compute t_i

	t=theta - sort(beta)  # expected value of z_i

	#2)compute alpha1_i, beta1_i, alpha2_i, beta2_i
	f=expression(exp(-u)/(1+exp(-u))^2) 
	df=D(f,"u")  #standart lojistik dagilim fonksiyonu t revi
	u=t
	F=plogis(t)

	g1=dlogis(t)/plogis(t) 
	g2=dlogis(t)/(1-plogis(t))  

	b1=(dlogis(t)^2-eval(df)*F)/(F^2)
	a1= g1+b1*t

	b2=(dlogis(t)^2+(1-F)*eval(df))/(1-F)^2
	a2=g2-b2*t

	#3&4)order statistics beta(i) and u[i] & compute m_i and delta_i 
	if(length(d)==1) {
	m=b1*data[1]+b2*(1-data[1])
	delta=a1*data[1]+a2*(1-data[1])
	} else{
	data=data[order(beta),] 
	m=b1*data[,1]+b2*(1-data[,1])
	delta=a1*data[,1]-a2*(1-data[,1])
	}

	#5) revise theta using weighted average of m, delta and beta;

	if(length(d)==1){
	theta_new=(sum(delta)/sum(m) + sum(m*data[2])/sum(m)) 
	} else {
	theta_new=(sum(delta)/sum(m) + sum(m*data[,2])/sum(m))
	}
	data=data=cbind(dat,beta)
}


#6)SE of theta
z=theta_new+data[,2]
Q=dlogis(z)^2/(plogis(z)*(1-plogis(z)))
sd_theta=sqrt(1/sum(Q))
}else{
theta_new=theta
}

#7) MLE i in hesaplnana sd

sd_theta_mle=sqrt(1/(sum(exp(theta-data[,2])/(1+exp(theta-data[,2]))^2)))

return(cbind(theta_new,sd_theta,sd_theta_mle))

}



######################
######################
###### FISHER INFO ###
######################
######################

#calculates information of a single item

info<-function(beta,theta){
logit=theta-beta
pr=1/(1+(exp(logit)^-1))
inf=pr*(1-pr)
return(inf)
}

