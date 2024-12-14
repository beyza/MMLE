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


