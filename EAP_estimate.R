
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

