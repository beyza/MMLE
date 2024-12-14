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
