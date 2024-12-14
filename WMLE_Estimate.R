##############################
##############################
######### Weighted MLE #######
##############################
##############################


# p = logit[exp((theta - b))]
p.br<-
function(b, theta){

  return(  1 / ( 1 + exp(-( theta - b ) )))

} 


# q = 1 - p = 1/(1+exp(theta-beta))
q.br <-
function(b, theta){
    
  return( 1 / ( 1 + exp(  ( theta - b ) ) ) )
    
} 

# p' =  exp(theta-beta)/(1+exp(theta-beta))^2
pder1.br <-
function(b, theta){
  
  ex <- exp(theta - b)
    
  return( ex / ( 1 + ex )^2 )
    
} 

# p''= (exp(theta-b)-exp(theta-b)^2)/(1+exp(theta-b))^3
pder2.br <-
function(b, theta){
  
  ex <- exp( theta - b )
    
  return( ( ex - ex^2 ) / ( 1 + ex )^3  )
    
} 

pder3.br<-
function(b,theta){
    ex <- exp(theta - b)
    return( ex *(ex^2-4*ex+1)/ ( 1 + ex )^4 )
}

#I= p'^2 / (p * q)
Ii <- 
function(b, theta){
  
  # Calculating the probability of response:
  p <- p.br(b, theta)
  q <- 1 - p
  
  # Calculating the first and second derivatives:
  pder1  <- pder1.br(b, theta)
  pder2 <- pder2.br(b, theta)
  I <- pder1^2 / (p * q)
  Ider1 <- pder1*(2*p*q*pder2-pder1^2*(q-p))/(p^2*q^2)
  return(list(I=I,Ider1=Ider1))
}


#Ji=(p'*p'')/(p*q)
Ji <- 
function(b, theta){
  
  # Calculating the probability of response:
  p <- p.br(b, theta)
  q <- 1 - p
  
  # Calculating the first and second derivatives:
  pder1  <- pder1.br(b, theta)
  pder2 <- pder2.br(b, theta)
  pder3 <- pder3.br(b,theta)
  J <- (pder1 * pder2)  / (p * q) 
  Jder1 <- (p*q*(pder2^2+pder1*pder3)-pder1^2*pder2*(q-p))/(p^2*q^2)
  return(list(J=J,Jder1=Jder1))
}
  
WLE_delta <-
  function(b,theta,u){
  # Calculating the probability of response:
  p <- p.br(b, theta)
  q <- 1 - p
  
  # Calculating the first and second derivatives:
  pder1  <- pder1.br(b, theta)
  pder2 <- pder2.br(b, theta)
  info <- sum(Ii(b,theta)$I)
  info.der <- sum(Ii(b,theta)$Ider1)
  J <- sum(Ji(b,theta)$J)
  J.der <- sum(Ji(b,theta)$Jder1)
  lder <- sum( (u - p) * pder1 / (p * q) )
  H <- ((info*J.der-info.der*J)/(2*info^2))+info
  delta <- ( lder + J / (2 * info)) / (H) 
  return(list(delta=delta,info=info,H=H))
}

WLE_est <-
  function(b,u,theta0=0,exac=0.001,max_iter=33){
  count <- 0
  ergv <- rep(NA,max_iter)
  repeat
	{
	del <- WLE_delta(b,theta=theta0,u)$delta
  if(is.na(del))
    {
  	theta0  <- NA
    SE_wle     <- NA
	SE_mle     <- NA
		break
    }
	if(abs(del)>2){del <- del/abs(del) * 2}
  theta.next <- theta0 + del
	theta0 <- theta.next
  count <- count + 1
  ergv[count] <- theta0
	if(abs(del) <= exac | count >= max_iter)
    {
    SE_wle <- 1/sqrt(WLE_delta(b,theta=theta0,u)$H)
	SE_mle <- 1/sqrt(WLE_delta(b,theta=theta0,u)$info)
    break
    }
	
	}

#return(list("resp"=u,"theta_wle"=theta0,"iterations"=count,"estproc"=ergv,"SE_wle"=SE_wle,"SE_mle"=SE_mle))
return(cbind(theta0,SE_wle))
}




  
