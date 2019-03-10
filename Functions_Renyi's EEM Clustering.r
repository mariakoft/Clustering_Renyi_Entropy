########## Functions defined
#Kernel Density Estimator of pdf
gausiankernel <- function(x1,x2){
  x<-x1-x2
  x
  #s <- cov(x1, x2)   #if used sigma^2=s
  s<-0.36
  #smatrix <- (s^2)*diag(length(x))
  #H <- base::norm(smatrix, "F")
  #A <- as.numeric(t(x)%*%matlib::Ginv(smatrix)%*%x)
  value <- ((2*pi*s)^(-length(x)/2))*exp(-norm(as.matrix(x))/(2*s^2))
  value
  return(value)
}


#calculate within cluster entropy for m-th cluster
withinV <- function(m){
  
    kernelmat<-matrix(, nrow=dim(cluster[[m]])[1], ncol=dim(cluster[[m]])[1])
    for (i in 1:dim(cluster[[m]])[1]){
      for (j in 1:dim(cluster[[m]])[1]){
        x1<-cluster[[m]][i,]
        x2<-cluster[[m]][j,]
        
        x1<-as.numeric(x1[-c(length(x1)-1, length(x1))])
        x2<-as.numeric(x2[-c(length(x2)-1, length(x2))])
        
        
        #made changes to kernel function
        #keeps only the last element, it is wrong
        kernelmat[i,j]<-gausiankernel(x1,x2)
        
      }
    }
    return(sum(kernelmat)/(dim(cluster[[m]])[1])^2)
 
}


withinClustEntropy <-function(m){
  return(-log(withinV(m)))
}
##os edw kala mallon


#M and delta are functions used to identify for which elements we should calculate 
#the

#compute between cluster information gain 
#only if the elements for which i am computing the gausian kernel
#do not belong in the same function
M <- function(x){
  if (x==FALSE){
    return(1)
  } else {
    return(0)
  }
}



#the next function fills the elements of a matrix
#when i,j do not belong to the same cluster
#more elements belong in the same cluster-->the matrix will produce a bigger sum

#B is Ξ (Nk) with Nk the cardinality of each cluster
#and it is a very large number

#bigger sum(A)-->leads to a smaller ratio  **probably
#leading to a bigger betweenEntropy??

#between cluster entropy
betweenV <- function(){
  #it is computed on kcentres
  #on every itteration the labels on kcentres change according to the
  #number of clusters
  #and the partioning number of every cluster ??
  
  A <- matrix(, ncol=dim(kcentres)[1], nrow=dim(kcentres)[1])
  for (i in 1:dim(kcentres)[1]){
    for (j in 1:dim(kcentres)[1]){
      x1<-kcentres[i,]
      x2<-kcentres[j,]
      x <- x1['clustvector']==x2['clustvector']
      
      
      #x1 <- as.numeric(x1[-c(length(x1)-1, length(x1))])
      #x2 <- as.numeric(x2[-c(length(x2)-1, length(x2))])
      
      x1 <- x1[-c(length(x1)-1, length(x1))]
      x2 <- x2[-c(length(x2)-1, length(x2))]
      
      x1<-as.numeric(as.character.numeric_version(x1))
      x2<-as.numeric(as.character.numeric_version(x2))
      
      A[i,j] <- M(x)*gausiankernel(x1,x2)
    }
  }
  B <-1
  for (m in 1:length(cluster)){
    B<-B*dim(cluster[[m]])[1]    #B is a huge number
  }
  return(sum(A)/(2*B))    #therefore this sum is very close to 0
}

 

#works with matrix W
betweenVV <- function(W){
  
  Wmat<-W
  #make a NxN matrix, where N=the number of elements
  for (i in 1:dim(kcentres)[1]){
    for (j in 1:dim(kcentres)[1]){
      x1<-kcentres[i,]
      x2<-kcentres[j,]
      
      if (x1['clustvector']==x2['clustvector']){
        Wmat[i,j]<-0
      }
    }
  }
  
  B <-1
  for (m in 1:length(cluster)){
    B<-B*dim(cluster[[m]])[1]
  }
  return(sum(Wmat)/(2*B))
}


betweenClusterEntropy <- function(W){
  return(-log(betweenVV(W)))
} 



#worst cluster criterion
#if anyone of x1 or x2 belong to cluster m then 
#delta returns 0
delta <- function(x1,x2,x,m){
  if  (x1['clustvector']==m |x2["clustvector"]==m){
    return(0)
  }
  else {
    return(1)
  }
}



worstClusterV <- function(m){
  #make a NxN matrix, where N=the number of elements
  W <- matrix(, ncol=dim(kcentres)[1], nrow=dim(kcentres)[1])
  for (j in 1:dim(kcentres)[1]){
    for (i in 1:dim(kcentres)[1]){
      x1<-kcentres[i,]
      x2<-kcentres[j,]
      
      x <- x1['clustvector']==x2['clustvector']
      d<-delta(x1,x2,x,m)
      
      #x1 <- as.numeric(x1[-c(length(x1)-1, length(x1))])
      #x2 <- as.numeric(x2[-c(length(x2)-1, length(x2))])
      
      x1 <- x1[-c(length(x1)-1, length(x1))]
      x2 <- x2[-c(length(x2)-1, length(x2))]
      
      x1<-as.numeric(as.character.numeric_version(x1))
      x2<-as.numeric(as.character.numeric_version(x2))
      
      W[i,j] <- d*gausiankernel(x1,x2)
      #W[i,j] <- gausiankernel(x1,x2)
    
    }
  }
  B <-1
  for (m in 1:length(cluster)){
    B<-B*dim(cluster[[m]])[1]
  }
  return(sum(W)/(2*B))
  #return(W)
}


#worstClusterEntropy <- function(m){
#  return(-log(worstClusterV(m)))
#}

#slightly different and maybe faster compute worstClusterV
#keep W and change if based on m
#so, we must compute in every iteration
#before calling the function
#W must be computed inside the iterrative algorithm and not inside the function
#name it worstClusterVV

#compute W outside


worstClusterVV <- function(W, m){
  Wmat<-W
  #make a NxN matrix, where N=the number of elements
  for (j in 1:dim(kcentres)[1]){
      x1<-kcentres[j,]
      
      if (x1['clustvector']==m){
        Wmat[j,]<-0
        Wmat[,j]<-0
      }
  }
    
  B <-1
  for (m in 1:length(cluster)){
    B<-B*dim(cluster[[m]])[1]
  }
  return(sum(Wmat)/(2*B))
  #return(Wmat)
}


worstClusterEntropy <- function(W,m){
  return(-log(worstClusterVV(W,m)))
}



withinVplusElem <- function(t){
  kernel <- matrix(, ncol=dim(t)[1], nrow=dim(t)[1])
  
    for (i in 1:dim(t)[1]){
      for (j in 1:dim(t)[1]){
        x1<-t[i,]
        x2<-t[j,]
        
        #x1 <- as.numeric(x1[-c(length(x1)-1, length(x1))])
        #x2 <- as.numeric(x2[-c(length(x2)-1, length(x2))])
        
        x1 <- x1[-c(length(x1)-1, length(x1))]
        x2 <- x2[-c(length(x2)-1, length(x2))]
        
        x1<-as.numeric(as.character.numeric_version(x1))
        x2<-as.numeric(as.character.numeric_version(x2))
        
        kernel[i,j] <- gausiankernel(x1,x2)
      }
    }
    return(sum(kernel)/(dim(t)[1])^2)
  
}

withinClustEntropyplusElem <-function(t){
  return(-log(withinVplusElem(t)))
}

