#set seed to get the same results
#a different algorithm on claster initialization
set.seed(123)
#call the data  
#setwd("C:/Users/Koft/Desktop/Δ4-Θεωρία Πληροφορίας")
wine.data2<- read.csv("wine2.txt", header = T)
wine.data<-wine.data2[,c(2:14)]

#wine.data<-wine.data[sample(1:dim(wine.data)[1],20),] 
library(matlib)

#step1

#initialy choose k random elements
X<-as.matrix(wine.data)
for (i in 1:dim(X)[1]){
  for (j in 1:dim(X)[2]){
    X[i,j] <- X[i,j]/max(X[,j])
  }
}
X<-cbind(X, wine.data2$Wine)
X<-as.matrix(X)
#rownames(iris)<-paste("L",1:dim(x)[1])
#x<-as.matrix(iris)
#x2<-sapply(x[,1:4], as.numeric)
#x3<-matrix(x2, nrow=dim(x)[1], byrow=FALSE)

#x<-cbind(scale(x3), x[5])



#set number of initial clusters 
#k<- round(dim(x)[1]/4,0)
k<-10

#randomly choose k "centers"
kinitial <- sample(1:dim(X)[1], k)

#extract them from the dataset
kcentres <- X[kinitial,]

#vector of kcentres labels
clustvector <- vector(mode="numeric", length = dim(kcentres)[1])
for ( i in 1:dim(kcentres)[1]){
  clustvector[i] <- i
}
kcentres <- cbind( kcentres, clustvector)

#set matrix of the remains
remains <- X[-kinitial,]
iterations<-dim(remains)[1]


for (m in 1:(iterations-1)){   #repeat until only one element in remains
  #distance matrix
  #sample the matrix every time
  distmat<-matrix(, nrow=dim(kcentres)[1], ncol=dim(remains)[1])
  for (i in 1:dim(kcentres)[1] ){
    for ( j in 1:dim(remains)[1]){
      distmat[i,j] <- dist(rbind(kcentres[i, -c((dim(kcentres)[2]-1),dim(kcentres)[2])], 
                                 remains[j, -dim(remains)[2]]))
    }
  }

  #distmat1<- distmat[,-as.vector(which (distmat1 == min(distmat1), arr.ind = TRUE))[2]]
  #remove <- which (distmat1 == min(distmat1), arr.ind = TRUE)
  mindistvalue <- as.vector( which (distmat == min(distmat), arr.ind = TRUE))
  
  if (length(mindistvalue)==2){
    remove <-as.vector(remains[mindistvalue[2],]) 
  } else {
    l <-length(mindistvalue)/2
    pairsmat<- matrix(, nrow=l, ncol=2)
    pairsmat[1,1] <- mindistvalue[1]
    pairsmat[1,2] <-mindistvalue[l+1]
    for ( rows in 2:l ){
        pairsmat[rows, 1] <- mindistvalue[rows]
        pairsmat[rows, 2] <- mindistvalue[l+rows]
    }
    
    pair <-pairsmat[sample(nrow(pairsmat),1),]
    remove <-as.vector(remains[pair[2],]) 
  }
  #add to kcenters
  remove[dim(kcentres)[2]] <- kcentres[mindistvalue[1],'clustvector']
  kcentres <- rbind(kcentres , remove)

  #remove from yet unclustered elements
  if (length(mindistvalue)==2){
    remains <-remains[-mindistvalue[2],] 
  } else {
    remains <-remains[-pair[2],]
  }
  
  # print(kcentres[dim(kcentres)[1],])
}


#for the last element in remains
distmat<-matrix(, nrow=dim(kcentres)[1], 1)

for (i in 1:dim(kcentres)[1] ){
  distmat[i,1] <- dist(rbind(kcentres[i, -dim(kcentres)[2]], remains))
}

#distmat1<- distmat[,-as.vector(which (distmat1 == min(distmat1), arr.ind = TRUE))[2]]
#remove <- which (distmat1 == min(distmat1), arr.ind = TRUE)
mindistvalue <- as.vector( which (distmat == min(distmat), arr.ind = TRUE))

#add to kcenters
remove <-remains 
remove[dim(kcentres)[2]] <- kcentres[mindistvalue[1],'clustvector']
kcentres <- rbind(kcentres , remove)

#print to see how many in each class
hist(as.numeric(kcentres[,'clustvector']),breaks = length(cluster))

#########################################################################
#Clustering algorithm

kcen<-as.data.frame(kcentres)
#first itteration
cluster <- list()
#assign every element to the data frame in the list
#containing each cluster label
for ( i in 1:k ){
  cluster[[i]] <- as.matrix(data.frame(kcen[kcen$clustvector==i,]))
}
#cluster[[i]] is the matrix containing every element belonging to cluster i


######################################################
######################################################
#irritative algorithm

#calculation of the following matrix is the most time consuming
#part of the algorithm
#it becomes more complicated, as we add features
W <- matrix(, ncol=dim(kcentres)[1], nrow=dim(kcentres)[1])
for (j in 1:dim(kcentres)[1]){
  for (i in 1:dim(kcentres)[1]){
    x1<-kcentres[i,]
    x2<-kcentres[j,]
    
    x1 <- x1[-c(length(x1)-1, length(x1))]
    x2 <- x2[-c(length(x2)-1, length(x2))]
    
    x1<-as.numeric(as.character.numeric_version(x1))
    x2<-as.numeric(as.character.numeric_version(x2))
    
    W[i,j] <- gausiankernel(x1,x2)
    diag(W)=0	#to get read of the e^0=1 that are huge numbers compared to the rest of the matrix values
  }
}

which(W==0, arr.ind=TRUE) #ok its none!!


for (step in 1:(k-2)){
  
  worstCluster <- c()
  
  
  for ( i in 1:(length(cluster))){
    worstCluster[i] <- worstClusterEntropy(W,i)
    #print(worstCluster[i])
  }
  worst <- as.vector( which (worstCluster == max(worstCluster), arr.ind = TRUE))
  print(worst)  #it is the cluster that we will eliminate
  
  
  #add one element at each cluster one by one 
  #and then compute again the new withinClustEntropyplusEleme(t)-...
  #d2<-list()
  #length(d2)<-dim(cluster[[worst]])[1]
  #a<-c()
  #for (i in 1:dim(cluster[[worst]])[1]){
  #i=1
  #d2<-c()
  #  for (m in 1:length(cluster)){
  #    if ( m !=worst){
  #      t<-rbind(cluster[[m]], cluster[[worst]][i,])
  #      d2[[i]][m]<-withinClustEntropyplusElem(t)-withinClustEntropy(m)
  #    }
  #  }
  #print(d2)
  #  a[i]<-as.vector( which (d2[[i]] == min(d2[[i]], na.rm=T), arr.ind = TRUE))
  #  cluster[[worst]][i,'clustvector'] <- a[i]
  #  cluster[[a[i]]]<-rbind(cluster[[a[i]]], cluster[[worst]][i,])
  #  #}
  #  #print(a)
  #}
  #a
  
  d1<-matrix(, nrow=dim(cluster[[worst]])[1], ncol=length(cluster))
  for (m in 1:length(cluster)){
    if (m!=worst){
      for (elemworst in 1:dim(cluster[[worst]])[1]){
        t<-rbind(cluster[[m]], cluster[[worst]][elemworst,])
        d1[elemworst,m]<-withinClustEntropyplusElem(t)-withinClustEntropy(m)
      }
    }
  }
  
  a<-c()
  for (i in 1:dim(d1)[1]){
    a[i]<-as.vector( which (d1[i,] == min(d1[i,], na.rm=T), arr.ind = TRUE))
    cluster[[worst]][i,'clustvector'] <- a[i]
    cluster[[a[i]]]<-rbind(cluster[[a[i]]], cluster[[worst]][i,])
  }
  print(a)	#a shows the cluster that each element of the Worst cluster assigned to
  
  
  cluster[[worst]]<-NULL	#destroy the worst cluster
  if (worst<=length(cluster)){
    for ( i in (worst):length(cluster)){
      for (j in 1:dim(cluster[[i]])[1]){
        cluster[[i]][j,'clustvector'] <- as.numeric(cluster[[i]][j,'clustvector'])-1
      }
    }
  }
  
  print(betweenClusterEntropy(W))
  
  kcentres <- do.call("rbind", lapply(cluster, as.data.frame))
  #kcentres
  kcentres$clustvector<-as.factor(kcentres$clustvector)
  kcentres$V5<-as.factor(kcentres$V5)
  #summary(kcentres)
  #summary(X)
  plot(kcentres$clustvector, col=rainbow(k))
  
  filename<-sprintf("iterIris%sCluster.csv", step)
  
  #cluster.df <- do.call("rbind", lapply(cluster, as.data.frame))
  write.csv(kcentres, file=filename)
}

#call the k-3 csv to create the comparisson table
cl<-read.csv("iterIris12Cluster.csv")
t<-table(cl$clustvector, iris$Species)
t

library(NMF)
purity(t)
entropy(t)
t
###  Calculate accuracy
lamda<-0
for (g in 1:3){
  ind<-which(t==max(t), arr.ind = T)
  lamda<-lamda+t[ind]
  t[ind[1],]<-0
  t[,ind[2]]<-0
}
accur<-lamda/150
accur