
#set seed to get the same results
#a different algorithm on claster initialization
set.seed(123)



#library(matlib)
#call the data            
#wine.data2<- read.csv("wine2.txt", header = T)
#wine.data<-wine.data2[,c(2:14)]

#xwine<-wine.data[,c(2,4,8,10)]

#wine.data<-wine.data[sample(1:dim(wine.data)[1],20),] 

######################################
##########normalize the data ###########OPTIONAL
X<-as.matrix(iris[,1:4])
for (i in 1:dim(X)[1]){
  for (j in 1:dim(X)[2]){
    X[i,j] <- X[i,j]/max(X[,j])
  }
}

X<-cbind(X, iris$Species)
X<-as.matrix(X)

#library(HSAR)
#X<-as.data.frame(X)
#with(X[1:13], pairs(X[1:13], col=c(1:8)[X$`wine.data2$Wine`]))
#X$V14<-as.factor(X$V14)
#summary(X)
#X<-as.matrix(X)


#####################################################
######set number of initial clusters 
#k<- round(dim(X)[1]/5,0)
k<-15

#we want the clusters to have equal number of elements

#so we set a maximum value of elements in each cluster as follows
Nmax<-round((dim(X)[1]-k)/k, digits=0)

#each cluster will have Nmax+1 elements, because cluster in the list starts 
#with already one element


clusters<-list()
length(clusters)<-k

#iterations<-dim(remains)[1]
#first distance matrix

distmat <- matrix(, nrow=dim(X)[1], ncol=dim(X)[1])
distmat <- dist(X[,-dim(X)[2]], upper=TRUE)
distmat <- as.matrix(distmat)
diag(distmat) <-max(distmat)+1

kinitial <- sample(1:dim(X)[1], 1)

clustlabels<-c()

#filling clusters proceess
for (j in 1:(length(clusters)-1)){
  #initialize the cluster
  clusters[[j]]<-X[kinitial,]
  
  #clusters[[j]]
  
  values<-c()
  values[[1]]<-kinitial
  values
  
  columnk<-NULL
  
  for ( i in 1:Nmax){
    
    columnk<-distmat[-values,values]
    #columnk
    mindistvalue <- as.vector( which (columnk == min(columnk), arr.ind = TRUE))
    mindistvalue
    #find the label that is on the mindistvalue[[1]]-th position
    removelabel<-as.numeric(rownames(as.data.frame(columnk))[mindistvalue[[1]]])
    #removelabel
    #bind the values together
    values<-cbind(values, removelabel)
    #values
    remove<-X[removelabel,]
    clusters[[j]]<-rbind(clusters[[j]], remove)
  }
  columnk<-distmat[-values,values]
  
  clustlabels<-c(clustlabels, as.numeric(colnames(columnk)))
  #print(clustlabels)
  #print(values)
  
  #clusters[j]
  #set the values we the elements we already used from the distmat as max(distmat)+1
  #so they wont be selevted again
  distmat[, as.numeric(colnames(columnk))] <- max(distmat)
  distmat[as.numeric(colnames(columnk)), ] <- max(distmat)
  kinitial <- sample(as.numeric(rownames(distmat[-clustlabels,])),1)
  
  #print(max(distmat))
  #sum(distmat==max(distmat))
}

#And for the last cluster, add all elements left
labels<-c()
m<-0
for (i in 1:dim(distmat)[1]){
  if (sum(distmat[,i])!=dim(distmat)[1]*max(distmat)){
    m<-m+1
    labels<-cbind(labels,i)
  }
}

labels
length(labels)
clusters[[k]]<-X[labels,]

for (i in 1:length(clusters)){
  clusters[[i]]<-cbind(clusters[[i]], 'clustvector'=i)
}

kcentres <- do.call("rbind", lapply(clusters, as.data.frame))
#kcentres
kcentres$clustvector<-as.factor(kcentres$clustvector)
kcentres$V5<-as.factor(kcentres$V5)
summary(kcentres)
#summary(X)
plot(kcentres$clustvector, col=rainbow(k))

sum(distmat[,1]==max(distmat))

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