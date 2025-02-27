install.packages("epimdr2")
install.packages("statnet")
install.packages("vioplot")
install.packages("GGally")
library(statnet)
library(GGally)

#Calculation scalefreeNetwork
set.seed(17)
barabasiAlbertModel=function(Nmax,NodeDegree) {          #Nmax: Total of Nodes| NodeDegree: Node Degree of new nodes| Method to create a skalefree Network
  adjMatrix=matrix(0,ncol=Nmax, nrow=Nmax)             #adjMatrix is 0 if there isn't a connection between 2 Nodes and 1 if there is a connection.
  adjMatrix[1,2]=adjMatrix[2,1]=1
  for(i in 3:Nmax){
    probs=apply(adjMatrix,1,sum)                      #calculate the probability for a link
    link=unique(sample(c(1:Nmax)[-i],size=min(c(NodeDegree,i-1)),prob=probs[-i]))       #Connact new Nodes with K links and Probability probab
    adjMatrix[i,link]=adjMatrix[link,i]=1              #Connect the new Node with K Nodes 
  }
  return(adjMatrix)
}
scalefreeNetwork=barabasiAlbertModel(1000,5)
testsum=summary(scalefreeNetwork)
#net=network(scalfreeNetwork, directed=FALSE)
#plot(net,displaylabels=TRUE, vertex.cex=2)             Theoretisch geeignet um Netz darzustellen aber nicht sinnvoll, da für zu viel Knoten nicht verständlich
#ggnet2(net, color="blue", edge.color="gray", size=3, label=TRUE)
degree_dist <- table(rowSums(scalefreeNetwork))  # Knotengrade berechnen
plot(as.numeric(names(degree_dist)), as.numeric(degree_dist),log="xy",xlab="Degree",ylab="Frequency")
plot(as.numeric(names(degree_dist)), as.numeric(degree_dist),xlab="Degree",ylab="Frequency")

#Funktion for the spread on the Network
i=0
Diseafunktion=function(N,S,I,R,I1,scalefreeNetwork,lambda,gamma){
  I[I1,1]=1                                       #I1 infectios
  S[I1,1]=0                                       #I1 no longer susceptible
  Srel=c(0)
  Irel=c(0)
  Rrel=c(0)
  Itotal=c(0)
  Itotal[1]=1/N
  Srel[1]=sum(S[,1])/N
  Irel[1]=sum(I[,1])/N
  Rrel[1]=sum(R[,1])/N
  t=1
  time=c(1)
  while(sum(I[,t-1], na.rm=TRUE) > 0|t==1){
    t=t+1
    y=scalefreeNetwork%*%I[,t-1]
    p=1-(1-lambda)^y                               #Total probability for infection of one node(Thesis page 5)
    newI = rbinom(N, size = S[,t-1], prob = p)# new infectious people
    Itotal[t-1]=sum(newI)/N
    newR = rbinom(N, size = I[,t-1], prob = gamma)  # new recovered perople
    nextS=S[,t-1]-newI                            #explaint in set theory:remove all new infected individuals from the Susceptibles
    nextI=I[,t-1]+newI-newR                       #explaint in set theory:add all new infected and remove all new removed/recovered from theset of infectious
    nextR=R[,t-1]+newR                            #explaint in set theory:add all new recovered to the recovered
    Srel[t-1]=sum(nextS)/N                        #percentage of S
    Irel[t-1]=sum(nextI)/N                        #              I
    Rrel[t-1]=sum(nextR)/N                        #              R
    I=cbind(I,nextI)                              #save the new S,I,R in the  matrix
    S=cbind(S,nextS)
    R=cbind(R,nextR)
    resultabs<<-list(S=S,I=I,R=R)
    resultrel=matrix(0,nrow=t,ncol=4)
    resultrel<<-cbind(Srel,Irel,Rrel,Itotal)
  }
}
#Simulating Disea on the Network
disea=function(scalefreeNetwork,lambda,gamma){
  N=dim(scalefreeNetwork)[1]
  S=matrix(rep(1, N),nrow=N,ncol=1)               #Susceptibles
  I=matrix(rep(0, N),nrow=N,ncol=1)               #Infektios
  R=matrix(rep(0, N),nrow=N,ncol=1)               #Recovered
  I1=sample(1:N, size=1)                          #Choose one random person which getting infected
  Diseafunktion(N,S,I,R,I1,scalefreeNetwork,lambda,gamma)
}

#Duration until the Disea dies out of the Network
duration=numeric(100)                                   
for (i in 1:100){                                   #Calculating the average duration for Dieseas with at least 7 Periods
  x=0
  while(x<=7){                                      #Diseas with a period of 2 aren't that relevant because there wasn't really a Disea
    disea(scalefreeNetwork,0.15,0.5)
    x=length(resultrel[,1])
  }
  duration[i]=x
}
durationavg=sum(duration)/100
round(durationavg)
#Average Values in each Period                      #Calculating the average relativ S,I,R for Diseas with average Duration
Savg=matrix(0,nrow=round(durationavg),ncol=100)
Iavg=matrix(0,nrow=round(durationavg),ncol=100)
Ravg=matrix(0,nrow=round(durationavg),ncol=100)
Itotalavg=matrix(0,nrow=round(durationavg),ncol=100)
for (i in 1:100){
  x=0
  while(x!=round(durationavg)){
    disea(scalefreeNetwork,0.15,0.5)
    x=length(resultrel[,1])
  }
  Savg[,i]=resultrel[,1]
  Iavg[,i]=resultrel[,2]
  Ravg[,i]=resultrel[,3]
  Itotalavg[,i]=resultrel[,4]
}
avg=matrix(0,nrow = round(durationavg),ncol=4)
colnames(avg) <- c("S", "I", "R","4") 
for (i in 1:round(durationavg)){
  avg[i,1]=sum(Savg[i,])/100
  avg[i,2]=sum(Iavg[i,])/100
  avg[i,3]=sum(Ravg[i,])/100
  avg[i,4]=sum(Itotalavg[i,])/100
}
Ikum=numeric(round(durationavg))
for(i in 1:round(durationavg)){
  Ikum[i] = sum(avg[1:i, 4])
}
avg[,1:3]
plot(x=(1:round(durationavg)),y=avg[,1],ylab="Fraction",xlab="Time",type="l",ylim=c(0,1))     #Plotting S,I,R
lines(x=(1:round(durationavg)),y=avg[,2],col="red")
lines(x=(1:round(durationavg)),y=Ikum,col="purple")
lines(x=(1:round(durationavg)),y=avg[,3],col="green")

#Random Vaccianation
diseawithrandvacc=function(scalefreeNetwork,lambda,gamma,percVacc){
  N=dim(scalefreeNetwork)[1]
  V=matrix(rep(0,N),nrow=N,ncol=1)
  S=matrix(rep(1, N),nrow=N,ncol=1)               #Susceptibles
  I=matrix(rep(0, N),nrow=N,ncol=1)               #Infektios
  R=matrix(rep(0, N),nrow=N,ncol=1)               #Recovered
  Ivacc=sample(1:N,size=round(percVacc*N),replace=FALSE)     #Create random nodes which getting vaccinated
  S[Ivacc,1]=0                                    #Take those from the group of susceptibles
  V[Ivacc,1]=1                                    #put them in the group of Vaccinated
  Nonvacc=which(V[,1]==0)
  I1=sample(Nonvacc, size=1)                      #Choose one random person of the susceptible group which getting infected
  Diseafunktion(N,S,I,R,I1,scalefreeNetwork,lambda,gamma)
}
#Duration until the Disea dies out of the Network
duration=numeric(100)
for (i in 1:100){                                   #Calculating the average duration for Dieseas with at least 7 Periods
  x=0
  while(x<=7){                                      #Diseas with a period of 2 arent that relevant because there wasn't really a Disea
    diseawithrandvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  duration[i]=x
}
durationavg=sum(duration)/100
round(durationavg)
#Average Values in each Period                      #Calculating the average relativ S,I,R for Diseas with average Duration
Savg=matrix(0,nrow=round(durationavg),ncol=100)
Iavg=matrix(0,nrow=round(durationavg),ncol=100)
Ravg=matrix(0,nrow=round(durationavg),ncol=100)
Itotalavg=matrix(0,nrow=round(durationavg),ncol=100)
for (i in 1:100){
  x=0
  while(x!=round(durationavg)){
    diseawithrandvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  Savg[,i]=resultrel[,1]
  Iavg[,i]=resultrel[,2]
  Ravg[,i]=resultrel[,3]
  Itotalavg[,i]=resultrel[,4]
}
avg=matrix(0,nrow = round(durationavg),ncol=4)
colnames(avg) <- c("S", "I", "R","4") 
for (i in 1:round(durationavg)){
  avg[i,1]=sum(Savg[i,])/100
  avg[i,2]=sum(Iavg[i,])/100
  avg[i,3]=sum(Ravg[i,])/100
  avg[i,4]=sum(Itotalavg[i,])/100
}
Ikum=numeric(round(durationavg))
for(i in 1:round(durationavg)){
  Ikum[i] = sum(avg[1:i, 4])
}
avg[,1:3]
plot(x=(1:round(durationavg)),y=avg[,1],ylab="Fraction",xlab="Time",type="l",ylim=c(0,1))     #Plotting S,I,R
lines(x=(1:round(durationavg)),y=avg[,2],col="red")
lines(x=(1:round(durationavg)),y=Ikum,col="purple")
lines(x=(1:round(durationavg)),y=avg[,3],col="green")

#Hub Vaccination
diseawithhubvacc=function(scalefreeNetwork,lambda,gamma,percVacc){
  N=dim(scalefreeNetwork)[1]
  V=matrix(rep(0,N),nrow=N,ncol=1)
  S=matrix(rep(1, N),nrow=N,ncol=1)               #Susceptibles
  I=matrix(rep(0, N),nrow=N,ncol=1)               #Infektios
  R=matrix(rep(0, N),nrow=N,ncol=1)               #Recovered
  maxcolsum=c(colSums(scalefreeNetwork))
  i=1
  while(i<=round(N*percVacc)){                    #Vaccine the Person with the highest sum in AdjMatrix because AdjMatrix is 1 if there is a contact and 0 else so the node with the highest sum is the node with the highest node degree
    y=which.max(maxcolsum)
    V[y,1]=1
    S[y,1]=0
    maxcolsum[y]=0
    i=i+1
  }
  Nonvacc=which(V[,1]==0)
  I1=sample(Nonvacc, size=1)                      #Choose one random person of the susceptible group which getting infected
  Diseafunktion(N,S,I,R,I1,scalefreeNetwork,lambda,gamma)
}
#Duration until the Disea dies out of the Network
duration=numeric(100)
for (i in 1:100){
  x=0
  while(x<=7){
    diseawithhubvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  duration[i]=x
}
durationavg=sum(duration)/100
round(durationavg)

#Average Values in each Period                      #Calculating the average relativ S,I,R for Diseas with average Duration
Savg=matrix(0,nrow=round(durationavg),ncol=100)
Iavg=matrix(0,nrow=round(durationavg),ncol=100)
Ravg=matrix(0,nrow=round(durationavg),ncol=100)
Itotalavg=matrix(0,nrow=round(durationavg),ncol=100)
for (i in 1:100){
  x=0
  while(x!=round(durationavg)){
    diseawithhubvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  Savg[,i]=resultrel[,1]
  Iavg[,i]=resultrel[,2]
  Ravg[,i]=resultrel[,3]
  Itotalavg[,i]=resultrel[,4]
}
avg=matrix(0,nrow = round(durationavg),ncol=4)
colnames(avg) <- c("S", "I", "R","4") 
for (i in 1:round(durationavg)){
  avg[i,1]=sum(Savg[i,])/100
  avg[i,2]=sum(Iavg[i,])/100
  avg[i,3]=sum(Ravg[i,])/100
  avg[i,4]=sum(Itotalavg[i,])/100
}
Ikum=numeric(round(durationavg))
for(i in 1:round(durationavg)){
  Ikum[i] = sum(avg[1:i, 4])
}
avg[,1:3]
plot(x=(1:round(durationavg)),y=avg[,1],ylab="Fraction",xlab="Time",type="l",ylim=c(0,1))     #Plotting S,I,R
lines(x=(1:round(durationavg)),y=avg[,2],col="red")
lines(x=(1:round(durationavg)),y=Ikum,col="purple")
lines(x=(1:round(durationavg)),y=avg[,3],col="green")

#Friendshipparadoxon
diseawithfriendvacc=function(scalefreeNetwork,lambda,gamma,percVacc){
  N=dim(scalefreeNetwork)[1]
  V=matrix(rep(0,N),nrow=N,ncol=1)
  S=matrix(rep(1, N),nrow=N,ncol=1)               #Susceptibles
  I=matrix(rep(0, N),nrow=N,ncol=1)               #Infektios
  R=matrix(rep(0, N),nrow=N,ncol=1)               #Recovered
  Ipers=sample(1:N,size=round(percVacc*N),replace=FALSE)    #Select a random group of individuals wich should be vaccined
  for(i in Ipers){
    partner=which(scalefreeNetwork[i,]==1)                  #Group of partner nodes of node i
    x=sample(partner,size=1)                                #Select a random neigbour
    while(V[x,1]==1){                                       #as long as the slected neighbour is Vaccinated
      x=sample(partner,size=1)
    }
    V[x,1]=1                                                #vaccinat the partner x of i which isn't vaccinaded yet
    S[x,1]=0
  }
  Nonvacc=which(V[,1]==0)
  I1=sample(Nonvacc, size=1)                      #Choose one random person of the susceptible group which getting infected
  Diseafunktion(N,S,I,R,I1,scalefreeNetwork,lambda,gamma)
}
#Duration until the Disea dies out of the Network
duration=numeric(100)
for (i in 1:100){
  x=0
  while(x<=7){                                              #Diseas with a period of 2 arent that relevant because there wasn't really a Disea
    diseawithfriendvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  duration[i]=x
}
durationavg=sum(duration)/100
round(durationavg)
#Average Values in each Period                      #Calculating the average relativ S,I,R for Diseas with average Duration
Savg=matrix(0,nrow=round(durationavg),ncol=100)
Iavg=matrix(0,nrow=round(durationavg),ncol=100)
Ravg=matrix(0,nrow=round(durationavg),ncol=100)
Itotalavg=matrix(0,nrow=round(durationavg),ncol=100)
for (i in 1:100){
  x=0
  while(x!=round(durationavg)){
    diseawithfriendvacc(scalefreeNetwork,0.15,0.5,0.05)
    x=length(resultrel[,1])
  }
  Savg[,i]=resultrel[,1]
  Iavg[,i]=resultrel[,2]
  Ravg[,i]=resultrel[,3]
  Itotalavg[,i]=resultrel[,4]
}
avg=matrix(0,nrow = round(durationavg),ncol=4)
colnames(avg) <- c("S", "I", "R","4") 
for (i in 1:round(durationavg)){
  avg[i,1]=sum(Savg[i,])/100
  avg[i,2]=sum(Iavg[i,])/100
  avg[i,3]=sum(Ravg[i,])/100
  avg[i,4]=sum(Itotalavg[i,])/100
}
Ikum=numeric(round(durationavg))
for(i in 1:round(durationavg)){
  Ikum[i] = sum(avg[1:i, 4])
}
avg[,1:3]
plot(x=(1:round(durationavg)),y=avg[,1],ylab="Fraction",xlab="Time",type="l",ylim=c(0,1))     #Plotting S,I,R
lines(x=(1:round(durationavg)),y=avg[,2],col="red")
lines(x=(1:round(durationavg)),y=Ikum,col="purple")
lines(x=(1:round(durationavg)),y=avg[,3],col="green")
