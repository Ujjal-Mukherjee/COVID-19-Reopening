#
# Simulation
#

##
## Initiate the process
##



# N total number of individuals
# I0 total number of initial infections 
# ContactMatrix is the matrix of pairs of index numbers who have contacted at a given period

SEIRModelInit = function(N, I0){
  S = seq(1,N)
	status = rep(0,N)
	index = sample.int(N, size = I0, replace = FALSE)
	status[index] = rep(1, length(index))
	ContactMatrix = matrix(0, nrow=0, ncol=2)
	return(list(S=S, status = status, ContactMatrix = ContactMatrix))
}

#SEIRModelInit(50,5)

# Testing is a fuction that updates status from 1 to -1
# P Master List structure of individuals. 2 vectors and 1 matrix. 
#     1st element vector s: id of individuals
#     2nd element vector status: infection status of individuals. 0  if susceptible, 1 if infected but not detected, -1 detected and removed 
#     3rd element matix with two columns: col 1 is id of individual, col 2 is individual. They were in contact in the last period. Contacts are genereted randomly. 
# c is a fraction of the tests that are used for contact tracing 

Testing = function(P, T, c){

	S = P$S
	status = P$status
	index = which(status != -1)

	# Test is a vector of the individuals that get tested in a period	
	if((1-c)*T<length(index)){
		ind = sample(index, size = ceiling((1-c)*T), replace = FALSE)	
		Test = index[ind]
	}else{
		Test = index
	}
	
	Detect = intersect(Test, which(status == 1))
	if(length(Detect)>0)status[Detect]=rep(-1,length(Detect))

	#Contact Tracing and Testing

	if(length(Detect)>0){

		Contacts = c()
		ContactMatrix = P$ContactMatrix

		for(i in 1:length(Detect)){

			ind1 = which(ContactMatrix[,2] == Detect[i])
			if(length(ind1)>0){
				Contacts = c(Contacts, ContactMatrix[ind1, 1])		
			}
		}
		
		Cont = min(floor(c*T), length(Contacts))
		if(length(Contacts)>0){
  		if(length(Cont)>0){
			if(length(Contacts)>length(Cont)){
  				ContactsSelected = sample(Contacts, size=length(Cont))
  				ind2 = which(status[ContactsSelected]==1)
			}else{
				ContactsSelected = Contacts
				ind2 = which(status[ContactsSelected]==1)
			}
  			if(length(ind2)>0)status[ContactsSelected[ind2]]=rep(-1, length(ind2))

				# Test Remining if any is a vector of the individuals that get tested in a period	
				if(c*T>length(ContactsSelected)){
					remainingTests = - length(ContactsSelected) + c*T
					index = which(status != -1)
					if(length(index)>remainingTests){
						ind = sample(index, size = ceiling(remainingTests), replace = FALSE)	
						Test1 = index[ind]
					}else{
						Test1 = index
					}
					Detect = intersect(Test1, which(status == 1))
					if(length(Detect)>0)status[Detect]=rep(-1,length(Detect))

				}  
  		}
		}
	}
	return(list(S=S, status = status, ContactMatrix = P$ContactMatrix))
}


# Infection is a fucntion that updates status and ContactMatrix 
# P Master List structure of individuals. 2 vectors and 1 matrix. 
#     1st element vector s: id of individuals
#     2nd element vector status: infection status of individuals. 0  if susceptible, 1 if infected but not detected, -1 detected and removed 
#     3rd element matix with two columns: col 1 is id of individual, col 2 is individual. They were in contact in the last period. Contacts are genereted randomly. 
# p probablility that a susceptible individual becomes infected if she meet an infected individual
# M meassure of mobility. Number of other unique individuals that one indivual comes into contact per unit of time  

Infection = function(P, p, M){

	S = P$S
	status = P$status
	ContactMatrix = matrix(0, nrow=0, ncol=2) #P$ContactMatrix
		
	index = which(status != -1)
	M = min(M, length(index)-1)
	
	if(M>0){
		ind1 = which(status == 1)
    for(i in 1:length(status)){
			if(status[i] == 0){
					ind = sample(index[-i], size = M)
					indint = intersect(ind, ind1)
					ln_int = length(indint)
					p1 = 1 - (1-p)^ln_int
					stat_new = sample(c(0,1), size=1, prob=c(1-p1, p1))
					status[i] = stat_new
					ContactMatrix = rbind(ContactMatrix, cbind(rep(i, length(ind)), ind))
			}
		}
	}
	return(list(S=S, status = status, ContactMatrix = ContactMatrix))
}


# Simulate is a function that runs an Agent Base infection simulation 
# T is the number of tests
# c is a fraction of the tests that are used for contact tracing 
# p probablility that a susceptible individual becomes infected if she meet an infected individual
# M meassure of mobility. Number of other unique individuals that one indivual comes into contact per unit of time  
# t is the number of time periods to simulate for
Simulate = function(N, I0, T, c, p, M, t){
  
  # Initialization
	P = SEIRModelInit(N, I0)
	nS = c()
	nI = c()
	nR = c()
	for(i in 1:t){
		nS = c(nS, length(which(P$status == 0)))	
		nI = c(nI, length(which(P$status == 1)))
		nR = c(nR, length(which(P$status == -1)))
		P = Infection(P, p, M)
		P = Testing(P, T, c)
	}
	return(list(Suseptibles = nS, Infected = nI, Removed = nR))
}

#
# N=500; I0 = 10; T = N*0.05; c = 0.5; p=0.1; M=5; t=10
#

#N = 5000; I0 = 5; T = N*0.1; c = 0.7; p=0.1; M=5; t=50

B = 100

N = 5000; t=50; I0=5

seq = 1

setwd("C:/Users/ukm/Box/COVID19/CodeAndData/R")

for(ii in c(0.05, 0.1, 0.2, 0.3, 0.5)){
for(c in c(0.1, 0.3, 0.5, 0.7, 0.9)){
for(p in c(0.05, 0.1, 0.15, 0.2, 0.25)){
for(M in c(1, 5, 10, 15, 20)){

#png(paste("Simulation_Runs_Continious_",seq,".png",sep=""), width = 960, height = 640)

T=N*ii

Suseptibles = matrix(nrow = t, ncol = B)
Infected = matrix(nrow = t, ncol = B)
Removed = matrix(nrow = t, ncol = B)

for(i in 1:B){

	Sim1 = Simulate(N=N, I0 = I0, T = T, c = c, p=p, M=M, t=t)

	Suseptibles[,i] = Sim1$Suseptibles
	Infected[,i] = Sim1$Infected
	Removed[,i] = Sim1$Removed

}

# Sim1= Simulate(N=50000, I0 = 10, T = 10000, c = 0.7, p=0.1, M=5, t=60)

#plot(Sim1$Suseptibles, type="b", col = 4, xlim=c(0,t), ylim=c(0,N), xlab="Time (Days)", ylab="Disease Dynamics",
#		main=paste("N:",N,"; Io:",I0,"; T:",T,"; c:",c,"; p:",p,"; M:",M,sep=""))
#lines(Sim1$Infected, type="b", col = 2)
#lines(Sim1$Removed, type="b", col = 3)
#legend("left", legend = c("Susceptibles", "Infected", "Removed"), lty=c(1,1,1), pch=c(1,1,1), col=c(4,2,3))

#library(ggplot2)

SusArea = apply(Suseptibles, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Susceptibles", dim(Suseptibles)[1]), 
				Median = Suseptibles[,order(SusArea)[floor(0.5*B)]], 
				Lower = Suseptibles[,which(SusArea == min(SusArea))[1]],
				Upper = Suseptibles[,which(SusArea == max(SusArea))[1]]))

plot(as.numeric(as.character(DataPlot$Days)), as.numeric(as.character(DataPlot$Median)), type="l", col=4, ylim=c(0, N), xlab="Days from Start", ylab="Number of Individuals",
		main=paste("N:",N,"; Io:",I0,"; T:",T,"; c:",c,"; p:",p,"; M:",M,sep=""))
lines(as.numeric(as.character(DataPlot$Lower)), col=4, lty=2)
lines(as.numeric(as.character(DataPlot$Upper)), col=4, lty=3)

#InfArea = apply(Infected, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Infected", dim(Suseptibles)[1]), 
					Median = Infected[,order(SusArea)[floor(0.5*B)]], 
					Lower = Infected[,which(SusArea == min(SusArea))[1]],
					Upper = Infected[,which(SusArea == max(SusArea))[1]]))

lines(as.numeric(as.character(DataPlot$Days)), as.numeric(as.character(DataPlot$Median)),  col=2)
lines(as.numeric(as.character(DataPlot$Lower)), col=2, lty=2)
lines(as.numeric(as.character(DataPlot$Upper)), col=2, lty=3)


#RemArea = apply(Removed, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Removed", dim(Suseptibles)[1]), 
					Median = Removed[,order(SusArea)[floor(0.5*B)]], 
					Lower = Removed[,which(SusArea == min(SusArea))[1]],
					Upper = Removed[,which(SusArea == max(SusArea))[1]]))


lines(as.numeric(as.character(DataPlot$Days)), as.numeric(as.character(DataPlot$Median)),  col=1)
lines(as.numeric(as.character(DataPlot$Lower)), col=1, lty=2)
lines(as.numeric(as.character(DataPlot$Upper)), col=1, lty=3)

legend("left", c("Not Infected", "Undetected Infected", "Detected Infected", "Optimistic Scenario", "Likely Scenario", "Pessimistic Scenario"),
			lty = c(1,1,1,3,1,2), col=c(4,2,1,1,1,1)) 


#dev.off()

seq = seq + 1

}}}}


#################################################

B = 100

N = 5000; t=50; I0=5

T=c(250,250,250,250,250,250,250,250,250,500,1000,2000,1000,1500,2500,250,1000,2500,1000,1000,1000,1000)
c=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.3,0.5,0.7)
p=c(0.05,0.05,0.05,0.05,0.20,0.15,0.10,0.05,0.1,0.1,0.1,0.1,0.05,0.1,0.15,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
M=c(20,10,5,1,5,5,5,5,10,10,10,10,5,5,5,1,5,15,10,10,10,10)

d=cbind(T=T, c=c, p=p, M=M)

setwd("C:/Users/ukm/Box/COVID19/CodeAndData/R")

for(j in 1:dim(d)[1]){

Suseptibles = matrix(nrow = t, ncol = B)
Infected = matrix(nrow = t, ncol = B)
Removed = matrix(nrow = t, ncol = B)

for(i in 1:B){

	Sim1 = Simulate(N=N, I0 = I0, T = T[j], c = c[j], p=p[j], M=M[j], t=t)

	Suseptibles[,i] = Sim1$Suseptibles
	Infected[,i] = Sim1$Infected
	Removed[,i] = Sim1$Removed

}

save(Suseptibles, file=paste("Raw_Susceptibles_",j,".RData",sep=""))
save(Infected, file=paste("Raw_Infected_",j,".RData",sep=""))
save(Removed, file=paste("Raw_Removed_",j,".RData",sep=""))


SusArea = apply(Suseptibles, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Susceptibles", dim(Suseptibles)[1]), 
				Median = Suseptibles[,order(SusArea)[floor(0.5*B)]], 
				Lower = Suseptibles[,which(SusArea == min(SusArea))[1]],
				Upper = Suseptibles[,which(SusArea == max(SusArea))[1]]))

save(DataPlot, file=paste("Susceptibles_",j,".RData",sep=""))

#InfArea = apply(Infected, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Infected", dim(Suseptibles)[1]), 
					Median = Infected[,order(SusArea)[floor(0.5*B)]], 
					Lower = Infected[,which(SusArea == min(SusArea))[1]],
					Upper = Infected[,which(SusArea == max(SusArea))[1]]))

save(DataPlot, file=paste("Infected_",j,".RData",sep=""))

#RemArea = apply(Removed, 2, FUN=sum)
DataPlot = as.data.frame(cbind(Days = seq(1, dim(Suseptibles)[1]), TimeSeries = rep("Removed", dim(Suseptibles)[1]), 
					Median = Removed[,order(SusArea)[floor(0.5*B)]], 
					Lower = Removed[,which(SusArea == min(SusArea))[1]],
					Upper = Removed[,which(SusArea == max(SusArea))[1]]))


save(DataPlot, file=paste("Recovered_",j,".RData",sep=""))

}


######################################################################
#	CREATE COMBINED IMAGES
######################################################################

rm(list=ls())

library(ggplot2)

#################################################
#Figure 1. Mobility 
#################################################

load("Susceptibles_1.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))

load("Susceptibles_2.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Susceptibles_3.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


load("Susceptibles_4.RData")

DataPlot4 = DataPlot

DataPlot4$Days=as.numeric(as.character(DataPlot4$Days))
DataPlot4$Median=as.numeric(as.character(DataPlot4$Median))
DataPlot4$Lower=as.numeric(as.character(DataPlot4$Lower))
DataPlot4$Upper=as.numeric(as.character(DataPlot4$Upper))


DataPlot1$Mobility = rep("1. M=20", dim(DataPlot1)[1])

DataPlot2$Mobility = rep("2. M=10", dim(DataPlot2)[1])

DataPlot3$Mobility = rep("3. M=5", dim(DataPlot3)[1])

DataPlot4$Mobility = rep("4. M=1", dim(DataPlot4)[1])


DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3, DataPlot4)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=Mobility)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)

####################################################################
#Figure 0. Example Output
####################################################################

load("Susceptibles_2.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))


load("Infected_2.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Recovered_2.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=TimeSeries)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)

########################################################################
#Mask Wearing
########################################################################

load("Susceptibles_5.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))

load("Susceptibles_6.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Susceptibles_7.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


load("Susceptibles_8.RData")

DataPlot4 = DataPlot

DataPlot4$Days=as.numeric(as.character(DataPlot4$Days))
DataPlot4$Median=as.numeric(as.character(DataPlot4$Median))
DataPlot4$Lower=as.numeric(as.character(DataPlot4$Lower))
DataPlot4$Upper=as.numeric(as.character(DataPlot4$Upper))


DataPlot1$Mask_Adherence = rep("1. m=10% (p=0.20)", dim(DataPlot1)[1])

DataPlot2$Mask_Adherence = rep("2. m=40% (p=0.15)", dim(DataPlot2)[1])

DataPlot3$Mask_Adherence = rep("3. m=70% (p=0.10)", dim(DataPlot3)[1])

DataPlot4$Mask_Adherence = rep("4. m=90% (p=0.05)", dim(DataPlot4)[1])


DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3, DataPlot4)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=Mask_Adherence)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)


########################################################################
#TESTING 1
########################################################################

load("Susceptibles_9.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))

load("Susceptibles_10.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Susceptibles_11.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


load("Susceptibles_12.RData")

DataPlot4 = DataPlot

DataPlot4$Days=as.numeric(as.character(DataPlot4$Days))
DataPlot4$Median=as.numeric(as.character(DataPlot4$Median))
DataPlot4$Lower=as.numeric(as.character(DataPlot4$Lower))
DataPlot4$Upper=as.numeric(as.character(DataPlot4$Upper))


load("Susceptibles_15.RData")

DataPlot5 = DataPlot

DataPlot5$Days=as.numeric(as.character(DataPlot5$Days))
DataPlot5$Median=as.numeric(as.character(DataPlot5$Median))
DataPlot5$Lower=as.numeric(as.character(DataPlot5$Lower))
DataPlot5$Upper=as.numeric(as.character(DataPlot5$Upper))


DataPlot1$Tests_Per_Day = rep("1. T=250", dim(DataPlot1)[1])

DataPlot2$Tests_Per_Day = rep("2. T=500", dim(DataPlot2)[1])

DataPlot3$Tests_Per_Day = rep("3. T=1000", dim(DataPlot3)[1])

DataPlot4$Tests_Per_Day = rep("4. T=2000", dim(DataPlot4)[1])

DataPlot5$Tests_Per_Day = rep("5. T=2500", dim(DataPlot5)[1])


DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3, DataPlot4, DataPlot5)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=Tests_Per_Day)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)



########################################################################
#TESTING 1
########################################################################

load("Susceptibles_13.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))

load("Susceptibles_14.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Susceptibles_15.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


load("Susceptibles_16.RData")

DataPlot4 = DataPlot

DataPlot4$Days=as.numeric(as.character(DataPlot4$Days))
DataPlot4$Median=as.numeric(as.character(DataPlot4$Median))
DataPlot4$Lower=as.numeric(as.character(DataPlot4$Lower))
DataPlot4$Upper=as.numeric(as.character(DataPlot4$Upper))


load("Susceptibles_17.RData")

DataPlot5 = DataPlot

DataPlot5$Days=as.numeric(as.character(DataPlot5$Days))
DataPlot5$Median=as.numeric(as.character(DataPlot5$Median))
DataPlot5$Lower=as.numeric(as.character(DataPlot5$Lower))
DataPlot5$Upper=as.numeric(as.character(DataPlot5$Upper))


load("Susceptibles_18.RData")

DataPlot6 = DataPlot

DataPlot6$Days=as.numeric(as.character(DataPlot6$Days))
DataPlot6$Median=as.numeric(as.character(DataPlot6$Median))
DataPlot6$Lower=as.numeric(as.character(DataPlot6$Lower))
DataPlot6$Upper=as.numeric(as.character(DataPlot6$Upper))


DataPlot1$Parameters = rep("1. T=1000, p=0.05, M=5", dim(DataPlot1)[1])

DataPlot2$Parameters = rep("2. T=1500, p=0.1, M=5", dim(DataPlot2)[1])

DataPlot3$Parameters = rep("3. T=2500, p=0.1, M=5", dim(DataPlot3)[1])

DataPlot4$Parameters = rep("4. T=250, p=0.05, M=1", dim(DataPlot4)[1])

DataPlot5$Parameters = rep("5. T=1000, p=0.05, M=5", dim(DataPlot5)[1])

DataPlot6$Parameters = rep("5. T=2500, p=0.05, M=15", dim(DataPlot6)[1])


DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3, DataPlot4, DataPlot5, DataPlot6)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=Parameters)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)




########################################################################
#CONTACT TRACING
########################################################################

load("Susceptibles_19.RData")

DataPlot1 = DataPlot

DataPlot1$Days=as.numeric(as.character(DataPlot1$Days))
DataPlot1$Median=as.numeric(as.character(DataPlot1$Median))
DataPlot1$Lower=as.numeric(as.character(DataPlot1$Lower))
DataPlot1$Upper=as.numeric(as.character(DataPlot1$Upper))

load("Susceptibles_20.RData")

DataPlot2 = DataPlot

DataPlot2$Days=as.numeric(as.character(DataPlot2$Days))
DataPlot2$Median=as.numeric(as.character(DataPlot2$Median))
DataPlot2$Lower=as.numeric(as.character(DataPlot2$Lower))
DataPlot2$Upper=as.numeric(as.character(DataPlot2$Upper))


load("Susceptibles_21.RData")

DataPlot3 = DataPlot

DataPlot3$Days=as.numeric(as.character(DataPlot3$Days))
DataPlot3$Median=as.numeric(as.character(DataPlot3$Median))
DataPlot3$Lower=as.numeric(as.character(DataPlot3$Lower))
DataPlot3$Upper=as.numeric(as.character(DataPlot3$Upper))


load("Susceptibles_22.RData")

DataPlot4 = DataPlot

DataPlot4$Days=as.numeric(as.character(DataPlot4$Days))
DataPlot4$Median=as.numeric(as.character(DataPlot4$Median))
DataPlot4$Lower=as.numeric(as.character(DataPlot4$Lower))
DataPlot4$Upper=as.numeric(as.character(DataPlot4$Upper))


DataPlot1$Contact_Tracing = rep("1. CT = 0.1", dim(DataPlot1)[1])

DataPlot2$Contact_Tracing = rep("2. CT = 0.3", dim(DataPlot2)[1])

DataPlot3$Contact_Tracing = rep("3. CT = 0.5", dim(DataPlot3)[1])

DataPlot4$Contact_Tracing = rep("4. CT = 0.7", dim(DataPlot4)[1])

DataPlot = rbind(DataPlot1, DataPlot2, DataPlot3, DataPlot4)

p<-ggplot(data=DataPlot, aes(x=Days, y=Median, colour=Contact_Tracing)) + geom_line(size=1)
p<-p+geom_ribbon(data=DataPlot, aes(ymin=Lower, ymax=Upper), linetype=2, alpha=0.1, size=1)
p<-p+xlab("Days after Reopening") + ylab("Number of Not-Infected Individuals")+ 
	theme(
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x= element_text(size = 15, face="bold"),
	axis.text.y= element_text(size=15, face="bold"),
	legend.text=element_text(size=15, face="bold"),
		legend.title=element_text(size=15, face="bold")
  )

plot(p)










