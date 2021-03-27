
library(mirt)
library(lavaan)
library(progress)

################################################################
# Evaluating Local Dependence (LD) Using Pairwise LD Indices in Psychological Data
## PSYC 859 Item Response Theory Project 
## Fall 2019
################################################################


# Simulation of IRT data using population parameters from 
# Lui et al., 2012 and Edwards, 2018
### Simulate data 1000 times

######################################3

#### DATA LOOP

data <- list()
set.seed(7)
f1 <- rnorm(6, mean = 1.7, sd = .3)
f2 <- rnorm(3, mean = 2.55, sd = .15)
a <- matrix(c(f1, rep(0,3), f2), ncol = 2) # Draw a parameters from dist
d <- rnorm(6, mean = 0, sd = 1) # Draw b parameters from N(0,1)
nd <-  1000 # Number of replications



###########################################################################################
# Use gnerated data and compute local dependence using 
# Pearson's ??2 statistic(P - ??2), jackknife slope index (JSI), and modification indicies (MI). 

###################################
### Pearson's ??2 statistic(P - ??2)
###################################

residuals <- list()
v <- rep(1,nd) # Empty Vector for LD Chi
v1 <- rep(1,nd)
v2 <- rep(1,nd)

for (i in 1:nd){  
  
  #MIRT Data Simulation
  fun <- simdata( a= a, d=d,  N = 1000, itemtype = "dich" )
  data <- c(data,list(fun))
  
  #replications of data frame
  f1 <- data.frame(data[[i]])
  perf <- mirt(f1, 1, itemtype = '2PL')
  
  # coef(perf, simplify=T, IRTpars=T)
  r1 <- residuals(perf, type = "LD")
  residuals <- c(residuals, list(r1))
  
  # Means of Pairs 
  
  
  for(j in 1:i){
    v[] <- lapply(residuals, "[", n = 2) # Mean of pair (1,2)
    v1[] <- lapply(residuals, "[", n = 16) # Mean of pair (3,4)
    v2[] <- lapply(residuals, "[", n = 30) #Mean of pair (5,6)
    mean <- mean(as.numeric(v[]))
    mean1 <- mean(as.numeric(v1[1:j]))
    mean2 <- mean(as.numeric(v2[1:j]))
    
  }}
## ETA
pb <- progress_bar$new(
  format = " downloading [:bar] :percent eta: :eta",
  total = 100, clear = FALSE, width= 60)
for (i in 1:100) {
  pb$tick()
  Sys.sleep(1 / 100)
}



##############################
# jackknife slope index (JSI)
##########################

jsi <- list()
j <- rep(1,nd) # Empty Vector for JSI
j1 <- rep(1,nd)
j2 <- rep(1,nd)

for (i in 1:nd){  
  
  #replications of data frame
  kj <- data.frame(data[[i]])
  perf <- mirt(kj, 1, itemtype = '2PL')
  
  
  r2 <- residuals(perf, type = "JSI", fold = T)
  jsi <- c(jsi, list(r2))
  
  for(k in 1:i){
    j[] <- lapply(residuals, "[", n = 2) # Mean of pair (1,2)
    j1[] <- lapply(residuals, "[", n = 16) # Mean of pair (3,4)
    j2[] <- lapply(residuals, "[", n = 30) #Mean of pair (5,6)
    Jmean <- mean(as.numeric(j[]))
    Jmean1 <- mean(as.numeric(j1[1:k]))
    Jmean2 <- mean(as.numeric(j2[1:k]))
  }}



#### Add elements of JSImatrix to calc threshold values to detect LD

k <- NULL

for (i in 1:nd){
  ki <- (c(jsi[[i]][2:6], jsi[[i]][3:6], jsi[[i]][5:6]))
  k <- cbind(k, ki)
}


info <- NULL
for (p in 1:nd){
  mean <- mean(k[,p])
  sd <- sd(k[,p])
  thresh <- mean +2*sd
  info <- rbind(c(mean, sd, thresh), info)
  colnames(info) <- c("mean", "sd", "tresh")
}


###################################################


#r1 <- data.frame(r1)
#print(paste0("LD Chi Mean of Pair (1,2)", mean))
#print(paste0("LD Chi Mean of Pair (3,4)", mean1 ))
#print(paste0("LD Chi Mean of Pair (5,6)", mean2 ))


# Tell me how many LD are greater than a value

#print(paste0("LD <0 Pair (1,2)", sum(v<0)))
#print(paste0("LD <0 (3,4)", sum(v1<0)))
#print(paste0("LD <0 (5,6)", sum(v2<0) ))  


#####################################################
# Modification Indicies (MI)
#####################################################

#fit structural equation model using generated data
t1 <- '
F1=~    Item_1 + Item_2 + Item_3 + Item_4 + Item_5 + Item_6

'
#### Extract all MI
mi <- list()
for (k in 1:nd){
  frame <- data.frame(data[[k]])
  tf1 <- cfa(t1, frame)
  
  if (!is.null(mod <- tryCatch(modindices(tf1), 
                               error=function(e){cat("ERROR", conditionMessage(e), "\n")}))) {
    
    mi[[k]] <- mod$mi }
}

nullToNA <- function(x) {
  x[sapply(x, is.null)] <- 0
  return(x)
}
miData <- nullToNA(mi)
#Create data frame of all MI across data
miDATA <- data.frame(miData)
#Rename variables in MI frame
for (i in 1:nd){
  names <- c(paste0("MI", seq(1:i)))
}
colnames(miDATA) <- names


### extract all MI of Interest
sum <- sapply(miDATA,  "[", n = 1)
sum1 <- sapply(miDATA,  "[", n = 10)
sum2 <- sapply(miDATA,  "[", n = 15)

meanMI <- mean(sum)
meanMI1 <- mean(sum1)
meanMI2 <- mean(sum2) 

#### Print number of MI with Criteria greater than value 
#print(paste0("MI >0  Pair (1,2)" sum(sum >0)))
#print(paste0("MI >0 Pair (3,4)", sum(sum1 >0) ))
#print(paste0("MI >0 Pair (5,6)", sum(sum2>0) ))

##########################################################################


###Subject to ALTERATION
c <-  3.84#Criteria to retain JSI 
c2 <- 3.84 # Criteria to retain MI



### Criteria for JSI - Mean + 2*SD for each sample
numb <- NULL
for (m in 1:nd){
  if (j[[m]] >= info[m,3]){
    numb[m] <- T
  }else{numb[m] <- F}
}
numb1 <- NULL
for (m in 1:nd){
  if (j1[[m]] >= info[m,3]){
    numb1[m] <- T
  }else{numb1[m] <- F}
}

numb2 <- NULL
for (m in 1:nd){
  if (j2[[m]] >= info[m,3]){
    numb2[m] <- T
  }else{numb2[m] <- F}
}


# Produce table 
table <- matrix(c( mean,    mean1,   mean2, 
                   meanMI,  meanMI1, meanMI2,
                   mean(as.numeric(j)), mean(as.numeric(j1)),mean(as.numeric(j2)),
                   sum(v >c),sum(v1 >c),  sum(v2 >c), 
                   sum (sum >c2),   sum (sum1 >c2),sum (sum2 >c2),
                   sum(numb), sum(numb1), sum(numb2)
), 3, 6)


colnames(table) <- c("MeanLD","MeanMI" , "MeanJSI", "LD","MI", "JSI")
rownames(table) <- c("Pair (1,2)","Pair (3,4)","Pair (5,6)")

table <- as.table(table)
table

#############################################################################################################












