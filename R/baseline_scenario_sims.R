# Stochastic parameters

# Our modelled height values

Hmax <- c(145.42035, 103.65811, 115.94892, 89.81047, 459.25922, 64.70510, 48.75292, 81.85831, 258.23389, 267.63473, 81.73318, 42.15995, 211.01529, 62.83146, 86.18343, 84.96240, 57.87685, 41.01645, 143.17702, 59.73425)
a <- c(4.226643, 4.543715, 3.849143, 3.226788, 4.580450, 4.489764, 4.770788, 
       3.483222, 6.622986, 6.644306, 4.597916, 4.892076, 4.385591, 4.410173, 
       4.462377, 4.700410, 4.362824, 5.073423, 4.689330, 3.371535)
b <- c(1.720544, 1.555347, 1.818385, 2.122570, 1.818912, 1.441581, 1.265560, 
       1.927080, 1.133050, 1.136036, 1.551051, 1.203531, 1.670703, 1.481089, 
       1.515627, 1.431935, 1.479525, 1.149516, 1.570982, 1.903793)
Time <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
row_names <- c('Acabra', 'Acamon', 'Acawil', 'Beyopa', 'Codon', 'Dodbur', 
                   'Erecrass', 'Eregla', 'EucBlue', 'EucGreen', 'Grehue', 
                   'Halcya', 'Mellan', 'Olemul', 'Olepim', 'Olesub', 
                   'Phebalium', 'ProstantheraGreen', 'Senart', 'Wesrig')
mu <- matrix(NA, nrow = length(Hmax), ncol=length(Time))
N.species <- dim(mu)[1]
row.names(mu) <- row_names
colnames(mu) <- TSFvalues
TSFvalues <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
detection_prob <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
occupancy <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

Occupancy_meta_matrix <- BuildMatrix(occupancy)
p_matrix <- BuildMatrix(detection_prob)

# Strategic sampling with a PriorityList of sites to visit
PriorityList <- matrix(NA, 1400, 2)
PriorityList[, 1]  <- rep(c(2, 8, 26, 41, 4, 15, 33, 86, 1, 6, 13, 3, 36, 28), 100)
PriorityList[, 2] <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 100)
colnames(PriorityList) <- c('TSF', 'Rep')

indexes <- list(7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 111, 150, 222) 
value_seqs <- list(c(rep(0,8), rep(1,3), rep(2,1),rep(3,1), rep(5,1), rep(10,1)),
                   c(rep(0,16), rep(1,5), rep(2,4), rep(3,2), rep(4,1), rep(5,1)),
                   c(rep(0,26), rep(1,8), rep(2,5), rep(3,2), rep(4,1), rep(5,1)),
                   c(rep(0,34), rep(1,12), rep(2,6),rep(3,2), rep(4,1), rep(5,1), rep(10,1)),
                   c(rep(0,43), rep(1,14), rep(2,7),rep(3,3), rep(4,1), rep(5,1), rep(6,1), rep(10,1)),
                   c(rep(0,51), rep(1,17), rep(2,9),rep(3,4), rep(4,1), rep(5,1), rep(6,1), rep(10,1)),
                   c(rep(0,60), rep(1,20), rep(2,10), rep(3,3), rep(4,1), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,70), rep(1,20), rep(2,11), rep(3,4), rep(4,2), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,77), rep(1,25), rep(2,13), rep(3,4), rep(4,3), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,85), rep(1,28), rep(2,14), rep(3,5), rep(4,4), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,93), rep(1,31), rep(2,16), rep(3,6), rep(4,3), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,102), rep(1,34), rep(2,17), rep(3,6), rep(4,4), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,110), rep(1,37), rep(2,18), rep(3,6), rep(4,4), rep(5,3), rep(6,2), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,118), rep(1,40), rep(2,20), rep(3,7), rep(4,5), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,127), rep(1,42), rep(2,21), rep(3,8), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,134), rep(1,45), rep(2,22), rep(3,9), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2)),
                   c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
                  )

######################################################################################
# Deterministic parameters

HrsPerDay <- 10 # long days
FieldMinsDay <- HrsPerDay * 60   
sd.obs <- 0.395
N.fieldtrips <- 100 
HometoSiteTravel <- (1 * 2) # 1 minute to get to sites and home
setuptime <- 1 # 1 minute to set up
mt <- 1 # 1 minute to measure 
TravelbwSites <- 1 # 1 minute to travel between sites
Maxtime <- 60 * 30   # 30 hours at one replicate
N.reps <- 1 # one replicate
N.TSF <- length(TSFvalues)

for (j in 1:length(Hmax)) { 
  for (x in 1:length(Time)) {
    mu[j, x] <- Hmax[j] / (1 + exp(-a[j] * (Time[x] - b[j])))  
  }  
}

#####################################################################################

BuildMatrix <- function(row_values) {
  matrx <- matrix(NA, nrow=N.species, ncol=length(TSFvalues))
  row.names(matrx) <- row_names
  colnames(matrx) <- TSFvalues
  for (w in 1:length(TSFvalues)) {
    matrx[,w] <- occupancy
  }
  return(matrx)
}

BuildFieldDays <- function(days, value_seq) {
  values <- matrix(0, nrow = 2, ncol = days * 2 + 1)
  values[2, ] <- value_seq 
  for (q in 2:length(values[1, ])) { 
    values[1, q] <- values[1, q - 1] + 0.5
  } 
  values[2,] <- values[2,]/sum(values[2,])
  for (q in 2:length(values[2, ])) {
    values[2, q] <- values[2, q - 1] + values[2, q]
  }
  return(values)
}

BuildScenarios <- function(days, values) {
  list(mu, values, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, 
       PriorityList, N.fieldtrips, TotalFieldDays=days, TSFvalues, N.reps, 
       setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, 
       N.TSF)
}

field_day_values <- Map(BuildFieldDays, indexes, value_seqs)
scenarios <- Map(BuildScenarios, indexes, field_day_values)


##################################

# collecting between 2 and 3 individuals per species

Fieldtrip <- function(Scenario) { 
  
  as.numeric(Sys.time())-> g; set.seed((g - floor(g)) * 1e8 -> seed); print(seed)
  
  N.species <- dim(mu)[1]
  
  species.ID <- c(rep(0,20))
  
  Rando <- runif(1, 0, 1)
  
  RealIndex <- min(which(Rando <= Scenario[[2]][2,])) 
  
  RealFieldDays <- Scenario[[2]][1, RealIndex] 
  
  RealFieldDaysMins <- RealFieldDays * Scenario[[3]] 
  
  TotalSiteTime <- RealFieldDaysMins - Scenario[[4]] * RealFieldDays  
  
  TotalSites <- round(TotalSiteTime/ (Scenario[[5]] + Scenario[[6]]))
  
  rm(Rando, RealIndex, RealFieldDays, RealFieldDaysMins, TotalSiteTime) 
  
  
  FittingData <- array(NA, c(40000, 4))
  colnames(FittingData) <- c('Sp', 'Height', 'TSF', 'Rep')
  
  N.plants <- 0
  
  Findsthrutime <- array(NA, c(40000, 6, N.reps, N.TSF))
  colnames(Findsthrutime) <- c('SpeciesID', ' ', ' ', ' ', 'Time(mins)','Height(cm)')
  
  SummaryData <- array(NA, c(N.species, 6, N.reps, N.TSF))
  colnames(SummaryData) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')
  
  
  for (t in 1:TotalSites) {
    
    TSF <- Scenario[[7]][t, 1]
    
    TSFIndex <- which(Scenario[[10]] == TSF)
    
    Rep <- Scenario[[7]][t, 2]
    
    matrix_ind_counts <- matrix(0, N.species, 1)
    matrix_test <- matrix(NA, 40000, 6)   
    colnames(matrix_test) <- c('ChosenSpecies','indivID',  'TSF',  'Rep', 'time_mins','H')
    
    Presence <- rbinom(N.species, 1, Scenario[[16]][,TSFIndex])
    
    Searchable <- Presence
    
    timeatsite <- 0
    
    setuptime <- Scenario[[12]] 
    timeatsite <- timeatsite + setuptime
    
    datarowcounter <- 0
    
    
    PossDetectiontime <- rexp(N.species, 1/Scenario[[15]][,TSFIndex])
    PossMinDetectiontime <- min(PossDetectiontime[Searchable==1]) 
    
    
    while((1 - prod(matrix_ind_counts[Presence==1] >= 2 )) & (timeatsite + PossMinDetectiontime + Scenario[[13]] < Scenario[[5]])) {  
      
      Detectiontime <- PossMinDetectiontime
      
      ChosenSpecies <- which(Detectiontime == PossDetectiontime)
      species.ID[ChosenSpecies] <- 1
      
      datarowcounter <- datarowcounter + 1
      
      N.plants <- N.plants + 1
      
      matrix_ind_counts[ChosenSpecies, 1] <- matrix_ind_counts[ChosenSpecies, 1] + 1
      
      matrix_test[datarowcounter, 1] <- ChosenSpecies
      matrix_test[datarowcounter, 2] <- N.plants
      matrix_test[datarowcounter, 3] <-  TSF
      matrix_test[datarowcounter, 4] <-  Rep
      
      k <- log(mu[ChosenSpecies, TSFIndex]) - Scenario[[17]]^2
      
      H <- rlnorm(1, k, Scenario[[17]])
      
      FittingData[N.plants, 1] <- ChosenSpecies
      FittingData[N.plants, 2] <- H
      FittingData[N.plants, 3] <- TSFvalues[TSFIndex]
      FittingData[N.plants, 4] <- Rep
      
      
      matrix_test[datarowcounter, 6] <- H
      
      timeatsite <- timeatsite + Detectiontime + Scenario[[13]]
      
      matrix_test[datarowcounter, 5] <- timeatsite
      
      
      if  (matrix_ind_counts[ChosenSpecies, 1] == 5 ) { 
        
        Searchable[ChosenSpecies] <- 0
        
      }
      
      PossDetectiontime <- rexp(N.species, 1/ Scenario[[15]][,TSFIndex])
      
      PossMinDetectiontime <- min(PossDetectiontime[Searchable == 1]) 
      
      rm(k,H, ChosenSpecies)
      
    }
    
    rm(timeatsite,setuptime,datarowcounter, PossDetectiontime, PossMinDetectiontime, Detectiontime)
    
    Findsthrutime[ ,  , Rep, TSFIndex] <- matrix_test
    
    SummaryData[ , 1, Rep, TSFIndex] <- matrix_ind_counts
    SummaryData[ , 2, Rep, TSFIndex] <- Presence     
    SummaryData[ , 3, Rep, TSFIndex] <- Scenario[[16]][,TSFIndex]
    SummaryData[ , 4, Rep, TSFIndex] <- Searchable
    SummaryData[ , 5, Rep, TSFIndex] <- Scenario[[1]][,6]
    SummaryData[ , 6, Rep, TSFIndex] <- TSF
    
    rm(TSF, TSFIndex, Rep, matrix_ind_counts, matrix_test, Presence, Searchable)
    
  }
  
  rm(TotalSites)
  
  return(list(FittingData, N.plants, Findsthrutime, SummaryData, species.ID))
  
  rm(FittingData, N.plants, Findsthrutime, SummaryData, species.ID)  
  
}




#########################################################################

# simulate 100 fieldtrips for dataset simulation - perhaps model less?

N.species <- dim(mu)[1]
N.scen <- 17
N.fieldtrips <- 100
All_FittingData <- array(NA, c(40000, 4, N.scen, N.fieldtrips))
#colnames(All_FittingData) <- c('Sp', 'Height', 'TSF', 'Reps')
Sample_sizes <- array(NA, c(N.scen, N.fieldtrips))
Sample_species <- array(NA, c(1,N.scen, N.fieldtrips))
All_Findsthrutime <- array(NA, c(40000, 6, N.reps, N.TSF, N.scen, N.fieldtrips))
#colnames(All_Findsthrutime) <- c('ChosenSpecies','N.plants',  'TSF',  'Rep', 'time_mins','H')
All_SummaryData <- array(NA, c(N.species, 6, N.reps, N.TSF, N.scen, N.fieldtrips))
#colnames(All_SummaryData) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')

All_species.ID <- array(NA, c(N.species, N.scen, N.fieldtrips))


##########################################################################

for (i in 1:N.scen) {
  
  for (f in 1:N.fieldtrips) {
    
    temp.list <- Fieldtrip(scenarios[[i]])
    All_FittingData[ , ,i,f] <-  temp.list[[1]]
    Sample_sizes[i,f] <- temp.list[[2]]
    All_Findsthrutime[, , , ,i,f] <- temp.list[[3]]
    All_SummaryData[, , , ,i,f] <- temp.list[[4]]
    All_species.ID[ ,i,f] <- temp.list[[5]]
    
    Sample_species[ , i,f] <- length(na.omit(unique(temp.list[[1]][,1])))
    
    temp.species <- as.numeric(length(unique(All_FittingData[ ,1 ,i,f])))
    
    #All_Results[1:(Sample_species[ ,i,f] * 3 + 7), ,i,f] <- Fitmodel(All_FittingData[1:Sample_sizes[i,f], ,i,f])
    
    rm(temp.list)
    
    
  }
  
  
}


############################################################################


#Sample_sizes
#All_SummaryData[1:20,,1,1:15,1,]
#All_FittingData[1:280,1:4,15,1]


#All_FittingData_1Indiv100f <- All_FittingData[1:280,1:4,15,]
#save(All_FittingData_1Indiv100f, file="~/Dropbox/PhD_1/Data/Optimisation/All_FittingData_1Indiv100f.rda")
#load("~/Dropbox/PhD_1/Data/Optimisation/All_FittingData_1Indiv100f.rda")


