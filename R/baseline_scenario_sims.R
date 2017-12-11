species <- read.csv("species.csv")

Time <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
tsf_values <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
priorities <- c(2, 8, 26, 41, 4, 15, 33, 86, 1, 6, 13, 3, 36, 28)

# Strategic sampling with a PriorityList of sites to visit

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
num_reps = 1 # one replicate
num_fieldtrips = 100 
sd_obs = 0.395
hours_per_day = 10 # long days
site_max_time = 60 * 30   # 30 hours at one replicate
setup_time = 1 # 1 minute to set up
measurement_time = 1 # 1 minute to measure 
home_site_travel = (1 * 2) # 1 minute to get to sites and home
bw_site_travel = 1 # 1 minute to travel between sites

#####################################################################################

mu <- matrix(NA, nrow = length(species$hmax), ncol=length(Time))
num_species <- dim(mu)[1]
row.names(mu) <- species$names
colnames(mu) <- tsf_values
occupancy_meta_matrix <- BuildMatrix(species$occupancy)
prob_matrix <- BuildMatrix(species$detection_prob)
num_tsf <- length(tsf_values)
field_mins_day <- hours_per_day * 60   

priority_list <- matrix(NA, 1400, 2)
priority_list[, 1] <- rep(priorities, 100)
priority_list[, 2] <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 100)
colnames(priority_list) <- c('TSF', 'Rep')

for (j in 1:length(species$hmax)) { 
  for (x in 1:length(Time)) {
    mu[j, x] <- species$hmax[j] / (1 + exp(-species$a[j] * (Time[x] - species$b[j])))  
  }  
}

BuildMatrix <- function(row_values) {
  matrx <- matrix(NA, nrow = num_species, ncol = length(tsf_values))
  row.names(matrx) <- species$name
  colnames(matrx) <- tsf_values
  for (w in 1:length(tsf_values)) {
    matrx[,w] <- row_values 
  }
  return(matrx)
}

BuildScenarios <- function(days, value_seq) {
  values <- matrix(0, nrow = 2, ncol = days * 2 + 1)
  values[2, ] <- value_seq 
  for (q in 2:length(values[1, ])) { 
    values[1, q] <- values[1, q - 1] + 0.5
  } 
  values[2,] <- values[2,]/sum(values[2,])
  for (q in 2:length(values[2, ])) {
    values[2, q] <- values[2, q - 1] + values[2, q]
  }
  list(
    mu = mu, 
    values = values, 
    field_mins_day = field_mins_day, 
    home_site_travel = home_site_travel, 
    site_max_time = site_max_time, 
    bw_site_travel = bw_site_travel, 
    priority_list = priority_list, 
    num_fieldtrips = num_fieldtrips, 
    total_field_days = days, 
    tsf_values = tsf_values, 
    num_reps = num_reps, 
    setup_time = setup_time, 
    measurement_time = measurement_time, 
    hours_per_day = hours_per_day, 
    prob_matrix = prob_matrix, 
    occupancy_meta_matrix = occupancy_meta_matrix, 
    sd_obs = sd_obs, 
    num_tsf = num_tsf
   )
}

scenarios <- Map(BuildScenarios, indexes, value_seqs)
scenario <- scenarios[[1]]
scenario$field_mins_day

##################################
# Collecting between 2 and 3 individuals per species

FieldTrip <- function(scenario) { 
  g <- as.numeric(Sys.time())
  set.seed((g - floor(g)) * 1e8 -> seed) 
  print(seed)

  species_id <- c(rep(0, 20))
  rando <- runif(1, 0, 1)
  realindex <- min(which(rando <= scenario$values[2, ])) 
  realfield_days <- scenario$values[1, realindex] 
  realfield_daysmins <- realfield_days * scenario$field_mins_day 
  total_site_time <- realfield_daysmins - scenario$home_site_travel * realfield_days  
  total_sites <- round(total_site_time / (scenario$site_max_time + scenario$bw_site_travel))
  num_plants <- 0

  fitting_data <- array(NA, c(40000, 4))
  colnames(fitting_data) <- c('Sp', 'Height', 'TSF', 'Rep')
  
  findsthrutime <- array(NA, c(40000, 6, num_reps, num_tsf))
  colnames(findsthrutime) <- c('SpeciesID', ' ', ' ', ' ', 'Time(mins)','Height(cm)')
  
  summary_data <- array(NA, c(num_species, 6, num_reps, num_tsf))
  colnames(summary_data) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')
  
  for (t in 1:total_sites) {
    tsf <- scenario$priority_list[t, 1]
    tsf_index <- which(scenario$tsf_values == tsf)
    reps <- scenario$priority_list[t, 2]
    matrix_ind_counts <- matrix(0, num_species, 1)
    matrix_test <- matrix(NA, 40000, 6)   
    colnames(matrix_test) <- c('ChosenSpecies', 'indivID',  'TSF',  'Rep', 'time_mins', 'H')
    
    presence <- rbinom(num_species, 1, scenario$occupancy_meta_matrix[, tsf_index])
    searchable <- presence
    site_time <- scenario$setup_time 
    row_counter <- 0
    
    poss_detection_time <- rexp(num_species, 1 / scenario$prob_matrix[, tsfindex])
    poss_min_detection_time <- min(poss_detection_time[searchable == 1]) 
    
    while((1 - prod(matrix_ind_counts[presence == 1] >= 2 )) & (site_time + poss_min_detection_time + scenario$measurement_time < scenario$site_max_time)) {  
      row_counter <- row_counter + 1
      num_plants <- num_plants + 1
      matrix_ind_counts[chosen_species, 1] <- matrix_ind_counts[chosen_species, 1] + 1
      
      detection_time <- poss_min_detection_time
      site_time <- site_time + detection_time + scenario$measurement_time
      
      chosen_species <- which(detection_time == poss_detection_time)
      species_id[chosen_species] <- 1
      
      k <- log(mu[chosenspecies, tsf_index]) - scenario$sd_obs^2
      H <- rlnorm(1, k, scenario$sd_obs)
      
      fitting_data[num_plants, 1] <- chosen_species
      fitting_data[num_plants, 2] <- H
      fitting_data[num_plants, 3] <- tsf_values[tsf_index]
      fitting_data[num_plants, 4] <- reps
      
      matrix_test[row_counter, 1] <- chosen_species
      matrix_test[row_counter, 2] <- num_plants
      matrix_test[row_counter, 3] <- tsf
      matrix_test[row_counter, 4] <- reps
      matrix_test[row_counter, 5] <- site_time
      matrix_test[row_counter, 6] <- H
      
      if  (matrix_ind_counts[chosen_species, 1] == 5 ) { 
        searchable[chosen_species] <- 0
      }
      
      poss_detection_time <- rexp(num_species, 1 / scenario$prob_matrix[, tsf_index])
      poss_min_detection_time <- min(poss_detection_time[searchable == 1]) 
    }
    
    findsthrutime[,  , reps, tsf_index] <- matrix_test
    
    summary_data[, 1, reps, tsf_index] <- matrix_ind_counts
    summary_data[, 2, reps, tsf_index] <- presence     
    summary_data[, 3, reps, tsf_index] <- scenario$occupancy_meta_matrix[, tsf_index]
    summary_data[, 4, reps, tsf_index] <- searchable
    summary_data[, 5, reps, tsf_index] <- scenario$mu[, 6]
    summary_data[, 6, reps, tsf_index] <- tsf
  }
  return(list(fitting_data, num_plants, findsthrutime, summary_data, species_id))
}

#########################################################################
# Simulate fieldtrips for dataset simulation - perhaps model less?

SimulateFieldTrips <- function(num_scen) {
  all_fitting_data <- array(NA, c(40000, 4, num_scen, num_fieldtrips))
  #colnames(all_fitting_data) <- c('Sp', 'Height', 'TSF', 'Reps')
  sample_sizes <- array(NA, c(num_scen, num_fieldtrips))
  sample_species <- array(NA, c(1,num_scen, num_fieldtrips))
  # all_findsthru_time <- array(NA, c(40000, 6, num_reps, num_TSF, num_scen, num_fieldtrips))
  #colnames(all_findsthru_time) <- c('ChosenSpecies','num_plants',  'TSF',  'Rep', 'time_mins','H')
  all_summary_data <- array(num_, c(num_species, 6, num_reps, num_TSF, num_scen, num_fieldtrips))
  colnames(all_summary_data) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')

  all_species_id <- array(NA, c(num_species, num_scen, num_fieldtrips))

  for (i in 1:num_scen) {
    for (f in 1:num_fieldtrips) {
      temp_list <- FieldTrip(scenarios[[i]])
      all_fitting_data[ , ,i,f] <-  temp_list$mu
      sample_sizes[i,f] <- temp_list$values
      all_findsthru_time[, , , ,i,f] <- temp_list$field_mins_day
      all_summary_data[, , , ,i,f] <- temp_list$home_site_travel
      all_species_id[ ,i,f] <- temp_list$site_max_time
      sample_species[ , i,f] <- length(na_omit(unique(temp_list$mu[,1])))
      temp_species <- as.numeric(length(unique(all_fitting_data[ ,1 ,i,f])))
      # all_results[1:(sample_species[ ,i,f] * 3 + 7), ,i,f] <- Fitmodel(all_fitting_data[1:sample_sizes[i,f], ,i,f])
    }
  }
}

############################################################################


#Sample_sizes
#All_SummaryData[1:20,,1,1:15,1,]
#All_FittingData[1:280,1:4,15,1]


#All_FittingData_1Indiv100f <- All_FittingData[1:280,1:4,15,]
#save(All_FittingData_1Indiv100f, file="~/Dropbox/PhD_1/Data/Optimisation/All_FittingData_1Indiv100f.rda")
#load("~/Dropbox/PhD_1/Data/Optimisation/All_FittingData_1Indiv100f.rda")


