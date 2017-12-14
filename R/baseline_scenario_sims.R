summary_colnames <- c("Indivs Counts", "Presence", "Occup", "Searchable", "Mean Height(cm)", "Time Since Fire")
findsthrutime_colnames <- c("Chosen Species", "Num Plants", "Time Since Fire",  "Replicates", "Time (mins)", "Height (cm)")
fitting_data_colnames <- c("Species", "Height", "Time Since Fire", "Replicates")
matrix_test_colnames <- c("Chosen Species", "Indiv ID",  "Time Since Fire",  "Replicates", "Time (mins)", "Height (cm)")
priority_list_colnames <- c("Time Since Fire", "Replicates")

build_mu <- function(species, time_since_fire) {
  mu <- matrix(NA, nrow = length(species$hmax), ncol = length(time_since_fire))
  row.names(mu) <- species$names
  colnames(mu) <- time_since_fire
  for (j in 1:length(species$hmax)) {
    for (x in 1:length(time_since_fire)) {
      mu[j, x] <- species$hmax[j] / (1 + exp(-species$a[j] * (time_since_fire[x] - species$b[j])))
    }
  }
  return(mu)
}

build_priority_list <- function(sites, settings) {
  num_fieldtrips <- settings$num_fieldtrips
  priority_list <- matrix(NA, length(sites$prioritised_tsf) * num_fieldtrips, 2)
  priority_list[, 1] <- rep(sites$prioritised_tsf, num_fieldtrips)
  priority_list[, 2] <- rep(sites$replication, num_fieldtrips)
  colnames(priority_list) <- priority_list_colnames 
  return(priority_list)
}

build_matrix <- function(row_values, time_since_fire, species) {
  matrx <- matrix(NA, nrow = length(species$hmax), ncol = length(time_since_fire))
  row.names(matrx) <- species$name
  colnames(matrx) <- time_since_fire
  for (w in 1:length(time_since_fire)) {
    matrx[, w] <- row_values
  }
  return(matrx)
}

build_arrays <- function(species, sites, settings) {
  time_since_fire <- sort(sites$prioritised_tsf) 
  occupancy_meta_matrix <- build_matrix(species$occupancy, time_since_fire, species)
  prob_matrix <- build_matrix(species$detection_prob, time_since_fire, species)
  mu <- build_mu(species, time_since_fire)
  priority_list <- build_priority_list(sites, settings)

  list(
    mu = mu,
    occupancy_meta_matrix = occupancy_meta_matrix,
    prob_matrix = prob_matrix,
    priority_list = priority_list,
    time_since_fire = time_since_fire
  )
}

build_scenario <- function(days, half_day_seq) {
  timesteps <- days * 2 + 1 # half days + possible zero 
  half_day_probs <- matrix(0, nrow = 2, ncol = timesteps)
  half_day_probs[1, ] <- seq(0, days, 0.5)
  half_day_probs[2, ] <- half_day_seq / sum(half_day_seq)
  for (q in 2:length(half_day_probs[2,])) {
    half_day_probs[2,q] <- half_day_probs[2,q-1] + half_day_probs[2,q]
  }

  list(
    half_day_probs = half_day_probs,
    total_field_days = days
  )
}

random_index <- function(scenario) {
  g <- as.numeric(Sys.time()) 
  set.seed(1e8 * (g - floor(g)) -> seed)

  rando <- runif(1, 0, 1)
  which(rando <= scenario$half_day_probs[2, ])
  index <- min(which(rando <= scenario$half_day_probs[2, ]))
  return(index)
}

total_sites <- function(index, scenario, settings) {
  field_days <- scenario$half_day_probs[1, index]
  field_mins <- field_days * settings$mins_per_day
  total_site_time <- field_mins - settings$home_site_travel * field_days
  total_sites <- trunc(total_site_time / (settings$max_site_time + settings$between_site_travel))
  return(total_sites)
}

fieldtrip <- function(scenario, settings, arrays) {
  index <- random_index(scenario)
  total_sites <- total_sites(index, scenario, settings)

  num_plants <- 0
  num_species <- dim(arrays$mu)[1]
  num_sites <- length(arrays$time_since_fire)

  species_sampled <- c(rep(0, num_species))

  fitting_data <- array(NA, c(settings$size, 4))
  colnames(fitting_data) <- fitting_data_colnames
  findsthrutime <- array(NA, c(settings$size, 6, settings$num_replicates, num_sites))
  colnames(findsthrutime) <- findsthrutime_colnames
  summary_data <- array(NA, c(num_species, 6, settings$num_replicates, num_sites))
  colnames(summary_data) <- summary_colnames 

  # TODO: Simplify these loops and make the logic clearer.
  for (t in 1:total_sites) {
    tsf <- arrays$priority_list[t, 1]
    tsf_index <- which(arrays$time_since_fire == tsf)
    replicates <- arrays$priority_list[t, 2]
    matrix_ind_counts <- matrix(0, num_species, 1)
    matrix_test <- matrix(NA, settings$size, 6)
    colnames(matrix_test) <- matrix_test_colnames

    poss_detection_time <- rexp(num_species, 1 / arrays$prob_matrix[, tsf_index])
    presence <- rbinom(num_species, 1, arrays$occupancy_meta_matrix[, tsf_index])
    searchable <- presence
    site_time <- settings$setup_time
    row_counter <- 0

    while (1 - prod(matrix_ind_counts[presence == 1] >= 2)) {
      detection_time <- min(poss_detection_time[searchable == 1])
      potential_site_time <- site_time + detection_time + settings$measurement_time
      if (potential_site_time > settings$max_site_time) {
        break
      }

      site_time <- potential_site_time
      row_counter <- row_counter + 1

      num_plants <- num_plants + 1

      chosen_species <- which(detection_time == poss_detection_time)
      species_sampled[chosen_species] <- 1
      matrix_ind_counts[chosen_species, 1] <- matrix_ind_counts[chosen_species, 1] + 1
      if  (matrix_ind_counts[chosen_species, 1] == settings$min_collection ) {
        searchable[chosen_species] <- 0
      }

      k <- log(arrays$mu[chosen_species, tsf_index]) - settings$sd_obs ^ 2
      height <- rlnorm(1, k, settings$sd_obs)

      fitting_data[num_plants, ] <- c(chosen_species, height, arrays$time_since_fire[tsf_index], replicates)
      matrix_test[row_counter, ] <- c(chosen_species, num_plants, tsf, replicates, site_time, height)
    }

    findsthrutime[, , replicates, tsf_index] <- matrix_test

    summary_data[, 1, replicates, tsf_index] <- matrix_ind_counts
    summary_data[, 2, replicates, tsf_index] <- presence
    summary_data[, 3, replicates, tsf_index] <- arrays$occupancy_meta_matrix[, tsf_index]
    summary_data[, 4, replicates, tsf_index] <- searchable
    summary_data[, 5, replicates, tsf_index] <- arrays$mu[, 6]
    summary_data[, 6, replicates, tsf_index] <- tsf
  }

  list(
    fitting_data = fitting_data,
    num_plants = num_plants,
    findsthrutime = findsthrutime,
    summary_data = summary_data,
    species_sampled = species_sampled
  )
}

simulate_fieldtrips <- function(scenarios, settings, arrays) {
  num_fieldtrips <- settings$num_fieldtrips
  num_replicates <- settings$num_replicates
  num_scenarios <- length(scenarios)
  num_species <- dim(arrays$mu)[1]
  num_sites <- length(arrays$time_since_fire)

  sample_sizes <- array(NA, c(num_scenarios, num_fieldtrips))
  sample_species <- array(NA, c(num_scenarios, num_fieldtrips))
  all_fitting_data <- array(NA, c(settings$size, 4, num_scenarios, num_fieldtrips))
  colnames(all_fitting_data) <- fitting_data_colnames
  all_findsthrutime <- array(NA, c(settings$size, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
  colnames(all_findsthrutime) <- findsthrutime_colnames
  all_summary_data <- array(NA, c(num_species, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
  colnames(all_summary_data) <- summary_colnames
  all_species_sampled <- array(NA, c(num_species, num_scenarios, num_fieldtrips))

  for (i in 1:num_scenarios) {
    for (f in 1:num_fieldtrips) {
      temp_list <- fieldtrip(scenarios[[i]], settings, arrays)
      temp_list$fitting_data
      all_fitting_data[, , i, f] <- temp_list$fitting_data
      sample_sizes[i, f] <- temp_list$num_plants
      all_findsthrutime[, , , , i, f] <- temp_list$findsthrutime
      all_summary_data[, , , , i, f] <- temp_list$summary_data
      all_species_sampled[, i, f] <- temp_list$species_sampled
      sample_species[i, f] <- length(na.omit(unique(temp_list$fitting_data[, 1])))
      temp_species <- as.numeric(length(unique(all_fitting_data[, 1, i, f])))
      # all_results[1:(sample_species[ ,i,f] * 3 + 7), ,i,f] <- Fitmodel(all_fitting_data[1:sample_sizes[i,f], ,i,f])
    }
  }

  return(list(
    all_fitting_data = all_fitting_data,
    sample_sizes = sample_sizes,
    all_findsthrutime = all_findsthrutime,
    all_summary_data = all_summary_data,
    all_species_sampled = all_species_sampled,
    sample_species = sample_species
  ))
}

# TODO: Separate and move this section into user-defined code (examples in readme/docs)
settings <- list(
  min_collection = 5,
  max_collection = 10,
  size = 2000,
  num_replicates = 1, # one replicate
  num_fieldtrips = 10, #
  sd_obs = 0.395,
  max_days = 222,
  mins_per_day = 10 * 60,
  max_site_time = 60 * 30,   # 30 hours at one replicate
  setup_time = 1, # 1 minute to set up
  measurement_time = 1, # 1 minute to measure
  home_site_travel = (1 * 2), # 1 minute to get to sites and home
  between_site_travel = 1 # 1 minute to travel between sites
)

indexes <- list(7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84) #, 91, 98, 105, 111, 150)#, 222) # seq(7, settings$max_days, 7) 
field_day_seqs <- list(c(rep(0,8), rep(1,3), rep(2,1),rep(3,1), rep(5,1), rep(10,1)),
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
                   c(rep(0,102), rep(1,34), rep(2,17), rep(3,6), rep(4,4), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1))#,
                   # c(rep(0,110), rep(1,37), rep(2,18), rep(3,6), rep(4,4), rep(5,3), rep(6,2), rep(7,1), rep(8,1), rep(10,1)),
                   # c(rep(0,118), rep(1,40), rep(2,20), rep(3,7), rep(4,5), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   # c(rep(0,127), rep(1,42), rep(2,21), rep(3,8), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   # c(rep(0,134), rep(1,45), rep(2,22), rep(3,9), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   # c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2)),
                   # c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
                  )

scenarios <- Map(build_scenario, indexes, field_day_seqs)
scenario <- scenarios[[1]]
species <- read.csv("mallee.csv")
sites <- read.csv("sites.csv")
arrays <- build_arrays(species, sites, settings)
trip <- fieldtrip(scenarios[[1]], settings, arrays)
simulations <- simulate_fieldtrips(scenarios, settings, arrays)
