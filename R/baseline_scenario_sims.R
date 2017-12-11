build_mu <- function(species, tsf_values, Time) {
  mu <- matrix(NA, nrow = length(species$hmax), ncol = length(Time))
  row.names(mu) <- species$names
  colnames(mu) <- tsf_values
  for (j in 1:length(species$hmax)) {
    for (x in 1:length(Time)) {
      mu[j, x] <- species$hmax[j] / (1 + exp(-species$a[j] * (Time[x] - species$b[j])))
    }
  }
  return(mu)
}

build_priority_list <- function(priorities) {
  priority_list <- matrix(NA, 1400, 2)
  priority_list[, 1] <- rep(priorities, 100)
  priority_list[, 2] <- rep(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 100)
  colnames(priority_list) <- c("TSF", "Rep")
  return(priority_list)
}

build_matrix <- function(row_values, tsf_values, species) {
  matrx <- matrix(NA, nrow = length(species$hmax), ncol = length(tsf_values))
  row.names(matrx) <- species$name
  colnames(matrx) <- tsf_values
  for (w in 1:length(tsf_values)) {
    matrx[, w] <- row_values
  }
  return(matrx)
}

build_arrays <- function(species, priorities, tsf_values, Time) {
  occupancy_meta_matrix <- build_matrix(species$occupancy, tsf_values, species)
  prob_matrix <- build_matrix(species$detection_prob, tsf_values, species)
  mu <- build_mu(species, tsf_values, Time)
  priority_list <- build_priority_list(priorities)

  list(
    mu = mu,
    occupancy_meta_matrix = occupancy_meta_matrix,
    prob_matrix = prob_matrix,
    priority_list = priority_list,
    tsf_values = tsf_values
  )
}

build_scenario <- function(days, field_day_seq) {
  field_day_values <- matrix(0, nrow = 2, ncol = days * 2 + 1)
  field_day_values[2, ] <- field_day_seq
  for (q in 2:length(field_day_values[1, ])) {
    field_day_values[1, q] <- field_day_values[1, q - 1] + 0.5
  }
  field_day_values[2, ] <- field_day_values[2, ] / sum(field_day_values[2, ])
  for (q in 2:length(field_day_values[2, ])) {
    field_day_values[2, q] <- field_day_values[2, q - 1] + field_day_values[2, q]
  }

  list(
    field_day_values = field_day_values,
    total_field_days = days
  )
}

random_index <- function(scenario) {
  g <- as.numeric(Sys.time())
  set.seed(1e8 * (g - floor(g)) -> seed)
  print(seed)

  rando <- runif(1, 0, 1)
  index <- min(which(rando <= scenario$field_day_values[2, ]))
  return(index)
}

total_sites <- function(index, scenario, settings) {
  field_days <- scenario$field_day_values[1, index]
  field_mins <- field_days * settings$mins_per_day
  total_site_time <- field_mins - settings$home_site_travel * field_days
  # TODO: Should this be trunc instead of round?
  total_sites <- round(total_site_time / (settings$max_site_time + settings$between_site_travel))
  return(total_sites)
}

fieldtrip <- function(scenario, settings, arrays) {
  index <- random_index(scenario)
  total_sites <- total_sites(index, scenario, settings)
  species_id <- c(rep(0, 20))

  num_plants <- 0
  num_species <- dim(arrays$mu)[1]
  num_tsf <- length(arrays$tsf_values)

  fitting_data <- array(NA, c(40000, 4))
  colnames(fitting_data) <- c("Sp", "Height", "TSF", "Rep")
  findsthrutime <- array(NA, c(40000, 6, settings$num_reps, num_tsf))
  colnames(findsthrutime) <- c("SpeciesID", " ", " ", " ", "Time(mins)","Height(cm)")
  summary_data <- array(NA, c(num_species, 6, settings$num_reps, num_tsf))
  colnames(summary_data) <- c("indivs_counts", "Presence", "Occup", "Searchable", "mean Height (cm)", "TSF")

  # TODO: Simplify these loops and make the logic clearer.
  for (t in 1:total_sites) {
    tsf <- arrays$priority_list[t, 1]
    tsf_index <- which(arrays$tsf_values == tsf)
    reps <- arrays$priority_list[t, 2]
    matrix_ind_counts <- matrix(0, num_species, 1)
    matrix_test <- matrix(NA, 40000, 6)
    colnames(matrix_test) <- c("ChosenSpecies", "indivID",  "TSF",  "Rep", "time_mins", "H")

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
      species_id[chosen_species] <- 1
      matrix_ind_counts[chosen_species, 1] <- matrix_ind_counts[chosen_species, 1] + 1
      if  (matrix_ind_counts[chosen_species, 1] == 5 ) {
        searchable[chosen_species] <- 0
      }

      k <- log(arrays$mu[chosen_species, tsf_index]) - settings$sd_obs ^ 2
      H <- rlnorm(1, k, settings$sd_obs)

      fitting_data[num_plants, 1] <- chosen_species
      fitting_data[num_plants, 2] <- H
      fitting_data[num_plants, 3] <- arrays$tsf_values[tsf_index]
      fitting_data[num_plants, 4] <- reps

      matrix_test[row_counter, 1] <- chosen_species
      matrix_test[row_counter, 2] <- num_plants
      matrix_test[row_counter, 3] <- tsf
      matrix_test[row_counter, 4] <- reps
      matrix_test[row_counter, 5] <- site_time
      matrix_test[row_counter, 6] <- H
    }

    findsthrutime[, , reps, tsf_index] <- matrix_test

    summary_data[, 1, reps, tsf_index] <- matrix_ind_counts
    summary_data[, 2, reps, tsf_index] <- presence
    summary_data[, 3, reps, tsf_index] <- arrays$occupancy_meta_matrix[, tsf_index]
    summary_data[, 4, reps, tsf_index] <- searchable
    summary_data[, 5, reps, tsf_index] <- arrays$mu[, 6]
    summary_data[, 6, reps, tsf_index] <- tsf
  }

  list(
    fitting_data = fitting_data,
    num_plants = num_plants,
    findsthrutime = findsthrutime,
    summary_data = summary_data,
    species_id = species_id
  )
}

############################################
# Simulate fieldtrips for dataset simulation

simulate_fieldtrips <- function(scenarios, settings, arrays) {
  num_scenarios <- settings$num_scenarios
  num_fieldtrips <- settings$num_fieldtrips
  num_reps <- settings$num_reps
  num_species <- dim(arrays$mu)[1]
  num_tsf <- length(arrays$tsf_values)

  sample_sizes <- array(NA, c(num_scenarios, num_fieldtrips))
  sample_species <- array(NA, c(1,num_scen, num_fieldtrips))
  all_fitting_data <- array(NA, c(40000, 4, num_scenarios, num_fieldtrips))
  colnames(all_fitting_data) <- c("Sp", "Height", "TSF", "Reps")
  all_findsthrutime <- array(NA, c(40000, 6, num_reps, num_tsf, num_scenarios, num_fieldtrips))
  colnames(all_findsthrutime) <- c("ChosenSpecies", "num_plants",  "TSF",  "Rep", "time_mins","H")
  all_summary_data <- array(NA, c(num_species, 6, num_reps,  num_tsf, num_scenarios, num_fieldtrips))
  colnames(all_summary_data) <- c("indivs_counts", "Presence", "Occup", "Searchable", "mean H(cm)", "TSF")

  all_species_id <- array(NA, c(num_species, num_scenarios, num_fieldtrips))

  for (i in 1:num_scen) {
    for (f in 1:num_fieldtrips) {
      temp_list <- FieldTrip(scenarios[[i]], settings, arrays)
      all_fitting_data[, , i, f] <- temp_list$fitting_data
      sample_sizes[i,f] <- temp_list$num_plants
      all_findsthrutime[, , , , i, f] <- temp_list$findsthrutime
      all_summary_data[, , , , i, f] <- temp_list$summary_data
      all_species_id[, i, f] <- temp_list$species_id
      sample_species[, i, f] <- length(na_omit(unique(temp_list$fitting_data[, 1])))
      temp_species <- as.numeric(length(unique(all_fitting_data[, 1, i, f])))
      # all_results[1:(sample_species[ ,i,f] * 3 + 7), ,i,f] <- Fitmodel(all_fitting_data[1:sample_sizes[i,f], ,i,f])
    }
  }
}

# TODO: Import these as data / calculate mathematically
indexes <- list(7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 111, 150, 222)
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
                   c(rep(0,102), rep(1,34), rep(2,17), rep(3,6), rep(4,4), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,110), rep(1,37), rep(2,18), rep(3,6), rep(4,4), rep(5,3), rep(6,2), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,118), rep(1,40), rep(2,20), rep(3,7), rep(4,5), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,127), rep(1,42), rep(2,21), rep(3,8), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,134), rep(1,45), rep(2,22), rep(3,9), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1)),
                   c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2)),
                   c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
                  )
Time <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
tsf_values <- c(0, 1, 2, 3, 4, 6, 8, 13, 15, 26, 28, 33, 36, 41, 86)
priorities <- c(2, 8, 26, 41, 4, 15, 33, 86, 1, 6, 13, 3, 36, 28)


# TODO: Separate and move this section into user-defined code (examples in readme/docs)
settings <- list(
  num_reps = 1, # one replicate
  num_fieldtrips = 100,
  sd_obs = 0.395,
  mins_per_day = 10 * 60,
  max_site_time = 60 * 30,   # 30 hours at one replicate
  setup_time = 1, # 1 minute to set up
  measurement_time = 1, # 1 minute to measure
  home_site_travel = (1 * 2), # 1 minute to get to sites and home
  between_site_travel = 1 # 1 minute to travel between sites
)

scenarios <- Map(build_scenario, indexes, field_day_seqs)
species <- read.csv("species.csv")
arrays <- build_arrays(species, priorities, tsf_values, Time)
trip <- fieldtrip(scenarios[[1]], settings, arrays)
