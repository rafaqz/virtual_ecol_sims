

# Initialisation functions ----------------------------------------------------

setup_simulation <- function(settings, sites, species, model_data) {
  # TODO: Use a formula to generate these lists and generalise this idea. 

  # This list could be replaced with a formula for a distribution of the 
  # number of days, optionally specified explicitly in a list?
  num_days <- list(7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 111, 150)#, 222) # seq(7, settings$max_days, 7)

  # And a formula to replace this list using num_days and parameters to change 
  # the probability distribution?
  field_day_scenarios <- list(c(rep(0,8), rep(1,3), rep(2,1),rep(3,1), rep(5,1), rep(10,1)),
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
                     c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2))#,
                     # c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
                    )

  scenarios = Map(build_scenario, num_days, field_day_scenarios)

  list(
    settings = settings,
    scenarios = scenarios,
    sites = sites,
    model_data = model_data,
    column_names = build_colnames(settings)
  )
}

# Build cumulative probability distribution.
build_scenario <- function(num_days, field_day_values) {
  timesteps <- num_days * 2 + 1 # half days + possible zero
  half_day_probs <- matrix(0, nrow = 2, ncol = timesteps)
  half_day_probs[1, ] <- seq(0, num_days, 0.5)
  half_day_probs[2, ] <- field_day_values / sum(field_day_values)
  for (q in 2:length(half_day_probs[2, ])) {
    half_day_probs[2,q] <- half_day_probs[2, q-1] + half_day_probs[2, q]
  }

  list(
    half_day_probs = half_day_probs,
    total_field_days = num_days
  )
}

# Allows customisation of output column names
build_colnames <- function(settings) {
  dependent_label <- settings$dependent_label 
  dependent_summary_label <- settings$dependent_summary_label
  site_name_label <- settings$site_name_label
  replicates <- "Replicates"
  chosen <- "Chosen Species"
  time_ <- "Time (mins)"

  return(list(
    summary = c("Indivs Counts", "Presence", "Occupancy", "Searchable", dependent_summary_label, site_name_label),
    findsthrutime = c(chosen, "Num Plants", site_name_label, replicates, time_, dependent_label),
    fitting_data = c(chosen, dependent_label, site_name_label, replicates),
    matrix_test = c(chosen, "Indiv ID", site_name_label, replicates, time_, dependent_label),
    priority_list = c(site_name_label, replicates)
  ))
}

generate_seed <- function() {
  g <- as.numeric(Sys.time())
  set.seed(1e8 * (g - floor(g)) -> seed)
  return(seed)
}

# Randomisation functions ---------------------------------------------

random_half_days <- function(scenario) {
  randomiser <- runif(1, 0, 1)
  half_days <- min(which(randomiser <= scenario$half_day_probs[2, ]))
  return(half_days)
}

random_presences <- function(sites, num_species, site_index) {
  rbinom(num_species, 1, sites$occupancy[, site_index])
}

random_detections <- function(sites, num_species, site_index) {
  rexp(num_species, 1 / sites$detection_probabilities[, site_index])
}

# Other math ----------------------------------------------------------

# Calculate the number of sites possible to visit in a scenario.
calculate_total_sites <- function(half_days, scenario, settings) {
  field_days <- scenario$half_day_probs[1, half_days]
  field_mins <- field_days * settings$mins_per_day
  total_site_time <- field_mins - settings$home_site_travel * field_days
  total_sites <- round(total_site_time / (settings$max_site_time + settings$between_site_travel))
  return(total_sites)
}

# Calculate the maximum possible finds for the feildtrip
calculate_max_finds <- function(settings, num_species, num_sites) {
  settings$max_collection * num_species * num_sites
}

# Core functions ---------------------------------------------

# Simulate a single fieldtrip in a given scenario.
fieldtrip <- function(model_data, scenario, settings, sites, column_names, 
                      output_test = FALSE) {
  half_days <- random_half_days(scenario)
  total_sites <- calculate_total_sites(half_days, scenario, settings)

  num_individuals <- 0
  num_species <- dim(model_data$mu)[1]
  num_sites <- length(sites$site_names)
  max_finds <- calculate_max_finds(settings, num_species, num_sites)

  species_sampled <- c(rep(0, num_species))

  fitting_data <- array(NA, c(settings$size, 4))
  colnames(fitting_data) <- column_names$fitting_data
  summary_data <- array(NA, c(num_species, 6, settings$num_replicates, num_sites))
  colnames(summary_data) <- column_names$summary
  if (output_test) {
    findsthrutime <- array(NA, c(max_finds, 6, settings$num_replicates, num_sites))
    colnames(findsthrutime) <- column_names$findsthrutime
  }

  for (s in 1:total_sites) {
    s_cyclic = (s - 1) %% num_sites + 1
    site_name <- sites$priorities[s_cyclic]
    site_index <- which(sites$site_names == site_name)
    replicates <- sites$replicates[s_cyclic]
    matrix_ind_counts <- matrix(0, num_species, 1)
    matrix_test <- matrix(NA, max_finds, 6)
    colnames(matrix_test) <- column_names$matrix_test

    poss_detection_time <- random_detections(sites, num_species, site_index)
    presence <- searchable <- random_presences(sites, num_species, site_index)
    site_time <- settings$setup_time
    row_counter <- 0
    pres_index <- presence == 1

    while (!prod(matrix_ind_counts[pres_index] >= settings$min_collection)) {
      detection_time <- min(poss_detection_time[searchable == 1])
      site_time <- site_time + detection_time + settings$measurement_time
      if (site_time > settings$max_site_time) {
        break
      }

      row_counter <- row_counter + 1
      num_individuals <- num_individuals + 1

      chosen_species <- which(detection_time == poss_detection_time)
      species_sampled[chosen_species] <- 1
      matrix_ind_counts[chosen_species, 1] <- matrix_ind_counts[chosen_species, 1] + 1
      if (matrix_ind_counts[chosen_species, 1] == settings$max_collection) {
        searchable[chosen_species] <- 0
      }

      measurement <- settings$measurement_func(settings, sites, model_data, chosen_species, site_index)
      fitting_data[num_individuals, ] <- c(chosen_species, measurement, sites$site_names[site_index], replicates)
      if (output_test) {
        matrix_test[row_counter, ] <- c(chosen_species, num_individuals, site, replicates, site_time, measurement)
      }
    }

    if (output_test) {
      findsthrutime[, , replicates, site_index] <- matrix_test
    }
    
    summary_data[, 1, replicates, site_index] <- matrix_ind_counts
    summary_data[, 2, replicates, site_index] <- presence
    summary_data[, 3, replicates, site_index] <- sites$occupancy[, site_index]
    summary_data[, 4, replicates, site_index] <- searchable
    summary_data[, 5, replicates, site_index] <- model_data$mu[, 6]
    summary_data[, 6, replicates, site_index] <- site_name
  }

  output <- list(
    fitting_data = fitting_data,
    num_individuals = num_individuals,
    summary_data = summary_data,
    species_sampled = species_sampled
  )
  if (output_test) {
    output_test$findsthrutime <- findsthrutime
  }

  return(output)
}


# Simulate multiple fieldtrips in multiple scenarios.
simulate_fieldtrips <- function(simulation, output_test = FALSE) {
  # Define local variables
  model_data <- simulation$model_data
  scenarios <- simulation$scenarios
  settings <- simulation$settings
  sites <- simulation$sites
  column_names <- simulation$column_names
  num_fieldtrips <- settings$num_fieldtrips
  num_replicates <- settings$num_replicates
  num_scenarios <- length(scenarios)
  num_species <- dim(sites$occupancy)[1]
  num_sites <- length(sites$site_names)
  max_finds <- calculate_max_finds(settings, num_species, num_sites)

  # Allocate arrays
  sample_sizes <- array(NA, c(num_scenarios, num_fieldtrips))
  sample_species <- array(NA, c(num_scenarios, num_fieldtrips))
  all_fitting_data <- array(NA, c(settings$size, 4, num_scenarios, num_fieldtrips))
  colnames(all_fitting_data) <- column_names$fitting_data
  all_summary_data <- array(NA, c(num_species, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
  colnames(all_summary_data) <- column_names$summary
  all_species_sampled <- array(NA, c(num_species, num_scenarios, num_fieldtrips))
  if (output_test) {
    all_findsthrutime <- array(NA, c(max_finds, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
    colnames(all_findsthrutime) <- column_names$findsthrutime
  }

  # Run simulations
  for (i in 1:num_scenarios) {
    for (j in 1:num_fieldtrips) {
      output_test = FALSE
      trip <- fieldtrip(model_data, scenarios[[i]], settings, sites, column_names, output_test)
      all_fitting_data[, , i, j] <- trip$fitting_data
      sample_sizes[i, j] <- trip$num_individuals
      all_summary_data[, , , , i, j] <- trip$summary_data
      all_species_sampled[, i, j] <- trip$species_sampled
      sample_species[i, j] <- length(na.omit(unique(trip$fitting_data[, 1])))
      if (output_test) {
        all_findsthrutime[, , , , i, j] <- trip$findsthrutime
      }
    }
  }

  # Compile output data
  output <- list(
    all_fitting_data = all_fitting_data,
    sample_sizes = sample_sizes,
    all_summary_data = all_summary_data,
    all_species_sampled = all_species_sampled,
    sample_species = sample_species
  )

  if (output_test) {
    output$all_findsthrutime <- all_findsthrutime
  }

  return(output)
}



# Project specific functions and lists ----------------------------------------
# TODO: Separate and move this section into user-defined code (examples in readme/docs)

# Conveniance function to fill a matrix with a vector when site specific
# priors are not available.
fill_matrix <- function(row_values, site_names, species) {
  matrx <- matrix(NA, nrow = length(species$hmax), ncol = length(site_names))
  row.names(matrx) <- species$name
  colnames(matrx) <- site_names
  for (w in 1:length(site_names)) {
    matrx[, w] <- row_values
  }
  return(matrx)
}

# Prepare data to be used in the measurment function during simulations.
prepare_model_data <- function(settings, species_data, site_data, site_names) {
  mu <- matrix(NA, nrow = length(species_data$hmax), ncol = length(site_names))
  row.names(mu) <- species_data$names
  colnames(mu) <- site_names
  for (j in 1:length(species_data$hmax)) {
    for (x in 1:length(site_names)) {
      mu[j, x] <- species_data$hmax[j] / (1 + exp(-1 * species_data$a[j] * (site_names[x] - species_data$b[j])))
    }
  }
  k <- log(mu) - settings$sd_obs^2

  list(mu = mu, k = k)
}

# A measurement simulation function for tree height.
random_measurement <- function(settings, sites, model_data, chosen_species, site_index) {
  rlnorm(1, model_data$k[chosen_species, site_index], settings$sd_obs)
}

# Settings for this study.
settings <- list(
  min_collection = 2,
  max_collection = 5,
  size = 5000,
  num_replicates = 1, # one replicate
  num_fieldtrips = 10,
  sd_obs = 0.395,
  max_days = 222,
  mins_per_day = 10 * 60,
  max_site_time = 60 * 30,   # 30 hours at one replicate
  setup_time = 1, # 1 minute to set up
  measurement_time = 1, # 1 minute to measure
  home_site_travel = 1 * 2, # 1 minute to get to sites and home
  between_site_travel = 1, # 1 minute to travel between sites
  dependent_label = "Height (cm)",
  dependent_summary_label = "Mean Height(cm)",
  site_name_label = "Time Since Fire",
  measurement_func = random_measurement
)

# Run Simulations -------------------------------------------------------------

# Import data from csv and prepare for simulation.
species_data <- read.csv("mallee.csv")
site_data <- read.csv("sites.csv")

site_names <- sort(site_data$name)
occupancy <- fill_matrix(species_data$occupancy, site_names, species_data) 
detection_probabilities <- fill_matrix(species_data$detection_prob, site_names, species_data)
sites = list(
  site_names = site_names,
  occupancy = occupancy,
  detection_probabilities = detection_probabilities,
  priorities = site_data$name,
  replicates = site_data$replicates
)
model_data <- prepare_model_data(settings, species_data, site_data, site_names)


# Run simulation.
simulation <- setup_simulation(settings, sites, species, model_data)
results <- simulate_fieldtrips(simulation)
