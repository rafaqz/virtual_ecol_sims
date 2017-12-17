
build_priority_list <- function(sites, settings) {
  num_fieldtrips <- settings$num_fieldtrips
  priority_list <- matrix(NA, length(sites$id) * num_fieldtrips, 2)
  priority_list[, 1] <- rep(sites$id, num_fieldtrips)
  priority_list[, 2] <- rep(sites$replication, num_fieldtrips)
  colnames(priority_list) <- priority_list_colnames
  return(priority_list)
}

build_matrix <- function(row_values, site_ids, species) {
  matrx <- matrix(NA, nrow = length(species$hmax), ncol = length(site_ids))
  row.names(matrx) <- species$name
  colnames(matrx) <- site_ids
  for (w in 1:length(site_ids)) {
    matrx[, w] <- row_values
  }
  return(matrx)
}

build_input <- function(species, sites, settings) {
  site_ids <- sort(sites$id)
  occupancy_matrix <- build_matrix(species$occupancy, site_ids, species)
  probability_matrix <- build_matrix(species$detection_prob, site_ids, species)
  priority_list <- build_priority_list(sites, settings)
  list(
    occupancy_matrix = occupancy_matrix,
    probability_matrix = probability_matrix,
    priority_list = priority_list,
    site_ids = site_ids
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

generate_seed <- function() {
  g <- as.numeric(Sys.time())
  set.seed(1e8 * (g - floor(g)) -> seed)
  return(seed)
}

random_half_days <- function(scenario) {
  randomiser <- runif(1, 0, 1)
  half_days <- min(which(randomiser <= scenario$half_day_probs[2, ]))
  return(half_days)
}

random_presences <- function(input, num_species, site_index) {
  rbinom(num_species, 1, input$occupancy_matrix[, site_index])
}

random_detections <- function(input, num_species, site_index) {
  rexp(num_species, 1 / input$probability_matrix[, site_index])
}

total_sites <- function(half_days, scenario, settings) {
  field_days <- scenario$half_day_probs[1, half_days]
  field_mins <- field_days * settings$mins_per_day
  total_site_time <- field_mins - settings$home_site_travel * field_days
  total_sites <- round(total_site_time / (settings$max_site_time + settings$between_site_travel))
  return(total_sites)
}

fieldtrip <- function(model_data, scenario, settings, input, output_test = FALSE) {
  half_days <- random_half_days(scenario)
  total_sites <- total_sites(half_days, scenario, settings)

  num_plants <- 0
  num_species <- dim(model_data$mu)[1]
  num_sites <- length(input$site_ids)

  species_sampled <- c(rep(0, num_species))

  fitting_data <- array(NA, c(settings$size, 4))
  colnames(fitting_data) <- fitting_data_colnames
  summary_data <- array(NA, c(num_species, 6, settings$num_replicates, num_sites))
  colnames(summary_data) <- summary_colnames
  if (output_test) {
    findsthrutime <- array(NA, c(settings$size, 6, settings$num_replicates, num_sites))
    colnames(findsthrutime) <- findsthrutime_colnames
  }

  # TODO: Simplify these loops and make the logic clearer.
  for (s in 1:total_sites) {
    s = 1
    site_id <- input$priority_list[s, 1]
    site_index <- which(input$site_ids == site_id)
    replicates <- input$priority_list[s, 2]
    matrix_ind_counts <- matrix(0, num_species, 1)
    matrix_test <- matrix(NA, settings$size, 6)
    colnames(matrix_test) <- matrix_test_colnames

    poss_detection_time <- random_detections(input, num_species, site_index)
    presence <- searchable <- random_presences(input, num_species, site_index)
    site_time <- settings$setup_time
    row_counter <- 0
    pres_index <- presence == 1

    while (!prod(matrix_ind_counts[pres_index] >= 2)) {
      detection_time <- min(poss_detection_time[searchable == 1])
      site_time <- site_time + detection_time + settings$measurement_time
      if (site_time > settings$max_site_time) {
        break
      }

      row_counter <- row_counter + 1
      num_plants <- num_plants + 1

      chosen_species <- which(detection_time == poss_detection_time)
      species_sampled[chosen_species] <- 1
      matrix_ind_counts[chosen_species, 1] <- matrix_ind_counts[chosen_species, 1] + 1
      if (matrix_ind_counts[chosen_species, 1] == settings$min_collection) {
        searchable[chosen_species] <- 0
      }

      measurement <- random_measurement(input, settings, chosen_species, site_index)
      fitting_data[num_plants, ] <- c(chosen_species, measurement, input$site_ids[site_index], replicates)
      if (output_test) {
        matrix_test[row_counter, ] <- c(chosen_species, num_plants, site, replicates, site_time, measurement)
      }
    }

    if (output_test) {
      findsthrutime[, , replicates, site_index] <- matrix_test
    }

    summary_data[, 1, replicates, site_index] <- matrix_ind_counts
    summary_data[, 2, replicates, site_index] <- presence
    summary_data[, 3, replicates, site_index] <- input$occupancy_matrix[, site_index]
    summary_data[, 4, replicates, site_index] <- searchable
    summary_data[, 5, replicates, site_index] <- model_data$mu[, 6]
    summary_data[, 6, replicates, site_index] <- site_id
  }

  output <- list(
    fitting_data = fitting_data,
    num_plants = num_plants,
    summary_data = summary_data,
    species_sampled = species_sampled
  )

  if (output_test) {
    output_test$findsthrutime <- findsthrutime
  }

  return(output)
}

simulate_fieldtrips <- function(model_data, scenarios, settings, input, output_test = FALSE) {
  num_fieldtrips <- settings$num_fieldtrips
  num_replicates <- settings$num_replicates
  num_scenarios <- length(scenarios)
  num_species <- dim(input$occupancy_matrix)[1]
  num_sites <- length(input$site_ids)

  sample_sizes <- array(NA, c(num_scenarios, num_fieldtrips))
  sample_species <- array(NA, c(num_scenarios, num_fieldtrips))
  all_fitting_data <- array(NA, c(settings$size, 4, num_scenarios, num_fieldtrips))
  colnames(all_fitting_data) <- fitting_data_colnames
  all_summary_data <- array(NA, c(num_species, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
  colnames(all_summary_data) <- summary_colnames
  all_species_sampled <- array(NA, c(num_species, num_scenarios, num_fieldtrips))
  if (output_test) {
    all_findsthrutime <- array(NA, c(settings$size, 6, num_replicates, num_sites, num_scenarios, num_fieldtrips))
    colnames(all_findsthrutime) <- findsthrutime_colnames
  }

  for (i in 1:num_scenarios) {
    for (j in 1:num_fieldtrips) {
      output_test = FALSE
      trip <- fieldtrip(model_data, scenarios[[i]], settings, input, output_test)
      all_fitting_data[, , i, j] <- trip$fitting_data
      sample_sizes[i, j] <- trip$num_plants
      all_summary_data[, , , , i, j] <- trip$summary_data
      all_species_sampled[, i, j] <- trip$species_sampled
      sample_species[i, j] <- length(na.omit(unique(trip$fitting_data[, 1])))
      if (output_test) {
        all_findsthrutime[, , , , i, j] <- trip$findsthrutime
      }
      # temp_species <- as.numeric(length(unique(all_fitting_data[, 1, i, j])))
      # all_results[1:(sample_species[ ,i,j] * 3 + 7), ,i,j] <- Fitmodel(all_fitting_data[1:sample_sizes[i,j], ,i,j])
    }
  }

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

dependent_label <- "Height (cm)"
site_label <- "Time Since Fire"
replicates <- "Replicates"
chosen <- "Chosen Species"
time_ <- "Time (mins)"
dependent_summary_label <- "Mean Height(cm)"

summary_colnames <- c("Indivs Counts", "Presence", "Occupancy", "Searchable", dependent_summary_label, site_label)
findsthrutime_colnames <- c(chosen, "Num Plants", site_label, replicates, time_, dependent_label)
fitting_data_colnames <- c(chosen, dependent_label, site_label, replicates)
matrix_test_colnames <- c(chosen, "Indiv ID", site_label, replicates, time_, dependent_label)
priority_list_colnames <- c(site_label, replicates)

build_model_data <- function(species, sites, settings) {
  site_ids <- sort(sites$id)
  mu <- build_mu(species, site_ids)
  k <- build_k(mu, settings)
  list(
    mu = mu,
    k = k
  )
}

build_mu <- function(species, site_ids) {
  mu <- matrix(NA, nrow = length(species$hmax), ncol = length(site_ids))
  row.names(mu) <- species$names
  colnames(mu) <- site_ids
  for (j in 1:length(species$hmax)) {
    for (x in 1:length(site_ids)) {
      mu[j, x] <- species$hmax[j] / (1 + exp(-species$a[j] * (site_ids[x] - species$b[j])))
    }
  }
  return(mu)
}

build_k <- function(mu, settings) {
  log(mu) - settings$sd_obs^2
}

random_measurement <- function(input, settings, chosen_species, site_index) {
  rlnorm(1, model_data$k[chosen_species, site_index], settings$sd_obs)
}

# TODO: Separate and move this section into user-defined code (examples in readme/docs)
settings <- list(
  min_collection = 5,
  max_collection = 10,
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
  between_site_travel = 1 # 1 minute to travel between sites
)

# TODO: We need a formula to generate these tables.
indexes <- list(7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 111, 150)#, 222) # seq(7, settings$max_days, 7)
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
                   c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2))#,
                   # c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
                  )

seed <- generate_seed()
set.seed(seed)
scenarios <- Map(build_scenario, indexes, field_day_seqs)
scenario <- scenarios[[1]]
species <- read.csv("mallee.csv")
sites <- read.csv("sites.csv")
model_data <- build_model_data(species, sites, settings)
input <- build_input(species, sites, settings)
simulations <- simulate_fieldtrips(model_data, scenarios, settings, input)
