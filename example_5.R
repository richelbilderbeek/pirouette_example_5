# Code of example 5
#
# Works under Linux and MacOS only
#
#
#

# Set the RNG seed
rng_seed <- 314

suppressMessages(library(ggplot2))
library(pirouette)
library(babette)
library(beautier)

root_folder <- getwd()
example_no <- 5
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())
phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)

alignment_params <- create_alignment_params(
  sim_true_alignment_fun = get_sim_true_alignment_with_std_site_model_fun(
    mutation_rate = 0.1
  ),
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed
)

# JC69, strict, Yule
generative_experiment <- create_gen_experiment()
check_experiment(generative_experiment)

# All non-Yule tree priors
candidate_experiments <- create_all_experiments(
  exclude_model = generative_experiment$inference_model
)
check_experiments(candidate_experiments)

experiments <- c(list(generative_experiment), candidate_experiments)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

check_experiments(experiments)

# Testing
if (1 == 2) {
  experiments <- experiments[1:2]
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc <- create_mcmc(chain_length = 10000, store_every = 1000)
    experiments[[i]]$est_evidence_mcmc <- create_mcmc_nested_sampling(
      chain_length = 10000,
      store_every = 1000,
      epsilon = 100.0
    )
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments
)

rm_pir_param_files(pir_params)

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))


pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
