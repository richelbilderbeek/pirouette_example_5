# Works under Linux and MacOS only

suppressMessages(library(ggplot2))
library(pirouette)
library(babette)
library(beautier)

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()
example_no <- 5
rng_seed <- 314
folder_name <- paste0("example_", example_no, "_", rng_seed)

################################################################################
# Create phylogeny
################################################################################
phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)

################################################################################
# Setup pirouette
################################################################################
pir_params <- create_std_pir_params(folder_name = folder_name)
pir_params$twinning_params <- NA

if (is_testing) {
  pir_params <- shorten_pir_params(pir_params)
}

################################################################################
# Run pirouette
################################################################################
pir_out <- pir_run(
  phylogeny,
  pir_params = pir_params
)

pir_save(
  phylogeny = phylogeny,
  pir_params = pir_params,
  pir_out = pir_out,
  folder_name = folder_name
)

