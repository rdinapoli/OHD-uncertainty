# 31_session_info.R - Capture R session information
# Purpose: Records R version, package versions, and system info for reproducibility
# Inputs: None
# Outputs: output/session_info.txt
# Runtime: <1 minute

library(tidyverse)
library(rstan)
library(bayesplot)
library(loo)
library(here)
library(Bchron)
library(rcarbon)

options(mc.cores = parallel::detectCores())

sink(here("output/session_info.txt"))
cat("Session information captured:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
cat("\n\nStan version:\n")
cat(rstan::stan_version(), "\n")
cat("\nCPU cores available:", parallel::detectCores(), "\n")
sink()

cat("Session info saved to output/session_info.txt\n")
