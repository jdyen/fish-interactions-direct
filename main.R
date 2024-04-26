# VEFMAP fish population modelling:
#   Analysis of fish data to determine strength of direct interactions
#   with invasive species and functional groups of species (incl. native 
#   species)
#
# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
#
# Last updated: 26 April 2024

# flags for slow steps
reload_data <- FALSE
sample_again <- TRUE

# need some packages
library(qs)
library(aae.db)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(brms)
library(ggplot2)
library(bayesplot)

# load helper scripts
source("R/data.R")
source("R/lookup.R")
source("R/utils.R")

# load data sets
cpue <- fetch_fish(recompile = reload_data)

# filter out comparison surveys
cpue <- cpue |> filter(category != "comparison")

# and some groups where some size categories do not occur
cpue_sum <- cpue |>
  group_by(species, category) |>
  summarise(catch = sum(catch)) |>
  ungroup() |>
  mutate(include = catch >= 10) |>
  select(-catch)
cpue <- cpue |>
  left_join(cpue_sum, by = c("species", "category")) |>
  filter(include) |>
  select(-include)

# collapse to a single survey per site and year
cpue <- cpue |>
  select(waterbody, id_site, survey_year, species, category, gear_type, reach_no, catch, effort_h) |>
  group_by(waterbody, id_site, survey_year, species, category, gear_type, reach_no) |>
  summarise(
    catch = sum(catch), 
    effort_h = sum(effort_h)
  ) |>
  ungroup()

# define predictors for each species, based on exotic species catch
#   and summed catch of all other species (large, small, or combined)
mc_predictors <- define_predictors(cpue, target = "Maccullochella peelii")
bf_predictors <- define_predictors(cpue, target = "Gadopsis marmoratus")
rb_predictors <- define_predictors(cpue, target = "Melanotaenia fluviatilis")

# set up response data sets for each species, adding the relevant predictors 
#   to each and filtering to relevant target systems
mc_cpue <- cpue |>
  filter(
    species == "Maccullochella peelii",
    waterbody %in% c(
      "Broken River", "Broken Creek", "Campaspe River",
      "Goulburn River", "Loddon River", "Ovens River"
    ),
    !(waterbody == "Loddon River" & reach_no %in% c(2, 3)),
  ) |>
  mutate(survey_year_minus1 = survey_year - 1) |>
  left_join(
    cpue |>
      group_by(waterbody, id_site, survey_year, species, gear_type, reach_no) |>
      summarise(log_cpue_ym1 = log(sum(catch) + 1) - log(median(effort_h))) |>
      ungroup(),
    by = c(
      "waterbody", "reach_no", "id_site", "species", "gear_type",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  left_join(
    mc_predictors |>
      mutate(
        across(contains("cpue_"), ~ log(. + 1), .names = "log_{.col}")
      ) |>
      select(waterbody, id_site, survey_year, gear_type, reach_no, contains("log_cpue_")),
    by = c(
      "waterbody", "id_site", "gear_type", "reach_no",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  filter(!is.na(log_cpue_ym1)) |>
  mutate(species_category = paste(species, category, sep = "_"))
bf_cpue <- cpue |>
  filter(
    species == "Gadopsis marmoratus",
    !waterbody %in% c(
      "Broken River", "Broken Creek", "Campaspe River",
      "Goulburn River", "Loddon River", "Ovens River"
    ),
    !(waterbody == "Loddon River" & reach_no > 2),
  ) |>
  mutate(survey_year_minus1 = survey_year - 1) |>
  left_join(
    cpue |>
      group_by(waterbody, id_site, survey_year, species, gear_type, reach_no) |>
      summarise(log_cpue_ym1 = log(sum(catch) + 1) - log(median(effort_h))) |>
      ungroup(),
    by = c(
      "waterbody", "reach_no", "id_site", "species", "gear_type",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  left_join(
    bf_predictors |>
      mutate(
        across(contains("cpue_"), ~ log(. + 1), .names = "log_{.col}")
      ) |>
      select(waterbody, id_site, survey_year, gear_type, reach_no, contains("log_cpue_")),
    by = c(
      "waterbody", "id_site", "gear_type", "reach_no",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  filter(!is.na(log_cpue_ym1)) |>
  mutate(species_category = paste(species, category, sep = "_"))
rb_cpue <- cpue |>
  filter(
    species == "Melanotaenia fluviatilis",
    waterbody %in% c(
      "Broken River", "Broken Creek", "Campaspe River",
      "Goulburn River", "Loddon River", "Ovens River"
    ),
    !(waterbody == "Loddon River" & reach_no %in% c(2, 3)),
  ) |>
  mutate(survey_year_minus1 = survey_year - 1) |>
  left_join(
    cpue |>
      group_by(waterbody, id_site, survey_year, species, gear_type, reach_no) |>
      summarise(log_cpue_ym1 = log(sum(catch) + 1) - log(median(effort_h))) |>
      ungroup(),
    by = c(
      "waterbody", "reach_no", "id_site", "species", "gear_type",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  left_join(
    rb_predictors |>
      mutate(
        across(contains("cpue_"), ~ log(. + 1), .names = "log_{.col}")
      ) |>
      select(waterbody, id_site, survey_year, gear_type, reach_no, contains("log_cpue_")),
    by = c(
      "waterbody", "id_site", "gear_type", "reach_no",
      "survey_year_minus1" = "survey_year"
    )
  ) |>
  filter(!is.na(log_cpue_ym1)) |>
  mutate(species_category = paste(species, category, sep = "_"))

# fit model
iter <- 200
chains <- 1
if (sample_again) {
  
  # sample from brms model
  stan_seed <- 2024-04-26
  
  # Murray Cod
  #   Predictors (|r| < 0.7): carp, LB, SB, redfin
  mod_mc <- brm(
    bf(
      catch ~
        log_cpue_ym1 +
        category * (
          log_cpue_lb + log_cpue_sb + 
            log_cpue_carp + log_cpue_redfin
        ) +
        (category | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           id_site:survey_year + gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | category +
          waterbody + waterbody:reach_no +
          gear_type + 
          category:waterbody +
          category:gear_type
      )
    ),
    data = mc_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # River Blackfish
  #   Predictors (|r| < 0.7): carp, SB, redfin
  mod_bf <- brm(
    bf(
      catch ~
        log_cpue_ym1 +
        category * (
          log_cpue_sb + log_cpue_carp + log_cpue_redfin
        ) +
        (category | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           id_site:survey_year + gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | category +
          waterbody + waterbody:reach_no +
          gear_type + 
          category:waterbody +
          category:gear_type
      )
    ),
    data = bf_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # Murray-Darling Rainbowfish
  #   Predictors (|r| < 0.7): carp, LB, SB, redfin
  mod_rb <- brm(
    bf(
      catch ~
        log_cpue_ym1 +
        category * (
          log_cpue_lb + log_cpue_sb + 
            log_cpue_carp + log_cpue_redfin
        ) +
        (category | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           id_site:survey_year + gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | category +
          waterbody + waterbody:reach_no +
          gear_type + 
          category:waterbody +
          category:gear_type
      )
    ),
    data = rb_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.8, max_treedepth = 10),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # save to file
  qsave(mod_mc, file = "outputs/fitted/draws-mc.qs")
  qsave(mod_bf, file = "outputs/fitted/draws-bf.qs")
  qsave(mod_rb, file = "outputs/fitted/draws-rb.qs")
  
} else {
  
  # load saved version
  mod_mc <- qread("outputs/fitted/draws-mc.qs")
  mod_bf <- qread("outputs/fitted/draws-bf.qs")
  mod_rb <- qread("outputs/fitted/draws-rb.qs")
  
}

# only repeat these steps if it's a newly fitted model
if (sample_again) {
  
  # basic diagnostics (returns NULL, saves outputs to files)
  diagnose_mcmc(mod_mc, suffix = "mc")
  diagnose_mcmc(mod_bf, suffix = "bf")
  diagnose_mcmc(mod_rb, suffix = "rb")
  
}
