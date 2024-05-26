# VEFMAP fish population modelling:
#   Analysis of fish data to determine strength of direct interactions
#   with invasive species and functional groups of species (incl. native 
#   species)
#
# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
#
# Last updated: 16 May 2024

# flags for slow steps
reload_data <- FALSE
sample_again <- FALSE

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
library(ggeffects)
library(patchwork)

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

# tidy up the categories for MC and BF (remove from analysis for RB)
cpue <- cpue |>
  mutate(
    category = ifelse(
      species == "Maccullochella peelii", 
      gsub("tiny", "small", category),
      category
    ),
    category = ifelse(
      species == "Maccullochella peelii", 
      gsub("large", "medium", category),
      category
    ),
    category = ifelse(
      species == "Maccullochella peelii", 
      gsub("huge", "large", category),
      category
    ),
    
    category = ifelse(
      species == "Gadopsis marmoratus", 
      gsub("small", "medium", category),
      category
    ),
    category = ifelse(
      species == "Gadopsis marmoratus", 
      gsub("tiny", "small", category),
      category
    ),
    category = ifelse(
      species == "Gadopsis marmoratus", 
      gsub("huge", "large", category),
      category
    ),
    
    category = ifelse(species == "Melanotaenia fluviatilis", "small", category)
  )

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
  filter(!is.na(log_cpue_ym1))
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
  filter(!is.na(log_cpue_ym1))
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
  filter(!is.na(log_cpue_ym1))

# scale all predictors
mc_cpue <- mc_cpue |>
  mutate(
    across(
      contains("log_cpue"),
      \(x) (x - mean(x)) / sd(x), 
      .names = "{.col}_std"
    )
  )
bf_cpue <- bf_cpue |>
  mutate(
    across(
      contains("log_cpue"),
      \(x) (x - mean(x)) / sd(x), 
      .names = "{.col}_std"
    )
  )
rb_cpue <- rb_cpue |>
  mutate(
    across(
      contains("log_cpue"),
      \(x) (x - mean(x)) / sd(x), 
      .names = "{.col}_std"
    )
  )

# fit model
iter <- 4000
chains <- 4
if (sample_again) {
  
  # sample from brms model
  stan_seed <- 2024-04-26
  
  # Murray Cod
  #   Predictors (|r| < 0.7): carp_large, carp_small, SB, redfin_large, redfin_small
  #   NOTE: carp_large associated with cpue_lb, so might be confounded
  mod_mc <- brm(
    bf(
      catch ~
        log_cpue_ym1_std +
        category * (
          log_cpue_sb_std + 
            log_cpue_carp_large_std + log_cpue_carp_small_std +
            log_cpue_redfin_large_std + log_cpue_redfin_small_std
        ) +
        (-1 + category | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | category +
          waterbody + waterbody:reach_no +
          gear_type
      )
    ),
    data = mc_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.85, max_treedepth = 12),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # River Blackfish
  #   Predictors (|r| < 0.7): carp_large, carp_small, SB, redfin_large, redfin_small
  #   NOTE: carp_large associated with cpue_lb, so might be confounded
  mod_bf <- brm(
    bf(
      catch ~
        log_cpue_ym1_std +
        category * (
          log_cpue_sb_std + 
            log_cpue_carp_large_std + log_cpue_carp_small_std +
            log_cpue_redfin_large_std + log_cpue_redfin_small_std
        ) +
        (-1 + category | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | category +
          waterbody + waterbody:reach_no +
          gear_type
      )
    ),
    data = bf_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.85, max_treedepth = 12),
    backend = "rstan",
    threads = threading(2),
    init_r = 1.0
  )
  
  # Murray-Darling Rainbowfish
  #   Predictors (|r| < 0.7): carp_large, carp_small, SB, redfin_large, redfin_small
  #   NOTE: carp_large associated with cpue_lb, so might be confounded
  mod_rb <- brm(
    bf(
      catch ~
        log_cpue_ym1_std +
        log_cpue_sb_std + 
        log_cpue_carp_large_std + log_cpue_carp_small_std +
        log_cpue_redfin_large_std + log_cpue_redfin_small_std +
        (1 | waterbody + waterbody:reach_no + id_site +
           survey_year + waterbody:survey_year +
           waterbody:reach_no:survey_year +
           gear_type) +
        offset(log_effort_h),
      shape ~ (
        1 | waterbody + 
          waterbody:reach_no +
          gear_type
      )
    ),
    data = rb_cpue |> mutate(log_effort_h = log(effort_h)),
    family = negbinomial(),
    chains = chains,
    cores = chains,
    seed = stan_seed,
    iter = iter, 
    warmup = floor(iter / 2),
    control = list(adapt_delta = 0.85, max_treedepth = 12),
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
  pp_mc <- diagnose_mcmc(mod_mc, suffix = "mc")
  pp_bf <- diagnose_mcmc(mod_bf, suffix = "bf")
  pp_rb <- diagnose_mcmc(mod_rb, suffix = "rb")
  
  pp_plots <- ((pp_mc[[1]] + theme(legend.position = "none")) |
      (pp_mc[[2]] + scale_x_log10() + theme(legend.position = "none")) | 
      (pp_mc[[3]] + theme(legend.position = "none"))) /
    ((pp_bf[[1]] + theme(legend.position = "none")) |
       (pp_bf[[2]] + scale_x_log10() + theme(legend.position = "none")) | 
       (pp_bf[[3]] + theme(legend.position = "none"))) /
    ((pp_rb[[1]] + theme(legend.position = "none")) |
       (pp_rb[[2]] + scale_x_log10() + theme(legend.position = "none")) | 
       (pp_rb[[3]] + theme(legend.position = "none"))) +
    plot_annotation(tag_levels = "a")
  standard_diag <- pp_mc[[4]] / pp_bf[[4]] / pp_rb[[4]] + 
    plot_annotation(tag_levels = "a")
      
  # save these
  ggsave(
    filename = "outputs/figures/ppcheck-all.png",
    plot = pp_plots,
    device = ragg::agg_png,
    width = 8, 
    height = 8,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = "outputs/figures/diagnostics-all.png",
    plot = standard_diag,
    device = ragg::agg_png,
    width = 6, 
    height = 8,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
}

# summarise the fitted models
vars <- c("log_cpue_sb_std [-1.9:2.3 by=0.4]",
          "log_cpue_carp_large_std [-0.5:3.3 by=0.4]",
          "log_cpue_carp_small_std [-0.4:4.0 by=0.4]",
          "log_cpue_redfin_large_std [-0.25:6.0 by=0.4]",
          "log_cpue_redfin_small_std [-0.35:6.0 by=0.4]")

# ranges deteremined from here (without seq bits of vars):
# mod_mc$data[, vars] |> apply(2, range)
# mod_bf$data[, vars] |> apply(2, range)
# mod_rb$data[, vars] |> apply(2, range)

pred_mc <- lapply(
  vars,
  \(x) ggpredict(mod_mc, terms = c(x, "category"), ci.lvl = 0.1)
)
pred_bf <- lapply(
  vars,
  \(x) ggpredict(mod_bf, terms = c(x, "category"), ci.lvl = 0.1)
)
pred_rb <- lapply(
  vars,
  \(x) ggpredict(mod_rb, terms = c(x), ci.lvl = 0.1)
)

pred_all <- bind_rows(
  mapply(
    \(.x, .y) .x |> as_tibble() |> mutate(var = .y),
    .x = pred_mc,
    .y = vars,
    SIMPLIFY = FALSE
  ) |>
    bind_rows() |> 
    mutate(species = "Murray Cod"),
  mapply(
    \(.x, .y) .x |> as_tibble() |> mutate(var = .y),
    .x = pred_bf,
    .y = vars,
    SIMPLIFY = FALSE
  ) |>
    bind_rows() |>
    mutate(species = "River Blackfish"),
  mapply(
    \(.x, .y) .x |> as_tibble() |> mutate(var = .y),
    .x = pred_rb,
    .y = vars,
    SIMPLIFY = FALSE
  ) |>
    bind_rows() |> 
    mutate(species = "Murray-Darling Rainbowfish")
)

# TODO: back-transform predictors to raw scales
rescale_values <- bind_rows(
  tibble(
    species = "Murray Cod",
    var = mc_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      colnames(),
    mid = mc_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, mean),
    scale = mc_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, sd)
  ),
  tibble(
    species = "River Blackfish",
    var = bf_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      colnames(),
    mid = bf_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, mean),
    scale = bf_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, sd)
  ),
  tibble(
    species = "Murray-Darling Rainbowfish",
    var = rb_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      colnames(),
    mid = rb_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, mean),
    scale = rb_cpue |> 
      select(contains("log_cpue") & !contains("_std")) |>
      apply(2, sd)
  )
)

var_lu <- c(
  "log_cpue_sb" = "Small-bodied species",
  "log_cpue_carp_large" = "Large Common Carp",
  "log_cpue_carp_small" = "Small Common Carp",
  "log_cpue_redfin_large" = "Large European Perch",
  "log_cpue_redfin_small" = "Small European Perch"
)
mean_field_effects <- pred_all |>
  mutate(
    var = sapply(var, \(x) strsplit(x, " \\[")[[1]][1]),
    var = gsub("_std", "", var)
  ) |>
  left_join(rescale_values, by = c("var", "species")) |>
  mutate(
    x = mid + scale * x,
    group = as.character(group),
    group = ifelse(group == "1", "small", group),
    group = factor(
      group,
      levels = c("small", "medium", "large"),
      labels = c("Small", "Medium", "Large")
    ),
    var = var_lu[var],
    var = factor(var, levels = var_lu[c(1, 3, 2, 5, 4)]),
    species = factor(
      species,
      levels = c("Murray Cod", "River Blackfish", "Murray-Darling Rainbowfish")
    )
  ) |>
  ggplot(
    aes(
      y = predicted,
      x = exp(x),
      ymin = conf.low,
      ymax = conf.high,
      col = group,
      fill = group
    ),
  ) +
  geom_line() +
  geom_ribbon(alpha = 0.25, col = NA) +
  scale_color_brewer(name = "", palette = "Set2") +
  scale_fill_brewer(name = "", palette = "Set2") +
  xlab("Predictor value") +
  ylab("Predicted CPUE") +
  facet_wrap(species ~ var, scales = "free", ncol = 5)

ggsave(
  filename = "outputs/figures/mean-field-effects.png",
  plot = mean_field_effects + theme(legend.position = "bottom"),
  device = ragg::agg_png,
  width = 11, 
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
