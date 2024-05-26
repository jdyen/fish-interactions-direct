diagnose_mcmc <- function(x, suffix = NULL) {
  
  # check that a suffix is provided for filenames
  if (is.null(suffix))
    stop("suffix must be provided for filenames", call. = FALSE)
  
  # posterior checks
  pp_plot <- pp_check(x, type = "dens_overlay") +
    scale_x_log10() +
    theme(plot.background = element_rect(fill = "white"))
  pp_max <- bayesplot::pp_check(
    x, type = "stat", stat = "max"
  ) +
    theme(plot.background = element_rect(fill = "white"))
  pp_pzero <- bayesplot::pp_check(
    x, 
    type = "stat", 
    stat = \(x) mean(x == 0)
  ) +
    theme(plot.background = element_rect(fill = "white"))
  rhat <- brms::rhat(x)
  neff <- brms::neff_ratio(x)
  standard_diag <- tibble(
    stat = c(rep("Rhat", length(rhat)), rep("Neff ratio", length(neff))),
    par = c(names(rhat), names(neff)),
    value = c(rhat, neff)
  ) |>
    ggplot(aes(x = value)) +
    geom_histogram() +
    xlab("Value") +
    ylab("Count") +
    facet_wrap( ~ stat, scales = "free")
  
  # model fit
  fitted_vals <- posterior_epred(x)
  fitted_vals <- t(fitted_vals)
  colnames(fitted_vals) <- paste("pred", seq_len(ncol(fitted_vals)), sep = "_")
  
  # include category if available
  if (!is.null(x$data$category)) {
    r2 <- x$data |>
      mutate(est = apply(fitted_vals, 1, median)) |>
      group_by(category) |>
      summarise(cor = cor(catch, est) ^ 2) |>
      ungroup()
  } else {
    r2 <- x$data |>
      mutate(est = apply(fitted_vals, 1, median)) |>
      summarise(cor = cor(catch, est) ^ 2) |>
      ungroup()
  }
  write.csv(r2, file = paste0("outputs/tables/r2-", suffix,".csv"))
  
  # return nothing
  list(pp_plot, pp_max, pp_pzero, standard_diag)
  
}
