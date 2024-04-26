diagnose_mcmc <- function(x, suffix = NULL) {
  
  # check that a suffix is provided for filenames
  if (is.null(suffix))
    stop("suffix must be provided for filenames", call. = FALSE)

  # posterior checks
  pp_plot <- pp_check(x, group = "species", type = "dens_overlay_grouped") +
    scale_x_log10() +
    theme(plot.background = element_rect(fill = "white"))
  pp_max <- bayesplot::pp_check(
    x, group = "species", type = "stat_grouped", stat = "max"
  ) +
    theme(plot.background = element_rect(fill = "white"))
  pp_pzero <- bayesplot::pp_check(
    x, 
    group = "species",
    type = "stat_grouped", 
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
  
  # save these
  ggsave(
    filename = paste0("outputs/figures/ppcheck-", suffix, ".png"),
    plot = pp_plot + theme(legend.position = "none"),
    device = ragg::agg_png,
    width = 6, 
    height = 6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = paste0("outputs/figures/ppmax-", suffix, ".png"),
    plot = pp_max + scale_x_log10() + theme(legend.position = "none"),
    device = ragg::agg_png,
    width = 6, 
    height = 6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = paste0("outputs/figures/ppzero-", suffix, ".png"),
    plot = pp_pzero + theme(legend.position = "none"),
    device = ragg::agg_png,
    width = 6, 
    height = 6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = paste0("outputs/figures/diagnostics-", suffix, ".png"),
    plot = standard_diag,
    device = ragg::agg_png,
    width = 6, 
    height = 6,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  # model fit
  fitted_vals <- posterior_epred(x)
  fitted_vals <- t(fitted_vals)
  colnames(fitted_vals) <- paste("pred", seq_len(ncol(fitted_vals)), sep = "_")
  r2 <- x$data |>
    mutate(est = apply(fitted_vals, 1, median)) |>
    group_by(category) |>
    summarise(cor = cor(catch, est) ^ 2) |>
    ungroup()
  write.csv(r2, file = paste0("outputs/tables/r2-", suffix,".csv"))
  
  # return nothing
  NULL
  
}
