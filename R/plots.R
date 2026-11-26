## =============================================================================
##  plots.R -- figures of the framework study (ggplot2)
## =============================================================================
engine_colours <- c("none" = "#7a7a7a", "one-way EXNEX" = "#0072B2",
                    "one-way EXNEX, scaled Beta" = "#56B4E9", "two-way EXNEX" = "#D55E00")

## Rejection rate per basket: null baskets (open symbols) against alpha, active baskets = power.
plot_baskets <- function(b, recruitment, alpha) {
  library(ggplot2)
  b <- b[b$recruitment == recruitment, ]
  b$type <- ifelse(b$null, "null basket (type I error)", "active basket (power)")
  ggplot(b, aes(basket, reject, colour = engine, shape = type)) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = pmax(reject - 2 * reject_se, 0), ymax = pmin(reject + 2 * reject_se, 1)),
                   position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    scale_shape_manual(values = c("null basket (type I error)" = 1, "active basket (power)" = 16)) +
    scale_colour_manual(values = engine_colours) +
    facet_wrap(~ scenario, ncol = 2) +
    labs(x = "basket (tumour type . biomarker group)", y = "rejection rate", colour = "borrowing", shape = NULL,
         title = sprintf("Rejection per basket (%s recruitment); dashed line: alpha per basket", recruitment)) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical",
                                     axis.text.x = element_text(angle = 45, hjust = 1))
}

## Power against expected total sample size, one panel per scenario with active baskets.
plot_tradeoff <- function(pl, recruitment) {
  library(ggplot2)
  pl <- pl[pl$recruitment == recruitment & !is.na(pl$power_basket_mean), ]
  ggplot(pl, aes(N_mean, power_basket_mean, colour = engine)) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = engine_colours) +
    facet_wrap(~ scenario, ncol = 3, scales = "free_x") +
    labs(x = "expected total sample size", y = "mean power over active baskets", colour = "borrowing",
         title = sprintf("Power against sample size (%s recruitment)", recruitment)) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
}

save_plot <- function(p, file, width, height) {
  ggplot2::ggsave(paste0(file, ".pdf"), p, width = width, height = height)
  ggplot2::ggsave(paste0(file, ".png"), p, width = width, height = height, dpi = 150)
}
