# Canonical plotting implementation extracted from the accepted figures.
# Only names and input-path resolution changed; no scientific recomputation.
f8_tissue_plot <- function (emp, summ, mode, p = NA_real_, title = "A  Tissue and model ploidy") 
{
    x <- emp$sites[emp$sites$shared_site, ]
    ord <- c("bladder", "lung", "ovary", "skin", "esoph./stom", "biliary", "breast", "bone/soft tiss.", "bowel", "cervix", 
        "uterus", "liver", "lymph", "thyroid", "cns", "kidney", "pancreas", "prostate")
    stopifnot(length(unique(x$site)) == 18, setequal(unique(x$site), ord))
    x$context <- ifelse(x$cohort == "PCAWG", "in vivo", "in vitro")
    x$source <- x$cohort
    pair <- summ[summ$condition_mode == mode & summ$initial_ploidy == 2 & ((summ$model_context == "in vivo" & summ$O2_pct == 
        0.5) | (summ$model_context == "in vitro" & summ$O2_pct %in% c(0.5, 20))), ]
    if (mode == "fixed_p_misseg") 
        pair <- pair[abs(pair$p_misseg - p) < 1e-12, ]
    stopifnot(nrow(pair) == 3, all(pair$n_seeds == 1000))
    pair$site <- paste0("Model: ", pair$model_context, "; O2 = ", pair$O2_pct, "%")
    model <- data.frame(site = pair$site, context = pair$model_context, source = "Model", median = pair$median, q25 = pair$q25, 
        q75 = pair$q75)
    dat <- rbind(x[, c("site", "context", "source", "median", "q25", "q75")], model)
    models <- c("Model: in vivo; O2 = 0.5%", "Model: in vitro; O2 = 0.5%", "Model: in vitro; O2 = 20%")
    dat$site <- factor(dat$site, levels = rev(c(ord, " ", models)))
    sub <- if (mode == "fitted_p_misseg") 
        "Model: own fitted p_misseg; 2N cohort; day 1000"
    else paste0("Model: fixed p_misseg = ", p, "; 2N cohort; day 1000")
    ggplot2::ggplot(dat, ggplot2::aes(site, median, color = context, shape = source)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, 
        ymax = q75), position = ggplot2::position_dodge(0.45), width = 0.2) + ggplot2::geom_point(position = ggplot2::position_dodge(0.45), 
        size = 1.8) + ggplot2::geom_vline(xintercept = 3.5, colour = "#999999", linetype = "dashed", linewidth = 0.3) + ggplot2::coord_flip() + 
        ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_color_manual(values = f78_context_colors) + ggplot2::scale_shape_manual(values = c(PCAWG = 16, 
        DepMap = 17, Model = 18)) + ggplot2::labs(title = title, subtitle = sub, x = NULL, y = "Ploidy (N): median and IQR") + 
        f78_plot_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}
f8_author_figure5 <- function (a) 
{
    z <- f78_read_tsv(f78ro_path(a, "author_figure5_points.tsv"))
    ggplot2::ggplot(z, ggplot2::aes(wgd_fraction, oxygen)) + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, 
        level = 0.95, colour = "#0072B2", fill = "#B8D4E5", linewidth = 0.5) + ggplot2::geom_errorbar(ggplot2::aes(ymin = oxygen_low, 
        ymax = oxygen_high), width = 0.012, linewidth = 0.35, na.rm = TRUE) + ggplot2::geom_point(size = 1.8) + ggrepel::geom_text_repel(ggplot2::aes(label = organ), 
        size = 2.5, seed = 20261005, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.3) + ggplot2::labs(title = "B  Organ oxygen and WGD prevalence", 
        x = "WGD fraction", y = "Oxygen partial pressure (mmHg)") + f78_plot_theme()
}
f8_ab_theme <- function () 
ggplot2::theme(aspect.ratio = 1, axis.text = ggplot2::element_text(size = 12, colour = "#222222"), axis.text.y = ggplot2::element_text(size = 12), 
    axis.title = ggplot2::element_text(size = 14), legend.text = ggplot2::element_text(size = 12), plot.title = ggplot2::element_text(size = 16, 
        face = "bold"))
f8_median_iqr <- function (a) 
{
    empirical <- list(sites = f78_read_tsv(f78ro_path(a, "empirical_tissue_ploidy_summary.tsv")))
    summaries <- f78_read_tsv(f78ro_path(a, "curve_summary.tsv"))
    p <- f8_tissue_plot(empirical, summaries, "fitted_p_misseg") + ggplot2::labs(subtitle = NULL) + f8_ab_theme()
    p$layers[[1]]$aes_params$linewidth <- 0.55
    p$layers[[2]]$aes_params$size <- 2.2
    stopifnot(nrow(p$data) == 39L, all(p$data$q25 <= p$data$median & p$data$median <= p$data$q75))
    list(plot = p, summary = p$data)
}
f8_labeled_author_figure5 <- function (a) 
{
    p <- f8_author_figure5(a) + f8_ab_theme()
    stopifnot(inherits(p$layers[[4]]$geom, "GeomTextRepel"))
    p$layers[[4]] <- ggrepel::geom_label_repel(ggplot2::aes(label = organ), size = 3.5, seed = 20261005, max.overlaps = Inf, 
        min.segment.length = 0, box.padding = 0.3, label.padding = grid::unit(0.08, "lines"), label.r = grid::unit(0, "lines"), 
        label.size = NA, fill = "white", max.time = 2, max.iter = 20000)
    p
}
f8_mixture_colors <- c("2N" = "#0072B2", "4N" = "#CC79A7", "2N + 4N (1:1)" = "#E69F00", "2N + 4N (1:2)" = "#D55E00", "2N + 4N (1:3)" = "#009E73", 
"2N + 4N (1:4)" = "#56B4E9", "2N + 4N (1:5)" = "#7A3E9D")
f8_mixture_conditions <- function () 
data.table::data.table(Q = c(0, Inf, 1:5), initial_condition = c("2N", "4N", paste0("2N + 4N (1:", 1:5, ")")), label = c("Q = 0 (2N)", 
    "Q = ∞ (4N)", paste0("Q = ", 1:5)))
f8_mixture_curve_plot <- function (z, title, ylim) 
{
    z <- data.table::copy(data.table::as.data.table(z))
    stopifnot(nrow(z) == 357L, setequal(z$Q, c(0, Inf, 1:5)), all(z$model_context == "in vivo"))
    z$initial_condition <- factor(z$initial_condition, levels = names(f8_mixture_colors))
    labels <- merge(z[O2_pct == 5], f8_mixture_conditions()[, .(Q, label)], by = "Q", sort = FALSE)
    ggplot2::ggplot(z, ggplot2::aes(O2_pct, median, colour = initial_condition, fill = initial_condition)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, 
        ymax = q75), alpha = 0.1, colour = NA) + ggplot2::geom_line(linewidth = 0.7, linetype = "solid") + ggrepel::geom_text_repel(data = labels, 
        ggplot2::aes(label = label), size = 3, direction = "y", nudge_x = 0.15, hjust = 0, seed = 20261006, force = 2, force_pull = 0.05, 
        box.padding = 0.3, point.padding = 0.08, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, show.legend = FALSE) + 
        ggplot2::scale_color_manual(values = f8_mixture_colors, drop = FALSE) + ggplot2::scale_fill_manual(values = f8_mixture_colors, 
        drop = FALSE) + ggplot2::scale_x_continuous(breaks = 0:5, limits = c(0, 6.05), expand = ggplot2::expansion(mult = c(0.02, 
        0.01))) + ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = c(0.06, 0.08))) + ggplot2::labs(title = paste0(title, 
        "\nQ = 4N:2N"), x = "Imposed model oxygen (%)", y = "Day-1000 ploidy (N): median and IQR") + f78_plot_theme() + ggplot2::theme(legend.position = "none")
}
f8_mixture_panels <- function (a) 
{
    library(data.table)
    s <- fread(f78ro_path(a, "fitted_family_curve_summary.tsv"))
    low <- fread(f78ro_path(a, "supp_fig8_3_curve_summary.tsv"))
    high <- fread(f78ro_path(a, "supp_fig8_4_curve_summary.tsv"))
    stopifnot(nrow(s) == 714L, all(s$n_seeds == 500L), all(low$n_seeds == 1000L), all(high$n_seeds == 1000L))
    main_limits <- range(s$q25, s$q75)
    supp_limits <- range(low$q25, low$q75, high$q25, high$q75)
    list(C = f8_mixture_curve_plot(s[pair_label == "C01"], "C  C01; in vivo; fitted p_misseg", main_limits), D = f8_mixture_curve_plot(s[pair_label == 
        "C02"], "D  C02; in vivo; fitted p_misseg", main_limits), supp3 = f8_mixture_curve_plot(low, "In vivo initial conditions; p_misseg = 0.01", 
        supp_limits), supp4 = f8_mixture_curve_plot(high, "In vivo initial conditions; p_misseg = 0.3", supp_limits))
}
f8_start_curve_plot <- function (z, title, ylim) 
{
    z <- data.table::copy(data.table::as.data.table(z))
    levels <- names(f8_mixture_colors)[1:3]
    stopifnot(nrow(z) == 153L, setequal(z$Q, c(0, 1, Inf)), all(z$n_seeds == 1000L))
    z$initial_condition <- factor(z$initial_condition, levels = levels)
    ggplot2::ggplot(z, ggplot2::aes(O2_pct, median, colour = initial_condition, fill = initial_condition)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, 
        ymax = q75), alpha = 0.16, colour = NA, show.legend = FALSE) + ggplot2::geom_line(linewidth = 0.8, linetype = "solid") + 
        ggplot2::scale_color_manual(values = f8_mixture_colors[levels], drop = FALSE) + ggplot2::scale_fill_manual(values = f8_mixture_colors[levels], 
        drop = FALSE, guide = "none") + ggplot2::scale_x_continuous(breaks = 0:5, limits = c(0, 5), expand = ggplot2::expansion(mult = c(0.02, 
        0.02))) + ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = c(0.06, 0.08))) + ggplot2::labs(title = title, 
        x = "Imposed model oxygen (%)", y = "Day-1000 ploidy (N)") + f78_plot_theme() + ggplot2::theme(axis.text = ggplot2::element_text(size = 12), 
        axis.title = ggplot2::element_text(size = 14), plot.title = ggplot2::element_text(size = 16, face = "bold"), legend.text = ggplot2::element_text(size = 12), 
        legend.position = "bottom")
}
f8_culture_family_plot <- function (z) 
{
    z <- data.table::copy(data.table::as.data.table(z))
    stopifnot(nrow(z) == 1428L, all(z$n_seeds == 500L), all(z$model_context == "in vitro"))
    z$initial_condition <- factor(z$initial_condition, levels = names(f8_mixture_colors))
    z$pair_label <- factor(z$pair_label, levels = c("C01", "C02"))
    z$p_setting <- factor(z$p_setting, levels = c("0.01", "0.3"))
    ggplot2::ggplot(z, ggplot2::aes(O2_pct, median, colour = initial_condition, fill = initial_condition)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, 
        ymax = q75), alpha = 0.12, colour = NA, show.legend = FALSE) + ggplot2::geom_line(linewidth = 0.75, linetype = "solid") + 
        ggplot2::facet_grid(rows = ggplot2::vars(p_setting), cols = ggplot2::vars(pair_label), switch = "y", labeller = ggplot2::labeller(p_setting = c(`0.01` = "p_misseg = 0.01", 
            `0.3` = "p_misseg = 0.3"))) + ggplot2::scale_color_manual(values = f8_mixture_colors, labels = f8_mixture_conditions()$label, 
        drop = FALSE) + ggplot2::scale_fill_manual(values = f8_mixture_colors, drop = FALSE, guide = "none") + ggplot2::scale_x_continuous(breaks = 0:5, 
        limits = c(0, 5), expand = ggplot2::expansion(mult = c(0.02, 0.02))) + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, 
        byrow = TRUE)) + ggplot2::labs(title = "In vitro initial-condition comparison", x = "Imposed model oxygen (%)", y = "Day-1000 ploidy (N)") + 
        f78_plot_theme() + ggplot2::theme(strip.placement = "outside", strip.text = ggplot2::element_text(size = 14, face = "bold"), 
        axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14), plot.title = ggplot2::element_text(size = 16, 
            face = "bold"), legend.text = ggplot2::element_text(size = 11), legend.position = "bottom", panel.spacing = grid::unit(1.2, 
            "lines"))
}
f8_time_matrix_plot <- function (z, context) 
{
    z <- data.table::copy(data.table::as.data.table(z))
    days <- seq(100L, 1000L, 100L)
    maxo <- if (context == "in vitro") 
        20
    else 5
    colors <- c(`2N` = "#0072B2", `4N` = "#CC79A7", `2N + 4N (1:1)` = "#E69F00")
    stopifnot(setequal(z$day, days), all(z$n_seeds == 1000L), all(z$n_c01 == 500L & z$n_c02 == 500L), setequal(z$Q, c(0, 
        1, Inf)), all(z$model_context == context), nrow(z) == length(days) * 3L * (maxo * 10L + 1L))
    z[, `:=`(initial_condition, factor(initial_condition, levels = names(colors)))]
    z[, `:=`(snapshot, factor(paste0(day, " days"), levels = paste0(days, " days")))]
    ggplot2::ggplot(z, ggplot2::aes(O2_pct, median, colour = initial_condition, fill = initial_condition)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, 
        ymax = q75), alpha = 0.16, colour = NA, show.legend = FALSE) + ggplot2::geom_line(linewidth = 0.75, linetype = "solid") + 
        ggplot2::facet_wrap(ggplot2::vars(snapshot), nrow = 2, ncol = 5, dir = "h", scales = "fixed") + ggplot2::scale_color_manual(values = colors, 
        drop = FALSE) + ggplot2::scale_fill_manual(values = colors, drop = FALSE, guide = "none") + ggplot2::scale_x_continuous(limits = c(0, 
        maxo), breaks = if (maxo == 20) 
        seq(0, 20, 5)
    else 0:5, expand = ggplot2::expansion(mult = c(0.015, 0.015))) + ggplot2::scale_y_continuous(limits = range(z$q25, z$q75), 
        expand = ggplot2::expansion(mult = c(0.04, 0.06))) + ggplot2::labs(title = paste0(context, " finite-time ploidy"), 
        x = "Imposed model oxygen (%)", y = "Mean ploidy (N)", colour = NULL, fill = NULL) + ggplot2::theme_classic(base_size = 13, 
        base_family = "Arial") + ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 15, 
        face = "bold"), axis.text = ggplot2::element_text(size = 11), axis.title = ggplot2::element_text(size = 15), plot.title = ggplot2::element_text(size = 18, 
        face = "bold"), legend.position = "bottom", legend.text = ggplot2::element_text(size = 13), panel.spacing = grid::unit(1.1, 
        "lines"), plot.margin = ggplot2::margin(10, 15, 10, 10))
}
f8_family_time_blocks <- function (z, context) 
{
    z <- data.table::copy(data.table::as.data.table(z))
    maxo <- if (context == "in vitro") 
        20
    else 5
    days <- seq(100L, 1000L, 100L)
    colors <- c(`2N` = "#0072B2", `4N` = "#CC79A7", `2N + 4N (1:1)` = "#E69F00")
    stopifnot(all(z$n_seeds == 500L & z$unique_seed_ids == 500L), setequal(z$pair_label, c("C01", "C02")), all(z$model_context == 
        context), setequal(z$day, days), nrow(z) == 2L * 10L * 3L * (maxo * 10L + 1L))
    z[, `:=`(initial_condition, factor(initial_condition, levels = names(colors)))]
    z[, `:=`(snapshot, factor(paste0(day, " days"), levels = paste0(days, " days")))]
    ylim <- range(z$q25, z$q75)
    lapply(c("C01", "C02"), function(family) {
        ggplot2::ggplot(z[pair_label == family], ggplot2::aes(O2_pct, median, colour = initial_condition, fill = initial_condition)) + 
            ggplot2::geom_ribbon(ggplot2::aes(ymin = q25, ymax = q75), alpha = 0.16, colour = NA, show.legend = FALSE) + 
            ggplot2::geom_line(linewidth = 0.75, linetype = "solid") + ggplot2::facet_wrap(ggplot2::vars(snapshot), nrow = 2, 
            ncol = 5, dir = "h", scales = "fixed") + ggplot2::scale_color_manual(values = colors, drop = FALSE) + ggplot2::scale_fill_manual(values = colors, 
            drop = FALSE, guide = "none") + ggplot2::scale_x_continuous(limits = c(0, maxo), breaks = if (maxo == 20) 
            seq(0, 20, 5)
        else 0:5, expand = ggplot2::expansion(mult = c(0.015, 0.015))) + ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = c(0.04, 
            0.06))) + ggplot2::labs(title = family, x = "Imposed model oxygen (%)", y = "Mean ploidy (N)", colour = NULL, 
            fill = NULL) + ggplot2::theme_classic(base_size = 13, base_family = "Arial") + ggplot2::theme(strip.background = ggplot2::element_blank(), 
            strip.text = ggplot2::element_text(size = 15, face = "bold"), axis.text = ggplot2::element_text(size = 11), axis.title = ggplot2::element_text(size = 15), 
            plot.title = ggplot2::element_text(size = 18, face = "bold"), legend.position = "bottom", legend.text = ggplot2::element_text(size = 13), 
            panel.spacing = grid::unit(1.1, "lines"), plot.margin = ggplot2::margin(10, 15, 10, 10))
    })
}
f8_family_time_plot <- function (z, context) 
{
    patchwork::wrap_plots(f8_family_time_blocks(z, context), ncol = 1, guides = "collect") + patchwork::plot_annotation(title = paste0(context, 
        " finite-time ploidy"), theme = ggplot2::theme(text = ggplot2::element_text(family = "Arial"), plot.title = ggplot2::element_text(size = 20, 
        face = "bold"))) & ggplot2::theme(legend.position = "bottom")
}
