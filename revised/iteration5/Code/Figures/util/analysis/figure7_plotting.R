# Canonical plotting implementation extracted from the accepted figures.
# Only names and input-path resolution changed; no scientific recomputation.
f7_beta_colors <- c("0" = "#222222", "0.5" = "#E69F00", "1" = "#D55E00")
f7_beta_labels <- c("0" = "β = 0", "0.5" = "β = 0.5", "1" = "β = 1")
f7_metric_labels <- function (version, current_growth = FALSE, direct_growth = FALSE) 
{
    if (version == "instantaneous") 
        function(metric) if (metric == "flux_percent") 
            "Immediate high-to-low inflow (%)"
        else if (direct_growth) 
            "Direct H-to-L-source share of low-cell proliferation (%)"
        else if (current_growth) 
            "H-to-L-source share of low-cell proliferation (%)"
        else "Immediate inflow contribution to low-cell production (%)"
    else function(metric) if (metric == "flux_percent") 
        "Historical H-to-L share of current low cells (%)"
    else "Historical-lineage contribution to low-cell proliferation (%)"
}
f7_line_plot <- function (z, metric = "flux_percent", axes = TRUE, y_limits = NULL) 
{
    library(ggplot2)
    z <- z[z$metric == metric, ]
    z$q <- factor(z$q, levels = c(0, 0.5, 1))
    p <- ggplot(z, aes(day, median, color = q, fill = q, linetype = q)) + geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.14, 
        color = NA, linetype = 0, na.rm = TRUE) + geom_line(linewidth = 0.55, na.rm = TRUE) + scale_color_manual(values = f7_beta_colors, 
        drop = FALSE, labels = f7_beta_labels) + scale_fill_manual(values = f7_beta_colors, drop = FALSE, labels = f7_beta_labels) + 
        scale_linetype_manual(values = c(`0` = "solid", `0.5` = "dashed", `1` = "dotdash"), drop = FALSE, labels = f7_beta_labels) + 
        scale_x_continuous(breaks = c(0, 500, 1000), limits = c(0, 1000), expand = expansion(mult = c(0.01, 0.02))) + labs(x = if (axes) 
        "Time (days)"
    else NULL, y = if (axes) 
        f7_metric_label(metric)
    else NULL) + f78_plot_theme()
    if (metric == "flux_percent") 
        p <- p + scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100))
    else p <- p + geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.2) + scale_y_continuous(limits = if (isTRUE(getOption("figure7.proliferation"))) 
        c(0, 100)
    else y_limits, breaks = if (isTRUE(getOption("figure7.proliferation"))) 
        c(0, 50, 100)
    else scales::breaks_pretty(n = 3))
    p
}
f7_panel_block <- function (z, oxygen, p, title, metric = "flux_percent") 
{
    library(ggplot2)
    library(cowplot)
    canvas <- ggdraw()
    cols <- c(C01 = "#C99700", C02 = "#6A3D9A")
    header <- function(txt, color, x, y, w, h, angle = 0, size = 9) {
        ggdraw() + draw_label(txt, color = "white", angle = angle, fontfamily = "Helvetica", fontface = "bold", size = size) + 
            theme(plot.background = element_rect(fill = color, colour = "#BEBEBE", linewidth = 0.2))
    }
    xpos <- c(0.18, 0.585)
    ww <- 0.395
    hh <- 0.205
    ytop <- c(0.72, 0.505, 0.265, 0.05)
    canvas <- canvas + draw_label(title, x = 0.18, y = 0.985, hjust = 0, fontface = "bold", size = 13)
    for (ci in 1:2) canvas <- canvas + draw_plot(header(names(cols)[ci], cols[ci], 0, 0, 1, 1), xpos[ci], 0.938, ww, 0.035)
    context <- c("in vivo", "in vivo", "in vitro", "in vitro")
    initial <- c(4, 2, 4, 2)
    for (r in 1:4) {
        row_data <- z[z$model_context == context[r] & z$initial_ploidy == initial[r] & abs(z$O2_pct - oxygen) < 1e-12 & abs(z$p_misseg - 
            p) < 1e-12 & z$threshold_N == 66 & z$metric == metric, ]
        y_limits <- if (metric == "growth_percent") {
            if (isTRUE(getOption("figure7.proliferation"))) 
                c(0, 100)
            else range(c(0, row_data$median, row_data$q25, row_data$q75), na.rm = TRUE)
        }
        else NULL
        if (!is.null(y_limits) && diff(y_limits) == 0) 
            y_limits <- y_limits + c(-1, 1)
        tag <- ggdraw() + draw_label(paste0(initial[r], "N"), angle = 90, size = 9) + theme(plot.background = element_rect(fill = "#F2F2F2", 
            colour = "#BEBEBE", linewidth = 0.2))
        canvas <- canvas + draw_plot(tag, 0.11, ytop[r], 0.028, hh)
        for (ci in 1:2) {
            zz <- z[z$model_context == context[r] & z$initial_ploidy == initial[r] & z$pair_label == names(cols)[ci] & abs(z$O2_pct - 
                oxygen) < 1e-12 & abs(z$p_misseg - p) < 1e-12 & z$threshold_N == 66 & z$metric == metric, ]
            stopifnot(nrow(zz) == 3 * 1001)
            plot <- f7_line_plot(zz, metric = metric, axes = FALSE, y_limits = y_limits) + theme(legend.position = "none", 
                plot.margin = margin(3, 3, 3, 3))
            if (r != 4) 
                plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
            if (ci != 1) 
                plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
            canvas <- canvas + draw_plot(plot, xpos[ci], ytop[r], ww, hh)
        }
    }
    for (rr in c(1, 3)) canvas <- canvas + draw_plot(header(context[rr], f78_context_colors[[context[rr]]], 0, 0, 1, 1, 90), 0.07, 
        ytop[rr + 1], 0.028, 2 * hh + 0.01)
    canvas + draw_label(f7_metric_label(metric), x = 0.025, y = 0.49, angle = 90, size = 10) + draw_label("Time (days)", 
        x = 0.59, y = 0.015, size = 10)
}
f7_matrix_pages <- function (a, z, metric, stem, publish = TRUE) 
{
    library(ggplot2)
    if (metric == "growth_percent" && !requireNamespace("ggh4x", quietly = TRUE)) 
        stop("Growth matrices require ggh4x for row/column strips with independent cell y scales")
    pages <- list()
    labels <- character()
    for (ctx in c("in vivo", "in vitro")) for (fam in c("C01", "C02")) for (initial in c(4, 2)) {
        oxygen <- if (ctx == "in vivo") 
            a$cfg$matrix_oxygen_invivo
        else a$cfg$matrix_oxygen_invitro
        zz <- z[z$model_context == ctx & z$pair_label == fam & z$initial_ploidy == initial & z$p_misseg == a$cfg$matrix_p_misseg & 
            z$O2_pct %in% oxygen & z$metric == metric, ]
        zz$q <- factor(zz$q, levels = c(0, 0.5, 1))
        zz$threshold <- factor(zz$threshold_N, levels = a$cfg$matrix_thresholds_N, labels = paste0(a$cfg$matrix_thresholds_N/22, 
            "N"))
        zz$oxygen <- factor(zz$O2_pct, levels = oxygen, labels = paste0(oxygen, "%"))
        ylab <- f7_metric_label(metric)
        p <- ggplot(zz, aes(day, median, color = q, fill = q, linetype = q)) + geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.14, 
            color = NA, linetype = 0, na.rm = TRUE) + geom_line(linewidth = 0.4, na.rm = TRUE) + (if (metric == "growth_percent") 
            ggh4x::facet_grid2(rows = vars(threshold), cols = vars(oxygen), scales = "free_y", independent = "y", axes = "y")
        else facet_grid(threshold ~ oxygen)) + scale_color_manual(values = f7_beta_colors, drop = FALSE, labels = f7_beta_labels) + 
            scale_fill_manual(values = f7_beta_colors, drop = FALSE, labels = f7_beta_labels) + scale_linetype_manual(values = c(`0` = "solid", 
            `0.5` = "dashed", `1` = "dotdash"), drop = FALSE, labels = f7_beta_labels) + scale_x_continuous(breaks = c(0, 
            500, 1000), limits = c(0, 1000)) + labs(title = paste(stem, ctx, fam, paste0(initial, "N start"), paste0("p_misseg = ", 
            a$cfg$matrix_p_misseg), sep = " | "), x = "Time (days)", y = ylab) + f78_plot_theme() + theme(strip.text = element_text(size = 9), 
            axis.text = element_text(size = 8))
        if (metric == "flux_percent") 
            p <- p + scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100))
        else p <- p + geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.2) + scale_y_continuous(limits = if (isTRUE(getOption("figure7.proliferation"))) 
            c(0, 100)
        else NULL, breaks = if (isTRUE(getOption("figure7.proliferation"))) 
            c(0, 50, 100)
        else scales::breaks_pretty(n = 3))
        pages[[length(pages) + 1L]] <- p
        label <- paste0(stem, "_", if (ctx == "in vivo") 
            "invivo"
        else "invitro", "_", tolower(fam), "_", initial, "N")
        labels <- c(labels, label)
        if (publish) 
            f78_publish(a, p, if (ctx == "in vivo") 
                15
            else 17, 10, label)
    }
    if (!publish) 
        return(list(plots = pages, labels = labels))
    stop("Use the canonical render transaction")
}
f7_draw_sources <- function (a, products = c("supp_fig7_1", "supp_fig7_2")) 
{
    library(ggplot2)
    library(patchwork)
    h <- f78_read_tsv(f78ro_path(a, "source_daughter_flux_binned.tsv"))
    b <- f78_read_tsv(f78ro_path(a, "winner_low_population_budget.tsv"))
    hp <- bp <- list()
    for (ctx in c("in vivo", "in vitro")) for (fam in c("C01", "C02")) {
        hh <- h[h$model_context == ctx & h$pair_label == fam, ]
        hh$oxygen <- paste0("O2 = ", hh$O2_pct, "%")
        hp[[length(hp) + 1L]] <- ggplot(hh, aes(source, daughter, fill = log10(flux + 1e-06))) + geom_tile() + geom_abline(slope = 1, 
            intercept = 0, color = "#999999", linewidth = 0.3) + facet_grid(branch ~ oxygen) + scale_fill_viridis_c(limits = range(log10(h$flux + 
            1e-06)), name = "log10 daughters /\n1000 live cells / day") + coord_equal() + labs(title = paste(ctx, fam), x = "Mother ploidy (N)", 
            y = "Daughter ploidy (N)") + f78_plot_theme() + theme(legend.title = element_text(size = 8))
        bb <- b[b$model_context == ctx & b$pair_label == fam, ]
        if (!is.null(a$revision)) {
            fields <- c("external_HL", "internal_LL", "maternal_removal", "direct_death", "low_net_change")
            for (field in fields) bb[[field]] <- ifelse(bb$low_fraction > a$revision$low_fraction_tolerance, bb[[field]]/bb$low_fraction, 
                NA_real_)
        }
        bb$oxygen <- paste0("O2 = ", bb$O2_pct, "%")
        ll <- data.table::melt(data.table::as.data.table(bb), id.vars = c("O2_pct", "oxygen", "low_net_change"), measure.vars = c("external_HL", 
            "internal_LL", "maternal_removal", "direct_death"))
        ll$value[ll$variable %in% c("maternal_removal", "direct_death")] <- -ll$value[ll$variable %in% c("maternal_removal", 
            "direct_death")]
        bp[[length(bp) + 1L]] <- ggplot(ll, aes(variable, 1000 * value, fill = variable)) + geom_col() + geom_hline(yintercept = 0, 
            linewidth = 0.25) + geom_hline(data = bb, aes(yintercept = 1000 * low_net_change), linetype = 2) + facet_wrap(~oxygen) + 
            scale_x_discrete(labels = c(external_HL = "High-to-low births", internal_LL = "Low-to-low births", maternal_removal = "Dividing mothers", 
                direct_death = "Direct death")) + scale_fill_manual(values = c(external_HL = "#D55E00", internal_LL = "#009E73", 
            maternal_removal = "#666666", direct_death = "#CC79A7"), guide = "none") + labs(title = paste(ctx, fam), x = NULL, 
            y = if (!is.null(a$revision)) 
                "Cells / 1000 low-ploidy cells / day"
            else "Cells / 1000 live cells / day") + f78_plot_theme() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
    if ("supp_fig7_1" %in% products) 
        f78_publish(a, wrap_plots(hp, ncol = 2), 14, 13, "supp_fig7_1")
    if ("supp_fig7_2" %in% products) 
        f78_publish(a, wrap_plots(bp, ncol = 2), 14, 10, "supp_fig7_2")
}
f7_draw_lineage_diagnostics <- function (a, products = c("supp_fig7_1_ver2", "supp_fig7_2_ver2")) 
{
    library(ggplot2)
    library(patchwork)
    library(data.table)
    h <- fread(f78ro_path(a, "history_source_daughter_flux_binned.tsv"))
    b <- fread(f78ro_path(a, "history_low_population_budget.tsv"))
    h[, `:=`(lineage, factor(lineage, levels = c("Historical H-to-L", "Never H-to-L")))]
    b[, `:=`(lineage, factor(lineage, levels = c("Historical H-to-L", "Never H-to-L")))]
    hp <- bp <- list()
    fields <- c("external_HL", "internal_LL", "maternal_removal", "direct_death", "low_net_change")
    for (ctx in c("in vivo", "in vitro")) for (fam in c("C01", "C02")) {
        hh <- h[model_context == ctx & pair_label == fam]
        hh[, `:=`(oxygen, factor(O2_pct, levels = c(0, 5), labels = c("O2 = 0%", "O2 = 5%")))]
        hp[[length(hp) + 1L]] <- ggplot(hh, aes(source, daughter, fill = log10(flux + 1e-06))) + geom_tile() + geom_abline(slope = 1, 
            intercept = 0, color = "#999999", linewidth = 0.3) + facet_grid(lineage + branch ~ oxygen) + scale_fill_viridis_c(limits = range(log10(h$flux + 
            1e-06)), name = "log10 daughters /\n1000 live cells / day") + coord_equal() + labs(title = paste(ctx, fam), x = "Mother ploidy (N)", 
            y = "Daughter ploidy (N)") + f78_plot_theme() + theme(strip.text = element_text(size = 9), legend.title = element_text(size = 9))
        bb <- copy(b[model_context == ctx & pair_label == fam])
        for (field in fields) bb[, `:=`((field), ifelse(low_fraction > 1e-14, get(field)/low_fraction, NA_real_))]
        bb[, `:=`(oxygen, factor(O2_pct, levels = c(0, 5), labels = c("O2 = 0%", "O2 = 5%")))]
        mm <- melt(bb, id.vars = c("oxygen", "lineage", "low_net_change"), measure.vars = fields[1:4])
        mm[variable %in% c("maternal_removal", "direct_death"), `:=`(value, -value)]
        bp[[length(bp) + 1L]] <- ggplot(mm, aes(variable, 1000 * value, fill = variable)) + geom_col() + geom_hline(yintercept = 0, 
            linewidth = 0.25) + geom_hline(data = bb, aes(yintercept = 1000 * low_net_change), linetype = 2) + facet_grid(lineage ~ 
            oxygen) + scale_x_discrete(labels = c(external_HL = "High-to-low births", internal_LL = "Low-to-low births", 
            maternal_removal = "Dividing mothers", direct_death = "Direct death")) + scale_fill_manual(values = c(external_HL = "#D55E00", 
            internal_LL = "#009E73", maternal_removal = "#666666", direct_death = "#CC79A7"), guide = "none") + labs(title = paste(ctx, 
            fam), x = NULL, y = "Cells / 1000 low-ploidy cells / day") + f78_plot_theme() + theme(axis.text.x = element_text(angle = 30, 
            hjust = 1, size = 9), strip.text = element_text(size = 9))
    }
    if ("supp_fig7_1_ver2" %in% products) 
        f78_publish(a, wrap_plots(hp, ncol = 2), 16, 23, "supp_fig7_1_ver2")
    if ("supp_fig7_2_ver2" %in% products) 
        f78_publish(a, wrap_plots(bp, ncol = 2), 16, 16, "supp_fig7_2_ver2")
}
