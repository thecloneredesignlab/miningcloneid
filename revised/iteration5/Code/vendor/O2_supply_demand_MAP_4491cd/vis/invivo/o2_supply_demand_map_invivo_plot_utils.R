# Pure plotting helpers for materialized in-vivo simulation tables.
# This file must not read fitted inputs, raw observations, or invoke simulation.

.invivo_plot_utils_dir <- local({
    frame_files <- Filter(
        nzchar,
        vapply(sys.frames(), function(env) {
            ofile <- env$ofile
            if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
        }, character(1))
    )
    own <- frame_files[
        basename(frame_files) == "o2_supply_demand_map_invivo_plot_utils.R"
    ]
    if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(
    file.path(.invivo_plot_utils_dir, "..", "o2_supply_demand_map_common_plot_utils.R"),
    local = TRUE,
    chdir = TRUE
)
rm(.invivo_plot_utils_dir)

ploidy_fraction_fill_scale <- o2sd_ploidy_fraction_fill_scale

weighted_mean_series_label <-
function(cfg) {
    mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
    if (identical(mode, "chr_number")) {
        "Weighted Mean Chromosome Number (N)"
    }
    else {
        "Weighted Mean Ploidy (P = N / N_UNIT)"
    }
}

functional_state_axis_label <-
function(cfg) {
    mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
    if (identical(mode, "chr_number"))
        "Chromosome Number (N)"
    else "Ploidy"
}

resource_death_language <-
function() {
    list(live_label = "Viable", component_label = "Hypoxia-Origin Dead", cin_component_label = "CIN-Associated Dead", adjective = "hypoxia-origin",
        figure_phrase = "hypoxia-origin dead", report_phrase = "hypoxia-origin dead")
}

log10_plot_floor <-
function(x, default = 1e-12) {
    x <- suppressWarnings(as.numeric(x))
    positive <- x[is.finite(x) & x > 0]
    if (length(positive) > 0L) {
        floor_use <- min(positive, na.rm = TRUE)/10
        if (is.finite(floor_use) && floor_use > 0)
            return(floor_use)
    }
    default <- suppressWarnings(as.numeric(default))
    if (!is.finite(default) || default <= 0)
        default <- 1e-12
    default
}

floor_for_log10_plot <-
function(x, floor) {
    x <- suppressWarnings(as.numeric(x))
    floor <- suppressWarnings(as.numeric(floor))
    if (!is.finite(floor) || floor <= 0)
        floor <- 1e-12
    out <- x
    finite <- is.finite(out)
    out[finite & out <= floor] <- floor
    out[!finite] <- NA_real_
    out
}

make_burden_decomp_long <-
function(burden_decomp, death_language) {
    burden_decomp %>% pivot_longer(cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"), names_to = "component",
        values_to = "value") %>% mutate(component = factor(component, levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
        labels = c(death_language$live_label, death_language$component_label, death_language$cin_component_label)))
}

make_burden_decomp_ribbon <-
function(burden_decomp, death_language, floor) {
    component_levels <- c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    burden_decomp %>% mutate(burden_live = pmax(as.numeric(burden_live), 0), burden_dead_hypoxia = pmax(as.numeric(burden_dead_hypoxia),
        0), burden_dead_buffer = pmax(as.numeric(burden_dead_buffer), 0)) %>% pivot_longer(cols = c("burden_live", "burden_dead_hypoxia",
        "burden_dead_buffer"), names_to = "component_raw", values_to = "value") %>% mutate(component = factor(component_raw,
        levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"), labels = component_levels), component_index = as.integer(component)) %>%
        arrange(cohort, day, component_index) %>% group_by(cohort, day) %>% mutate(ymax_raw = cumsum(value), ymin_raw = ymax_raw -
        value, ymin = floor_for_log10_plot(ymin_raw, floor), ymax = floor_for_log10_plot(ymax_raw, floor)) %>% ungroup() %>%
        mutate(component = factor(component, levels = component_levels))
}

format_original_scale_labels <-
function(x) {
    vapply(as.numeric(x), function(z) {
        if (!is.finite(z))
            return("")
        if (z != 0 && (abs(z) >= 10000 || abs(z) < 0.01)) {
            return(formatC(z, format = "e", digits = 1))
        }
        format(signif(z, 3), trim = TRUE, scientific = FALSE)
    }, character(1))
}

format_log10_axis_labels <-
function(x) {
    vapply(log10(as.numeric(x)), function(z) {
        if (!is.finite(z))
            return("")
        if (abs(z - round(z)) < 1e-08)
            return(as.character(as.integer(round(z))))
        format(round(z, 2), trim = TRUE, scientific = FALSE)
    }, character(1))
}

log10_burden_y_scale <-
function() {
    scale_y_log10(labels = format_log10_axis_labels)
}

log10_original_breaks <-
function(x, floor, n = 4) {
    x <- suppressWarnings(as.numeric(x))
    floor <- suppressWarnings(as.numeric(floor))
    positive <- x[is.finite(x) & x > 0]
    if (!length(positive))
        return(floor)
    span <- c(floor, max(positive, na.rm = TRUE))
    span <- span[is.finite(span) & span > 0]
    if (length(span) < 2L)
        span <- c(min(span), min(span) * 10)
    exponents <- pretty(log10(span), n = n)
    out <- 10^exponents
    out[is.finite(out) & out > 0]
}

cohort_strip_layers <-
function(horizon_day, y, height, text_size = 3) {
    horizon_day <- suppressWarnings(as.numeric(horizon_day))
    if (!is.finite(horizon_day) || horizon_day <= 0)
        horizon_day <- 100
    strip_width <- max(horizon_day * 0.018, 1)
    labels <- c("2N", "4N")
    y_raw <- suppressWarnings(as.numeric(y))
    if (!is.null(names(y))) {
        y_use <- y_raw[match(labels, names(y))]
    }
    else {
        y_use <- rep(y_raw, length.out = length(labels))
    }
    y_use[!is.finite(y_use)] <- 1
    strip_df <- data.frame(cohort = factor(labels, levels = labels), x = horizon_day + strip_width/2, y = y_use, label = labels,
        stringsAsFactors = FALSE)
    list(geom_tile(data = strip_df[strip_df$label == "2N", , drop = FALSE], aes(x = x, y = y), inherit.aes = FALSE, width = strip_width,
        height = height, fill = "#9ecae1", color = "grey70", linewidth = 0.3), geom_tile(data = strip_df[strip_df$label ==
        "4N", , drop = FALSE], aes(x = x, y = y), inherit.aes = FALSE, width = strip_width, height = height, fill = "lightpink",
        color = "grey70", linewidth = 0.3), geom_text(data = strip_df, aes(x = x, y = y, label = label), inherit.aes = FALSE,
        angle = 270, size = text_size, color = "black"))
}

cohort_facet_grid <-
function() {
    if (!requireNamespace("ggh4x", quietly = TRUE)) {
        stop("Package 'ggh4x' is required for colored 2N/4N facet strips in Predicted plots.", call. = FALSE)
    }
    ggh4x::facet_grid2(rows = vars(cohort), strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = c("#9ecae1",
        "lightpink"), colour = "grey70"), text_y = ggh4x::elem_list_text(colour = "black", size = 8)))
}

make_predict_annotation_track_plot <-
function(df, value_col, y_label, legend_title, day_width, horizon_day, x_breaks, colors, transform = "identity", limits = NULL,
    breaks = NULL, labels = NULL, title = NULL, subtitle = NULL, show_legend = TRUE) {
    plot_df <- data.frame(cohort = factor(as.character(df$cohort), levels = c("2N", "4N")), day = as.numeric(df$day), value = suppressWarnings(as.numeric(df[[value_col]])),
        stringsAsFactors = FALSE)
    plot_df <- plot_df[is.finite(plot_df$day), , drop = FALSE]
    plot_df$cohort_y <- ifelse(as.character(plot_df$cohort) == "2N", 2, 1)
    if (identical(transform, "log10")) {
        floor <- log10_plot_floor(plot_df$value, default = 1e-12)
        original_breaks <- log10_original_breaks(plot_df$value, floor = floor, n = 4)
        if (length(original_breaks) > 2L) {
            original_breaks <- original_breaks[c(1L, length(original_breaks))]
        }
        plot_df$value_fill <- log10(floor_for_log10_plot(plot_df$value, floor))
        fill_scale <- scale_fill_gradientn(colors = colors, breaks = log10(original_breaks), labels = format_original_scale_labels(original_breaks),
            na.value = "grey90", name = legend_title, guide = guide_colorbar(barheight = grid::unit(0.2, "in"), barwidth = grid::unit(0.12,
                "in")))
    }
    else {
        plot_df$value_fill <- plot_df$value
        fill_scale <- scale_fill_gradientn(colors = colors, limits = limits, breaks = breaks, labels = labels, oob = scales::squish,
            na.value = "grey90", name = legend_title, guide = guide_colorbar(barheight = grid::unit(0.2, "in"), barwidth = grid::unit(0.12,
                "in")))
    }
    ggplot(plot_df, aes(x = day, y = cohort_y, fill = value_fill)) + geom_tile(width = day_width, height = 0.9) + cohort_strip_layers(horizon_day = horizon_day,
        y = c(`2N` = 2, `4N` = 1), height = 0.9, text_size = 2.8) + scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
        scale_y_continuous(breaks = NULL, limits = c(0.5, 2.5), expand = c(0, 0)) + fill_scale + coord_cartesian(xlim = c(0,
        horizon_day), expand = FALSE, clip = "off") + labs(title = title, subtitle = subtitle, x = NULL, y = y_label) + theme_bw(base_size = 10) +
        theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8), legend.position = if (isTRUE(show_legend))
                "right"
            else "none", legend.title = element_text(size = 6.5), legend.text = element_text(size = 5.5), plot.title = element_text(size = 12),
            plot.subtitle = element_text(size = 8), plot.margin = margin(1, 10, 1, 1))
}

make_predicted_annotation_legend <-
function(annotation_summary, o2_plot_min, o2_plot_max) {
    format_log10_labels <- function(x) {
        vapply(as.numeric(x), function(z) {
            if (!is.finite(z))
                return("")
            if (abs(z - round(z)) < 1e-08)
                return(as.character(as.integer(round(z))))
            format(signif(z, 3), trim = TRUE, scientific = FALSE)
        }, character(1))
    }
    log10_legend_spec <- function(x, default = 1e-12, n = 5) {
        floor <- log10_plot_floor(x, default = default)
        values <- log10(floor_for_log10_plot(x, floor))
        values <- values[is.finite(values)]
        if (!length(values))
            values <- log10(c(floor, floor * 10))
        range_use <- range(values, na.rm = TRUE)
        if (!all(is.finite(range_use)) || diff(range_use) <= 0) {
            range_use <- range_use[1] + c(0, 1)
        }
        breaks <- pretty(range_use, n = n)
        breaks <- breaks[is.finite(breaks) & breaks >= range_use[1] & breaks <= range_use[2]]
        if (length(breaks) < 2L)
            breaks <- range_use
        list(range = range_use, breaks = breaks, labels = format_log10_labels(breaks))
    }
    linear_legend_spec <- function(range_use, n = 5) {
        range_use <- suppressWarnings(as.numeric(range_use))
        range_use <- range_use[is.finite(range_use)]
        if (length(range_use) < 2L)
            range_use <- c(0, 1)
        range_use <- range(range_use)
        if (diff(range_use) <= 0)
            range_use <- range_use[1] + c(0, 1)
        breaks <- pretty(range_use, n = n)
        breaks <- breaks[is.finite(breaks) & breaks >= range_use[1] & breaks <= range_use[2]]
        if (length(breaks) < 2L)
            breaks <- range_use
        list(range = range_use, breaks = breaks, labels = format(signif(breaks, 3), trim = TRUE))
    }
    make_bar_df <- function(x0, x1, colors, n = 80L) {
        data.frame(x = seq(x0, x1, length.out = n), y = 0.42, fill_col = (grDevices::colorRampPalette(colors))(n), stringsAsFactors = FALSE)
    }
    make_tick_df <- function(x0, x1, spec) {
        rng <- spec$range
        pos <- if (diff(rng) > 0) {
            x0 + (spec$breaks - rng[1])/diff(rng) * (x1 - x0)
        }
        else {
            rep(mean(c(x0, x1)), length(spec$breaks))
        }
        data.frame(x = pos, label = spec$labels, stringsAsFactors = FALSE)
    }
    burden_spec <- log10_legend_spec(annotation_summary$burden, default = 1e-12, n = 5)
    live_spec <- log10_legend_spec(annotation_summary$live_cells, default = 1, n = 5)
    o2_spec <- linear_legend_spec(c(o2_plot_min, o2_plot_max), n = 5)
    burden_ticks <- make_tick_df(0.02, 0.28, burden_spec)
    live_ticks <- make_tick_df(0.37, 0.63, live_spec)
    o2_ticks <- make_tick_df(0.72, 0.98, o2_spec)
    tick_df <- dplyr::bind_rows(dplyr::mutate(burden_ticks, track = "burden"), dplyr::mutate(live_ticks, track = "live"),
        dplyr::mutate(o2_ticks, track = "o2"))
    bars <- dplyr::bind_rows(make_bar_df(0.02, 0.28, c("#ffffbf", "#542788")), make_bar_df(0.37, 0.63, c("#ffffbf", "#2c7bb6")),
        make_bar_df(0.72, 0.98, c("#f7f7f7", "#9ecae1", "#08519c")))
    ggplot(bars, aes(x = x, y = y, fill = fill_col)) + geom_tile(width = 0.004, height = 0.22) + scale_fill_identity() +
        geom_segment(data = tick_df, aes(x = x, xend = x, y = 0.27, yend = 0.33), inherit.aes = FALSE, linewidth = 0.18,
            color = "grey20") + geom_text(data = tick_df, aes(x = x, y = 0.12, label = label), inherit.aes = FALSE, size = 1.9,
        color = "black") + annotate("text", x = 0.02, y = 0.82, label = "log10 Burden (mm^3)", hjust = 0, size = 2.6) + annotate("text",
        x = 0.37, y = 0.82, label = "log10 Viable cells", hjust = 0, size = 2.6) + annotate("text", x = 0.72, y = 0.82, label = "Effective O2 (%)",
        hjust = 0, size = 2.6) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void() + theme(plot.margin = margin(0,
        0, 0, 0))
}

plot_missegregation_probability_timecourse <-
function(ms_timecourse, out_dir, fit_dir = NULL, report_dt = NULL) {
    if (is.null(ms_timecourse) || !is.data.frame(ms_timecourse) || !nrow(ms_timecourse)) {
        return(FALSE)
    }
    plot_df <- ms_timecourse %>% mutate(day = as.numeric(day), mean_p_misseg = as.numeric(mean_p_misseg), cohort = as.character(cohort),
        sample_id = as.character(sample_id)) %>% filter(is.finite(day), is.finite(mean_p_misseg))
    if (!nrow(plot_df))
        return(FALSE)
    subtitle_parts <- character()
    if (!is.null(fit_dir))
        subtitle_parts <- c(subtitle_parts, paste0("fit_dir=", basename(fit_dir)))
    if (!is.null(report_dt) && is.finite(as.numeric(report_dt))) {
        subtitle_parts <- c(subtitle_parts, paste0("report_dt=", signif(as.numeric(report_dt), 4)))
    }
    subtitle <- if (length(subtitle_parts))
        paste(subtitle_parts, collapse = " | ")
    else NULL
    p <- ggplot(plot_df, aes(x = day, y = mean_p_misseg, color = cohort, group = sample_id)) + geom_line(linewidth = 0.8,
        alpha = 0.9) + facet_wrap(~harvest, ncol = 2) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE,
        digits = 3)) + labs(title = "Resource-Stress Model: Mean Per-Chromosome Missegregation Probability Over Time", subtitle = subtitle,
        x = "Day", y = "Viable-population-weighted mean per-chromosome missegregation probability", color = "Cohort") + theme_bw(base_size = 11)
    pdf_path <- file.path(out_dir, "missegregation_probability_over_time.pdf")
    png_path <- file.path(out_dir, "missegregation_probability_over_time.png")
    ggsave(pdf_path, p, width = 13, height = 9)
    ggsave(png_path, p, width = 13, height = 9, dpi = 180, bg = "white")
    TRUE
}

plot_terminal_ploidy_violin_compare <-
function(compare_df, fit_dir, out_dir) {
    if (nrow(compare_df) == 0L)
        return(NULL)
    endpoint_label <- unique(compare_df$endpoint_mode)
    endpoint_label <- if (length(endpoint_label) == 1L && identical(endpoint_label[[1]], "chr_number")) {
        "Chromosome Number (N)"
    }
    else {
        "Ploidy (2N scale)"
    }
    weighted_quantile_local <- function(x, w, probs) {
        x <- as.numeric(x)
        w <- as.numeric(w)
        keep <- is.finite(x) & is.finite(w) & (w > 0)
        x <- x[keep]
        w <- w[keep]
        if (!length(x))
            return(rep(NA_real_, length(probs)))
        ord <- order(x)
        x <- x[ord]
        w <- w[ord]
        cw <- cumsum(w)/sum(w)
        vapply(probs, function(p) {
            idx <- which(cw >= p)[1]
            if (!length(idx) || is.na(idx))
                x[[length(x)]]
            else x[[idx]]
        }, numeric(1))
    }
    box_df <- compare_df %>% group_by(cohort, source, fill_group) %>% summarise(q1 = weighted_quantile_local(endpoint_value,
        weight, 0.25), median = weighted_quantile_local(endpoint_value, weight, 0.5), q3 = weighted_quantile_local(endpoint_value,
        weight, 0.75), ymin_raw = min(endpoint_value[weight > 0], na.rm = TRUE), ymax_raw = max(endpoint_value[weight > 0],
        na.rm = TRUE), .groups = "drop") %>% mutate(iqr = q3 - q1, ymin = pmax(ymin_raw, q1 - 1.5 * iqr), ymax = pmin(ymax_raw,
        q3 + 1.5 * iqr))
    fill_values <- c(`2N Observed` = "#8FBF7A", `2N Predicted` = "#D7EFCC", `4N Observed` = "#F2A9BC", `4N Predicted` = "#F9D9E3")
    p <- ggplot(compare_df, aes(x = source, y = endpoint_value, weight = weight, fill = fill_group)) + geom_violin(trim = FALSE,
        scale = "width", quantiles = NULL, color = "grey35", linewidth = 0.35, alpha = 0.95) + geom_boxplot(data = box_df,
        aes(x = source, ymin = ymin, lower = q1, middle = median, upper = q3, ymax = ymax), stat = "identity", inherit.aes = FALSE,
        width = 0.18, fill = "white", color = "black", linewidth = 0.45, outlier.shape = NA) + facet_wrap(~cohort, nrow = 1,
        ncol = 2) + scale_fill_manual(values = fill_values, drop = FALSE) + labs(title = paste0("Observed vs Predicted ",
        endpoint_label, " Distributions Used in Endpoint Objective"), subtitle = if (identical(endpoint_label, "Chromosome Number (N)")) {
        paste0("fit_dir=", basename(fit_dir))
    }
    else {
        paste0("Observed ploidy is mapped to the chromosome-number grid used by the objective | fit_dir=", basename(fit_dir))
    }, x = NULL, y = endpoint_label, fill = NULL) + theme_bw(base_size = 11) + theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "grey80"))
    ggsave(file.path(out_dir, "terminal_ploidy_observed_vs_predicted_violin.pdf"), p, width = 6.5, height = 6.5)
    p
}

extract_ggplot_legend_grob <-
function(plot) {
    if (!inherits(plot, "ggplot"))
        return(NULL)
    plot_grob <- ggplot2::ggplotGrob(plot + theme(legend.position = "right"))
    guide_idx <- which(vapply(plot_grob$grobs, function(g) g$name %||% "", character(1)) == "guide-box")
    if (!length(guide_idx))
        return(NULL)
    plot_grob$grobs[[guide_idx[[1]]]]
}

save_aligned_plot_row_with_shared_legend <-
function(plots, path, width = 18, height = 6.5) {
    plots <- Filter(function(p) inherits(p, "ggplot"), plots)
    if (!length(plots))
        return(invisible(NULL))
    legend_grob <- extract_ggplot_legend_grob(plots[[1]])
    plot_grobs <- lapply(plots, function(p) ggplot2::ggplotGrob(p + theme(legend.position = "none")))
    height_lengths <- vapply(plot_grobs, function(g) length(g$heights), integer(1))
    if (length(unique(height_lengths)) == 1L) {
        max_heights <- do.call(grid::unit.pmax, lapply(plot_grobs, function(g) g$heights))
        plot_grobs <- lapply(plot_grobs, function(g) {
            g$heights <- max_heights
            g
        })
    }
    grDevices::pdf(path, width = width, height = height, onefile = TRUE)
    grid::grid.newpage()
    if (is.null(legend_grob)) {
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = length(plot_grobs))))
        for (i in seq_along(plot_grobs)) {
            grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
            grid::grid.draw(plot_grobs[[i]])
            grid::upViewport()
        }
        grid::upViewport()
    }
    else {
        lay <- grid::grid.layout(nrow = 1, ncol = length(plot_grobs) + 1L, widths = grid::unit.c(rep(grid::unit(1, "null"),
            length(plot_grobs)), grid::grobWidth(legend_grob) + grid::unit(0.25, "in")))
        grid::pushViewport(grid::viewport(layout = lay))
        for (i in seq_along(plot_grobs)) {
            grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = i))
            grid::grid.draw(plot_grobs[[i]])
            grid::upViewport()
        }
        grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = length(plot_grobs) + 1L))
        grid::grid.draw(legend_grob)
        grid::upViewport()
        grid::upViewport()
    }
    grDevices::dev.off()
    invisible(path)
}

plot_predict_burden_live_dead_decomposition_combined <-
function(predict_results, out_dir, death_language) {
    predict_results <- Filter(function(x) {
        is.list(x) && is.data.frame(x$burden_decomp_predict) && is.data.frame(x$burden_decomp_predict_long) && is.finite(as.numeric(x$horizon_day))
    }, predict_results)
    if (!length(predict_results))
        return(invisible(NULL))
    predict_results <- predict_results[order(vapply(predict_results, function(x) as.numeric(x$horizon_day), numeric(1)))]
    fill_values <- stats::setNames(c("#1f77b4", "#d62728", "#2ca02c"), c(death_language$live_label, death_language$component_label,
        death_language$cin_component_label))
    plots <- lapply(predict_results, function(res) {
        horizon_day <- as.numeric(res$horizon_day)
        decomp_floor <- log10_plot_floor(c(res$burden_decomp_predict$burden_live, res$burden_decomp_predict$burden_dead_hypoxia,
            res$burden_decomp_predict$burden_dead_buffer, res$burden_decomp_predict$burden_total), default = 1e-12)
        burden_decomp_predict <- res$burden_decomp_predict %>% mutate(burden_total_log_plot = floor_for_log10_plot(burden_total,
            decomp_floor))
        burden_decomp_predict_ribbon <- make_burden_decomp_ribbon(burden_decomp_predict, death_language, decomp_floor)
        ggplot(burden_decomp_predict_ribbon, aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)) +
            geom_ribbon(alpha = 0.55) + geom_line(data = burden_decomp_predict, aes(x = day, y = burden_total_log_plot),
            inherit.aes = FALSE, color = "black", linewidth = 0.65) + facet_wrap(~cohort, ncol = 1, scales = "free_y") +
            scale_fill_manual(values = fill_values, drop = FALSE) + log10_burden_y_scale() + coord_cartesian(xlim = c(0,
            horizon_day)) + labs(title = paste0("0-", as.integer(round(horizon_day)), " days"), subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
            x = "Day", y = "log10 Tumor burden (mm^3)", fill = "Component") + theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
            legend.position = "right")
    })
    save_aligned_plot_row_with_shared_legend(plots, file.path(out_dir, "predict_burden_live_dead_decomposition_combined.pdf"),
        width = 18, height = 6.5)
}

read_invivo_simulation_table <- function(simulation_dir, filename) {
    path <- file.path(simulation_dir, filename)
    if (!file.exists(path)) stop("Missing required simulation table: ", path)
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

plot_functional_response_curves <- function(simulation_dir, cfg, out_dir) {
    o2_curve <- read_invivo_simulation_table(simulation_dir, "functional_curve_oxygen.tsv")
    o2_curve_multi <- read_invivo_simulation_table(simulation_dir, "functional_curve_oxygen_multi_ploidy.tsv")
    viability_curve <- read_invivo_simulation_table(simulation_dir, "functional_curve_ploidy.tsv")
    read_invivo_simulation_table(simulation_dir, "functional_curve_ploidy_by_o2.tsv")

    state_axis_label <- functional_state_axis_label(cfg)
    ref_states <- unique(o2_curve_multi[, c("cohort", "ploidy_multiple", "N_ref"), drop = FALSE])
    ref_states <- ref_states[order(ref_states$ploidy_multiple), , drop = FALSE]
    ref_state_subtitle <- paste0("Reference states: ", paste(ref_states$cohort, collapse = ", "))
    o2_curve_multi$cohort <- factor(o2_curve_multi$cohort, levels = unique(o2_curve_multi$cohort))

    save_o2_curve_plot <- function(value_col, title, y_label, filename) {
        plot_df <- o2_curve_multi
        plot_df$value <- suppressWarnings(as.numeric(plot_df[[value_col]]))
        plot_df <- plot_df[is.finite(plot_df$oxygen_pct) & is.finite(plot_df$value), , drop = FALSE]
        if (!nrow(plot_df)) return(invisible(FALSE))
        p <- ggplot(plot_df, aes(x = oxygen_pct, y = value, color = cohort)) +
            geom_line(linewidth = 1) +
            labs(
                title = title,
                subtitle = "In vivo fitted rate function across oxygen levels",
                x = "Effective oxygen (%)",
                y = y_label,
                color = "Reference state"
            ) +
            theme_bw(base_size = 11)
        ggsave(file.path(out_dir, filename), p, width = 10, height = 7)
        invisible(TRUE)
    }
    save_o2_curve_plot("ms_rate", "Effective Oxygen vs Missegregation Rate", "Missegregation rate",
                       "oxygen_vs_missegregation_rate_multi_ploidy.pdf")
    save_o2_curve_plot("proliferation_rate", "Effective Oxygen vs Proliferation Rate", "Proliferation rate",
                       "oxygen_vs_proliferation_rate.pdf")
    save_o2_curve_plot("death_rate", "Effective Oxygen vs Stress-Associated Death Rate", "Stress-associated death rate",
                       "oxygen_vs_death_rate.pdf")

    multi_colors <- stats::setNames(grDevices::hcl.colors(nrow(ref_states), palette = "Dark 3"), ref_states$cohort)
    p_death_msr <- ggplot(o2_curve, aes(x = death_rate, y = ms_rate, color = cohort, group = cohort)) +
        geom_path(linewidth = 1) +
        scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
        labs(
            title = "Stress-Associated Death Rate vs Missegregation Rate",
            subtitle = "Same effective-oxygen sweep and reference chromosome-number states as Effective Oxygen vs Missegregation Rate",
            x = "Stress-associated death rate", y = "Missegregation rate", color = "Cohort"
        ) +
        theme_bw(base_size = 11)
    p_msr_nonviable_daughter_fraction <- ggplot(
        o2_curve_multi,
        aes(x = ms_rate, y = misseg_nonviable_daughter_fraction, color = factor(cohort, levels = ref_states$cohort))
    ) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = multi_colors, drop = FALSE) +
        labs(
            title = "Nonviable Daughter Fraction vs Missegregation Rate Across Reference Chromosome-Number States",
            subtitle = paste0(
                "Missegregation-linked nonviable daughters / all daughters per division; excludes boundary-drop losses | ",
                ref_state_subtitle
            ),
            x = "Missegregation rate", y = "Nonviable daughters / all daughters", color = "Reference state"
        ) +
        theme_bw(base_size = 11)
    p_viability <- ggplot(viability_curve, aes(x = endpoint_value, y = viability_after_ms)) +
        geom_line(color = "#2ca02c", linewidth = 1) +
        labs(
            title = paste0(state_axis_label, " vs Post-Missegregation Survival"),
            subtitle = "Ploidy-dependent survival for a one-copy missegregation event",
            x = state_axis_label, y = "Post-missegregation survival"
        ) +
        theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "death_rate_vs_missegregation_rate.pdf"), p_death_msr, width = 10, height = 7)
    ggsave(file.path(out_dir, "ms_rate_vs_nonviable_daughter_fraction.pdf"), p_msr_nonviable_daughter_fraction, width = 10, height = 7)
    ggsave(file.path(out_dir, "ploidy_vs_viability_after_ms.pdf"), p_viability, width = 10, height = 7)
    invisible(list(p_death_msr = p_death_msr,
                   p_msr_nonviable_daughter_fraction = p_msr_nonviable_daughter_fraction,
                   p_viability = p_viability))
}

plot_prediction_horizon_population_average_cin <- function(cin_table, out_dir, target_days) {
    if (!is.data.frame(cin_table) || !nrow(cin_table)) return(FALSE)
    cohort_levels <- c("2N", "4N")
    cohort_labels <- c("2N" = "2N-derived", "4N" = "4N-derived")
    cohort_colors <- c("2N" = "#1f77b4", "4N" = "#d62728")
    cohort_linetypes <- c("2N" = "solid", "4N" = "solid")
    any_saved <- FALSE
    for (target_day in sort(unique(as.numeric(target_days)))) {
        plot_df <- cin_table %>%
            filter(abs(as.numeric(.data$target_day) - target_day) < 1e-8) %>%
            mutate(initial_cohort = factor(as.character(initial_cohort), levels = cohort_levels)) %>%
            arrange(cohort_order, day)
        if (!nrow(plot_df)) next
        target_label <- format(target_day, trim = TRUE, scientific = FALSE)
        target_tag <- gsub("\\.", "p", target_label)
        p <- ggplot(plot_df, aes(x = day, y = population_average_cin, color = initial_cohort,
                                linetype = initial_cohort, group = initial_cohort)) +
            geom_line(linewidth = 0.9, alpha = 0.95) +
            scale_color_manual(values = cohort_colors, labels = cohort_labels, drop = FALSE) +
            scale_linetype_manual(values = cohort_linetypes, labels = cohort_labels, drop = FALSE) +
            scale_x_continuous(limits = c(0, target_day), breaks = pretty(c(0, target_day), n = 5)) +
            scale_y_continuous(labels = function(x) formatC(x, format = "f", digits = 4)) +
            labs(
                title = paste0("0-", target_label, " Day Population-average CIN rate over time"),
                subtitle = "Canonical 2N-derived and 4N-derived trajectories simulated from the fitted parameters.",
                x = "Day",
                y = "Population-average CIN rate\n(mean per-chromosome missegregation probability)",
                color = "Initial cohort", linetype = "Initial cohort"
            ) +
            theme_bw(base_size = 11)
        out_base <- paste0("population_average_cin_by_initial_cohort_day", target_tag)
        ggsave(file.path(out_dir, paste0(out_base, ".pdf")), p, width = 10, height = 7)
        ggsave(file.path(out_dir, paste0(out_base, ".png")), p, width = 10, height = 7, dpi = 180, bg = "white")
        any_saved <- TRUE
    }
    any_saved
}

plot_predict_horizon <- function(simulation_dir, cfg, out_dir, horizon_day, report_dt = 1.0) {
  save_aligned_plot_stack <- function(plots, path, width = 14, height = 7, row_heights = NULL) {
    plots <- Filter(function(p) inherits(p, "ggplot"), plots)
    if (!length(plots)) return(invisible(NULL))
    if (is.null(row_heights)) {
      row_heights <- rep(1, length(plots))
    }
    row_heights <- as.numeric(row_heights)
    if (length(row_heights) != length(plots) || any(!is.finite(row_heights)) || any(row_heights <= 0)) {
      row_heights <- rep(1, length(plots))
    }

    if (requireNamespace("cowplot", quietly = TRUE)) {
      aligned <- cowplot::align_plots(plotlist = plots, align = "v", axis = "lr")
      combined <- cowplot::plot_grid(plotlist = aligned, ncol = 1, rel_heights = row_heights)
      ggsave(path, combined, width = width, height = height, device = grDevices::pdf)
      return(invisible(path))
    }

    grobs <- lapply(plots, ggplot2::ggplotGrob)
    width_lengths <- vapply(grobs, function(g) length(g$widths), integer(1))
    if (length(unique(width_lengths)) == 1L) {
      max_widths <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
      grobs <- lapply(grobs, function(g) {
        g$widths <- max_widths
        g
      })
    }

    grDevices::pdf(path, width = width, height = height, onefile = TRUE)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = length(grobs),
      ncol = 1,
      heights = grid::unit(row_heights, "null")
    )))
    for (i in seq_along(grobs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
      grid::grid.draw(grobs[[i]])
      grid::upViewport()
    }
    grid::upViewport()
    grDevices::dev.off()
    invisible(path)
  }
  horizon_tag <- paste0("0_", as.integer(round(horizon_day)), "day")
  .fit_label <- cfg$fit_label %||% basename(cfg$fit_dir %||% simulation_dir)
  burden_all <- read_invivo_simulation_table(simulation_dir, paste0("predict_burden_", horizon_tag, ".tsv"))
  ploidy_all <- read_invivo_simulation_table(simulation_dir, paste0("predict_ploidy_", horizon_tag, ".tsv"))
  ploidy_mean <- read_invivo_simulation_table(simulation_dir, paste0("predict_ploidy_weighted_mean_", horizon_tag, ".tsv"))
  chr_density_df <- read_invivo_simulation_table(simulation_dir, paste0("predict_chromosome_density_", horizon_tag, ".tsv"))
  death_ratio_predict <- read_invivo_simulation_table(simulation_dir, paste0("predict_death_ratio_", horizon_tag, ".tsv"))
  resource_death_fraction_predict <- read_invivo_simulation_table(simulation_dir, paste0("predict_resource_death_fraction_", horizon_tag, ".tsv"))
  burden_decomp_predict <- death_ratio_predict[, c("cohort", "harvest", "dose", "day", "burden_live", "burden_dead_hypoxia", "burden_dead_buffer", "burden_dead_total", "burden_total"), drop = FALSE]
  burden_decomp_predict$cohort <- factor(as.character(burden_decomp_predict$cohort), levels = c("2N", "4N"))
  death_language <- resource_death_language()
  cin_table <- read_invivo_simulation_table(simulation_dir, "population_average_cin_by_initial_cohort_horizons.tsv")
  plot_prediction_horizon_population_average_cin(cin_table, out_dir, target_days = horizon_day)
  horizon_tag <- paste0("0_", as.integer(round(horizon_day)), "day")
  # Remove deprecated multi-file prediction plot outputs for this horizon to avoid stale files.
  unlink(file.path(out_dir, c(
    paste0("predict_burden_normalized_", horizon_tag, ".pdf"),
    paste0("predict_burden_absolute_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("predict_ploidy_weighted_mean_", horizon_tag, ".pdf"),
    paste0("forecast_burden_normalized_", horizon_tag, ".pdf"),
    paste0("forecast_burden_absolute_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_heatmap_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_top_states_", horizon_tag, ".pdf"),
    paste0("forecast_ploidy_weighted_mean_", horizon_tag, ".pdf"),
    paste0("predict_g_timecourse_", horizon_tag, ".pdf"),
    paste0("predict_live_resource_death_fraction_", horizon_tag, ".pdf"),
    paste0("predict_live_weighted_pms_", horizon_tag, ".pdf"),
    paste0("predict_death_ratio_", horizon_tag, ".pdf"),
    paste0("predict_o2_timecourse_", horizon_tag, ".pdf"),
    paste0("predicted_", horizon_tag, ".pdf")
  )), force = TRUE)
  burden_plot_df <- burden_all %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value = as.numeric(pred_burden),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
    )
  burden_plot_floor <- log10_plot_floor(burden_plot_df$value, default = 1e-12)
  burden_plot_df <- burden_plot_df %>%
    mutate(value_log_plot = floor_for_log10_plot(value, burden_plot_floor))

  endpoint_plot_df <- ploidy_mean %>%
    transmute(
      harvest = as.character(harvest),
      cohort = as.character(cohort),
      dose = as.numeric(dose),
      day = as.numeric(day),
      value_chr = as.numeric(weighted_mean_N),
      sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
    )

  predict_x_breaks <- unique(as.numeric(seq(0, as.numeric(horizon_day), length.out = 5)))
  predict_x_scale <- function() {
    scale_x_continuous(
      breaks = predict_x_breaks,
      expand = c(0, 0)
    )
  }
  predict_curve_theme <- theme(
    axis.title = element_text(size = 8),
    axis.title.y.right = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.text.y.right = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 8)
  )
  ploidy_n_range <- range(endpoint_plot_df$value_chr, na.rm = TRUE)
  if (!all(is.finite(ploidy_n_range))) {
    ploidy_n_range <- c(as.numeric(cfg$N_MIN), as.numeric(cfg$N_MAX))
  }
  ploidy_n_pad <- diff(ploidy_n_range) * 0.04
  if (!is.finite(ploidy_n_pad) || ploidy_n_pad <= 0) ploidy_n_pad <- 1
  ploidy_n_limits <- c(
    floor(max(as.numeric(cfg$N_MIN), ploidy_n_range[[1]] - ploidy_n_pad)),
    ceiling(min(as.numeric(cfg$N_MAX), ploidy_n_range[[2]] + ploidy_n_pad))
  )
  if (!all(is.finite(ploidy_n_limits)) || ploidy_n_limits[[2]] <= ploidy_n_limits[[1]]) {
    ploidy_n_limits <- c(as.numeric(cfg$N_MIN), as.numeric(cfg$N_MAX))
  }
  ploidy_n_breaks <- pretty(ploidy_n_limits, n = 4)
  ploidy_n_breaks <- ploidy_n_breaks[ploidy_n_breaks >= ploidy_n_limits[[1]] & ploidy_n_breaks <= ploidy_n_limits[[2]]]

  o2_plot_min <- 0
  o2_plot_max <- 5
  o2_plot_df <- if ("pred_o2_pct" %in% names(burden_all)) {
    burden_all %>%
      filter(is.finite(pred_o2_pct)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        value = as.numeric(clip(pred_o2_pct, o2_plot_min, o2_plot_max)),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
      )
  } else {
    data.frame()
  }

  live_cells_plot_df <- if ("pred_burden_live_cells" %in% names(burden_all)) {
    burden_all %>%
      filter(is.finite(pred_burden_live_cells)) %>%
      transmute(
        harvest = as.character(harvest),
        cohort = as.character(cohort),
        dose = as.numeric(dose),
        day = as.numeric(day),
        value = as.numeric(pred_burden_live_cells),
        sample_id = paste(harvest, cohort, format(dose, trim = TRUE, scientific = FALSE), sep = "__")
      )
  } else {
    data.frame()
  }
  live_cells_plot_floor <- if (nrow(live_cells_plot_df) > 0L) {
    log10_plot_floor(live_cells_plot_df$value, default = 1)
  } else {
    1
  }
  if (nrow(live_cells_plot_df) > 0L) {
    live_cells_plot_df <- live_cells_plot_df %>%
      mutate(value_log_plot = floor_for_log10_plot(value, live_cells_plot_floor))
  }

  chr_density_day_width <- {
    day_vals <- sort(unique(as.numeric(ploidy_all$day[is.finite(ploidy_all$day)])))
    day_step <- diff(day_vals)
    day_step <- day_step[is.finite(day_step) & day_step > 0]
    width_use <- if (length(day_step) > 0L) stats::median(day_step, na.rm = TRUE) else as.numeric(report_dt)
    if (!is.finite(width_use) || width_use <= 0) width_use <- 1
    width_use
  }
  chr_density_n_min <- as.integer(cfg$N_MIN)
  chr_density_n_max <- as.integer(cfg$N_MAX)
  chr_density_bin_width <- 5L
  p_predict_burden <- ggplot(
    burden_plot_df,
    aes(x = day, y = value_log_plot, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
    scale_y_log10() +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      title = paste0("Predict Curves: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = paste0("Single summary plot (all scenarios overlaid) | fit_dir=", .fit_label, " | report_dt=", report_dt),
      x = "Day",
      y = "Burden (log10 scale)",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  p_predict_endpoint <- ggplot(
    endpoint_plot_df,
    aes(x = day, y = value_chr, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
    scale_y_continuous(
      name = "Mean chr. number",
      breaks = ploidy_n_breaks,
      sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Mean ploidy")
    ) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      x = "Day",
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  p_predict_chr_density <- if (nrow(chr_density_df) > 0L) {
    chr_density_fill_max <- max(chr_density_df$density, na.rm = TRUE)
    if (!is.finite(chr_density_fill_max) || chr_density_fill_max <= 0) {
      chr_density_fill_max <- 1
    }
    ggplot(
      chr_density_df,
      aes(x = day, y = chr_bin_mid, fill = density)
    ) +
      geom_tile(width = chr_density_day_width, height = chr_density_bin_width) +
      facet_grid(cohort ~ .) +
      predict_x_scale() +
      ploidy_fraction_fill_scale(chr_density_fill_max, name = "Probability\ndensity") +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      labs(
        x = "Day",
        y = "Chromosome count (N)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "grey80")
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      labs(
        x = "Day",
        y = "Chromosome count (N)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  }

  p_predict_o2 <- if (nrow(o2_plot_df) > 0L) {
    ggplot(
      o2_plot_df,
      aes(x = day, y = value, group = sample_id, color = cohort)
    ) +
      geom_line(linewidth = 0.65, alpha = 0.8) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_continuous(
        limits = c(o2_plot_min, o2_plot_max),
        breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = "Day",
        y = "Effective O2 (%)",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_continuous(
        limits = c(o2_plot_min, o2_plot_max),
        breaks = seq(o2_plot_min, o2_plot_max, by = 1),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      labs(
        x = "Day",
        y = "Effective O2 (%)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank()
      ) +
      predict_curve_theme
  }

  p_predict_live_cells <- if (nrow(live_cells_plot_df) > 0L) {
    ggplot(
      live_cells_plot_df,
      aes(x = day, y = value_log_plot, group = sample_id, color = cohort)
    ) +
      geom_line(linewidth = 0.65, alpha = 0.8) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_log10() +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = NULL,
        y = "Viable cells (log10 scale)",
        color = "Cohort"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()
      ) +
      predict_curve_theme
  } else {
    ggplot() +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), expand = FALSE) +
      scale_y_log10() +
      labs(
        x = NULL,
        y = "Viable cells (log10 scale)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()
      ) +
      predict_curve_theme
  }
  p_resource_death_fraction_predict <- ggplot(
    resource_death_fraction_predict,
    aes(x = day, y = resource_death_fraction, group = sample_id, color = cohort)
  ) +
    geom_line(linewidth = 0.65, alpha = 0.8) +
    predict_x_scale() +
    coord_cartesian(xlim = c(0, horizon_day), ylim = c(0, 1), expand = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    labs(
      x = "Day",
      y = paste0(death_language$adjective, " death / all deaths"),
      color = "Cohort"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    predict_curve_theme

  if (identical(as.integer(round(horizon_day)), 1000L)) {
    annotation_summary <- burden_all %>%
      filter(as.character(cohort) %in% c("2N", "4N"), is.finite(day)) %>%
      transmute(
        cohort = factor(as.character(cohort), levels = c("2N", "4N")),
        day = as.numeric(day),
        burden = as.numeric(pred_burden),
        live_cells = as.numeric(.first_non_null_local(pred_burden_live_cells, NA_real_)),
        o2_pct = as.numeric(clip(.first_non_null_local(pred_o2_pct, NA_real_), o2_plot_min, o2_plot_max))
      ) %>%
      group_by(cohort, day) %>%
      summarise(
        burden = mean(burden, na.rm = TRUE),
        live_cells = mean(live_cells, na.rm = TRUE),
        o2_pct = mean(o2_pct, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        burden = ifelse(is.nan(burden), NA_real_, burden),
        live_cells = ifelse(is.nan(live_cells), NA_real_, live_cells),
        o2_pct = ifelse(is.nan(o2_pct), NA_real_, o2_pct)
      )

    p_annotation_title <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0.95,
        hjust = 0,
        vjust = 1,
        label = paste0("Predicted (0-", as.integer(round(horizon_day)), " day)"),
        size = 3.8
      ) +
      annotate(
        "text",
        x = 0,
        y = 0.25,
        hjust = 0,
        vjust = 1,
        label = paste0("Column annotations are cohort-level means; fit_dir=", .fit_label, " | report_dt=", report_dt),
        size = 2.4
      ) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    p_annotation_legend <- make_predicted_annotation_legend(
      annotation_summary = annotation_summary,
      o2_plot_min = o2_plot_min,
      o2_plot_max = o2_plot_max
    )

    p_annotation_burden <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "burden",
      y_label = "Burden",
      legend_title = "Burden\n(mm^3)",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#ffffbf", "#542788"),
      transform = "log10",
      show_legend = FALSE
    )
    p_annotation_live <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "live_cells",
      y_label = "Viable cells",
      legend_title = "Viable cells",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#ffffbf", "#2c7bb6"),
      transform = "log10",
      show_legend = FALSE
    )
    p_annotation_o2 <- make_predict_annotation_track_plot(
      df = annotation_summary,
      value_col = "o2_pct",
      y_label = "Effective O2 (%)",
      legend_title = "Effective O2 (%)",
      day_width = chr_density_day_width,
      horizon_day = horizon_day,
      x_breaks = predict_x_breaks,
      colors = c("#f7f7f7", "#9ecae1", "#08519c"),
      transform = "identity",
      limits = c(o2_plot_min, o2_plot_max),
      breaks = unique(c(o2_plot_min, o2_plot_max)),
      labels = format(unique(c(o2_plot_min, o2_plot_max)), trim = TRUE),
      show_legend = FALSE
    )

    chr_density_fill_max_annot <- max(chr_density_df$density, na.rm = TRUE)
    if (!is.finite(chr_density_fill_max_annot) || chr_density_fill_max_annot <= 0) {
      chr_density_fill_max_annot <- 1
    }
    p_annotation_chr_density <- ggplot(
      chr_density_df,
      aes(x = day, y = chr_bin_mid, fill = density)
    ) +
      geom_tile(width = chr_density_day_width, height = chr_density_bin_width) +
      cohort_facet_grid() +
      predict_x_scale() +
      ploidy_fraction_fill_scale(chr_density_fill_max_annot, name = "Chromosome\nprobability") +
      scale_y_continuous(
        breaks = ploidy_n_breaks,
        expand = c(0, 0),
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Ploidy")
      ) +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      labs(
        x = NULL,
        y = "Chromosome count (N)"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    mean_chr_summary <- endpoint_plot_df %>%
      filter(as.character(cohort) %in% c("2N", "4N"), is.finite(day), is.finite(value_chr)) %>%
      mutate(cohort = factor(as.character(cohort), levels = c("2N", "4N"))) %>%
      group_by(cohort, day) %>%
      summarise(value_chr = mean(value_chr, na.rm = TRUE), .groups = "drop")
    p_annotation_mean_chr <- ggplot(
      mean_chr_summary,
      aes(x = day, y = value_chr, group = cohort, color = cohort)
    ) +
      geom_line(linewidth = 0.85, alpha = 0.95) +
      predict_x_scale() +
      coord_cartesian(xlim = c(0, horizon_day), ylim = ploidy_n_limits, expand = FALSE) +
      scale_y_continuous(
        name = "Mean chr. number",
        breaks = ploidy_n_breaks,
        sec.axis = sec_axis(~ . / as.numeric(cfg$N_UNIT), name = "Mean ploidy")
      ) +
      scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
      labs(
        x = "Day",
        color = "Cohort"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
      )

    save_aligned_plot_stack(
      list(
        p_annotation_title,
        p_annotation_legend,
        p_annotation_burden,
        p_annotation_live,
        p_annotation_o2,
        p_annotation_chr_density,
        p_annotation_mean_chr
      ),
      file.path(out_dir, paste0("predicted_", horizon_tag, ".pdf")),
      width = 12,
      height = 4.8,
      row_heights = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.25, 0.25)
    )
  }

  save_aligned_plot_stack(
    list(
      p_predict_burden,
      p_predict_live_cells,
      p_resource_death_fraction_predict,
      p_predict_endpoint,
      p_predict_chr_density,
      p_predict_o2
    ),
    file.path(out_dir, paste0("predict_curves_", horizon_tag, ".pdf")),
    width = 12,
    height = 10.8,
    row_heights = c(1, 1, 1, 1, 1.5, 1)
  )

  burden_decomp_predict <- burden_decomp_predict %>%
    filter(!is.na(cohort)) %>%
    group_by(cohort, day) %>%
    summarise(
      burden_live = mean(burden_live, na.rm = TRUE),
      burden_dead_hypoxia = mean(burden_dead_hypoxia, na.rm = TRUE),
      burden_dead_buffer = mean(burden_dead_buffer, na.rm = TRUE),
      burden_total = mean(burden_total, na.rm = TRUE),
      .groups = "drop"
    )

  burden_decomp_predict_floor <- log10_plot_floor(
    c(
      burden_decomp_predict$burden_live,
      burden_decomp_predict$burden_dead_hypoxia,
      burden_decomp_predict$burden_dead_buffer,
      burden_decomp_predict$burden_total
    ),
    default = 1e-12
  )
  burden_decomp_predict <- burden_decomp_predict %>%
    mutate(burden_total_log_plot = floor_for_log10_plot(burden_total, burden_decomp_predict_floor))
  burden_decomp_predict_long <- make_burden_decomp_long(burden_decomp_predict, death_language)
  burden_decomp_predict_ribbon <- make_burden_decomp_ribbon(
    burden_decomp_predict,
    death_language,
    burden_decomp_predict_floor
  )

  p_burden_decomp_predict <- ggplot(
    burden_decomp_predict_ribbon,
    aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)
  ) +
    geom_ribbon(alpha = 0.55) +
    geom_line(
      data = burden_decomp_predict,
      aes(x = day, y = burden_total_log_plot),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.65
    ) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = stats::setNames(
      c("#1f77b4", "#d62728", "#2ca02c"),
      c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    )) +
    log10_burden_y_scale() +
    coord_cartesian(xlim = c(0, horizon_day)) +
    labs(
      title = paste0("Predicted Total Burden Viable/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days"),
      subtitle = "Cohort-level mean across scenarios (2N top, 4N bottom)",
      x = "Day",
      y = "log10 Tumor burden (mm^3)",
      fill = "Component"
    ) +
    theme_bw(base_size = 11)

  burden_decomp_plot_for_cohort <- function(cohort_use, show_fill_legend = TRUE) {
    row_df <- burden_decomp_predict %>%
      filter(cohort == cohort_use)
    row_ribbon <- burden_decomp_predict_ribbon %>%
      filter(cohort == cohort_use)
    ggplot(
      row_ribbon,
      aes(x = day, ymin = ymin, ymax = ymax, fill = component, group = component)
    ) +
      geom_ribbon(alpha = 0.55) +
      geom_line(
        data = row_df,
        aes(x = day, y = burden_total_log_plot),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.65
      ) +
      scale_fill_manual(
        values = stats::setNames(
          c("#1f77b4", "#d62728", "#2ca02c"),
          c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
        )
      ) +
      log10_burden_y_scale() +
      coord_cartesian(xlim = c(0, horizon_day)) +
      labs(
        title = if (identical(as.character(cohort_use), "2N")) {
          paste0("Predicted Total Burden Viable/Dead Decomposition: 0-", as.integer(round(horizon_day)), " days")
        } else {
          NULL
        },
        subtitle = paste0(as.character(cohort_use), " cohort mean across scenarios"),
        x = "Day",
        y = "log10 Tumor burden (mm^3)",
        fill = "Component"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = if (show_fill_legend) "right" else "none"
      )
  }

  save_aligned_plot_stack(
    list(
      burden_decomp_plot_for_cohort("2N", show_fill_legend = TRUE),
      burden_decomp_plot_for_cohort("4N", show_fill_legend = TRUE)
    ),
    file.path(out_dir, paste0("predict_burden_live_dead_decomposition_", horizon_tag, ".pdf")),
    width = 14,
    height = 7
  )

  invisible(list(
    burden_decomp_predict = burden_decomp_predict,
    burden_decomp_predict_long = burden_decomp_predict_long,
    horizon_day = as.numeric(horizon_day)
  ))
}

render_invivo_observed_plots <- function(simulation_dir, cfg, out_dir) {
  fit_dir <- cfg$fit_dir %||% simulation_dir
  report_dt <- as.numeric(cfg$report_dt %||% 1)
  death_language <- resource_death_language()
  burden_all <- read_invivo_simulation_table(simulation_dir, "burden_timecourse.tsv")
  burden_decomp <- read_invivo_simulation_table(simulation_dir, "burden_live_dead_decomposition.tsv")
  ploidy_mean <- read_invivo_simulation_table(simulation_dir, "ploidy_weighted_mean_timecourse.tsv")
  misseg_timecourse <- read_invivo_simulation_table(simulation_dir, "missegregation_probability_timecourse.tsv")
  terminal_ploidy_compare <- read_invivo_simulation_table(simulation_dir, "terminal_ploidy_observed_vs_predicted.tsv")
  o2_lag_df <- read_invivo_simulation_table(simulation_dir, "o2_lag_timecourse.tsv")
  o2_burden_df <- read_invivo_simulation_table(simulation_dir, "predict_burden_vs_o2.tsv")

  if (nrow(misseg_timecourse)) {
    plot_missegregation_probability_timecourse(misseg_timecourse, out_dir, fit_dir, report_dt)
  }
  if (nrow(terminal_ploidy_compare)) {
    plot_terminal_ploidy_violin_compare(terminal_ploidy_compare, fit_dir, out_dir)
  }

  o2_plot_min <- 0
  o2_plot_max <- 5
  if (nrow(o2_lag_df)) {
    o2_lag_long <- o2_lag_df %>%
      select(harvest, cohort, dose, day, sample_id, o2_target_pct, o2_eff_pct) %>%
      pivot_longer(cols = c("o2_target_pct", "o2_eff_pct"), names_to = "o2_series", values_to = "o2_pct") %>%
      mutate(o2_series = factor(o2_series, levels = c("o2_target_pct", "o2_eff_pct"),
                                labels = c("O2 target", "Effective O2")))
    p_o2_lag <- ggplot(o2_lag_long,
                       aes(x = day, y = o2_pct, color = o2_series, linetype = o2_series,
                           group = interaction(sample_id, o2_series))) +
      geom_line(linewidth = 0.7, alpha = 0.85) +
      facet_wrap(~ harvest, ncol = 2) +
      scale_color_manual(values = c("O2 target" = "#ff7f0e", "Effective O2" = "#1f77b4")) +
      coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
      labs(
        title = "Resource-Stress Model: Effective Oxygen Relaxation Over Time",
        subtitle = "Oxygen supply-demand target vs lagged effective oxygen state",
        x = "Day", y = "Effective oxygen (%)", color = NULL, linetype = NULL
      ) +
      theme_bw(base_size = 11)
    ggsave(file.path(out_dir, "o2_target_vs_eff_timecourse.pdf"), p_o2_lag, width = 13, height = 9)
  }

  rho_2N_min <- as.numeric(cfg$rho_2N_min %||% 3.2e4)
  rho_2N_max <- as.numeric(cfg$rho_2N_max %||% 5.6e4)
  if (!is.finite(rho_2N_min) || rho_2N_min <= 0) rho_2N_min <- 3.2e4
  if (!is.finite(rho_2N_max) || rho_2N_max <= 0) rho_2N_max <- 5.6e4
  if (rho_2N_min > rho_2N_max) {
    tmp <- rho_2N_min
    rho_2N_min <- rho_2N_max
    rho_2N_max <- tmp
  }
  rho_2N_mid <- sqrt(rho_2N_min * rho_2N_max)
  pred_cell_col <- if ("pred_burden_cells" %in% names(burden_all)) "pred_burden_cells" else "pred_burden"
  burden_all_real <- burden_all %>%
    mutate(
      pred_burden_cell_number = as.numeric(.data[[pred_cell_col]]),
      obs_burden_cell_number_low = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_min, NA_real_),
      obs_burden_cell_number_mid = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_mid, NA_real_),
      obs_burden_cell_number_high = ifelse(is.finite(obs_burden), as.numeric(obs_burden) * rho_2N_max, NA_real_)
    )
  p_burden_abs_real <- ggplot(burden_all_real, aes(x = day, y = pred_burden_cell_number)) +
    geom_line(color = "#1f77b4", linewidth = 0.7) +
    geom_ribbon(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_low) & !is.na(obs_burden_cell_number_high)),
      aes(x = day, ymin = obs_burden_cell_number_low, ymax = obs_burden_cell_number_high),
      inherit.aes = FALSE, fill = "grey50", alpha = 0.18
    ) +
    geom_line(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid), color = "black", linewidth = 0.45, linetype = "dashed"
    ) +
    geom_point(
      data = burden_all_real %>% filter(!is.na(obs_burden_cell_number_mid)),
      aes(y = obs_burden_cell_number_mid), color = "black", size = 1
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    labs(
      title = "Resource-Stress Model: In Vivo Tumor Burden Trajectory",
      subtitle = paste0(
        "fit_dir=", basename(fit_dir), " | report_dt=", report_dt,
        " | obs burden -> CellNumber using rho_2N range=[", signif(rho_2N_min, 4), ", ",
        signif(rho_2N_max, 4), "] cells/mm^3 (mid=", signif(rho_2N_mid, 4), ")"
      ),
      x = "Day", y = "CellNumber (2N-equivalent range)"
    ) +
    theme_bw(base_size = 11)

  burden_decomp_long <- burden_decomp %>%
    pivot_longer(cols = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
                 names_to = "component", values_to = "value") %>%
    mutate(component = factor(
      component,
      levels = c("burden_live", "burden_dead_hypoxia", "burden_dead_buffer"),
      labels = c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    ))
  p_burden_decomp <- ggplot(
    burden_decomp_long,
    aes(x = day, y = value, fill = component, group = interaction(component, harvest, cohort, dose))
  ) +
    geom_area(alpha = 0.55, position = "stack") +
    geom_line(
      data = burden_decomp,
      aes(x = day, y = burden_total, group = interaction(harvest, cohort, dose)),
      inherit.aes = FALSE, color = "black", linewidth = 0.6
    ) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = stats::setNames(
      c("#1f77b4", "#d62728", "#2ca02c"),
      c(death_language$live_label, death_language$component_label, death_language$cin_component_label)
    )) +
    labs(
      title = "Resource-Stress Model: Total Tumor Burden Viable/Dead Decomposition",
      subtitle = paste0("Total burden (black) = viable + ", death_language$figure_phrase,
                        " + CIN-associated dead"),
      x = "Day", y = "Tumor burden (mm^3)", fill = "Component"
    ) +
    theme_bw(base_size = 11)

  p_burden_vs_o2 <- ggplot(o2_burden_df,
                           aes(x = burden_mm3, y = o2_pct, color = cohort, group = sample_id)) +
    geom_path(linewidth = 0.75, alpha = 0.9) +
    facet_wrap(~ harvest, ncol = 2, scales = "free_x") +
    scale_color_manual(values = c("2N" = "#1f77b4", "4N" = "#d62728")) +
    coord_cartesian(ylim = c(o2_plot_min, o2_plot_max)) +
    labs(
      title = "Resource-Stress Model: Predicted Effective Oxygen vs Tumor Burden",
      subtitle = paste0("fit_dir=", basename(fit_dir), " | report_dt=", report_dt),
      x = "Tumor burden (mm^3)", y = "Effective oxygen (%)", color = "Cohort"
    ) +
    theme_bw(base_size = 11)
  p_ploidy_weighted_mean <- ggplot(ploidy_mean, aes(x = day, y = weighted_mean_ploidy)) +
    geom_line(color = "#d62728", linewidth = 0.9) +
    facet_wrap(~ harvest, ncol = 2) +
    coord_cartesian(ylim = c(min(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE),
                             max(ploidy_mean$weighted_mean_ploidy, na.rm = TRUE))) +
    labs(
      title = paste0("Resource-Stress Model: ", weighted_mean_series_label(cfg), " Over Time"),
      subtitle = "Weighted by predicted viable chromosome-number state fractions",
      x = "Day", y = weighted_mean_series_label(cfg)
    ) +
    theme_bw(base_size = 11)

  ggsave(file.path(out_dir, "burden_trend_absolute(real_scale).pdf"), p_burden_abs_real, width = 13, height = 9)
  ggsave(file.path(out_dir, "burden_live_dead_decomposition.pdf"), p_burden_decomp, width = 13, height = 9)
  ggsave(file.path(out_dir, "predict_burden_vs_o2.pdf"), p_burden_vs_o2, width = 13, height = 9)
  ggsave(file.path(out_dir, "ploidy_weighted_mean_over_time.pdf"), p_ploidy_weighted_mean, width = 13, height = 9)
  invisible(TRUE)
}
