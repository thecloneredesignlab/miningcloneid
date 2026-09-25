# In-vivo ploidy simulation products. This module is data-only.

compute_ploidy_weighted_mean <-
function(ploidy_all, cfg) {
    start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
    ploidy_all %>% group_by(harvest, cohort, dose, day) %>% summarise(weighted_mean_N = sum(N * fraction, na.rm = TRUE)/pmax(sum(fraction,
        na.rm = TRUE), 1e-12), .groups = "drop") %>% mutate(weighted_mean_ploidy = if (identical(start_with_mode, "chr_number")) {
        weighted_mean_N
    }
    else {
        weighted_mean_N/cfg$N_UNIT
    }, weighted_mean_endpoint = weighted_mean_ploidy, start_with = start_with_mode)
}

build_terminal_ploidy_compare_df <-
function(scenarios, ploidy_all, cfg) {
    start_with_mode <- assert_canonical_start_with_mode(.first_non_null_local(cfg$start_with, "ploidy"))
    meta_rows <- vector("list", length(scenarios))
    obs_rows <- vector("list", length(scenarios))
    for (i in seq_along(scenarios)) {
        sc <- scenarios[[i]]
        obs_value <- if (identical(start_with_mode, "chr_number")) {
            as.numeric(sc$chr_number_obs)
        }
        else {
            obs_z_raw <- as.numeric(sc$ploidy_obs_z)
            obs_N <- map_ploidy_to_N_by_chrlen(ploidy_values = obs_z_raw/as.numeric(cfg$N_UNIT), N_grid = cfg$N_MIN:cfg$N_MAX,
                chr_lengths_bp = cfg$chr_lengths_bp)
            obs_N <- as.integer(clip(obs_N, cfg$N_MIN, cfg$N_MAX))
            as.numeric(obs_N)/as.numeric(cfg$N_UNIT)
        }
        obs_value <- obs_value[is.finite(obs_value)]
        if (length(obs_value) == 0L)
            next
        target_day <- as.numeric(sc$sim_end_day)
        if (!is.finite(target_day))
            next
        meta_rows[[i]] <- data.frame(harvest = as.character(sc$harvest), cohort = as.character(sc$cohort), dose = as.numeric(sc$dose),
            target_day = target_day, stringsAsFactors = FALSE)
        obs_rows[[i]] <- data.frame(harvest = as.character(sc$harvest), cohort = as.character(sc$cohort), dose = as.numeric(sc$dose),
            target_day = target_day, source = "Observed", ploidy = as.numeric(obs_value), endpoint_value = as.numeric(obs_value),
            endpoint_mode = start_with_mode, weight = 1, stringsAsFactors = FALSE)
    }
    meta_df <- bind_rows(meta_rows)
    obs_df <- bind_rows(obs_rows)
    if (nrow(meta_df) == 0L || nrow(obs_df) == 0L || nrow(ploidy_all) == 0L) {
        return(data.frame())
    }
    pred_df <- ploidy_all %>% inner_join(meta_df, by = c("harvest", "cohort", "dose")) %>% mutate(day_gap = abs(as.numeric(day) -
        as.numeric(target_day))) %>% group_by(harvest, cohort, dose, target_day) %>% filter(day_gap == min(day_gap, na.rm = TRUE)) %>%
        ungroup() %>% transmute(harvest = as.character(harvest), cohort = as.character(cohort), dose = as.numeric(dose),
        target_day = as.numeric(target_day), source = "Predicted", ploidy = if (identical(start_with_mode, "chr_number"))
            as.numeric(N)
        else as.numeric(N)/as.numeric(cfg$N_UNIT), endpoint_value = if (identical(start_with_mode, "chr_number"))
            as.numeric(N)
        else as.numeric(N)/as.numeric(cfg$N_UNIT), endpoint_mode = start_with_mode, weight = as.numeric(fraction)) %>% filter(is.finite(ploidy),
        is.finite(weight), weight > 0)
    bind_rows(obs_df, pred_df) %>% mutate(cohort = factor(as.character(cohort), levels = c("2N", "4N")), source = factor(as.character(source),
        levels = c("Observed", "Predicted")), fill_group = factor(paste(as.character(cohort), as.character(source)), levels = c("2N Observed",
        "2N Predicted", "4N Observed", "4N Predicted")))
}
