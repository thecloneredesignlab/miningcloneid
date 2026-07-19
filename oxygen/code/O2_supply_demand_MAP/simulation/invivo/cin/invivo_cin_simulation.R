# In-vivo cin simulation products. This module is data-only.

compute_missegregation_probability_timecourse <-
function(ploidy_all, burden_all, run_params) {
    empty <- data.frame(harvest = character(), cohort = character(), dose = numeric(), day = numeric(), sample_id = character(),
        o2_eff_pct = numeric(), mean_p_misseg = numeric(), weighted_mean_N = numeric(), total_fraction = numeric(), stringsAsFactors = FALSE)
    required_ploidy <- c("harvest", "cohort", "dose", "day", "N", "fraction")
    required_burden <- c("harvest", "cohort", "dose", "day", "pred_o2_pct")
    if (!all(required_ploidy %in% names(ploidy_all)) || !all(required_burden %in% names(burden_all))) {
        return(empty)
    }
    ploidy_df <- ploidy_all %>% transmute(harvest = as.character(harvest), cohort = as.character(cohort), dose = as.numeric(dose),
        day = as.numeric(day), day_key = round(as.numeric(day), 8), N = as.numeric(N), fraction = pmax(as.numeric(fraction),
            0)) %>% filter(is.finite(day_key), is.finite(N), is.finite(fraction), fraction > 0)
    o2_df <- burden_all %>% transmute(harvest = as.character(harvest), cohort = as.character(cohort), dose = as.numeric(dose),
        day = as.numeric(day), day_key = round(as.numeric(day), 8), o2_eff_pct = as.numeric(pred_o2_pct)) %>% filter(is.finite(day_key),
        is.finite(o2_eff_pct)) %>% distinct(harvest, cohort, dose, day_key, .keep_all = TRUE)
    joined <- inner_join(ploidy_df, o2_df[, c("harvest", "cohort", "dose", "day_key", "o2_eff_pct"), drop = FALSE], by = c("harvest",
        "cohort", "dose", "day_key"))
    if (!nrow(joined))
        return(empty)
    joined$p_misseg <- pmax(pmin(as.numeric(.pmisseg_of_O2(O2 = joined$o2_eff_pct, run_params = run_params, N = joined$N)),
        1), 0)
    joined <- joined[is.finite(joined$p_misseg), , drop = FALSE]
    if (!nrow(joined))
        return(empty)
    out <- joined %>% group_by(harvest, cohort, dose, day, o2_eff_pct) %>% summarise(mean_p_misseg = sum(fraction * p_misseg,
        na.rm = TRUE)/pmax(sum(fraction, na.rm = TRUE), 1e-12), weighted_mean_N = sum(fraction * N, na.rm = TRUE)/pmax(sum(fraction,
        na.rm = TRUE), 1e-12), total_fraction = sum(fraction, na.rm = TRUE), .groups = "drop") %>% mutate(sample_id = paste(as.character(harvest),
        as.character(cohort), format(as.numeric(dose), trim = TRUE, scientific = FALSE), sep = "__")) %>% select(harvest,
        cohort, dose, day, sample_id, o2_eff_pct, mean_p_misseg, weighted_mean_N, total_fraction) %>% arrange(harvest, cohort,
        dose, day)
    as.data.frame(out, stringsAsFactors = FALSE)
}

generate_population_average_cin_rows <-
function(ms_timecourse, target_day) {
    cohort_levels <- c("2N", "4N")
    cohort_labels <- c(`2N` = "2N-derived", `4N` = "4N-derived")
    ms_timecourse %>% transmute(initial_cohort = as.character(cohort), day = as.numeric(day), sample_id = as.character(sample_id),
        population_average_cin = as.numeric(mean_p_misseg)) %>% filter(initial_cohort %in% cohort_levels, is.finite(day),
        day >= -1e-06, day <= target_day + 1e-06, is.finite(population_average_cin)) %>% distinct(initial_cohort, day, sample_id,
        population_average_cin, .keep_all = TRUE) %>% arrange(initial_cohort, sample_id, day) %>% group_by(initial_cohort) %>%
        mutate(.first_sample_id = first(sample_id)) %>% ungroup() %>% filter(sample_id == .first_sample_id) %>% transmute(target_day = as.numeric(target_day),
        day = as.numeric(day), initial_cohort = as.character(initial_cohort), cohort_label = as.character(cohort_labels[initial_cohort]),
        cohort_order = as.integer(match(initial_cohort, cohort_levels)), sample_id = as.character(sample_id), population_average_cin = as.numeric(population_average_cin)) %>%
        arrange(target_day, cohort_order, day)
}
