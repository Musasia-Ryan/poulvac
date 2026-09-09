# POULVAC E. COLI CLUSTER-RANDOMISED CROSSOVER TRIAL

# load the packages
library(pacman)

pacman::p_load(
  readxl, 
  dplyr,
  tidyr, 
  tibble, 
  stringr, 
  forcats, 
  purrr, 
  ggplot2,
  geepack,      
  MASS,         
  lme4,         
  lmerTest,     
  sandwich,     
  lmtest,       
  ape,          
  patchwork, 
  cowplot, 
  sf
)

# arm colours
# control group (orange)
col_ctrl    <- "#E0A33E"   
# vaccinated group (blue)
col_vacc    <- "#2F6DB4"   
# rejected group
col_reject  <- "#C0392B"   
col_grey    <- "#7F7F7F"
theme_set(theme_classic(base_size = 12))
#dir.create("output", showWarnings = FALSE)


# loading data ------------------------------------------------------------
vaccine_trial <- read_excel("data/Poulvac_Master_44farm_Analytical_27_June_2026_Final.xlsx",
                             sheet = "Poulvac_Vaccine_Trial_Tool_F...")
bird_weighing <- read_excel("data/Poulvac_Master_44farm_Analytical_27_June_2026_Final.xlsx",
                             sheet = "bird_weighing_repeat")
farms <- read_excel("data/Final Exit Farmer Questionnaire _Poulvac Study_June 2026_final August.xlsx",
                             sheet = "44 Farms")
excluded_farms <- read_excel("data/Final Exit Farmer Questionnaire _Poulvac Study_June 2026_final August.xlsx",
                             sheet = "Excluded farms")


# data cleaning -----------------------------------------------------------

# 2.1 antimicrobial classification
# Active-ingredient stems 
# An antibiotic-USE event requires 
# (a) antibiotics_used_past_7days == "yes" AND
# (b) at least one antibacterial named. Records naming ONLY an anticoccidial (amprolium,toltrazuril) or a mucolytic (ambroxonil) are reclassified OUT of the outcome.

abx_stem <- "oxytet|tetracyc|chlortetra|\\botc\\b|doxycyc|doxyclin|doxynor|neomyc|gentam|enroflox|enrolflox|ciproflox|tylos|tyrosin|tylo|colistin|amoxic|amoxyc|sulfa|sulph|trimethoprim|cotrima"
hpcia_stem <- "enroflox|enrolflox|ciproflox|tylos|tyrosin|tylo|colistin"  # WHO HPCIA: fluoroquinolone, macrolide, polymyxin
fq_stem <- "enroflox|enrolflox|ciproflox"
macro_stem <- "tylos|tyrosin|tylo"
poly_stem <- "colistin"
excl_stem <- "amprol|toltraz|ambrox"

# 2.2 visit-level enrichment
day_map <- c(day_0 = 0, 
             day_3 = 3, 
             day_5 = 5, 
             day_14 = 14,
             day_21 = 21, 
             day_28 = 28, 
             day_35 = 35)

# two visit rows carry no record id (_id)
visit <- vaccine_trial |>
  dplyr::mutate(
    farm_id = ifelse(!is.na(Farm_identification_number_Arusha) &
                       Farm_identification_number_Arusha != "",
                     as.character(Farm_identification_number_Arusha),
                     as.character(Farm_identification_number_Moshi)),
    location = str_to_title(location),
    study_arm = factor(ifelse(study_arm == "treatment", "Vaccinated", "Control"),
                     levels = c("Control", "Vaccinated")),
    period  = factor(Which_phase_of_the_trial),
    day = day_map[study_day],
    lat = as.numeric(`_gps_location_latitude`),
    lon = as.numeric(`_gps_location_longitude`),
    abx_text = tolower(paste(replace_na(antibiotic_type, ""),
                             replace_na(treatment1_ingredient, ""),
                             replace_na(treatment2_ingredient, ""))),
    yes_abx = antibiotics_used_past_7days == "yes",
    has_abx = str_detect(abx_text, abx_stem),
    amu_event = yes_abx & has_abx,
    hpcia_event = amu_event & str_detect(abx_text, hpcia_stem),
    fq_event = amu_event & str_detect(abx_text, fq_stem),
    macro_event = amu_event & str_detect(abx_text, macro_stem),
    poly_event = amu_event & str_detect(abx_text, poly_stem),
    amu_raw = yes_abx, 
    window = ifelse(day <= 14, "restriction", "free")
  )

# 2.3 Treatment sequence per farm (VC vs CV)
# using the phase-1 arm from any visit 
# one farm-cycle has no day-0 record
sequence_tbl <- visit |> 
  dplyr::filter(period == "phase_1") |>
  distinct(farm_id, study_arm) |>
  dplyr::mutate(sequence = ifelse(study_arm == "Vaccinated", "VC", "CV")) |>
  dplyr::select(farm_id, sequence)


# 2.4 Farm-cycle analytic dataset (88 rows)
# birds placed is computed inside the farm-cycle, one cycle -mo12 phase-1- has no day-0 row, so birds_placed falls back to the earliest recorded live-bird count.
fc <- visit |>
  group_by(farm_id, location, study_arm, period) |>
  summarise(
    n_visits = n(),
    bp_max  = suppressWarnings(max(birds_purchased, na.rm = TRUE)),
    cc_early = suppressWarnings(max(current_bird_count[day <= 5], na.rm = TRUE)),
    deaths21 = sum(dead_birds_count[day <= 21], na.rm = TRUE),
    deaths28 = sum(dead_birds_count[day <= 28], na.rm = TRUE),
    amu_full = sum(amu_event, na.rm = TRUE),
    amu_raw_full= sum(amu_raw,   na.rm = TRUE),
    amu_restr = as.integer(any(amu_event & window == "restriction", na.rm = TRUE)),
    amu_free = as.integer(any(amu_event & window == "free", na.rm = TRUE)),
    amu_free_raw= as.integer(any(amu_raw   & window == "free", na.rm = TRUE)),
    hpcia_full = sum(hpcia_event, na.rm = TRUE),
    hpcia_any = as.integer(any(hpcia_event, na.rm = TRUE)),
    hpcia_free = sum(hpcia_event[window == "free"], na.rm = TRUE),
    fq_full = sum(fq_event, na.rm = TRUE),
    macro_full = sum(macro_event, na.rm = TRUE),
    poly_full = sum(poly_event, na.rm = TRUE),
    sick_visits = sum(sick_birds_present == "yes", na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    .groups = "drop") |>
  mutate(birds_placed = ifelse(is.finite(bp_max), bp_max, cc_early)) |>
  left_join(sequence_tbl, by = "farm_id") |>
  mutate(mort21 = 100 * deaths21 / birds_placed,
         mort28 = 100 * deaths28 / birds_placed,
         farm   = factor(farm_id),
         sequence = factor(sequence)) |>
  arrange(farm_id, study_arm)

# birds-placed lookup for downstream joins
placed <- fc |> 
  dplyr::select(farm_id, study_arm, period, birds_placed)

# 2.6 weighing dataset (mean body weight per farm-cycle-day)
wt_cycle_day <- bird_weighing |>
  dplyr::filter(`bird_weighing_repeat/bird_weight_g` > 0) |>
  inner_join(visit |> 
               dplyr::filter(!is.na(`_id`)) |> 
               dplyr::select(`_id`, farm_id, study_arm, period, day),
             by = c("_submission__id" = "_id"), relationship = "many-to-one",
             na_matches = "never") |>
  dplyr::group_by(farm_id, study_arm, period, day) |>
  dplyr::summarise(mean_wt = mean(`bird_weighing_repeat/bird_weight_g`),
                   sd_wt   = sd(`bird_weighing_repeat/bird_weight_g`),
                   n_birds = n(),
                   cv_wt   = 100 * sd_wt / mean_wt,
            .groups = "drop")

## average daily weight gain (g/bird/day) at farm-cycle level
adwg <- wt_cycle_day |>
  dplyr::filter(day %in% c(3, 14, 21, 28)) |>
  dplyr::select(farm_id, study_arm, period, day, mean_wt) |>
  pivot_wider(names_from = day, values_from = mean_wt, names_prefix = "d") |>
  dplyr::mutate(
    adwg_3_21  = (d21 - d3) / 18,
    adwg_3_28  = (d28 - d3) / 25,
    adwg_3_14  = (d14 - d3) / 11,
    adwg_14_21 = (d21 - d14) / 7)

# prints
cat("\n[CHECK] farm-cycles:", nrow(fc), "| farms:", n_distinct(fc$farm_id),
    "| ADWG(3-21) cycles:", sum(!is.na(adwg$adwg_3_21)), "\n")


# descriptive counts by arm 
# cycles
# birds placed
# deaths
# mortality
# events
fc |> group_by(study_arm) |>
  summarise(cycles = n(), birds_placed = sum(birds_placed),
            deaths21 = sum(deaths21), mortality21_pct = round(100*sum(deaths21)/sum(birds_placed),2),
            amu_events = sum(amu_full), hpcia_events = sum(hpcia_full),
            restr_use_cycles = sum(amu_restr), free_use_cycles = sum(amu_free)) |>
  as.data.frame() |> 
  print()

# full-cycle antibiotic use
m_amu_gee <- geeglm(amu_full ~ study_arm + period, 
                    id = farm,
                    data = fc,
                    family = poisson, 
                    corstr = "exchangeable", 
                    offset = log(n_visits))


irr <- exp(coef(m_amu_gee)["study_armVaccinated"])
irci <- exp(confint.default(m_amu_gee)["study_armVaccinated", ])
cat(sprintf("AMU IRR (Vaccinated vs Control) = %.2f (95%% CI %.2f-%.2f), p=%.1e\n",
            irr, irci[1], irci[2], 
            summary(m_amu_gee)$coefficients["study_armVaccinated","Pr(>|W|)"]))


# within farm fixed-effects Poisson (supportive)
m_amu_fe <- glm(amu_full ~ study_arm + period + farm,
                data = fc,
                family = poisson, 
                offset = log(n_visits))

cat(sprintf("Within-farm FE Poisson IRR = %.2f\n", exp(coef(m_amu_fe)["study_armVaccinated"])))
cat("Overdispersion (Pearson chi2/df):",
    round(sum(residuals(m_amu_gee, "pearson")^2)/(nrow(fc)-length(coef(m_amu_gee))),2), "\n")

# restriction window (days 0-14): McNemar's test -------------
restr_wide <- fc |> 
  dplyr::select(farm_id, study_arm, amu_restr) |>
  pivot_wider(names_from = study_arm, values_from = amu_restr)

cat("Vaccinated used:", sum(restr_wide$Vaccinated), "/44 ; Control used:",
    sum(restr_wide$Control), "/44\n")

# view the mcnemar test result
print(mcnemar.test(table(restr_wide$Vaccinated, restr_wide$Control)))

if (has_logistf) {
  fdat <- fc |> mutate(vacc = as.integer(study_arm == "Vaccinated"))
  mf <- logistf::logistf(amu_restr ~ vacc, data = fdat)
  cat("Firth OR (restriction) =", round(exp(coef(mf)["vacc"]),3), "\n")
}

# free-choice window (days 15-35): logistic GEE -------
m_free <- geeglm(amu_free ~ study_arm + period, 
                 id = farm, 
                 data = fc,
                 family = binomial,
                 corstr = "exchangeable")

orf  <- exp(coef(m_free)["study_armVaccinated"]); orci <- exp(confint.default(m_free)["study_armVaccinated",])


cat(sprintf("Free-window OR (reclassified) = %.2f (%.2f-%.2f), p=%.2f\n",
            orf, orci[1], orci[2], summary(m_free)$coefficients["study_armVaccinated","Pr(>|W|)"]))

# sensitivity: unreclassified counts
m_free_raw <- geeglm(amu_free_raw ~ study_arm + period,
                     id = farm, 
                     data = fc,
                     family = binomial,
                     corstr = "exchangeable")

cat(sprintf("Free-window OR (UN-reclassified) = %.2f\n",
            exp(coef(m_free_raw)["study_armVaccinated"])))

# HPCIA use --------------------------------------------
m_hpcia <- geeglm(hpcia_full ~ study_arm + period,
                  id = farm, 
                  data = fc,
                  family = poisson, 
                  corstr = "exchangeable", 
                  offset = log(n_visits))

cat(sprintf("HPCIA IRR = %.2f (%.2f-%.2f)\n", exp(coef(m_hpcia)["study_armVaccinated"]),
            exp(confint.default(m_hpcia)["study_armVaccinated",])[1],
            exp(confint.default(m_hpcia)["study_armVaccinated",])[2]))
cat("HPCIA events: Control", sum(fc$hpcia_full[fc$study_arm=="Control"]),
    "vs Vaccinated", sum(fc$hpcia_full[fc$study_arm=="Vaccinated"]), "\n")



# Mortality ---------------------------------------------------------------
mort_spec <- data.frame()
for (ep in c(21, 28)) {
  dv <- if (ep == 21) "deaths21" else "deaths28"
  fc$y  <- fc[[dv]]
  off   <- log(fc$birds_placed)
  
  # farm-weighted cluster-level RR
  # each farm its own control
  w <- fc |> 
    dplyr::select(farm_id, 
                  study_arm, 
                  y,
                  birds_placed) |>
    pivot_wider(names_from = study_arm, values_from = c(y, birds_placed))
  lrr <- log((w$y_Vaccinated + 0.5)/w$birds_placed_Vaccinated) -
    log((w$y_Control    + 0.5)/w$birds_placed_Control)
  K <- length(lrr); mrr <- mean(lrr); se <- sd(lrr)/sqrt(K); tc <- qt(.975, K-1)
  mort_spec <- rbind(mort_spec, 
                     data.frame(
                       endpoint = ep,
                       model = "Farm-weighted cluster-level RR (primary)",
                       RR = exp(mrr), 
                       lo = exp(mrr - tc*se), 
                       hi = exp(mrr + tc*se), 
                       p = 2*pt(-abs(mrr/se), K-1)))
  
  cat(sprintf("Day %d: Wilcoxon signed-rank p on farm log-RR = %.2f\n", ep, wilcox.test(lrr)$p.value))
  
  # crude bird-weighted RR (unadjusted)
  rr_crude <- (sum(fc$y[fc$study_arm=="Vaccinated"])/sum(fc$birds_placed[fc$study_arm=="Vaccinated"])) /
    (sum(fc$y[fc$study_arm=="Control"])   /sum(fc$birds_placed[fc$study_arm=="Control"]))
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Crude bird-weighted RR (unadjusted)",
                                           RR=rr_crude, lo=NA, hi=NA, p=NA))
  
  # Naive Poisson
  # ignores clustering AND overdispersion
  m1 <- glm(y ~ study_arm, data = fc, family = poisson, offset = off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep,
                                           model = if (ep==21) "Naive Poisson (rejected)" else "Naive Poisson",
                                           RR=exp(coef(m1)["study_armVaccinated"]), lo=exp(confint.default(m1)["study_armVaccinated",1]),
                                           hi=exp(confint.default(m1)["study_armVaccinated",2]), p=coef(summary(m1))["study_armVaccinated",4]))
  
  # Naive Poisson + period
  m2 <- glm(y ~ study_arm + period, data = fc, family = poisson, offset = off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Naive Poisson + period",
                                           RR=exp(coef(m2)["study_armVaccinated"]), lo=exp(confint.default(m2)["study_armVaccinated",1]),
                                           hi=exp(confint.default(m2)["study_armVaccinated",2]), p=coef(summary(m2))["study_armVaccinated",4]))
  
  # Poisson GEE, exchangeable, robust SE
  g_ex <- geeglm(y ~ study_arm, id=farm, data=fc, family=poisson, corstr="exchangeable", offset=off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Poisson GEE, exchangeable, robust SE",
                                           RR=exp(coef(g_ex)["study_armVaccinated"]), lo=exp(confint.default(g_ex)["study_armVaccinated",1]),
                                           hi=exp(confint.default(g_ex)["study_armVaccinated",2]), p=summary(g_ex)$coefficients["study_armVaccinated","Pr(>|W|)"]))
  
  # Poisson GEE, independence, robust SE
  g_in <- geeglm(y ~ study_arm, id=farm, data=fc, family=poisson, corstr="independence", offset=off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Poisson GEE, independence, robust SE",
                                           RR=exp(coef(g_in)["study_armVaccinated"]), lo=exp(confint.default(g_in)["study_armVaccinated",1]),
                                           hi=exp(confint.default(g_in)["study_armVaccinated",2]), p=summary(g_in)$coefficients["study_armVaccinated","Pr(>|W|)"]))
  
  # Poisson GEE, exchangeable + period
  g_exp <- geeglm(y ~ study_arm + period, id=farm, data=fc, family=poisson, corstr="exchangeable", offset=off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Poisson GEE, exchangeable + period",
                                           RR=exp(coef(g_exp)["study_armVaccinated"]), lo=exp(confint.default(g_exp)["study_armVaccinated",1]),
                                           hi=exp(confint.default(g_exp)["study_armVaccinated",2]), p=summary(g_exp)$coefficients["study_armVaccinated","Pr(>|W|)"]))
  
  # Poisson GEE, exchangeable + period + sequence
  g_exps <- geeglm(y ~ study_arm + period + sequence, id=farm, data=fc, family=poisson, corstr="exchangeable", offset=off)
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Poisson GEE, exchangeable + period + sequence",
                                           RR=exp(coef(g_exps)["study_armVaccinated"]), lo=exp(confint.default(g_exps)["study_armVaccinated",1]),
                                           hi=exp(confint.default(g_exps)["study_armVaccinated",2]), p=summary(g_exps)$coefficients["study_armVaccinated","Pr(>|W|)"]))
  
  # NB-GLM, cluster-robust SE
  nb  <- MASS::glm.nb(y ~ study_arm + offset(off), data = fc)
  vc  <- sandwich::vcovCL(nb, cluster = fc$farm)
  cinb<- coef(nb)["study_armVaccinated"] + c(-1.96,1.96)*sqrt(vc["study_armVaccinated","study_armVaccinated"])
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="NB-GLM, cluster-robust SE",
                                           RR=exp(coef(nb)["study_armVaccinated"]), lo=exp(cinb[1]), hi=exp(cinb[2]), p=NA))
  
  # NB-GLM + period, cluster-robust SE
  nbp  <- MASS::glm.nb(y ~ study_arm + period + offset(off), data = fc)
  vcp  <- sandwich::vcovCL(nbp, cluster = fc$farm)
  cinbp<- coef(nbp)["study_armVaccinated"] + c(-1.96,1.96)*sqrt(vcp["study_armVaccinated","study_armVaccinated"])
  mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="NB-GLM + period, cluster-robust SE",
                                           RR=exp(coef(nbp)["study_armVaccinated"]), lo=exp(cinbp[1]), hi=exp(cinbp[2]), p=NA))
  
  # Day-21-only extra specifications
  if (ep == 21) {
    # farm random intercept
    nbmm <- tryCatch(lme4::glmer.nb(y ~ study_arm + (1|farm) + offset(off), data = fc), error=function(e) NULL)
    if (!is.null(nbmm)) {
      se_mm <- sqrt(diag(vcov(nbmm)))["study_armVaccinated"]
      mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="NB-GLMM, farm random intercept",
                                               RR=exp(fixef(nbmm)["study_armVaccinated"]), lo=exp(fixef(nbmm)["study_armVaccinated"]-1.96*se_mm),
                                               hi=exp(fixef(nbmm)["study_armVaccinated"]+1.96*se_mm), p=NA))
    }
    # farm fixed effects
    nbfe <- tryCatch(MASS::glm.nb(y ~ study_arm + farm + offset(off), data = fc), error=function(e) NULL)
    if (!is.null(nbfe)) mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="NB, farm fixed effects",
                                                                 RR=exp(coef(nbfe)["study_armVaccinated"]), lo=exp(confint.default(nbfe)["study_armVaccinated",1]),
                                                                 hi=exp(confint.default(nbfe)["study_armVaccinated",2]), p=NA))
    # longitudinal NB with time-varying population-at-risk offset
    long_dat <- visit |> filter(day %in% c(3,5,14,21), !is.na(dead_birds_count),
                                 !is.na(current_bird_count), current_bird_count > 0)
    nblt <- tryCatch(MASS::glm.nb(dead_birds_count ~ study_arm + period + offset(log(current_bird_count)),
                                  data = long_dat), error=function(e) NULL)
    if (!is.null(nblt)) {
      vlt <- sandwich::vcovCL(nblt, cluster = long_dat$farm_id)
      cilt<- coef(nblt)["study_armVaccinated"] + c(-1.96,1.96)*sqrt(vlt["study_armVaccinated","study_armVaccinated"])
      mort_spec <- rbind(mort_spec, data.frame(endpoint=ep, model="Longitudinal NB, time-varying offset",
                                               RR=exp(coef(nblt)["study_armVaccinated"]), lo=exp(cilt[1]), hi=exp(cilt[2]), p=NA))
    }
  }
}
mort_spec$endpoint <- paste0("Day ", mort_spec$endpoint)
print(mort_spec |> mutate(across(c(RR,lo,hi), ~round(.,2)), p=round(p,3)), row.names = FALSE)

# growth adjusted weight gain: linear mixed models --------------------------
for (v in c("adwg_3_21","adwg_3_28","adwg_3_14","adwg_14_21")) {
  dd <- adwg[!is.na(adwg[[v]]), ]
  m  <- lmer(as.formula(paste(v, "~ study_arm + period + (1|farm_id)")), data = dd)
  co <- summary(m)$coefficients["study_armVaccinated", ]
  icc <- as.numeric(VarCorr(m)$farm_id) / (as.numeric(VarCorr(m)$farm_id) + sigma(m)^2)
  cat(sprintf("%-11s n=%d  Vacc-Ctrl = %+.2f (%.2f, %.2f)  p=%.3f  ICC=%.2f\n",
              v, nrow(dd), co[1], co[1]-1.96*co[2], co[1]+1.96*co[2], co[5], icc))
}

# vaccine dose adherence -------------------------------
# target reconstituted dose = 1 mL per 500 live birds present at vaccination.
dose <- visit |> filter(study_arm == "Vaccinated", day == 5) |>
  transmute(farm_id, period,
            vacc_ml = as.numeric(vaccine_quantity_ml),
            live_birds = as.numeric(current_bird_count),
            target_ml = live_birds / 500,
            dose_ratio = vacc_ml / target_ml) |>
  filter(is.finite(dose_ratio), dose_ratio > 0)
dose <- dose |> mutate(in_range = dose_ratio >= 0.8 & dose_ratio <= 1.2)
cat(sprintf("Within +/-20%% of target: %d/%d (%.0f%%); median dose ratio = %.2f\n",
            sum(dose$in_range), nrow(dose), 100*mean(dose$in_range), median(dose$dose_ratio)))

# clinical signs (field-recorded) ----------------------
m_sick <- geeglm(sick_visits ~ study_arm + period,
                 id = farm, 
                 data = fc,
                 family = poisson, 
                 corstr = "exchangeable", 
                 offset = log(n_visits))

cat(sprintf("Sick-visit IRR = %.2f (%.2f-%.2f); Vacc %d vs Ctrl %d visits\n",
            exp(coef(m_sick)["study_armVaccinated"]),
            exp(confint.default(m_sick)["study_armVaccinated",])[1],
            exp(confint.default(m_sick)["study_armVaccinated",])[2],
            sum(fc$sick_visits[fc$study_arm=="Vaccinated"]), sum(fc$sick_visits[fc$study_arm=="Control"])))



# distribution of the three primary outcomes -----
d_death <- fc$deaths21
p2a <- ggplot(data.frame(x = d_death), aes(x)) +
  geom_histogram(binwidth = 2, fill = "#AEC7D6", colour = "white", boundary = 0) +
  stat_function(fun = function(z) dnorm(z, mean(d_death), sd(d_death)) * length(d_death) * 5,
                colour = col_reject, linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Variance/mean = %.1f", var(d_death)/mean(d_death)), colour = col_grey) +
  labs(title = "(A) Day-21 death counts", x = "Deaths per farm-cycle", y = "Farm-cycles")

p2b <- ggplot(fc, aes(amu_full)) +
  geom_histogram(binwidth = 1, fill = "#AEC7D6", colour = "white", boundary = -0.5) +
  stat_function(fun = function(z) dnorm(z, mean(fc$amu_full), sd(fc$amu_full)) * nrow(fc),
                colour = col_reject, linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Variance/mean = %.1f", var(fc$amu_full)/mean(fc$amu_full)), colour = col_grey) +
  labs(title = "(B) Antibiotic-use events", x = "Events per farm-cycle", y = "Farm-cycles")

adwg_v <- adwg$adwg_3_21[!is.na(adwg$adwg_3_21)]

p2c <- ggplot(data.frame(x = adwg_v), aes(x)) +
  geom_histogram(bins = 10, fill = "#AEC7D6", colour = "white") +
  stat_function(fun = function(z) dnorm(z, mean(adwg_v), sd(adwg_v)) * length(adwg_v) *
                  (diff(range(adwg_v))/10), colour = col_reject, linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           label = sprintf("Shapiro-Wilk W = %.2f, p = %.2f",
                           shapiro.test(adwg_v)$statistic, shapiro.test(adwg_v)$p.value), colour = col_grey) +
  labs(title = "(C) ADWG, days 3-21", x = "g/bird/day", y = "Farm-cycles")

fig2 <- cowplot::plot_grid(p2a, p2b, p2c, nrow = 1, labels = NULL)

ggsave("output/Figure2_distributions.png", fig2, width = 12, height = 3.6, dpi = 150)

# antibiotic use across the cycle by arm ---------
amu_by_day <- visit |> group_by(study_arm, day) |>
  summarise(pct = 100 * mean(amu_event), .groups = "drop")

fig3 <- ggplot(amu_by_day, aes(day, pct, colour = study_arm, shape = study_arm)) +
  annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 80, alpha = .08, fill = col_vacc) +
  annotate("text", x = 7, y = 76, label = "Antibiotic-restriction window\n(days 0-14)", colour = col_vacc) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_colour_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc)) +
  scale_x_continuous(breaks = c(0,3,5,14,21,28,35)) +
  labs(x = "Study day", y = "Monitoring visits with antibiotic use (%)", colour = NULL, shape = NULL) +
  theme(legend.position = c(.85,.85))

ggsave("output/Figure3_AMU_timecourse.png", fig3, width = 9, height = 4.5, dpi = 150)

# mortality rate-ratio forest (day 21) -----------
forest <- mort_spec |> filter(endpoint == "Day 21",
                               model %in% c("Farm-weighted cluster-level RR (primary)","Poisson GEE, exchangeable, robust SE",
                                            "NB-GLM, cluster-robust SE","NB, farm fixed effects","Longitudinal NB, time-varying offset",
                                            "Crude bird-weighted RR (unadjusted)","Naive Poisson (rejected)"))
forest$model <- factor(forest$model, levels = rev(c(
  "Farm-weighted cluster-level RR (primary)","Poisson GEE, exchangeable, robust SE",
  "NB-GLM, cluster-robust SE","NB, farm fixed effects","Longitudinal NB, time-varying offset",
  "Crude bird-weighted RR (unadjusted)","Naive Poisson (rejected)")))
forest$col <- ifelse(grepl("primary", forest$model), "primary",
                     ifelse(grepl("rejected|Crude", forest$model), "reject", "support"))
fig4 <- ggplot(forest, aes(RR, model, colour = col)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = col_grey) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 1, na.rm = TRUE) +
  geom_point(size = 4) +
  geom_text(aes(x = 1.6, label = ifelse(is.na(lo), sprintf("%.2f", RR),
                                        sprintf("%.2f (%.2f-%.2f)", RR, lo, hi))), colour = "black", hjust = 0) +
  scale_x_log10(limits = c(0.4, 2.2)) +
  scale_colour_manual(values = c(primary = col_vacc, support = col_grey, reject = col_reject), guide = "none") +
  labs(x = "Rate ratio, vaccinated versus control (log scale)", y = NULL)

ggsave("output/Figure4_mortality_forest.png", fig4, width = 10, height = 4.2, dpi = 150)

# interval-specific ADWG difference 
iv <- data.frame()
for (v in c("adwg_3_14","adwg_14_21","adwg_3_21")) {
  dd <- adwg[!is.na(adwg[[v]]), ]
  m <- lmer(as.formula(paste(v, "~ study_arm + period + (1|farm_id)")), data = dd)
  co <- summary(m)$coefficients["study_armVaccinated", ]
  iv <- rbind(iv, data.frame(window = v, diff = co[1], lo = co[1]-1.96*co[2], hi = co[1]+1.96*co[2]))
}
iv$window <- factor(iv$window, levels = c("adwg_3_14","adwg_14_21","adwg_3_21"),
                    labels = c("Days 3-14\n(vaccination & restriction)","Days 14-21\n(post-restriction)","Days 3-21\n(primary endpoint)"))
fig5 <- ggplot(iv, aes(window, diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = col_grey) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = window == levels(window)[1]), linewidth = 1.2) +
  geom_point(aes(colour = window == levels(window)[1]), shape = 18, size = 5) +
  geom_text(aes(label = sprintf("%+.2f", diff)), hjust = -0.3) +
  scale_colour_manual(values = c(`TRUE` = col_vacc, `FALSE` = col_grey), guide = "none") +
  labs(x = NULL, y = "Difference in ADWG (g/bird/day)\nvaccinated - control")

ggsave("output/Figure5_interval_ADWG.png", fig5, width = 9, height = 4.5, dpi = 150)

# vaccine dose adherence -------------------------
dose_ord <- dose |> arrange(dose_ratio) |> mutate(idx = row_number())
fig6 <- ggplot(dose_ord, aes(idx, dose_ratio, fill = in_range)) +
  annotate("rect", xmin = 0, xmax = nrow(dose_ord)+1, ymin = 0.8, ymax = 1.2,
           alpha = .12, fill = col_vacc) +
  annotate("text", x = 1, y = 1.26, hjust = 0,                
           label = "Acceptable range (\u00B120%)", colour = col_vacc) +
  geom_col(width = 0.75) + geom_hline(yintercept = 1, colour = col_grey) +
  scale_fill_manual(values = c(`TRUE` = col_vacc, `FALSE` = col_reject), guide = "none") +
  annotate("text", x = nrow(dose_ord)*0.5, y = 1.55,
           label = sprintf("%d/%d (%.0f%%) within range; median %.2f",
                           sum(dose$in_range), nrow(dose), 100*mean(dose$in_range),
                           median(dose$dose_ratio))) +
  labs(x = "Vaccinated farm-cycles, ordered by dose ratio",
       y = "Delivered dose / protocol target")

ggsave("output/Figure6_dose_adherence.png", fig6, width = 9, height = 4.5, dpi = 150)

# cumulative mortality, bird- vs farm-weighted --
cum <- visit |> arrange(farm_id, study_arm, period, day) |>
  group_by(farm_id, study_arm, period) |>
  mutate(cum_deaths = cumsum(replace_na(dead_birds_count, 0))) |> ungroup() |>
  left_join(placed, by = c("farm_id","study_arm","period")) |>
  filter(day %in% c(3,5,14,21,28)) |> mutate(cmr = 100 * cum_deaths / birds_placed)

s1_bird <- cum |> 
  group_by(study_arm, day) |>
  summarise(m = 100*sum(cum_deaths)/sum(birds_placed), .groups="drop") |> 
  mutate(w="Bird-weighted (descriptive)")

s1_farm <- cum |> 
  group_by(study_arm, day) |> 
  summarise(m = mean(cmr), .groups="drop") |>
  mutate(w="Farm-weighted (design-consistent)")

figS1 <- ggplot(bind_rows(s1_bird, s1_farm), aes(day, m, colour = study_arm, shape = study_arm)) +
  annotate("rect", xmin=0, xmax=14, ymin=0, ymax=Inf, alpha=.07, fill=col_vacc) +
  annotate("text", x = 7, y = Inf, vjust = 1.4, label = "Restriction\nwindow",   
           colour = col_vacc, size = 3.3, lineheight = 0.9) +
  geom_vline(xintercept = 21, linetype = "dashed", colour = col_grey) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) + facet_wrap(~w) +
  scale_colour_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc)) +
  scale_x_continuous(breaks = c(3,5,14,21,28)) +
  labs(x = "Study day", y = "Cumulative mortality (% of birds placed)", colour=NULL, shape=NULL)

ggsave("output/FigureS1_cumulative_mortality.png", figS1, width = 11, height = 4.2, dpi = 150)

# paired within-farm day-21 mortality -----------
s2 <- fc |>
  dplyr::select(farm_id, study_arm, mort21) |>
  pivot_wider(names_from = study_arm, values_from = mort21) |>
  arrange(desc(Control)) |> mutate(idx = row_number())

# farms to label: the few highest by either cycle
s2_lab <- s2 |>
  mutate(ymax = pmax(Control, Vaccinated, na.rm = TRUE)) |>
  slice_max(ymax, n = 4)            

figS2 <- ggplot(s2) +
  geom_segment(aes(x = idx, xend = idx, y = Control, yend = Vaccinated), colour = "grey75") +
  geom_point(aes(idx, Control, colour = "Control cycle"), size = 2.5) +
  geom_point(aes(idx, Vaccinated, colour = "Vaccinated cycle"), size = 2.5) +
  geom_text(data = s2_lab, aes(idx, ymax, label = farm_id),
            vjust = -0.6, size = 3.2, colour = "grey30") +
  scale_colour_manual(values = c("Control cycle" = col_ctrl, "Vaccinated cycle" = col_vacc), name = NULL) +
  labs(x = "Farms (n = 44), ordered by control-cycle mortality", y = "Cumulative mortality to day 21 (%)") +
  theme(legend.position = c(.85,.85))

ggsave("output/FigureS2_paired_mortality.png", figS2, width = 11, height = 4.6, dpi = 150)

# leave-one-farm-out crude RR -------------------
overall_v <- sum(fc$deaths21[fc$study_arm=="Vaccinated"]); overall_pv <- sum(fc$birds_placed[fc$study_arm=="Vaccinated"])
overall_c <- sum(fc$deaths21[fc$study_arm=="Control"]);    overall_pc <- sum(fc$birds_placed[fc$study_arm=="Control"])
loo <- fc |> group_by(farm_id) |>
  summarise(dv = sum(deaths21[study_arm=="Vaccinated"]), pv = sum(birds_placed[study_arm=="Vaccinated"]),
            dc = sum(deaths21[study_arm=="Control"]),    pc = sum(birds_placed[study_arm=="Control"]), .groups="drop") |>
  mutate(rr_without = ((overall_v-dv)/(overall_pv-pv)) / ((overall_c-dc)/(overall_pc-pc))) |>
  arrange(rr_without) |> mutate(idx = row_number())

rr_all <- (overall_v/overall_pv)/(overall_c/overall_pc)

loo_lab <- loo |> filter(abs(rr_without - rr_all) > 0.02) 

figS3 <- ggplot(loo, aes(idx, rr_without, fill = abs(rr_without - rr_all) > 0.02)) +
  geom_col() +
  geom_hline(yintercept = rr_all, linetype = "dashed") +
  geom_text(data = loo_lab, aes(idx, rr_without, label = farm_id),  
            vjust = -0.5, colour = col_reject, size = 3.2) +
  scale_fill_manual(values = c(`FALSE` = col_grey, `TRUE` = col_reject), guide = "none") +
  annotate("text", x = 8, y = rr_all + 0.004, label = sprintf("All farms: RR = %.3f", rr_all)) +
  coord_cartesian(ylim = c(0.65, 0.855)) +                              
  scale_y_continuous(breaks = seq(0.65, 0.85, 0.025)) +
  labs(x = "Farm omitted (n = 44)", y = "Crude bird-weighted RR with that farm omitted")

ggsave("output/FigureS3_leave_one_out.png", figS3, width = 11, height = 4.4, dpi = 150)

# mean body weight by day & arm -----------------
s4 <- wt_cycle_day |> filter(day %in% c(3,14,21,28)) |> group_by(study_arm, day) |>
  summarise(m = mean(mean_wt), se = sd(mean_wt)/sqrt(n()), .groups="drop")

figS4 <- ggplot(s4, aes(day, m, colour = study_arm, shape = study_arm)) +
  annotate("rect", xmin=0, xmax=14, ymin=-Inf, ymax=Inf, alpha=.07, fill=col_vacc) +
  annotate("text", x = 10, y = Inf, vjust = 1.4, label = "Restriction\nwindow",   
           colour = col_vacc, size = 3.3, lineheight = 0.9) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = m-1.96*se, ymax = m+1.96*se), width = 0.2) +
  scale_colour_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc)) +
  scale_x_continuous(breaks = c(3,14,21,28)) +
  labs(x = "Study day", y = "Mean body weight (g)", colour=NULL, shape=NULL) +
  theme(legend.position = c(.15,.85))

ggsave("output/FigureS4_body_weight.png", figS4, width = 9, height = 4.5, dpi = 150)


# 5. NEW additions 

# 1: Antibiotic use by active ingredient & WHO class -------------
amu_class <- visit |> filter(amu_event) |>
  mutate(class = case_when(
    str_detect(abx_text, poly_stem)  ~ "Polymyxin (colistin) [HPCIA]",
    str_detect(abx_text, fq_stem) ~ "Fluoroquinolone [HPCIA]",
    str_detect(abx_text, macro_stem) ~ "Macrolide [HPCIA]",
    str_detect(abx_text, "sulfa|sulph|trimethoprim|cotrima") ~ "Sulfonamide/trimethoprim",
    str_detect(abx_text, "amoxic|amoxyc") ~ "Aminopenicillin",
    str_detect(abx_text, "neomyc|gentam") ~ "Aminoglycoside",
    str_detect(abx_text, "oxytet|tetracyc|chlortetra|doxycyc|doxyclin|\\botc\\b") ~ "Tetracycline",
    TRUE ~ "Other/unclassified")) |>
  count(study_arm, class) |> pivot_wider(names_from = study_arm, values_from = n, values_fill = 0)

print(as.data.frame(amu_class), row.names = FALSE)

n1 <- visit |> filter(amu_event) |>
  mutate(class = case_when(
    str_detect(abx_text, poly_stem) ~ "Polymyxin",
    str_detect(abx_text, fq_stem) ~ "Fluoroquinolone",
    str_detect(abx_text, macro_stem) ~ "Macrolide",
    str_detect(abx_text, "sulfa|sulph|trimethoprim|cotrima") ~ "Sulfonamide/TMP",
    str_detect(abx_text, "amoxic|amoxyc") ~ "Aminopenicillin",
    str_detect(abx_text, "neomyc|gentam") ~ "Aminoglycoside",
    TRUE ~ "Tetracycline")) |> count(study_arm, class)
figN1 <- ggplot(n1, aes(reorder(class, n), n, fill = study_arm)) +
  geom_col(position = "dodge") + coord_flip() +
  scale_fill_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc), name = NULL) +
  labs(x = NULL, y = "Antibiotic-use events", title = "Antibiotic events by WHO class and arm")
ggsave("output/am_by_class.png", figN1, width = 8.5, height = 4.5, dpi = 150)

# 2: feed intake (g/bird/day) by arm -----------------------------
feed <- visit |>
  mutate(feed_g_bird = 1000 * as.numeric(feed_quantity_24hrs_kg) / as.numeric(current_bird_count)) |>
  filter(is.finite(feed_g_bird), feed_g_bird > 0, day > 0)

feed_cycle <- feed |> group_by(farm_id, study_arm, period, window) |>
  summarise(feed_g_bird = mean(feed_g_bird), .groups = "drop")

m_feed <- lmer(feed_g_bird ~ study_arm + period + window + (1|farm_id), data = feed_cycle)

print(round(summary(m_feed)$coefficients["study_armVaccinated", ], 3))

figN2 <- ggplot(feed |> filter(day <= 28), aes(factor(day), feed_g_bird, fill = study_arm)) +
  geom_boxplot(outlier.size = .6) +
  scale_fill_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc), name = NULL) +
  labs(x = "Study day", y = "Feed intake (g/bird/day)", title = "Daily feed intake by arm")

ggsave("output/feed_intake.png", figN2, width = 9, height = 4.5, dpi = 150)

# 3: Harvest weight & offtake (production endpoint) ---------------
harvest <- visit |> filter(harvest_started == "yes") |>
  group_by(farm_id, study_arm, period) |>
  summarise(
    harvest_wt = if (all(is.na(average_harvest_weight))) NA_real_
    else max(as.numeric(average_harvest_weight), na.rm = TRUE),
    birds_sold = if (all(is.na(birds_sold))) NA_real_
    else max(as.numeric(birds_sold), na.rm = TRUE),
    .groups = "drop") |>
  filter(is.finite(harvest_wt), harvest_wt > 0)

m_harv <- lmer(harvest_wt ~ study_arm + period + (1|farm_id), data = harvest)
print(round(summary(m_harv)$coefficients["study_armVaccinated", ], 2))

# 4: Flock weight uniformity (CV) by day & arm --------------------
# CV of individual bird weights per farm-cycle-day: a subclinical-stress marker the farm-cycle-mean growth analysis cannot see.
unif <- wt_cycle_day |> filter(day %in% c(3,14,21,28), n_birds >= 5)
m_cv <- lmer(cv_wt ~ study_arm + period + factor(day) + (1|farm_id), 
             data = unif)

print(round(summary(m_cv)$coefficients["study_armVaccinated", ], 3))

figN4 <- ggplot(unif, aes(factor(day), cv_wt, fill = study_arm)) +
  geom_boxplot(outlier.size = .6) +
  scale_fill_manual(values = c(Control = col_ctrl, Vaccinated = col_vacc), name = NULL) +
  labs(x = "Study day", y = "Within-flock weight CV (%)", title = "Flock uniformity by arm")

ggsave("output/uniformity.png", figN4, width = 9, height = 4.5, dpi = 150)

# 5: symptom-specific clinical signs (within-farm McNemar) --------
sign_cols <- c(bloody_diarrhoea = "Je_Dalili_za_magonj_huona_kwa/kuharisha_damu",
               respiratory      = "Je_Dalili_za_magonj_huona_kwa/kupumua_kwa_shida",
               stunting         = "Je_Dalili_za_magonj_huona_kwa/kudumaa",
               weakness         = "Je_Dalili_za_magonj_huona_kwa/udhaifu",
               ruffled_feathers = "Je_Dalili_za_magonj_huona_kwa/kususha_mbwa__kuvaa_koti",
               swollen_eyes     = "Je_Dalili_za_magonj_huona_kwa/kuvimba_kwa_macho",
               sudden_death     = "Je_Dalili_za_magonj_huona_kwa/vifo_vya_gafla",
               white_diarrhoea  = "Je_Dalili_za_magonj_huona_kwa/kuharisha_nyeupe",
               green_diarrhoea  = "Je_Dalili_za_magonj_huona_kwa/kuharisha_kijani")
sign_res <- data.frame()
for (nm in names(sign_cols)) {
  col <- sign_cols[[nm]]
  tmp <- visit |> mutate(flag = as.integer(.data[[col]] == 1 | .data[[col]] == "1" | .data[[col]] == "yes")) |>
    mutate(flag = replace_na(flag, 0L)) |>
    group_by(farm_id, study_arm) |> summarise(any_sign = as.integer(any(flag == 1)), .groups = "drop") |>
    pivot_wider(names_from = study_arm, values_from = any_sign)
  if (all(c("Control","Vaccinated") %in% names(tmp))) {
    b <- sum(tmp$Vaccinated == 1 & tmp$Control == 0)   
    c_ <- sum(tmp$Vaccinated == 0 & tmp$Control == 1)
    p  <- tryCatch(mcnemar.test(matrix(c(0,b,c_,0),2))$p.value, error=function(e) NA)
    p_exact <- tryCatch(binom.test(b, b + c_)$p.value, error = function(e) NA)
    sign_res <- rbind(sign_res, data.frame(sign = nm, vacc_only = b, control_only = c_,
                                           p_mcnemar = round(p,3), p_exact = round(p_exact,3)))
  }
}
print(sign_res, row.names = FALSE)

# 6: Spatial autocorrelation of farm mortality (Moran's I) --------
visit <- visit |> mutate(lat = ifelse(lat < -3.6 | lat > -3.1, NA, lat),
                         lon = ifelse(lon < 36.6 | lon > 37.5, NA, lon))

cat("Stray GPS fixes removed:", sum(is.na(visit$lat)), "\n")

mort_farm <- fc |> group_by(farm_id) |>
  summarise(mort = 100*sum(deaths21)/sum(birds_placed), .groups = "drop")

farm_geo <- visit |> group_by(farm_id) |>
  summarise(lat = median(lat, na.rm = TRUE),
            lon = median(lon, na.rm = TRUE), .groups = "drop") |>
  inner_join(mort_farm, by = "farm_id") |>
  filter(is.finite(lat), is.finite(lon))

# shapefile
tz <- st_read("data/gadm41_TZA_2.shp", quiet = TRUE)

# filter to Moshi and Arusha  
study  <- tz |> 
  dplyr::filter(NAME_1 %in% c("Arusha", "Kilimanjaro"))

# sf object
farm_sf <- st_as_sf(farm_geo, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(st_crs(tz))

# panel A: locator
pA <- ggplot() +
  geom_sf(data = tz, fill = "grey97", colour = "black", linewidth = 0.15) +
  geom_sf(data = study, fill = "grey88", colour = "black", linewidth = 0.3) +
  geom_sf(data = farm_sf, fill = "#C0392B", shape = 21, colour = "white",
          size = 1.6, stroke = 0.2) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.6) +
  coord_sf(expand = FALSE) +
  labs(tag = "A") + theme_void() +
  theme(panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4))

# panel B: zoomed to the farms
pB <- ggplot() +
  geom_sf(data = study, fill = "white", colour = "black", linewidth = 0.3) +
  geom_sf_text(data = study, aes(label = NAME_2), size = 2.6,
               colour = "black", fontface = "italic", check_overlap = TRUE) +
  geom_sf(data = farm_sf, aes(fill = mort, size = mort),
          shape = 21, colour = "white", stroke = 0.4, alpha = 0.9) +
  scale_fill_viridis_c(option = "C", name = "Day-21\nmortality (%)") +
  scale_size(range = c(2, 8), guide = "none") +
  annotation_scale(location = "bl", width_hint = 0.28,
                   line_width = 0.5, text_cex = 0.7) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1.0, "cm"), width = unit(0.9, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) +    
  labs(title = "Farm-level day-21 mortality in women-owned smallholder broiler farms",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        panel.grid       = element_line(colour = "white", linewidth = 0.25),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.title       = element_text(face = "bold", size = 12))



pacman::p_load(cowplot)
figN6 <- plot_grid(pA, pB, labels = c("A", "B"), rel_widths = c(1, 1.4))

ggsave("output/spatial_mortality.png", 
       figN6, width = 11,
       height = 5.5,
       dpi = 300, 
       bg = "white")


# moran's I on the coordinates
farm_m <- st_transform(farm_sf, 32737)
dmat   <- as.matrix(dist(st_coordinates(farm_m)))     
winv   <- 1 / dmat; diag(winv) <- 0
winv[!is.finite(winv)] <- max(winv[is.finite(winv)])
moran  <- ape::Moran.I(farm_geo$mort, winv)
cat(sprintf("Moran's I = %.3f (expected %.3f), p = %.3f  [n = %d farms]\n",
            moran$observed, moran$expected, moran$p.value, nrow(farm_geo)))


# 7: Farmer perception vs measured within-farm effect ------------
# Do farmers who PERCEIVED fewer deaths / better growth in the vaccinated cycle actually have a measured within-farm advantage? (Crossover -> paired.)
measured <- fc |> 
  dplyr::select(farm_id, study_arm, mort21) |> 
  pivot_wider(names_from = study_arm, values_from = mort21) |>
  dplyr::mutate(mort_diff = Vaccinated - Control) |>          
  left_join(adwg |> 
              dplyr::select(farm_id, study_arm, adwg_3_21) |>
              pivot_wider(names_from = study_arm, values_from = adwg_3_21) |>
              mutate(adwg_diff = Vaccinated - Control) |> 
              dplyr::select(farm_id, adwg_diff),
            by = "farm_id")

perc <- farms |> transmute(farm_id = tolower(as.character(farm_id)),
                            perc_fewer_deaths = f1_lower_mortality,
                            perc_better_growth = f2_birds_grew_better)
pm <- measured |> inner_join(perc, by = "farm_id")

if (nrow(pm) > 10) {
  cat("Mean measured mortality difference (vacc-ctrl) by perceived-fewer-deaths group:\n")
  print(pm |> group_by(perc_fewer_deaths) |>
          summarise(n = n(), mean_mort_diff = round(mean(mort_diff, na.rm=TRUE),2),
                    mean_adwg_diff = round(mean(adwg_diff, na.rm=TRUE),2)) |> as.data.frame())
  pm2 <- pm |> mutate(perc = as.integer(perc_fewer_deaths %in% c("yes","Yes",1,"1")))
  if (length(unique(pm2$perc)) > 1)
    cat("Wilcoxon (measured mort diff ~ perceived fewer deaths) p =",
        round(wilcox.test(mort_diff ~ perc, data = pm2)$p.value, 3), "\n")
}

# 8: RA recording effect -----------
ra_cycle <- visit |> group_by(farm_id, study_arm, period) |>
  summarise(amu = sum(amu_event), nv = n(),
            ra = names(sort(table(ra_name), decreasing = TRUE))[1], .groups = "drop")
m_ra <- glmer(amu ~ study_arm + period + (1|farm_id) + (1|ra), data = ra_cycle,
              family = poisson, offset = log(nv))
cat("AMU IRR with RA random effect =", round(exp(fixef(m_ra)["study_armVaccinated"]), 2), "\n")
vc_ra <- as.data.frame(VarCorr(m_ra))[, c("grp", "vcov")]
vc_ra$vcov <- round(vc_ra$vcov, 3)
