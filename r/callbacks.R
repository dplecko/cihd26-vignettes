
#' * callbacks for Vignette II: Heterogeneity *

oxygen_miiv_cb <- function(...) {
  
  res <- merge(list(...)[[1]], list(...)[[2]], all = TRUE)
  res <- merge(res, list(...)[[3]], all = TRUE)
  
  # inferring oxygenation
  res[fio2 > 21, oxy_ind := 1]
  res[!is.na(o2flow), oxy_ind := 1]
  res[o2device != "None", oxy_ind := 1]
  res[o2device == "None", oxy_ind := 0]
  
  res[, c(meta_vars(res), "oxy_ind"), with=FALSE]
}

oxygen_impute_miiv_cb <- function(oxy_ind, ...) {
  
  res <- fill_gaps(oxy_ind)
  res[, lag_oxy_ind := oxy_ind]
  for (tlag in seq_len(4)) {
    
    res[, lag_oxy_ind := shift(lag_oxy_ind), by = c(id_vars(res))]
    res[is.na(oxy_ind) & !is.na(lag_oxy_ind), oxy_ind := lag_oxy_ind]
  }
  
  res[, oxy_imp := oxy_ind]
  res[is.na(oxy_imp), oxy_imp := 0]
  res[, c(meta_vars(res), "oxy_imp"), with=FALSE]
}

#' * callbacks for Vignette III: Race & Mortality in ICU *
acute_dayone <- function(sofa, ...) {
  
  ind_var <- index_var(sofa)
  sofa <- sofa[get(ind_var) == hours(24L)]
  sofa[, acu_24 := sofa]
  sofa[, c(ind_var, "sofa") := NULL]
  sofa
}

lact_dayone <- function(...) {
  
  tbl <- list(...)[[1]]
  
  ind_var <- index_var(tbl)
  val_var <- names(list(...)[1])
  tbl <- tbl[get(ind_var) >- hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(lact_24 = max(lact)), by = c(id_vars(tbl))]
}

ast_dayone <- function(...) {
  
  tbl <- list(...)[[1]]
  
  ind_var <- index_var(tbl)
  val_var <- names(list(...)[1])
  tbl <- tbl[get(ind_var) >- hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(ast_24 = max(ast)), by = c(id_vars(tbl))]
}

pafi_dayone <- function(...) {
  
  tbl <- list(...)[[1]]
  
  ind_var <- index_var(tbl)
  val_var <- names(list(...)[1])
  tbl <- tbl[get(ind_var) >- hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(pafi_24 = min(pafi)), by = c(id_vars(tbl))]
}

safi_dayone <- function(...) {
  
  tbl <- list(...)[[1]]
  
  ind_var <- index_var(tbl)
  val_var <- names(list(...)[1])
  tbl <- tbl[get(ind_var) >- hours(0L) & get(ind_var) <= hours(24L)]
  tbl[, list(safi_24 = min(safi)), by = c(id_vars(tbl))]
}

miiv_adm_epi_cb <- function(x, val_var, extra_var, env) {
  
  x <- merge(x, env$icustays[, c("stay_id", "intime")],
             by = "stay_id")
  x <- as_id_tbl(x, id_vars = "subject_id")
  x <- data.table::setorderv(x, cols = c("subject_id", "intime"))
  x[, adm_episode := seq_along(stay_id), by = "subject_id"]
  x <- as_id_tbl(x, id_vars = "stay_id")
  x[, subject_id := adm_episode]
}

race_mim_callback <- function(x, val_var, env) {
  
  groups <- list(
    Caucasian = c("WHITE", "WHITE - BRAZILIAN", "WHITE - EASTERN EUROPEAN",
                  "WHITE - OTHER EUROPEAN", "WHITE - RUSSIAN"),
    Asian = c("ASIAN", "ASIAN - ASIAN INDIAN", "ASIAN - CAMBODIAN",
              "ASIAN - CHINESE", "ASIAN - FILIPINO", "ASIAN - JAPANESE",
              "ASIAN - KOREAN", "ASIAN - OTHER", "ASIAN - THAI",
              "ASIAN - VIETNAMESE"),
    Hispanic = c("HISPANIC/LATINO - CENTRAL AMERICAN (OTHER)",
                 "HISPANIC/LATINO - COLOMBIAN",
                 "HISPANIC/LATINO - CUBAN", "HISPANIC/LATINO - DOMINICAN",
                 "HISPANIC/LATINO - GUATEMALAN",
                 "HISPANIC/LATINO - HONDURAN", "HISPANIC/LATINO - MEXICAN",
                 "HISPANIC/LATINO - PUERTO RICAN",
                 "HISPANIC/LATINO - SALVADORAN", "HISPANIC OR LATINO"),
    `African American` = c("BLACK/AFRICAN AMERICAN"),
    Other = c("AMERICAN INDIAN/ALASKA NATIVE FEDERALLY RECOGNIZED TRIBE",
              "UNABLE TO OBTAIN", "UNKNOWN/NOT SPECIFIED", "MIDDLE EASTERN",
              "MULTI RACE ETHNICITY",
              "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "OTHER",
              "PATIENT DECLINED TO ANSWER",
              "PORTUGUESE", "SOUTH AMERICAN",
              "AMERICAN INDIAN/ALASKA NATIVE",
              "AMERICAN INDIAN/ALASKA NATIVE FEDERALLY RECOGNIZED TRIBE",
              "CARIBBEAN ISLAND", "BLACK/AFRICAN", "BLACK/CAPE VERDEAN",
              "BLACK/HAITIAN")
  )
  map <- unlist(groups)
  names(map) <- rep(names(groups), times = lapply(groups, length))

  x[, ethnicity := names(map)[match(ethnicity, map)]]
}

race_miiv_callback <- function(x, val_var, env) {
  
  groups <- list(
    Caucasian = c("WHITE", "WHITE - BRAZILIAN", "WHITE - EASTERN EUROPEAN",
                  "WHITE - OTHER EUROPEAN", "WHITE - RUSSIAN"),
    Asian = c("ASIAN", "ASIAN - ASIAN INDIAN", "ASIAN - CHINESE",
              "ASIAN - KOREAN", "ASIAN - SOUTH EAST ASIAN"),
    Hispanic = c("HISPANIC OR LATINO", "HISPANIC/LATINO",
                 "HISPANIC/LATINO - CENTRAL AMERICAN",
                 "HISPANIC/LATINO - COLUMBIAN", "HISPANIC/LATINO - CUBAN",
                 "HISPANIC/LATINO - DOMINICAN", "HISPANIC/LATINO - GUATEMALAN",
                 "HISPANIC/LATINO - HONDURAN", "HISPANIC/LATINO - MEXICAN",
                 "HISPANIC/LATINO - PUERTO RICAN",
                 "HISPANIC/LATINO - SALVADORAN", "SOUTH AMERICAN"),
    `African American` = c("BLACK/AFRICAN", "BLACK/AFRICAN AMERICAN",
                           "BLACK/CAPE VERDEAN", "BLACK/CARIBBEAN ISLAND"),
    Other = c("OTHER", "UNKNOWN", "UNABLE TO OBTAIN", "MULTIPLE RACE/ETHNICITY",
              "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
              "PATIENT DECLINED TO ANSWER", "PORTUGUESE")
  )
  
  map <- unlist(groups)
  names(map) <- rep(names(groups), times = lapply(groups, length))
  
  x[, race := names(map)[match(race, map)]]
}

mimic_adm_epi_cb <- function(x, ...) {
  
  x <- merge(x, list(...)$env$icustays[, c("icustay_id", "intime")], 
             by = "icustay_id")
  x <- setorderv(x, cols = c("subject_id", "intime"))
  x[, adm_episode := seq_len(.N), by = "subject_id"]
  
  x <- x[, c(id_vars(x), "adm_episode"), with=FALSE]
  rename_cols(x, "subject_id", "adm_episode")
}

mimic_adm_diag <- function(x, val_var, ...) {
  
  mapp <- list(
    OTH = c("DENT", "PSYCH", "NBB", "NB", "ENT", "GU", "PSURG"),
    GYN = c("GYN", "OBS")
  )
  for (i in seq_along(mapp)) {
    x[get(val_var) %in% mapp[[i]], c(val_var) := names(mapp)[i]]
  }
  x
}

miiv_elixhauser_cb <- function(x, ...) {
  
  weights <- TRUE
  
  ch9 <- my_icd9_comorbid_elix(x[icd_version == 9])
  ch10 <- my_icd10_comorbid_elix(x[icd_version == 10])
  
  w_all <- rbind(ch9, ch10)
  w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
  
  w_all <- data.table::melt(w_all, id.vars = c(id_vars(w_all)))
  w_all <- as_id_tbl(w_all)
  if (weights) {
    
    wgh_vec <- c(
      CHF           = 7,
      Arrhythmia    = 5,
      Valvular      = 0,
      PHTN          = 0,
      PVD           = 2,
      HTN           = -1,
      HTNcx         = 1,
      Paralysis     = 5,
      NeuroOther    = 6,
      Pulmonary     = 3,
      DM            = 0,
      DMcx          = 2,
      Hypothyroid   = -2,
      Renal         = 5,
      Liver         = 11,
      PUD           = 0,
      HIV           = 0,
      Lymphoma      = 9,
      Mets          = 14,
      Tumor         = 4,
      Rheumatic     = 0,
      Coagulopathy  = 3,
      Obesity       = -5,
      WeightLoss    = 6,
      FluidsLytes   = 5,
      BloodLoss     = 0,
      Anemia        = -2,
      Alcohol       = 0,
      Drugs         = 0,
      Psychoses     = 0,
      Depression    = -3
    )
    
    elix_weights <- data.table(
      variable = names(wgh_vec),
      weight = as.integer(wgh_vec)
    )
    
    
    w_all <- merge(w_all, elix_weights, by = "variable")
    
  } else w_all[, weight := 1]
  
  # save(w_all, file = file.path("miiv_elix.rda"))
  w_all[, list(icd_code = sum(value * weight)), by = c(id_vars(w_all))]
}

miiv_charlson_cb <- function(x, ...) {
  
  weights <- TRUE
  
  ch9 <- my_icd9_comorbid_charlson(x[icd_version == 9])
  ch10 <- my_icd10_comorbid_charlson(x[icd_version == 10])
  
  w_all <- rbind(ch9, ch10)
  w_all <- w_all[, lapply(.SD, max), by = c(id_vars(w_all))]
  
  w_all <- data.table::melt(w_all, id.vars = c(id_vars(w_all)))
  w_all <- as_id_tbl(w_all)
  if (weights) {
    
    charlson_weights <- data.table::data.table(
      variable = c(
        "MI", "CHF", "PVD", "Stroke", "Dementia", "Pulmonary", "Rheumatic",
        "PUD", "LiverMild", "LiverSevere", "DM", "DMcx", "Paralysis",
        "Renal", "Cancer", "Mets", "HIV"
      ),
      weight = c(
        1, 1, 1, 1, 1, 1, 1,
        1, 1, 2, 1, 1, 2,
        2, 2, 4, 6
      )
    )
    
    w_all <- merge(w_all, charlson_weights, by = "variable")
    
  } else w_all[, weight := 1]
  
  
  save(w_all, file = file.path("miiv_charlson.rda"))
  w_all[, list(icd_code = sum(value * weight)), by = c(id_vars(w_all))]
}

miiv_elective <- function(x, ...) {
  
  x[, elective := fifelse(
    admission_type %in% c("AMBULATORY OBSERVATION", "ELECTIVE") |
      (admission_type == "DIRECT OBSERVATION" & admission_location %in% c("PROCEDURE SITE", "PACU", "AMBULATORY SURGERY TRANSFER")) |
      (admission_type == "OBSERVATION ADMIT" & admission_location %in% c("PROCEDURE_SITE", "PACU")) |
      (admission_type == "SURGICAL SAME DAY ADMISSION" & admission_location %in% c("PROCEDURE_SITE", "PACU", "PHYSICIAN REFERRAL")) ,
    1, # Elective (coded as 1)
    0  # Non-elective (coded as 0)
  )]
  
  x[, admission_type := elective]
}

mimic_elective <- function(x, ...) {
  
  x[, elective := admission_type == "ELECTIVE"]
  x[, admission_type := elective]
}

#' * callbacks for Vignette V: Obesity Paradox Measurement Error *
miiv_omr_bmi_cb <- function(x, ...) {

  x[, result_value := as.numeric(result_value)]
  x[, bmi_omr := mean(result_value, na.rm = TRUE), by = c(id_vars(x))]
  
  x
}

bmi_all_cb <- function(bmi, bmi_omr, ...) {
  
  bmi <- rename_cols(bmi, "bmi_all", "bmi")
  bmi_omr <- rename_cols(bmi_omr, "bmi_all", "bmi_omr")
  
  x <- rbind(bmi, bmi_omr)
  x[, list(bmi_all = mean(bmi_all)), by = c(id_vars(x))]
}

miiv_elix_by <- function(...) {
  
  idvar <- id_var(list(...)[[1]])
  x <- load_src("diagnoses_icd", "miiv")
  x <- as_id_tbl(x[, c(idvar, "icd_code", "icd_version"), with=FALSE])
  ch9 <- my_icd9_comorbid_elix(x[icd_version == 9])
  ch10 <- my_icd10_comorbid_elix(x[icd_version == 10])
  res <- rbind(ch9, ch10)
  res[, elix_components := rowSums(res[, -1])]
  res
}

#' * callbacks for Vignette VI: Opioids -> Delirium *
ts_to_win_12hours <- function(x, dur_var, ...) {
  
  x[, c(list(...)$val_var) := NULL]
  x[, c(list(...)$val_var) := TRUE]
  x[, c(dur_var) := NULL]
  x[, c(dur_var) := mins(720L)]
  
  as_win_tbl(x, dur_var = dur_var, by_ref = TRUE)
}

mimic_notes_cb <- function(x, grp1_var, grp2_var, index2_var, val_var, ...) {
  
  # for Echo, ECG, charttime inherits chartdate + hours(24L)
  x <- x[category %in% c("ECG", "Echo", "Discharge summary"),
         c(index2_var) := get(index_var(x)) + hours(24L)]
  
  # chartdate gets charttime
  x[, c(index_var(x)) := get(index2_var)]
  
  # x <- x[category != "Discharge summary"] # remove all discharge summaries
  
  x <- x[!(category %in% c("ECG", "Echo", "Radiology"))]
  
  # collapse notes to avoid losing information
  x[, list(text = paste(text, collapse = "\n")), 
    by = c(meta_vars(x))]
}
