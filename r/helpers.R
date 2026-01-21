
mean_or_na <- function (x) {
  if (all(is.na(x))) 
    return(x[1L])
  res <- mean(x, na.rm = TRUE)
  if (is.na(res)) 
    x[1L]
  else res
}

rjson <- function(path) {
  
  jsonlite::read_json(path, simplifyVector = TRUE,
                      simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
}

wjson <- function(x, path) {
  
  jsonlite::write_json(x, path, simplifyVector = TRUE,
                       simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
}

norm_code <- function(dt) {
  dt[, icd_code := toupper(gsub("\\.", "", icd_code))]
}

# Compute comorbidity flags by stay_id
compute_flags <- function(dt, map) {
  out_list <- vector("list", length(map))
  names(out_list) <- names(map)
  
  for (cat in names(map)) {
    codes <- map[[cat]]
    tmp <- dt[icd_code %in% codes, .(flag = 1L), by = c(id_vars(dt))]
    if (nrow(tmp) == 0L) {
      out_list[[cat]] <- unique(dt[, c(id_vars(dt)), with=FALSE])[, (cat) := 0L]
    } else {
      setnames(tmp, "flag", cat)
      out_list[[cat]] <- tmp
    }
  }
  
  res <- Reduce(function(a, b) merge(a, b, by = id_vars(dt), all = TRUE), out_list)
  for (j in setdiff(names(res), id_vars(dt))) {
    set(res, which(is.na(res[[j]])), j, 0L)
  }
  
  res
}

# ICD-9 wrapper
my_icd9_comorbid_elix <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  root <- rprojroot::find_root(rprojroot::has_dir(".git"))
  compute_flags(dt, rjson(file.path(root, "config/icd9-elix.json")))
}

# ICD-10 wrapper
my_icd10_comorbid_elix <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  root <- rprojroot::find_root(rprojroot::has_dir(".git"))
  compute_flags(dt, rjson(file.path(root, "config/icd10-elix.json")))
}

my_icd9_comorbid_charlson <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  compute_flags(dt, rjson("config/icd9-charlson.json"))
}

# ICD-10 wrapper
my_icd10_comorbid_charlson <- function(dt) {
  dt <- copy(dt)
  norm_code(dt)
  compute_flags(dt, rjson("config/icd10-charlson.json"))
}
