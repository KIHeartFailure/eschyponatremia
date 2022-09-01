

# Impute missing values ---------------------------------------------------

edataforimp <- edata %>%
  select(patientid, natremia, !!!syms(coxvars), !!!syms(logvars), !!!syms(predvars), starts_with("out"), d_hsSod_cat, d_dcSod_cat, num_dcVital)

noimpvars <- names(edataforimp)[!names(edataforimp) %in% c(coxvars, logvars, predvars)]

# Nelson-Aalen estimator
na <- basehaz(coxph(Surv(outtime_hosphf, out_deathhf == 1) ~ 1,
  data = edataforimp, method = "breslow"
))
edataforimp <- left_join(edataforimp, na, by = c("outtime_hosphf" = "time"))


ini <- mice(edataforimp, maxit = 0, print = F, m = 1)

pred <- ini$pred
pred[, noimpvars] <- 0
pred[noimpvars, ] <- 0 # redundant

# change method used in imputation to prop odds model
meth <- ini$method
meth[noimpvars] <- ""

## check no cores
cores_2_use <- detectCores() - 1
if (cores_2_use >= 10) {
  cores_2_use <- 10
  m_2_use <- 1
} else if (cores_2_use >= 5) {
  cores_2_use <- 5
  m_2_use <- 2
} else {
  stop("Need >= 5 cores for this computation")
}

cl <- makeCluster(cores_2_use)
clusterSetRNGStream(cl, 49956)
registerDoParallel(cl)

imp <-
  foreach(
    no = 1:cores_2_use,
    .combine = ibind,
    .export = c("meth", "pred", "edataforimp"),
    .packages = "mice"
  ) %dopar% {
    mice(edataforimp,
      m = m_2_use, maxit = 10, method = meth,
      predictorMatrix = pred,
      printFlag = FALSE
    )
  }
stopImplicitCluster()


# keep original imputed data (rsdataimp) just in case
imp.org <- imp

# Convert to Long
long <- mice::complete(imp, action = "long", include = TRUE)

long <- long %>%
  mutate(
    d_age_cat = factor(case_when(
      num_age < 65 ~ 1,
      num_age >= 65 ~ 2
    ), levels = 1:2, labels = c("<65", ">=65"))
  )

# Convert back to mids object
imp <- as.mids(long)
