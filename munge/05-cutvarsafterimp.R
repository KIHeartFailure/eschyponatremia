
# Creating variables for competing risk analysis in imputed dataset ------------

# keep org imp
imp.org <- imp

# Convert to Long
long <- mice::complete(imp, action = "long", include = TRUE)

## Create numeric variables needed for comp risk model
long <- long %>%
  mutate(
    d_dcEf_cat = factor(case_when(
      num_dcEf < 40 ~ 1,
      num_dcEf <= 49 ~ 2,
      num_dcEf >= 50 ~ 3
    ),
    levels = 1:3,
    labels = c("<40", "40-49", ">=50")
    ),
    d_age_cat = factor(case_when(
      num_age < 65 ~ 1,
      num_age >= 65 ~ 2
    ), levels = 1:2, labels = c("<65", ">=65")),
    d_dcCKDEPI_cat = factor(case_when(
      d_dcCKDEPI < 60 ~ 2,
      d_dcCKDEPI >= 60 ~ 1
    ),
    levels = 1:2,
    labels = c(">=60", "<60")
    )
  )

# Convert back to Mids
imput.short <- as.mids(long)
imp <- imput.short
