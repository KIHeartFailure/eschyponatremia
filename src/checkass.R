
ProjectTemplate::reload.project()

dataass <- mice::complete(imp, 3)
dataass <- mice::complete(imp, 6)

# check assumptions for cox models ----------------------------------------

mod <- coxph(formula(paste0(
  "Surv(outtime_hosphf, out_deathhosphf == 1) ~ relevel(natremia, ref = 'No/No') + ",
  paste0(coxvars, collapse = "+")
)),
data = dataass
)

testpat <- cox.zph(mod)
print(sig <- testpat$table[testpat$table[, 3] < 0.05, ])

x11()
plot(testpat[1], resid = F, ylim = c(-4, 4))

mod <- coxph(Surv(outtime_hosphf, out_deathhosphf == 1) ~ relevel(natremia, ref = "No/No") +
  ns(num_age, 4) + num_dmgender + ns(num_dmBmi, 4) + num_dmEtio_c1 +
  num_dmStroke + num_dmPvd + num_dmDiab_c1 + num_dmCopd + num_dmHepa +
  d_dmThy + num_dmDis + num_dmDepr + num_dmSmoking_c1 + ns(num_dcEf, 4) +
  num_dcAorSte + num_dcMitReg + num_dcTriCus + d_reasonforhosp +
  num_dcRyth + num_hsIntr + num_dcPci + d_residual_congestion + d_dcNyha_cat +
  ns(num_dcBpm, 4) + ns(num_dcBp1, 4) + ns(d_dcCKDEPI, 4) + ns(num_dcPot, 4) + ns(num_dcHb, 4) +
  d_loopDiurd + num_mdALd + d_arb_or_ace_or_arnid + num_mdBBd + num_mdAmid,
data = dataass
)

termplot(mod, terms = 2, ask = F, rug = T) # age
termplot(mod, terms = 4, ask = F, rug = T) # BMI
termplot(mod, terms = 15, ask = F, rug = T) # EF
termplot(mod, terms = 25, ask = F, rug = T) # bpm
termplot(mod, terms = 26, ask = F, rug = T) # bp1
termplot(mod, terms = 27, ask = F, rug = T) # egfr ok
termplot(mod, terms = 28, ask = F, rug = T) # pot
termplot(mod, terms = 29, ask = F, rug = T) # hb


# For log reg models
checkassfunc <- function(event, val, xvars = c("d_hsSod_cat", logvars), data = dataass) {
  logvarsuse <- xvars[xvars != event]
  logvarsuse2 <- logvarsuse
  logvarscont <- c(
    "num_age", "num_dmBmi", "num_dcEf", "num_dmBpm",
    "num_dmBp1", "d_hsCKDEPI", "num_hsPot", "num_hsHb"
  )
  logvarsuse2[logvarsuse2 %in% logvarscont] <- paste0("ns(", logvarsuse2[logvarsuse2 %in% logvarscont], ", 4)")
  ormod <- glm(formula(paste0(event, "=='", val, "' ~ ", paste(logvarsuse2, collapse = " + "))),
    family = binomial(link = "logit"), data = data
  )

  par(mfrow = c(2, 4))
  for (i in (1:length(logvarsuse2))[logvarsuse %in% logvarscont]) {
    termplot(ormod, terms = i, ask = F, rug = T)
  }
  
  ormod <- glm(formula(paste0(event, "=='", val, "' ~ ", paste(logvarsuse, collapse = " + "))),
    family = binomial(link = "logit"), data = data
  )

  print(car::vif(ormod))
  print(ResourceSelection::hoslem.test(ormod$y, ormod$fitted))
}

checkassfunc(
  event = "d_dcNyha_cat",
  val = "III-IV"
)
checkassfunc(
  event = "num_dcVital",
  val = "Dead"
)
checkassfunc(
  event = "d_lengtofstay_cat",
  val = ">7days" 
) # p stat sig
checkassfunc(
  event = "d_dcIccu_cat",
  val = ">2days"
) # p stat sig
checkassfunc(
  event = "d_change_weight_cat",
  val = "<-2kg"
)

checkassfunc(
  event = "d_hsSod_cat",
  val = "<135",
  xvars = logvars
)

checkassfunc(
  event = "d_dcSod_cat",
  val = ">=135",
  xvars = predvarsmult,
  data = dataass %>% filter(d_hsSod_cat == "<135")
)



# Outliers
kontcoxvars <- edata %>%
  select(!!!syms(coxvars)) %>%
  select(where(is.numeric))
nams <- colnames(kontcoxvars)
for (i in 1:ncol(kontcoxvars)) {
  x11()
  plot(kontcoxvars[, i], xlab = nams[i])
}


kontvars <- edata %>%
  select(!!!syms(logvars)) %>%
  select(where(is.numeric))
nams <- colnames(kontvars)
for (i in 1:ncol(kontvars)) {
  x11()
  plot(kontvars[, i], xlab = nams[i])
}


kontvars <- edata %>%
  select(!!!syms(predvarsmult)) %>%
  select(where(is.numeric))
nams <- colnames(kontvars)
for (i in 1:ncol(kontvars)) {
  x11()
  plot(kontvars[, i], xlab = nams[i])
}


ormod <- glm(d_dcSod_cat == ">=135" ~ num_age + num_dmgender + num_dmBmi +
  num_dmEtio_c1 + num_dmStroke + num_dmPvd + num_dmDiab_c1 + num_dmCopd +
  num_dmHepa + d_dmThy + num_dmDis + num_dmDepr + num_dmSmoking_c1 +
  num_dcEf + d_reasonforhosp + num_dcRyth +
  num_hsIntr + d_changepercent_weight + d_residual_congestion + d_dcNyha_cat + num_dcBp1 +
  d_dcCKDEPI + d_thiazideDiurh + d_loopDiurhmod + num_mdALh + d_arb_or_ace_or_arnih +
  num_mdBBh + num_mdAmih + num_mdAntiarh + num_mdDigoh,
family = binomial(link = "logit"), data = dataass
)

car::vif(ormod)
ResourceSelection::hoslem.test(ormod$y, ormod$fitted)
