
ProjectTemplate::reload.project()

dataass <- mice::complete(imp, 3)
dataass <- mice::complete(imp, 6)

# dataass <- dataass[edataforimp$num_nation != "LATVIA", ]

# check assumptions for cox models ----------------------------------------

mod <- coxph(formula(paste0("Surv(outtime_hosphf, out_deathhosphf == 1) ~ relevel(natremia, ref = 'Normo/Normo') + ", 
                            paste0(coxvars, collapse = "+"))),
  data = dataass
)

testpat <- cox.zph(mod)
print(sig <- testpat$table[testpat$table[, 3] < 0.05, ])

x11()
plot(testpat[1], resid = F, ylim = c(-4, 4))

mod <- coxph(Surv(outtime_hosphf, out_deathhosphf == 1) ~ relevel(natremia, ref = 'Normo/Normo') + 
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


# 
ormod <- glm(d_dcNyha_cat == "III-IV" ~ relevel(d_hsSod_cat, ref = '>=135') + 
                   ns(num_age, 4) + num_dmgender + ns(num_dmBmi, 4) + num_dmEtio_c1 +
                   num_dmStroke + num_dmPvd + num_dmDiab_c1 + num_dmCopd + num_dmHepa + d_dmThy +
                   num_dmDis + num_dmDepr + num_dmSmoking_c1 + ns(num_dcEf, 4) + num_dcAorSte + 
                   num_dcMitReg + num_dcTriCus + d_reasonforhosp + num_dcRyth + d_hsNyha_cat + 
                   ns(num_dmBpm, 4) + ns(num_dmBp1, 4) + ns(d_hsCKDEPI, 4) + ns(num_hsPot, 4) +
                   ns(num_hsHb, 4) + d_thiazideDiurp + d_loopDiurpmod + num_mdALp +
                   d_arb_or_ace_or_arnip + num_mdBBp + num_mdAmip + num_mdAntiarp + 
                   num_mdDigop + num_mdCcbp + num_mdAdepp + d_xanthinep,
                       family = binomial(link = "logit"), data = dataass
)
termplot(ormod, terms = 2, ask = F, rug = T) # age
termplot(ormod, terms = 4, ask = F, rug = T) # BMI 
termplot(ormod, terms = 15, ask = F, rug = T) # EF
termplot(ormod, terms = 22, ask = F, rug = T) # bpm
termplot(ormod, terms = 23, ask = F, rug = T) # bp1
termplot(ormod, terms = 24, ask = F, rug = T) # egfr ok
termplot(ormod, terms = 25, ask = F, rug = T) # pot
termplot(ormod, terms = 26, ask = F, rug = T) # hb


ormod <- glm(d_dcSod_cat == '>=135' ~ ns(num_age, 4) + num_dmgender + ns(num_dmBmi, 4) + 
num_dmEtio_c1 + num_dmStroke + num_dmPvd + num_dmDiab_c1 + num_dmCopd + 
  num_dmHepa + d_dmThy + num_dmDis + num_dmDepr + num_dmSmoking_c1 + 
  ns(num_dcEf, 4) + d_reasonforhosp + num_dcRyth + 
  num_hsIntr + ns(d_changepercent_weight, 4) + d_residual_congestion + d_dcNyha_cat + ns(num_dcBp1, 4) + 
  ns(d_dcCKDEPI, 4) + d_thiazideDiurh + d_loopDiurhmod + num_mdALh + d_arb_or_ace_or_arnih + 
  num_mdBBh + num_mdAmih + num_mdAntiarh + num_mdDigoh, 
                           family = binomial(link = "logit"), data = dataass)

termplot(ormod, terms = 1, ask = F, rug = T) # age
termplot(ormod, terms = 3, ask = F, rug = T) # BMI
termplot(ormod, terms = 14, ask = F, rug = T) # EF
termplot(ormod, terms = 18, ask = F, rug = T) # weight
termplot(ormod, terms = 21, ask = F, rug = T) # bp1
termplot(ormod, terms = 22, ask = F, rug = T) # egfr ok


# Outliers
kontcoxvars <- edata %>% 
  select(!!!syms(coxvars)) %>%
  select(where(is.numeric))
nams <- colnames(kontcoxvars)
for(i in 1:ncol(kontcoxvars)){
  x11()
  plot(kontcoxvars[, i], xlab = nams[i])
}


kontvars <- edata %>% 
  select(!!!syms(logvars)) %>%
  select(where(is.numeric))
nams <- colnames(kontvars)
for(i in 1:ncol(kontvars)){
  x11()
  plot(kontvars[, i], xlab = nams[i])
}


kontvars <- edata %>% 
  select(!!!syms(predvarsmult)) %>%
  select(where(is.numeric))
nams <- colnames(kontvars)
for(i in 1:ncol(kontvars)){
  x11()
  plot(kontvars[, i], xlab = nams[i])
}