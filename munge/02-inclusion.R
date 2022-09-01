

# Primary criteria --------------------------------------------------------

flow <- c(paste0("Number of patients in ESC "), nrow(esc))

edata <- esc %>%
  filter(num_dmPtype == "Hospital")
flow <- rbind(flow, c("Hospitalized", nrow(edata)))

#edata <- edata %>%
#  filter(num_dcVital == "Alive" & !is.na(num_dcVital))
#flow <- rbind(flow, c("Discharged alive", nrow(edata)))

#edata <- edata %>%
#  filter(!is.na(num_hsSod) & !is.na(num_dcSod))
#flow <- rbind(flow, c("Non-missing sodium at admission and discharge", nrow(edata)))

edata <- edata %>%
  filter(!is.na(num_hsSod))
flow <- rbind(flow, c("Non-missing sodium at admission", nrow(edata)))

flow <- rbind(flow, c(
  ". Non-missing sodium at discharge",
  nrow(edata %>% filter(!is.na(num_dcSod)))
))

edata <- edata %>%
  mutate(
    enddtm = coalesce(num_f1DeathDt, num_f1contDt),
    startdtm = coalesce(num_dcDischdt, num_dmVisitdt),
    outtime_death = as.numeric(enddtm - startdtm),
    survpop = num_f1lost == "No" & outtime_death >= 0 & !is.na(outtime_death) & num_dcVital == "Alive" & !is.na(num_dcVital)
  )
flow <- rbind(flow, c(
  ". Non-missing sodium at discharge, discharged alive and not lost to follow-up and not negative follow-up times (long-term outcome population)",
  nrow(edata %>% filter(!is.na(num_dcSod) & survpop))
))

colnames(flow) <- c("Criteria", "N")
