

# Create vars PRIOR to imputation (used in imp model) ------

edata <- edata %>%
  mutate(
    d_dcEf_cat = factor(case_when(
      num_dcEf < 40 ~ 1,
      num_dcEf <= 49 ~ 2,
      num_dcEf >= 50 ~ 3
    ),
    levels = 1:3,
    labels = c("<40", "40-49", ">=50")
    ),
    d_age_cat = case_when(
      num_age < 65 ~ "<65",
      num_age >= 65 ~ ">=65"
    ),
    d_dmBmi_cat = case_when(
      is.na(num_dmBmi) ~ NA_character_,
      num_dmBmi < 25 ~ "1.<25",
      num_dmBmi >= 25 ~ "2.>=25"
    ),
    prev_mi_cabg_pci = factor(case_when(
      is.na(num_dmMi) | is.na(num_dmPci) | is.na(num_dmCabg) ~ NA_real_,
      num_dmMi == "Yes" | num_dmPci == "Yes" | num_dmCabg == "Yes" ~ 1,
      TRUE ~ 0
    ),
    levels = 0:1,
    labels = c("No", "Yes")
    ),
    d_dmBpm_cat = factor(case_when(
      num_dmBpm <= 80 ~ 1,
      num_dmBpm > 80 ~ 2
    ),
    levels = 1:2,
    labels = c("<=80", ">80")
    ),
    d_dmBp1_cat = factor(case_when(
      num_dmBp1 <= 140 ~ 1,
      num_dmBp1 > 140 ~ 2
    ),
    levels = 1:2,
    labels = c("<=140", ">140")
    ),
    d_dcBp1_cat = factor(case_when(
      num_dcBp1 < 110 ~ 2,
      num_dcBp1 >= 110 ~ 1
    ),
    levels = 1:2,
    labels = c(">=110", "<110")
    ),
    d_hsSod_cat = factor(case_when(
      num_hsSod < 135 ~ 2,
      num_hsSod >= 135 ~ 1
    ),
    levels = 1:2,
    labels = c(">=135", "<135")
    ),
    d_dcSod_cat = factor(case_when(
      num_dcSod < 135 ~ 2,
      num_dcSod >= 135 ~ 1
    ),
    levels = 1:2,
    labels = c(">=135", "<135")
    ),
    natremia = factor(case_when(
      d_hsSod_cat == "<135" & d_dcSod_cat == "<135" ~ 1,
      d_hsSod_cat == "<135" & d_dcSod_cat == ">=135" ~ 2,
      d_hsSod_cat == ">=135" & d_dcSod_cat == "<135" ~ 3,
      d_hsSod_cat == ">=135" & d_dcSod_cat == ">=135" ~ 4,
    ),
    levels = 1:4,
    labels = c("Hypo/Hypo", "Hypo/Normo", "Normo/Hypo", "Normo/Normo")
    ),
    d_change_Sod = num_dcSod - num_hsSod,
    d_changepercent_Sod = (num_dcSod - num_hsSod) / num_hsSod * 100,
    d_change_Pot = num_dcPot - num_hsPot,
    d_changepercent_Pot = (num_dcPot - num_hsPot) / num_hsPot * 100,
    d_change_Cre = num_dcCre - num_hsCre,
    d_changepercent_Cre = (num_dcCre - num_hsCre) / num_hsCre * 100,
    d_changepercent_Cre_cat = factor(case_when(
      d_changepercent_Cre <= 30 ~ 0,
      d_changepercent_Cre > 30 ~ 1
    ),
    levels = 0:1, labels = c("<=30", ">30")
    ),

    # eGFR according to CKD-EPI 2021 https://www.nejm.org/doi/full/10.1056/NEJMoa2102953
    tmp_k = if_else(num_dmgender == "Female", 0.7, 0.9),
    tmp_a = if_else(num_dmgender == "Female", -0.241, -0.302),
    tmp_add = if_else(num_dmgender == "Female", 1.012, 1),
    d_hsCKDEPI = 142 * pmin(num_hsCre / tmp_k, 1)^tmp_a * pmax(num_hsCre / tmp_k, 1)^-1.200 * 0.9938^num_age * tmp_add,
    d_hsCKDEPI = if_else(d_hsCKDEPI == Inf, NA_real_, d_hsCKDEPI),
    d_dcCKDEPI = 142 * pmin(num_dcCre / tmp_k, 1)^tmp_a * pmax(num_dcCre / tmp_k, 1)^-1.200 * 0.9938^num_age * tmp_add,
    d_dcCKDEPI = if_else(d_dcCKDEPI == Inf, NA_real_, d_dcCKDEPI),
    d_change_CKDEPI = d_dcCKDEPI - d_hsCKDEPI,
    d_changepercent_CKDEPI = (d_dcCKDEPI - d_hsCKDEPI) / d_hsCKDEPI * 100,
    d_dcCKDEPI_cat = factor(case_when(
      d_dcCKDEPI < 60 ~ 2,
      d_dcCKDEPI >= 60 ~ 1
    ),
    levels = 1:2,
    labels = c(">=60", "<60")
    ),
    d_hsCKDEPI_cat = factor(case_when(
      d_hsCKDEPI < 60 ~ 2,
      d_hsCKDEPI >= 60 ~ 1
    ),
    levels = 1:2,
    labels = c(">=60", "<60")
    ),
    d_dcHb_cat = factor(case_when(
      is.na(num_dcHb) | is.na(num_dmgender) ~ NA_real_,
      num_dcHb < 12 & num_dmgender == "Female" | num_dcHb < 13 & num_dmgender == "Male" ~ 2,
      TRUE ~ 1
    ),
    levels = 1:2,
    labels = c(">=12/13(women/men)", "<12/13(women/men)")
    ),
    d_change_Hb = num_dcHb - num_hsHb,
    d_changepercent_Hb = (num_dcHb - num_hsHb) / num_hsHb * 100,
    d_dcBpm_cat = factor(case_when(
      is.na(num_dcBpm) | is.na(num_dcRyth) ~ NA_real_,
      num_dcBpm >= 70 & num_dcRyth %in% c("Sinus", "Other") |
        num_dcBpm >= 80 & num_dcRyth == "Atrial Fibrillation/Flutter" ~ 2,
      TRUE ~ 1
    ),
    levels = 1:2,
    labels = c("<70/80(sinus/af)", ">=70/80(sinus/af)")
    ),
    d_dmHF_history = case_when(
      is.na(num_dmHF) ~ NA_character_,
      num_dmHF %in% c("Yes with previous hospitalisation", "Yes without previous hospitalisation") ~ "Yes",
      TRUE ~ "No"
    ),
    d_dmHF_cat = case_when(
      is.na(num_dmHF) ~ NA_character_,
      num_dmHF == "Yes with previous hospitalisation" ~ "Previous HF hosp",
      TRUE ~ "No previous HF hosp"
    ),
    d_HFdiagnosis = case_when(
      num_dmHF == "No" ~ "<12mo",
      num_dmMonth %in% c("< 6 months", "6 - 12 months") ~ "<12mo",
      num_dmMonth %in% c("> 12 months") ~ ">12mo"
    ),
    num_dmEtio_c1 = relevel(num_dmEtio_c1, ref = "Non-ischemic heart disease"),
    num_dmEtio = factor(case_when(
      num_dmEtio == "Ischemic heart disease documented by coronary angiography" ~ "IHD doc by ca",
      num_dmEtio == "Ischemic heart disease not documented by coronary angiography" ~ "IHD not documented by ca",
      TRUE ~ as.character(num_dmEtio)
    )),
    d_dmDev_cat = case_when(
      is.na(num_dmDev) ~ NA_character_,
      num_dmDev %in% c("No") ~ "1.No",
      num_dmDev %in% c("PM") ~ "2.PM",
      num_dmDev %in% c("CRT-P", "CRT-D", "ICD") ~ "3.CRT/ICD"
    ),

    # no of non-cardiac comorbs
    d_dmThy = case_when(
      num_dmThy == "No" ~ "No",
      num_dmThy %in% c("Hypothyroidism", "Hyperthyroidism") ~ "Yes"
    ),
    d_anemia = case_when(
      is.na(num_hsHb) | is.na(num_dmgender) ~ NA_character_,
      num_hsHb < 13 & num_dmgender == "Male" ~ "Yes",
      num_hsHb < 12 & num_dmgender == "Female" ~ "Yes",
      TRUE ~ "No"
    ),
    d_X_pulmc_alvoedema = case_when(
      num_dcXrn == "Yes" ~ "No",
      is.na(num_dcXpu) | is.na(num_dcXal) ~ NA_character_,
      num_dcXpu == "No" & num_dcXal == "No" ~ "No",
      num_dcXpu == "Yes" | num_dcXal == "Yes" ~ "Yes",
      TRUE ~ NA_character_
    ),
    d_dcNyha_cat = case_when(
      num_dcNyha %in% c("NYHA I", "NYHA II") ~ "I-II",
      num_dcNyha %in% c("NYHA III", "NYHA IV") ~ "III-IV"
    ),
    d_hsNyha_cat = case_when(
      num_hsNyha %in% c("NYHA I", "NYHA II") ~ "I-II",
      num_hsNyha %in% c("NYHA III", "NYHA IV") ~ "III-IV"
    ),
    tmp_dcnyha = case_when(
      num_dcNyha == "NYHA I" ~ 1,
      num_dcNyha == "NYHA II" ~ 2,
      num_dcNyha == "NYHA III" ~ 3,
      num_dcNyha == "NYHA IV" ~ 4
    ),
    tmp_hsnyha = case_when(
      num_hsNyha == "NYHA I" ~ 1,
      num_hsNyha == "NYHA II" ~ 2,
      num_hsNyha == "NYHA III" ~ 3,
      num_hsNyha == "NYHA IV" ~ 4
    ),
    improvment1class_dcNyha = if_else(tmp_dcnyha < tmp_hsnyha, "Yes", "No"),
    d_qtcfridericia = ifelse(!num_dcRyth_c1 %in% c("Paced", "Other"), num_dcQt / ((60 / num_dcHr2)^0.33), NA),
    d_bsa = sqrt(num_dmHeight * num_dmWeight / 3600),
    d_LAVI = num_dcLavol / d_bsa,
    d_change_weight = num_dcWeight - num_dmWeight,
    d_changepercent_weight = (num_dcWeight - num_dmWeight) / num_dmWeight * 100,
    d_days_in_hosp = as.numeric(num_dcDischdt - num_dmVisitdt),
    d_days_in_hosp_cat = factor(case_when(
      d_days_in_hosp <= 7 ~ 1,
      d_days_in_hosp > 7 ~ 2
    ),
    levels = 1:2,
    labels = c("<= 7 days", ">7 days")
    ),
    d_loopDiurp = case_when(
      is.na(num_mdDiurp_c2) ~ NA_character_,
      num_mdDiurp_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") |
        num_mdDiur2p_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ "Yes",
      TRUE ~ "No"
    ),
    tmp_dosedp1 = case_when(
      num_mdDiurp_c2 == "Flurosemide" ~ num_mdDiurpdo / 40 * 40,
      num_mdDiurp_c2 == "Torasemide" ~ num_mdDiurpdo / 10 * 40,
      num_mdDiurp_c2 == "Bumetanide" ~ num_mdDiurpdo / 1 * 40
    ),
    tmp_dosedp2 = case_when(
      num_mdDiur2p_c2 == "Flurosemide" ~ num_mdDiur2pdo / 40 * 40,
      num_mdDiur2p_c2 == "Torasemide" ~ num_mdDiur2pdo / 10 * 40,
      num_mdDiur2p_c2 == "Bumetanide" ~ num_mdDiur2pdo / 1 * 40
    ),
    d_loopDiurpdose_eqFurosemide = case_when(
      num_mdDiurp_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") &
        num_mdDiur2p_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedp1 + tmp_dosedp2,
      num_mdDiurp_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedp1,
      num_mdDiur2p_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedp2
    ),
    d_loopDiurh = case_when(
      is.na(num_mdDiurh_c2) ~ NA_character_,
      num_mdDiurh_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") |
        num_mdDiur2h_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ "Yes",
      TRUE ~ "No"
    ),
    tmp_dosedh1 = case_when(
      num_mdDiurh_c2 == "Flurosemide" ~ num_mdDiurhdo / 40 * 40,
      num_mdDiurh_c2 == "Torasemide" ~ num_mdDiurhdo / 10 * 40,
      num_mdDiurh_c2 == "Bumetanide" ~ num_mdDiurhdo / 1 * 40
    ),
    tmp_dosedh2 = case_when(
      num_mdDiur2h_c2 == "Flurosemide" ~ num_mdDiur2hdo / 40 * 40,
      num_mdDiur2h_c2 == "Torasemide" ~ num_mdDiur2hdo / 10 * 40,
      num_mdDiur2h_c2 == "Bumetanide" ~ num_mdDiur2hdo / 1 * 40
    ),
    d_loopDiurhdose_eqFurosemide = case_when(
      num_mdDiurh_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") &
        num_mdDiur2h_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedh1 + tmp_dosedh2,
      num_mdDiurh_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedh1,
      num_mdDiur2h_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosedh2
    ),
    d_loopDiurd = case_when(
      is.na(num_mdDiurd_c2) ~ NA_character_,
      num_mdDiurd_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") |
        num_mdDiur2d_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ "Yes",
      TRUE ~ "No"
    ),
    tmp_dosed1 = case_when(
      num_mdDiurd_c2 == "Flurosemide" ~ num_mdDiurddo / 40 * 40,
      num_mdDiurd_c2 == "Torasemide" ~ num_mdDiurddo / 10 * 40,
      num_mdDiurd_c2 == "Bumetanide" ~ num_mdDiurddo / 1 * 40
    ),
    tmp_dosed2 = case_when(
      num_mdDiur2d_c2 == "Flurosemide" ~ num_mdDiur2ddo / 40 * 40,
      num_mdDiur2d_c2 == "Torasemide" ~ num_mdDiur2ddo / 10 * 40,
      num_mdDiur2d_c2 == "Bumetanide" ~ num_mdDiur2ddo / 1 * 40
    ),
    d_loopDiurddose_eqFurosemide = case_when(
      num_mdDiurd_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") &
        num_mdDiur2d_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosed1 + tmp_dosed2,
      num_mdDiurd_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosed1,
      num_mdDiur2d_c2 %in% c("Flurosemide", "Torasemide", "Bumetanide") ~ tmp_dosed2
    ),
    d_arb_or_ace_or_arnip = case_when(
      is.na(num_mdACEp) | is.na(num_mdATp) ~ NA_character_,
      num_mdACEp == "Yes" | num_mdATp == "Yes" | num_mdARNIp == "Yes" ~ "Yes",
      TRUE ~ "No"
    ),
    d_arb_or_ace_or_arnid = case_when(
      is.na(num_mdACEd) | is.na(num_mdATd) ~ NA_character_,
      num_mdACEd == "Yes" | num_mdATd == "Yes" | num_mdARNId == "Yes" ~ "Yes",
      TRUE ~ "No"
    ),
    d_change_f1Nt = num_f1Nt - num_hsNt,
    d_changepercent_f1Nt = (num_f1Nt - num_hsNt) / num_hsNt * 100,
    d_change_f1Bnp = num_f1Bnp - num_hsBnp,
    d_changepercent_f1Bnp = (num_f1Bnp - num_hsBnp) / num_hsBnp * 100,
    d_change_dcNt = num_dcNt - num_hsNt,
    d_changepercent_dcNt = (num_dcNt - num_hsNt) / num_hsNt * 100,
    d_change_dcBnp = num_dcBnp - num_hsBnp,
    d_changepercent_dcBnp = (num_dcBnp - num_hsBnp) / num_hsBnp * 100,
    tmpbothntbnp <- coalesce(d_changepercent_dcNt, d_change_dcBnp),
    d_changepercent_either_dcNtBnp = factor(case_when(
      is.na(tmpbothntbnp) ~ NA_real_,
      tmpbothntbnp <= 30 ~ 0,
      tmpbothntbnp > 30 ~ 1
    ),
    levels = 0:1, labels = c("<=30", ">30")
    ),

    # Outcomes
    enddtm = coalesce(num_f1DeathDt, num_f1contDt),
    startdtm = coalesce(num_dcDischdt, num_dmVisitdt),
    outtime_death = as.numeric(enddtm - startdtm),
    out_death = case_when(
      num_f1vital == "Alive" ~ 0,
      num_f1vital == "Dead" ~ 1
    ),
    out_deathcv = case_when(
      is.na(out_death) ~ NA_real_,
      num_f1DeathCs %in% c("Cardiac", "Vascular") ~ 1,
      TRUE ~ 0
    ), # pats with missing info are NOT included in CV

    out_deathcv_exclhf = case_when(
      is.na(out_death) ~ NA_real_,
      num_f1DeathCs %in% c("Cardiac", "Vascular") &
        (num_f1DthCa != "Heart Failure" | is.na(num_f1DthCa)) ~ 1,
      TRUE ~ 0
    ),
    out_deathnoncv = case_when(
      is.na(out_death) ~ NA_real_,
      num_f1DeathCs == c("Non cardiovascular") ~ 1,
      TRUE ~ 0
    ), # pats with missing info are NOT included in nonCV

    out_deathhf = case_when(
      is.na(out_death) ~ NA_real_,
      num_f1DthCa == "Heart Failure" ~ 1,
      TRUE ~ 0
    ),
    out_deathscd = case_when(
      is.na(out_death) ~ NA_real_,
      num_f1DeathCs == "Cardiac" & num_f1DthCaMd == "Sudden" ~ 1,
      TRUE ~ 0
    ),
    out_deathunknown = case_when(
      is.na(out_death) ~ NA_real_,
      is.na(num_f1DeathCs) & num_f1vital == "Dead" ~ 1,
      TRUE ~ 0
    ),

    # All-cause hosp
    out_hosp = case_when(
      num_f1lost != "No" ~ NA_real_,
      num_f1hosp1 == "Yes" |
        num_f1hosp2 == "Yes" |
        num_f1hosp3 == "Yes" |
        num_f1hosp4 == "Yes" |
        num_f1hosp5 == "Yes" ~ 1,
      TRUE ~ 0
    ),
    out_hospdtm = coalesce(
      num_f1hosp1dt, num_f1hosp2dt, num_f1hosp3dt,
      num_f1hosp4dt, num_f1hosp5dt
    ),
    outtime_hosp = as.numeric(out_hospdtm - startdtm),
    outtime_hospmissing = case_when(
      out_hosp == 1 & is.na(outtime_hosp) ~ 1,
      out_hosp == 1 ~ 0
    ),
    outtime_hosp = ifelse(out_hosp == 1 & is.na(outtime_hosp), outtime_death / 2, outtime_hosp),
    outtime_hosp = pmin(outtime_hosp, outtime_death, na.rm = TRUE),

    # CV
    out_hospcv = case_when(
      num_f1lost != "No" ~ NA_real_,
      num_f1hosp1cs %in% c("Cardiac, non HF", "HF", "Vascular") |
        num_f1hosp2cs %in% c("Cardiac, non HF", "HF", "Vascular") |
        num_f1hosp3cs %in% c("Cardiac, non HF", "HF", "Vascular") |
        num_f1hosp4cs %in% c("Cardiac, non HF", "HF", "Vascular") |
        num_f1hosp5cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ 1,
      TRUE ~ 0
    ),
    out_hospcvdtm = case_when(
      num_f1hosp1cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ num_f1hosp1dt,
      num_f1hosp2cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ num_f1hosp2dt,
      num_f1hosp3cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ num_f1hosp3dt,
      num_f1hosp4cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ num_f1hosp4dt,
      num_f1hosp5cs %in% c("Cardiac, non HF", "HF", "Vascular") ~ num_f1hosp5dt
    ),
    outtime_hospcv = as.numeric(out_hospcvdtm - startdtm),
    outtime_hospcv = ifelse(out_hospcv == 1 & is.na(outtime_hospcv), outtime_death / 2, outtime_hospcv),
    outtime_hospcv = pmin(outtime_hospcv, outtime_death, na.rm = TRUE),

    # CV excl HF
    out_hospcv_exclhf = case_when(
      num_f1lost != "No" ~ NA_real_,
      num_f1hosp1cs %in% c("Cardiac, non HF", "Vascular") |
        num_f1hosp2cs %in% c("Cardiac, non HF", "Vascular") |
        num_f1hosp3cs %in% c("Cardiac, non HF", "Vascular") |
        num_f1hosp4cs %in% c("Cardiac, non HF", "Vascular") |
        num_f1hosp5cs %in% c("Cardiac, non HF", "Vascular") ~ 1,
      TRUE ~ 0
    ),
    out_hospcv_exclhfdtm = case_when(
      num_f1hosp1cs %in% c("Cardiac, non HF", "Vascular") ~ num_f1hosp1dt,
      num_f1hosp2cs %in% c("Cardiac, non HF", "Vascular") ~ num_f1hosp2dt,
      num_f1hosp3cs %in% c("Cardiac, non HF", "Vascular") ~ num_f1hosp3dt,
      num_f1hosp4cs %in% c("Cardiac, non HF", "Vascular") ~ num_f1hosp4dt,
      num_f1hosp5cs %in% c("Cardiac, non HF", "Vascular") ~ num_f1hosp5dt
    ),
    outtime_hospcv_exclhf = as.numeric(out_hospcv_exclhfdtm - startdtm),
    outtime_hospcv_exclhf = ifelse(out_hospcv_exclhf == 1 & is.na(outtime_hospcv_exclhf), outtime_death / 2, outtime_hospcv_exclhf),
    outtime_hospcv_exclhf = pmin(outtime_hospcv_exclhf, outtime_death, na.rm = TRUE),

    # HF hosp
    out_hosphf = case_when(
      num_f1lost != "No" ~ NA_real_,
      num_f1hosp1cs == "HF" |
        num_f1hosp2cs == "HF" |
        num_f1hosp3cs == "HF" |
        num_f1hosp4cs == "HF" |
        num_f1hosp5cs == "HF" ~ 1,
      TRUE ~ 0
    ),
    out_hosphfdtm = case_when(
      num_f1hosp1cs == "HF" ~ num_f1hosp1dt,
      num_f1hosp2cs == "HF" ~ num_f1hosp2dt,
      num_f1hosp3cs == "HF" ~ num_f1hosp3dt,
      num_f1hosp4cs == "HF" ~ num_f1hosp4dt,
      num_f1hosp5cs == "HF" ~ num_f1hosp5dt
    ),
    outtime_hosphf = as.numeric(out_hosphfdtm - startdtm),
    outtime_hosphf = ifelse(out_hosphf == 1 & is.na(outtime_hosphf), outtime_death / 2, outtime_hosphf),
    outtime_hosphf = pmin(outtime_hosphf, outtime_death, na.rm = TRUE),

    # Non-CV
    out_hospnoncv = case_when(
      num_f1lost != "No" ~ NA_real_,
      num_f1hosp1cs %in% c("Non CV", "Renal dysfunction") |
        num_f1hosp2cs %in% c("Non CV", "Renal dysfunction") |
        num_f1hosp3cs %in% c("Non CV", "Renal dysfunction") |
        num_f1hosp4cs %in% c("Non CV", "Renal dysfunction") |
        num_f1hosp5cs %in% c("Non CV", "Renal dysfunction") ~ 1,
      TRUE ~ 0
    ),
    out_hospnoncvdtm = case_when(
      num_f1hosp1cs %in% c("Non CV", "Renal dysfunction") ~ num_f1hosp1dt,
      num_f1hosp2cs %in% c("Non CV", "Renal dysfunction") ~ num_f1hosp2dt,
      num_f1hosp3cs %in% c("Non CV", "Renal dysfunction") ~ num_f1hosp3dt,
      num_f1hosp4cs %in% c("Non CV", "Renal dysfunction") ~ num_f1hosp4dt,
      num_f1hosp5cs %in% c("Non CV", "Renal dysfunction") ~ num_f1hosp5dt
    ),
    outtime_hospnoncv = as.numeric(out_hospnoncvdtm - startdtm),
    outtime_hospnoncv = ifelse(out_hospnoncv == 1 & is.na(outtime_hospnoncv), outtime_death / 2, outtime_hospnoncv),
    outtime_hospnoncv = pmin(outtime_hospnoncv, outtime_death, na.rm = TRUE),

    # all-cause death or hf hosp
    out_deathhosphf = ifelse(out_hosphf == 1, 1, out_death),
    # cv death or hf hosp
    out_deathcvhosphf = ifelse(out_hosphf == 1, 1, out_deathcv),

    # competing risk outcomes
    out_deathcv_cr = if_else(out_deathcv == 0 & out_death == 1, 2, out_deathcv),
    out_deathcv_exclhf_cr = if_else(out_deathcv_exclhf == 0 & out_death == 1, 2, out_deathcv_exclhf),
    out_deathnoncv_cr = if_else(out_deathnoncv == 0 & out_death == 1, 2, out_deathnoncv),
    out_deathhf_cr = if_else(out_deathhf == 0 & out_death == 1, 2, out_deathhf),
    out_deathscd_cr = if_else(out_deathscd == 0 & out_death == 1, 2, out_deathscd),
    out_deathunknown_cr = if_else(out_deathunknown == 0 & out_death == 1, 2, out_deathunknown),
    out_hosp_cr = if_else(out_hosp == 0 & out_death == 1, 2, out_hosp),
    out_hospcv_cr = if_else(out_hospcv == 0 & out_death == 1, 2, out_hospcv),
    out_hospcv_exclhf_cr = if_else(out_hospcv_exclhf == 0 & out_death == 1, 2, out_hospcv_exclhf),
    out_hospnoncv_cr = if_else(out_hospnoncv == 0 & out_death == 1, 2, out_hospnoncv),
    out_hosphf_cr = if_else(out_hosphf == 0 & out_death == 1, 2, out_hosphf),
    out_deathcvhosphf_cr = if_else(out_deathcvhosphf == 0 & out_death == 1, 2, out_deathcvhosphf)
  ) %>%
  mutate(d_no_noncardiac_comorbs = rowSums(select(
    ., num_dmStroke, num_dmPvd, num_dmVte, num_dmDiab_c1,
    num_dmHyChol, num_dmCopd, num_dmApn,
    num_dmHepa, d_dmThy,
    num_dmDis, num_dmDepr,
    num_dmPark, num_dmRheu,
    d_anemia
  ) == "Yes")) %>%
  mutate_if(is.character, as.factor) %>%
  select(-starts_with("tmp_"))
