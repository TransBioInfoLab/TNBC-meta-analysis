# -------------------------------------------------------------------------------------------------- #
# Description: Extract pCR, survival and relevant clinical features
# -------------------------------------------------------------------------------------------------- #
# Load packages and create directory
# -------------------------------------------------------------------------------------------------- #

setwd(".")
dir.base <- getwd()
source(file.path(dir.base, "code/utility/utility_function.R"))
dir.exp <- file.path(dir.base, "/Data/Expr/")
dir.cData <- file.path(dir.base, "/Data/cData/")
dir.clinical <- file.path(dir.base, "/Data/clinical/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)


# -------------------------------------------------------------------------------------------------- #
# GSE25055
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE25055c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "GSE25055_GSE25065_pinfo.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age, 
  TStage = factor(clinical_t_stage),
  NStage = factor(clinical_nodal_status),
  Stage = factor(clinical_ajcc_stage),
  Grade = factor(stringr::str_extract(grade,"[1-9]")),
  pCR = factor(pcr_rd, levels = c("RD", "pCR")),
  Event = drfs_event,
  Time = drfs_time,
  .keep = "none"
)
pinfo_TNBC$pCR <- relevel(pinfo_TNBC$pCR, ref = "RD")
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE25055.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE25065
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE25065c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE25065c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age_years.ch1, 
  TStage = factor(clinical_t_stage.ch1),
  NStage = factor(clinical_nodal_status.ch1),
  Stage = factor(clinical_ajcc_stage.ch1),
  Grade = factor(stringr::str_extract(grade.ch1,"[1-9]")),
  pCR = factor(pathologic_response_pcr_rd.ch1),
  Event = drfs_1_event_0_censored.ch1,
  Time = drfs_even_time_years.ch1,
  .keep = "none"
)
pinfo_TNBC$pCR <- relevel(pinfo_TNBC$pCR, ref = "RD")
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE25065.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE41998
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE41998c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE41998c.csv"),
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  pCR = factor(
    ifelse(PCR == "YES", "pCR", ifelse(PCR == "NO", "RD", NA))
  ),
  .keep = "none"
)
pinfo_TNBC$pCR <- relevel(pinfo_TNBC$pCR, ref = "RD")
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE41998.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE16446
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE16446c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE16446.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = agebin.ch1, 
  TStage = factor(t.ch1),
  NStage = factor(n.ch1),
  Grade = factor(grade.ch1),
  pCR = factor(ifelse(pcr.ch1 == 1, "pCR",ifelse( pcr.ch1 == 0, "RD", NA)), levels = c("RD","pCR")),
  Event = dmfs_event.ch1,
  Time = dmfs_time.ch1/365,
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE16446.csv")
)


# -------------------------------------------------------------------------------------------------- #
# GSE18728
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE18728c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE18728.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  pCR = factor(
    ifelse(response.category.ch1 == "NR", "RD", ifelse(response.category.ch1 == "R", "pCR", NA)),
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE18728.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE20194
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE20194c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE20194c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1, 
  pCR = factor(pcr_vs_rd.ch1 , levels = c("RD","pCR")),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE20194.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE20271
# -------------------------------------------------------------------------------------------------- #

#### Load pre-prossessed data

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE20271c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE20271c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1, 
  Grade = factor(ifelse(bmn.grade.ch1 == "N/A", NA, bmn.grade.ch1)),
  pCR = factor(pcr.or.rd.ch1 , levels = c("RD","pCR")),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE20271.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE18864
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE18864c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE18864c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1, 
  pCR = factor(
    ifelse(`miller.payne.response.ch1` %in% c("0", "1", "2", "3"), "RD", 
                ifelse(`miller.payne.response.ch1` %in% c("4", "5"), "pCR", NA)),
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE18864.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE22093
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE22093c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE22093c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1, 
  Grade = factor(bmn.grade.ch1),
  pCR = factor(
    pcr.v.rd.ch1,
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE22093.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE32646
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE32646c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE32646c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1, 
  Stage = factor(clinical.stage.ch1),
  TStage = factor(clinical.t.stage.ch1),
  Grade = factor(histological.grade.ch1),
  pCR = factor(
    ifelse(pathologic.response.pcr.ncr.ch1 == 'pCR', 'pCR', ifelse(pathologic.response.pcr.ncr.ch1 == "nCR", "RD", NA)),
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE32646.csv")
)

# -------------------------------------------------------------------------------------------------- #
# VUMC
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read_csv(
  file.path(dir.exp, '/VUMC_fpkm_log2.csv'),
)[,-1]

pinfo <- read_csv(
  file.path(dir.cData, "Vanderbilt_meta.csv")
) %>% 
  filter( Timepoint == "PRE") %>%
  column_to_rownames("RNA_seq_ID")

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = AGE, 
  Stage = factor(`Clinical stage`),
  Grade = factor(`Tumor grade`),
  Arm = factor(Arm),
  pCR = factor(
    ifelse(pCR == 'NO', 'RD', 'pCR'),
    levels = c("RD", "pCR")
  ),
  Event = ifelse(is.na(`Disease free interval (months)`), 0, 1),
  Time = as.numeric(`Time to event (death or censoring) in months`)/12, 
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_VUMC.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE22226
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE22226.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE22226c.csv"), 
  row.names = 1
)

#### match TNBC samples

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch2, 
  Grade = factor(histologic.grade..1.grade.i..low....2..grade.ii..intermediate...3..grade.iii..high...4.indeterminate..ch2),
  pCR = factor(
    ifelse(pathological.complete.response..pcr..ch2 == "Yes", "pCR", ifelse(pathological.complete.response..pcr..ch2 == "No", "RD" , NA)),
    levels = c("RD", "pCR")
  ),
  Event = relapse.free.survival.indicator..1.event..local.or.distant.progression.or.death..0.censor.at.last.follow.up..ch2,
  Time = as.numeric(gsub("[^[:digit:].]", "", relapse.free.survival.time...time.from.chemo.start.date.until.earliest.ch2))/365,
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE22226.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE149322
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE149322.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE149322c.csv"), 
  row.names = 1
)

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Arm = factor(arm.ch1),
  pCR = factor(
    ifelse(pcr.ch1 == 1, "pCR", ifelse(pcr.ch1 == 0, 'RD', NA)),
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE149322.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE154524
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE154524_log2.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE154524c.csv"), 
  row.names = 1
)

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Arm = factor(IndRx),
  Age = factor(AGE_C),
  Grade = factor(GRADE),
  TStage = factor(ClinicalT),
  NStage = factor(ClinicalN),
  pCR = factor(
    ifelse(pCR_Breast == 1, "pCR", ifelse(pCR_Breast == 0, 'RD', NA)),
    levels = c("RD", "pCR")
  ),
  Event = efs_stat,
  Time = efs_time/365, 
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE154524.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE192341
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE192341_tpms_log2.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE192341c.csv"), 
  row.names = 1
)

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1,
  TStage = factor(t.stage.ch1),
  NStage = factor(n.stage..0..neg..1..pos..ch1),
  pCR = factor(
    ifelse(pathological.complete.response..pcr..ch1 == "pCR", "pCR", ifelse(pathological.complete.response..pcr..ch1 == "No pCR", 'RD', NA)),
    levels = c("RD", "pCR")
  ),
  Event =  ifelse(recurrence.ch1 == 'Yes', 1, ifelse(
    recurrence.ch1 == 'No', 0, NA
  )),
  Time = rfs..years..ch1, 
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE192341.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE163882
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE163882_tpms_log2.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE163882c.csv"), 
  row.names = 1
)

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Age = age.ch1,
  Stage = factor(breast.cancer.stage.ch1),
  Grade = tumor.grage.ch1,
  pCR = factor(
    response.to.nac.ch1,
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE163882.csv")
)

# -------------------------------------------------------------------------------------------------- #
# GSE164458
# -------------------------------------------------------------------------------------------------- #

exp_TNBC <- read.csv(
  file.path(dir.exp, '/exp_GSE164458c.csv')
)[,-1]

pinfo <- read.csv(
  file.path(dir.cData, "pinfo_GSE164458c.csv"), 
  row.names = 1
)

pinfo_TNBC <- pinfo[match(colnames(exp_TNBC), rownames(pinfo)),]
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Arm = planned_arm_code.,
  Stage = factor(pretreatment_lymphnode_stage.),
  pCR = factor(
    pathologic_complete_response.,
    levels = c("RD", "pCR")
  ),
  .keep = "none"
)
identical(colnames(exp_TNBC), rownames(pinfo_TNBC))
pinfo_TNBC <- pinfo_TNBC %>% mutate(
  Sample_ID = rownames(pinfo_TNBC),
  .before = 1
)

write_csv(
  pinfo_TNBC,
  file.path(dir.clinical, "clinical_GSE164458.csv")
)
