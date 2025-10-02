#################################
## Title: Setup data for Spitzer et al figures
## Date: 2025.02.24
## Author: Avishay Spitzer
## Description: This script sets up relevant data tables for plotting the Spitzer et al figure panels
#################################

####################################################################################################################################
# Global definitions
####################################################################################################################################

pt_pairs <- pt_pairs %>%
  filter(n >= 50) %>%
  group_by(Patient) %>%
  filter("T1" %in% Timepoint & "T2" %in% Timepoint) %>%
  ungroup()
dim(pt_pairs)

feature_scores_all_per_sample <- feature_scores_all_per_sample %>%
  left_join(gbm_subtypes_tbl %>%
              select(Sample, PvsR), by = "Sample")

comp_cluster_stats <- comp_cluster_stats %>%
  left_join(gbm_subtypes_tbl %>%
              select(ID, PvsR), by = "ID")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Add survival data to GBM subtypes table
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

surg_int <- as.numeric(sample_data$SurgicalIntervalmonth)
names(surg_int) <- sample_data$Patient
surg_int <- surg_int[!is.na(surg_int)]

gbm_subtypes_tbl$SurgicalInterval <- surg_int[gbm_subtypes_tbl$Patient]

os <- as.numeric(sample_data$OSmonth)
names(os) <- sample_data$Patient
os <- os[!is.na(os)]

gbm_subtypes_tbl$OSmonth <- os[gbm_subtypes_tbl$Patient]

vt <- sample_data$PatientStatus
vt[vt %ni% c("Dead", "Alive", "Censored")] <- NA
names(vt) <- sample_data$Patient
vt <- vt[!is.na(vt)]
vt <- ifelse(vt == "Dead", 1 , 0)

gbm_subtypes_tbl$VT <- vt[gbm_subtypes_tbl$Patient]

gbm_subtypes_tbl$SurgicalInterval2 <- gbm_subtypes_tbl$OSmonth - gbm_subtypes_tbl$SurgicalInterval

n_rec_vec <- sample_data %>%
  group_by(Patient) %>%
  summarise(Nrec = case_when("3rd recurrence" %in% PvsR ~ 3,
                             "2nd recurrence" %in% PvsR ~ 2,
                             TRUE ~ 1))
n_rec_vec <- setNames(n_rec_vec$Nrec, n_rec_vec$Patient)

gbm_subtypes_tbl$Nrecurrences <- n_rec_vec[gbm_subtypes_tbl$Patient]

age_vec <- setNames(sample_data$Age.at.initial.diagnosis, sample_data$Patient)
age_vec <- age_vec[!is.na(age_vec)]

gender_vec <- setNames(sample_data$Gender, sample_data$Patient)
gender_vec <- gender_vec[!is.na(gender_vec)]

mgmt_vec <- setNames(sample_data$MGMTstatus[sample_data$Timepoint == "T1"], sample_data$Patient[sample_data$Timepoint == "T1"])
mgmt_vec[mgmt_vec == "UA"] <- "NOS"

gbm_subtypes_tbl$Age <- age_vec[gbm_subtypes_tbl$Patient]
gbm_subtypes_tbl$Gender <- gender_vec[gbm_subtypes_tbl$Patient]
gbm_subtypes_tbl$MGMT <- mgmt_vec[gbm_subtypes_tbl$Patient]

gbm_subtypes_tbl <- gbm_subtypes_tbl %>%
  left_join(sample_data %>%
              select(Sample, PvsR), by = "Sample")

gbm_subtypes_tbl <- gbm_subtypes_tbl %>%
  mutate(PvsR = as.character(PvsR)) %>%
  mutate(PvsR = ifelse(PvsR == "Primary", "Primary", "Recurrence"))
