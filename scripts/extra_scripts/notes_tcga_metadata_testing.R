# Load the package
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)

# Query LUAD clinical data
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
# Download and parse the clinical data
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
colnames(clinical)
head(clinical, n=2)

anc_spectra %>% filter(IID %in% clinical$bcr_patient_barcode)
luad_afr <- filter(anc_spectra, tumor_type == "LUAD", anc3 == "afr") %>%
  pull(IID)
luad_eas <- filter(anc_spectra, tumor_type == "LUAD", anc3 == "eas") %>%
  pull(IID)
luad_eur <- filter(anc_spectra, tumor_type == "LUAD", anc3 == "eur") %>%
  pull(IID)
pys_afr <- filter(clinical, bcr_patient_barcode %in% luad_afr) %>%
  pull(number_pack_years_smoked)
pys_eas <- filter(clinical, bcr_patient_barcode %in% luad_eas) %>%
  pull(number_pack_years_smoked)
pys_eur <- filter(clinical, bcr_patient_barcode %in% luad_eur) %>%
  pull(number_pack_years_smoked)

pys_afr[is.na(pys_afr)] <- 0
pys_eas[is.na(pys_eas)] <- 0
pys_eur[is.na(pys_eur)] <- 0

mean(pys_afr, na.rm = TRUE)
mean(pys_eas, na.rm = TRUE)
mean(pys_eur, na.rm = TRUE)
t.test(pys_afr,pys_eur)

# Query ESCA clinical data
query2 <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
# Download and parse the clinical data
GDCdownload(query2)
clinical2 <- GDCprepare_clinic(query2, clinical.info = "patient")
colnames(clinical2)
clinical3 <- GDCprepare_clinic(query2, clinical.info = "drug")
colnames(clinical3)


select(clinical2, alcohol_history_documented,
       frequency_of_alcohol_consumption,amount_of_alcohol_consumption_per_day) 

esca_afr <- filter(anc_spectra, tumor_type == "ESCA", anc3 == "afr") %>%
  pull(IID)
esca_eas <- filter(anc_spectra, tumor_type == "ESCA", anc3 == "eas") %>%
  pull(IID)
esca_eur <- filter(anc_spectra, tumor_type == "ESCA", anc3 == "eur") %>%
  pull(IID)

apd_eas <- filter(clinical2, bcr_patient_barcode %in% esca_eas) %>%
  pull(amount_of_alcohol_consumption_per_day)
apd_eur <- filter(clinical2, bcr_patient_barcode %in% esca_eur) %>%
  pull(amount_of_alcohol_consumption_per_day)
apd_eas[is.na(apd_eas)] <- 0
apd_eur[is.na(apd_eur)] <- 0
mean(apd_eas, na.rm = TRUE)
mean(apd_eur, na.rm = TRUE)
t.test(apd_eas,apd_eur)
apd_eas <- filter(clinical2, bcr_patient_barcode %in% esca_eas) %>%
  pull(frequency_of_alcohol_consumption)
apd_eur <- filter(clinical2, bcr_patient_barcode %in% esca_eur) %>%
  pull(frequency_of_alcohol_consumption)
apd_eas[is.na(apd_eas)] <- 0
apd_eur[is.na(apd_eur)] <- 0
mean(apd_eas, na.rm = TRUE)
mean(apd_eur, na.rm = TRUE)
t.test(apd_eas,apd_eur)


apd_eas <- filter(clinical2, bcr_patient_barcode %in% esca_eas) %>%
  pull(alcohol_history_documented)
apd_eur <- filter(clinical2, bcr_patient_barcode %in% esca_eur) %>%
  pull(alcohol_history_documented)
table(apd_eas)
table(apd_eur)
t.test(apd_eas,apd_eur)

mat <- matrix(
  c(11, 33,   # EAS: NO, YES
    39, 84),  # EUR: NO, YES
  nrow = 2,
  byrow = TRUE
)

chisq.test(mat)

####################
#ESCA chemotherapy
head(clinical3)


filter(clinical3, bcr_patient_barcode %in% esca_eas) %>%
  pull(drug_name) %>%
  table()
filter(clinical3, bcr_patient_barcode %in% esca_eur) %>%
  pull(drug_name) %>%
  table()
length(esca_eas)
length(esca_eur)

clinical3$drug_name %>%
  table()

###################

# Query LIHC clinical data
query2 <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)
# Download and parse the clinical data
GDCdownload(query2)
clinical4 <- GDCprepare_clinic(query2, clinical.info = "patient")
colnames(clinical4)


lihc_eas <- filter(anc_spectra, tumor_type == "LIHC", anc3 == "eas") %>%
  pull(IID)
lihc_eur <- filter(anc_spectra, tumor_type == "LIHC", anc3 == "eur") %>%
  pull(IID)


clinical4$alcohol_related <- grepl("alcohol", clinical4$history_hepato_carcinoma_risk_factors, ignore.case = TRUE)
table(clinical4$alcohol_related, useNA = "ifany")

apd_eas <- filter(clinical4, bcr_patient_barcode %in% lihc_eas) %>%
  pull(alcohol_related)
apd_eur <- filter(clinical4, bcr_patient_barcode %in% lihc_eur) %>%
  pull(alcohol_related)
table(apd_eas)
table(apd_eur)
t.test(apd_eas,apd_eur)

mat <- matrix(
  c(112, 46,   # EAS: NO, YES
    102, 70),  # EUR: NO, YES
  nrow = 2,
  byrow = TRUE
)

chisq.test(mat)

# both esca and lihc
mat <- matrix(
  c(123, 79,   # EAS: NO, YES
    141, 84),  # EUR: NO, YES
  nrow = 2,
  byrow = TRUE
)

chisq.test(mat)