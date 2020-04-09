# Purpose: to include the novel and known miRNA directly from the literature.

setwd("/Users/helending/Documents/BCB330/miRNA-seq_literature")

# ====================== LIVER CANCER ===========================
liver_cancer <- read.csv("29587854_liver.csv")
colnames(liver_cancer)[1] <- "literature_miR"
liver_cancer <- liver_cancer[-1,]
literature_miR <- liver_cancer$literature_miR[1:104]
literature_miR <- as.character(literature_miR)
length(unique(literature_miR)) # all unique

literature_liver_cancer <- data.frame(matrix(nrow = length(literature_miR), ncol = 2))
colnames(literature_liver_cancer)[1] <- "miR"
colnames(literature_liver_cancer)[2] <- "liver.cancer"
literature_liver_cancer$miR <- literature_miR
literature_liver_cancer <- literature_liver_cancer[-95, ]
literature_liver_cancer$liver_cancer <- "1"

# save
write.table(literature_liver_cancer, "literature_liver_cancer_29587854.csv")
View(literature_liver_cancer)

# ====================== HEAD, NECK ======================
head_and_neck <- read.table("head_neck.csv")
colnames(head_and_neck)[1] <- "miR"
head_and_neck <- cbind(head_and_neck, "head.neck"="1")

# save
write.table(head_and_neck, "literature_head_and_neck.csv")
View(head_and_neck)

# ============= Gastric Adenocarcinoma ==============
gastric_adenocarcinoma <- read.table("gastric_adenocarcinoma.csv")
colnames(gastric_adenocarcinoma)[1] <- "miR"
gastric_adenocarcinoma <- cbind(gastric_adenocarcinoma, "gastric.adenocarcinoma"="1")

# save
write.table(gastric_adenocarcinoma, "literature_gastric_adenocarcinoma.csv")
View(gastric_adenocarcinoma)

# ============= papillary thyroid carcinoma ==============
novel_mirs <- replicate(234, "Tnov-mir-")
novel_mirs <- paste0(novel_mirs, 1:234)
final_table <- data.frame(matrix(nrow=234, ncol=2))
colnames(final_table) <- c("miR", "papillary.thyroid.carcinoma")
final_table$miR <- novel_mirs
final_table$papillary.thyroid.carcinoma <- "1"

# save
write.table(final_table, "literature_papillary_thyroid_carcinoma.csv")
View(final_table)

# ============= lung cancer ==============
lung <- read.csv("lung_cancer.csv")
lung$SUM <- rowSums(lung[,2:6])
expressed_mir <- lung$miRNA[which(lung$SUM > 1)] # 495
final_table <- data.frame(matrix(nrow = length(expressed_mir), ncol = 2))
colnames(final_table) <- c("miR", "lung.cancer")
final_table$miR <- expressed_mir
final_table$lung.cancer <- "1"

# save
write.table(final_table, "literature_lung_cancer.csv")