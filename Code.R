## devtools Package
#install.packages("devtools")
library(devtools)

## ICBioMark Package
# devtools::install_github("cobrbra/ICBioMark", force = TRUE)
library(ICBioMark)
load_all("../../ICBioMark/")

## ecTMB Package
#
library(ecTMB)

## Other R packages needed
#install.packages("cowplot")
library(cowplot)
#install.packages("magrittr")
library(magrittr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("readr")
library(readr)
#install.packages("tidyr")
library(tidyr)
#install.package("latex2exp")
library(latex2exp)
#install.packages("purrr")
library(purrr)
#install.packages("tibble")
library(tibble)

## Figures path
fig_path <- "figures/"



### Main workflow
message("Getting tables")
nsclc_tables <- get_mutation_tables(maf = nsclc_maf, include_synonymous = FALSE,
                                    acceptable_genes = ensembl_gene_lengths$Hugo_Symbol)

message("Getting generative model")
# nsclc_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = nsclc_tables$train,
#                                  progress = TRUE)
# write_rds(x = nsclc_gen_model, file = "data/results/nsclc_gen_model")

nsclc_gen_model <- read_rds("data/results/nsclc_gen_model")

message("Getting first-fit")
# nsclc_pred_first_tmb <- pred_first_fit(gen_model = nsclc_gen_model, lambda = exp(seq(-18, -26, length.out = 100)),
                                         # gene_lengths = ensembl_gene_lengths, training_matrix = nsclc_tables$train$matrix)
# write_rds(x = nsclc_pred_first_tmb, file = "data/results/nsclc_pred_first_tmb")

nsclc_pred_first_tmb <- read_rds("data/results/nsclc_pred_first_tmb")

message("Getting refit")
nsclc_pred_refit_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb, gene_lengths = ensembl_gene_lengths)

nsclc_tmb_values <- get_biomarker_tables(nsclc_maf, biomarker = "TMB")

# nsclc_pred_first_tib <- pred_first_fit(gen_model = nsclc_gen_model, lambda = exp(seq(-18, -30, length.out = 100)),
#                                        gene_lengths = ensembl_gene_lengths, training_matrix = nsclc_tables$train$matrix,
#                                        biomarker = "TIB")
# write_rds(x = nsclc_pred_first_tib, file = "data/results/nsclc_pred_first_tib")

nsclc_pred_first_tib <- read_rds("data/results/nsclc_pred_first_tib")
nsclc_pred_refit_tib <- pred_refit_range(pred_first = nsclc_pred_first_tib, gene_lengths = ensembl_gene_lengths,
                                         biomarker = "TIB")
nsclc_tib_values <- get_biomarker_tables(nsclc_maf, biomarker = "TIB")

message("Getting count estimators")
nsclc_pred_count_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb, gene_lengths = ensembl_gene_lengths, model = "Count", biomarker = "TMB", training_data = nsclc_tables$train, training_values = nsclc_tmb_values$train)
nsclc_pred_count_tib <- pred_refit_range(pred_first = nsclc_pred_first_tib, gene_lengths = ensembl_gene_lengths, model = "Count", biomarker = "TIB", training_data = nsclc_tables$train, training_values = nsclc_tib_values$train)

message("Getting OLM estimators")
# To give linear estimators a fighting chance, we train them on data than isn't
# separated into indel/non-indel (this helps with overfitting).

nsclc_linear_tables <- get_mutation_tables(maf = nsclc_maf, acceptable_genes = ensembl_gene_lengths$Hugo_Symbol, for_biomarker = "TMB", include_synonymous = FALSE)

nsclc_pred_linear_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb, gene_lengths = ensembl_gene_lengths,
                                          model = "OLM", biomarker = "TMB", training_data = nsclc_linear_tables$train, training_values = nsclc_tmb_values$train, max_panel_length = 1000000)
nsclc_pred_linear_tib <- pred_refit_range(pred_first = nsclc_pred_first_tib, gene_lengths = ensembl_gene_lengths,
                                          model = "OLM", biomarker = "TIB", training_data = nsclc_linear_tables$train, training_values = nsclc_tib_values$train, max_panel_length = 1000000)



### ecTMB Workflow
# extdataDir = "./data/ecTMB_data/references"
# exomef                 = file.path(extdataDir, "exome_hg38_vep.Rdata" )  #### hg38 exome file
# covarf                 = file.path(extdataDir,"gene.covar.txt")   ### gene properties
# mutContextf            = file.path(extdataDir,"mutation_context_96.txt" )  ### 96 mutation contexts
# ref                    = file.path(extdataDir,"GRCh38.d1.vd1.fa" )

# nsclc_maf_grch38 <- read_tsv("data/nsclc_maf_grch38.tsv", comment = "#")

# URL = "https://github.com/bioinform/ecTMB/releases/download/v0.1.0/ecTMB_data.tar.gz"
# download.file(URL,destfile = "data/ecTMB.example.tar.gz")
# untar("./data/ecTMB.example.tar.gz")

# URL_ref = "https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834"
# download.file(URL_ref,destfile = "./data/ecTMB_data/GRCh38.d1.vd1.fa.tar.gz")
# untar("./data/ecTMB_data/GRCh38.d1.vd1.fa.tar.gz")

# trainset = nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$train$sample_list) %>%
#   filter(!(Tumor_Sample_Barcode %in% "TCGA-86-A4P8-01")) %>% #hacky solution: There were no non-silent mutations in this observation, which was breaking ecTMB
#   readData(exomef = exomef, covarf = covarf, mutContextf = mutContextf, ref = ref)
# write_rds(x = trainset, file = "data/temporary_storage/trainset")
#
# trainset <- read_rds("data/temporary_storage/trainset")
#
# MRtriProb = getBgMRtri(trainset)

# trainedModel = fit_model(trainset, MRtriProb, cores = 1)
# write_rds(trainedModel, "data/temporary_storage/trainedModel")
#
# trainedModel = read_rds("data/temporary_storage/trainedModel")
#
# TSO_500_panel <- "data/tso_500_bed.bed"
# F1_panel <- "data/foundation_bed.bed"
# MSK_panel <- "data/msk_impact_bed.bed"
# TST_170_panel <- "data/tst_170_bed.bed"

# valset_WES <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$val$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref)
# write_rds(valset_WES, "data/temporary_storage/valset_WES")
# valset_WES <- read_rds("data/temporary_storage/valset_WES")


# testset_WES <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$test$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref)
# write_rds(testset_WES, "data/temporary_storage/testset_WES")
# testset_WES <- read_rds("data/temporary_storage/testset_WES")

## TST-170 Panel
# sample_tst_170_val = data.frame(SampleID = nsclc_tables$val$sample_list, BED = TST_170_panel, stringsAsFactors = FALSE)
# 
# valset_tst_170 <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$val$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref, samplef = sample_tst_170_val)
# write_rds(valset_tst_170, "data/temporary_storage/valset_tst_170")
# valset_tst_170 <- read_rds("data/temporary_storage/valset_tst_170")
#
# val_pred_tst_170 <- pred_TMB(valset_tst_170, WES = valset_WES, cores = 1, params = trainedModel, mut.nonsil = T,
#                              gid_nonsil_p = trainset$get_nonsil_passengers(0.95)) %>%
#   mutate(estimated_values = mean(nsclc_tmb_values$val$TMB) *ecTMB_panel_TMB/ mean(WES_TMB),
#          true_values = mean(nsclc_tmb_values$val$TMB)*WES_TMB / mean(WES_TMB),
#          Tumor_Sample_Barcode = sample) %>%
#   select(Tumor_Sample_Barcode, estimated_values, true_values) %>%
#   mutate(panel = "TST-170")
#
# write_tsv(val_pred_tst_170, "data/results/val_pred_tst_170.tsv")
#
# ## Foundation One Panel
# sample_f1_val = data.frame(SampleID = nsclc_tables$val$sample_list, BED = F1_panel, stringsAsFactors = FALSE)
#
# valset_f1 <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$val$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref, samplef = sample_f1_val)
# write_rds(valset_f1, "data/temporary_storage/valset_f1")
# valset_f1 <- read_rds("data/temporary_storage/valset_f1")
#
# val_pred_f1 <- pred_TMB(valset_f1, WES = valset_WES, cores = 1, params = trainedModel, mut.nonsil = T,
#                              gid_nonsil_p = trainset$get_nonsil_passengers(0.95)) %>%
#   mutate(estimated_values = mean(nsclc_tmb_values$val$TMB) *ecTMB_panel_TMB/ mean(WES_TMB),
#          true_values = mean(nsclc_tmb_values$val$TMB)*WES_TMB / mean(WES_TMB),
#          Tumor_Sample_Barcode = sample) %>%
#   select(Tumor_Sample_Barcode, estimated_values, true_values) %>%
#   mutate(panel = "F1")
#
# write_tsv(val_pred_f1, "data/results/val_pred_f1.tsv")
#
# ## TSO-500 Panel
# sample_tso_500_val = data.frame(SampleID = nsclc_tables$val$sample_list, BED = TSO_500_panel, stringsAsFactors = FALSE)
#
# valset_tso_500 <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$val$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref, samplef = sample_tso_500_val)
# write_rds(valset_tso_500, "data/temporary_storage/valset_tso_500")
#
# val_pred_tso_500 <- pred_TMB(valset_tso_500, WES = valset_WES, cores = 1, params = trainedModel, mut.nonsil = T,
#                         gid_nonsil_p = trainset$get_nonsil_passengers(0.95)) %>%
#   mutate(estimated_values = mean(nsclc_tmb_values$val$TMB) *ecTMB_panel_TMB/ mean(WES_TMB),
#          true_values = mean(nsclc_tmb_values$val$TMB)*WES_TMB / mean(WES_TMB),
#          Tumor_Sample_Barcode = sample) %>%
#   select(Tumor_Sample_Barcode, estimated_values, true_values) %>%
#   mutate(panel = "TSO-500")
#
# write_tsv(val_pred_tso_500, "data/results/val_pred_tso_500.tsv")
#
# ## MSK-IMPACT Panel
# sample_msk_val = data.frame(SampleID = nsclc_tables$val$sample_list, BED = MSK_panel, stringsAsFactors = FALSE)
#
# valset_msk <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$val$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref, samplef = sample_msk_val)
# write_rds(valset_msk, "data/temporary_storage/valset_msk")
#
# val_pred_msk <- pred_TMB(valset_msk, WES = valset_WES, cores = 1, params = trainedModel, mut.nonsil = T,
#                              gid_nonsil_p = trainset$get_nonsil_passengers(0.95)) %>%
#   mutate(estimated_values = mean(nsclc_tmb_values$val$TMB) *ecTMB_panel_TMB/ mean(WES_TMB),
#          true_values = mean(nsclc_tmb_values$val$TMB)*WES_TMB / mean(WES_TMB),
#          Tumor_Sample_Barcode = sample) %>%
#   select(Tumor_Sample_Barcode, estimated_values, true_values) %>%
#   mutate(panel = "MSK-I")
#
# write_tsv(val_pred_msk, "data/results/val_pred_msk.tsv")



### Figure 1
message("Creating Figure 1")



## Subfigure 1
fig1p1 <- nsclc_survival %>%
  add_count(SEX, name = "n_sex") %>%
  mutate(Sex = SEX, Age = as.numeric(AGE)) %>%
  filter(Sex %in% c("Male", "Female")) %>%
  ggplot(aes(x = Sex, y = Age)) + geom_violin(fill = "black", alpha = 0.5, colour = "black") +
  geom_text(aes(x = Sex, y = 67, label = paste("n = ", n_sex))) + theme_minimal() +
  theme(axis.title.x=element_blank())

stage_groups <- c("I", "I", "I", "II", "II", "II", "III", "III", "III", "IV")
names(stage_groups) <- c("I", "IA", "IB", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IV")

smoking_groups <- c("Lifelong Non-Smoker", "Current Reformed \n Smoker (>15 Years)",
                    "Current Reformed \n Smoker (<= 15 Years)", "Current Reformed \n Smoker (Duration UnSpecified)",
                    "Current Smoker")
names(smoking_groups) <- c("Lifelong Non-Smoker", "Current Reformed Smoker For > 15 Years",
                           "Current Reformed Smoker For < Or = 15 Years", "Current Reformed Smoker, Duration Not Specified",
                           "Current Smoker")

fig1p2 <- nsclc_survival %>%
  group_by(SMOKING_HISTORY, STAGE) %>%
  count() %>%
  ungroup() %>%
  mutate(Stage = stage_groups[STAGE]) %>%
  filter(!is.na(Stage)) %>%
  mutate(Stage = factor(Stage, levels = c("IV","III","II","I"))) %>%
  mutate(Smoking_History = smoking_groups[SMOKING_HISTORY])  %>%
  filter(!is.na(Smoking_History)) %>%
  select(Stage, Smoking_History, n) %>%
  group_by(Stage, Smoking_History) %>%
  mutate(n = sum(n)) %>%
  distinct() %>%
  ggplot(aes(x = Smoking_History, y = n, fill = Stage, order = Stage)) +
  geom_col(position = position_stack(), alpha = 0.8, colour = "black") + scale_fill_brewer(palette = "Greys") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  labs(y = "Frequency") + theme(axis.title.x=element_blank())

fig1 <- plot_grid(fig1p1, fig1p2, labels = "AUTO", align = "h")
ggsave(paste0(fig_path, "fig1.png"), fig1, width = 10, height = 4)



### Figure 2
message("Creating Figure 2")



tmb_tib_train <- inner_join(get_biomarker_tables(maf = nsclc_maf, biomarker = "TMB")$train, get_biomarker_tables(maf = nsclc_maf, biomarker = "TIB")$train, by = "Tumor_Sample_Barcode")

## Subfigure 1
fig2p1 <- tmb_tib_train %>%
  gather(key = "Biomarker", value = "Count", -Tumor_Sample_Barcode) %>%
  mutate(Biomarker = factor(Biomarker, levels = c("TMB", "TIB"))) %>%
  ggplot(aes(x = Biomarker, y = Count)) +
  geom_violin(fill = "black", alpha = 0.5, colour = "black") +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Mutation Count") + theme(axis.title.x = element_blank())


## Subfigure 2
variant_names <- c("Missense", "Nonsense", "Splice site", "Translation start",
                   "Nonstop", "In-frame indel", "In-frame indel", "Frame-shift indel",
                   "Frame-shift indel")
names(variant_names) <- names(get_mutation_dictionary(include_synonymous = FALSE))

fig2p2 <- nsclc_maf %>%
  filter(Tumor_Sample_Barcode %in% nsclc_tables$train$sample_list) %>%
  filter(Variant_Classification %in% names(variant_names)) %>%
  mutate(Variant_Classification = variant_names[Variant_Classification]) %>%
  group_by(Variant_Classification) %>%
  count() %>%
  ungroup() %>%
  mutate(n = n/sum(n)) %>%
  ggplot(aes(x = reorder(Variant_Classification, - n), y = n)) +
  geom_col(fill = "black", alpha = 0.5, colour = "black") + labs(y = "Proportion of Nonsynonymous Mutations") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.8)) +
  theme(axis.title.x=element_blank())


fig2 <- plot_grid(fig2p1, fig2p2, labels = "AUTO", align = "h")
ggsave(paste0(fig_path, "fig2.png"), fig2, width = 10, height = 4)




### Section 3 Intro Stats
n_train_samples <- length(nsclc_tables$train$sample_list)
n_val_samples <- length(nsclc_tables$va$sample_list)
n_test_samples <- length(nsclc_tables$test$sample_list)
print(paste("Training Samples:", n_train_samples))
print(paste("Validation Samples:", n_val_samples))
print(paste("Test Samples:", n_test_samples))

n_genes <- length(unique(nsclc_maf$Hugo_Symbol))
print(paste("Genes Analysed:", n_genes))

tmb_mean <- mean(tmb_tib_train$TMB)
tib_mean <- mean(tmb_tib_train$TIB)
print(paste("TMB Mean:", signif(tmb_mean, 3)))
print(paste("TIB Mean:", signif(tib_mean, 3)))

s3.intro.stats <- data.frame(n_train_samples = n_train_samples,
                              n_val_samples = n_val_samples,
                              n_test_samples = n_test_samples,
                              n_genes = n_genes,
                              tmb_mean = tmb_mean,
                              tib_mean = tib_mean)
write_tsv(x = s3.intro.stats, file = 'data/results/s3.intro.stats.tsv')


### Figure 3
message("Creating Figure 3")



fig3 <- vis_model_fit(nsclc_gen_model)
ggsave(paste0(fig_path, "fig3.png"), fig3, width = 10, height = 5)



### Figure 4
message("Creating Figure 4")



gene_locations <- nsclc_maf %>%
  mutate(Chromosome = if_else(Chromosome == "23", "X", Chromosome)) %>%
  mutate(Chromosome = if_else(Chromosome == "24", "Y", Chromosome)) %>%
  mutate(Start_Position = as.numeric(Start_Position)) %>%
  filter(Hugo_Symbol %in% nsclc_gen_model$names$gene_list) %>%
  select(Hugo_Symbol, Chromosome, Start_Position) %>%
  distinct() %>%
  group_by(Hugo_Symbol) %>%
  mutate(Position = mean(Start_Position)) %>%
  ungroup() %>%
  select(Hugo_Symbol, Chromosome, Position) %>%
  distinct() %>%
  arrange(Hugo_Symbol) %>%
  filter(Hugo_Symbol != lead(Hugo_Symbol))

gene_locations_and_coefficients <- data.frame(Hugo_Symbol = nsclc_gen_model$names$gene_list,
                                              ns_coefficient = as.vector(nsclc_gen_model$fit$beta[seq(1,length(nsclc_gen_model$names$col_names), 2), nsclc_gen_model$s_min]),
                                              i_coefficient = as.vector(nsclc_gen_model$fit$beta[seq(2,length(nsclc_gen_model$names$col_names), 2), nsclc_gen_model$s_min])) %>%
  inner_join(gene_locations, by = "Hugo_Symbol")

chromosomes <- c(1:22, "Y","X")
chromosomes <- factor(chromosomes, chromosomes)
names(chromosomes) <- chromosomes
chrom_info <- gene_locations_and_coefficients %>%
  group_by(Chromosome) %>%
  mutate(Chromosome_length = max(Position)) %>%
  select(Chromosome, Chromosome_length) %>%
  distinct() %>%
  mutate(Chromosome = chromosomes[Chromosome]) %>%
  arrange(Chromosome) %>%
  filter(!is.na(Chromosome)) %>%
  ungroup()

chrom_info[23,'Chromosome_length'] <- chrom_info[24, 'Chromosome_length']

chrom_info  <- chrom_info %>%
  mutate(cum_chrom_length = cumsum(Chromosome_length))

manhat_data <- gene_locations_and_coefficients %>%
  full_join(chrom_info, by = "Chromosome")

extreme_data_ns <- manhat_data %>%
  arrange(ns_coefficient) %>%
  filter(!is.na(Chromosome)) %>%
  {.[(nrow(.) - 4):nrow(.),]}

chrom_alternate <- rep(c(0,1), 12)
names(chrom_alternate) <- names(chromosomes)

fig4 <- manhat_data %>%
  filter(!is.na(Chromosome)) %>%
  mutate(full_position = cum_chrom_length + Position) %>%
  mutate(Chromosome = factor(chrom_alternate[Chromosome])) %>%
  ggplot(aes(x = full_position, y = ns_coefficient,colour = Chromosome, size = ns_coefficient)) + geom_point(alpha = 0.75) +
  geom_text(data = extreme_data_ns, aes(x = cum_chrom_length + Position + 2.1*nchar(Hugo_Symbol)*10^7, y = ns_coefficient, label = Hugo_Symbol, colour= factor(chrom_alternate[Chromosome])), size = 3, position=position_jitter()) +
  scale_x_continuous(label = names(chromosomes), breaks = chrom_info$cum_chrom_length + chrom_info$Chromosome_length/2 ) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), 12)) +
  scale_size_continuous(range = c(0.5,1)) +
  labs(x = NULL,
       y = TeX("$\\hat{\\lambda}_g$")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )
ggsave(filename = "figures/fig4.png", plot = fig4, width = 8, height = 4)



### Figure 5
message("Creating Figure 5")


extreme_data_i <- manhat_data %>%
  arrange(i_coefficient) %>%
  filter(!is.na(Chromosome)) %>%
  {.[c(1:5, (nrow(.) - 4):nrow(.)),]}

fig5 <- manhat_data %>%
  filter(!is.na(Chromosome)) %>%
  mutate(full_position = cum_chrom_length + Position) %>%
  mutate(Chromosome = factor(chrom_alternate[Chromosome])) %>%
  ggplot(aes(x = full_position, y = i_coefficient,colour = Chromosome, size = i_coefficient)) + geom_point(alpha = 0.75) +
  geom_text(data = extreme_data_i, aes(x = cum_chrom_length + Position + 2.1*nchar(Hugo_Symbol)*10^7, y = i_coefficient, label = Hugo_Symbol, colour= factor(chrom_alternate[Chromosome])), size = 3) +
  scale_x_continuous(label = names(chromosomes), breaks = chrom_info$cum_chrom_length + chrom_info$Chromosome_length/2 ) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), 12)) +
  scale_size_continuous(range = c(0.5,1)) +
  labs(x = NULL,
       y = TeX("$\\hat{\\eta}_{g, indel$")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )
ggsave(filename = "figures/fig5.png", plot = fig5, width = 8, height = 4)



### Section 3.1 Stats
ns_sparsity = 1 - mean(manhat_data$ns_coefficient != 0)
i_sparsity = 1 - mean(manhat_data$i_coefficient != 0)
print(paste0("Non-synonymous sparsity:", signif(100*ns_sparsity, 3), "%"))
print(paste0("Indel sparsity: ", signif(100*i_sparsity,3), "%"))

s3.1.stats <- data.frame(ns_sparsity = ns_sparsity,
                         i_sparsity = i_sparsity)
write_tsv(x = s3.1.stats, file = "data/results/s3.1.stats.tsv")



### Figure 6
message("Creating Figure 6")



first_stats_tmb <- nsclc_pred_first_tmb %>%
  get_predictions(new_data = nsclc_tables$val) %>%
  get_stats(biomarker_values = nsclc_tmb_values$val, model = "First-fit T", threshold = 300)

refit_stats_tmb <- nsclc_pred_refit_tmb %>%
  get_predictions(new_data = nsclc_tables$val) %>%
  get_stats(biomarker_values = nsclc_tmb_values$val, model = "Refitted T", threshold = 300)

fig6 <- bind_rows(refit_stats_tmb, first_stats_tmb) %>%
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  filter(panel_length <= 2000000) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line(size = 1) + ylim(0, 1) +
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") +
  theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))

ggsave(filename = "figures/fig6.png", plot = fig6, height = 4, width = 8)



### Figure 7
message("Creating Figure 7")



panels <- list("TST-170" = read_tsv("data/tst_170_genes.tsv")$Hugo_Symbol,
               "F1" = read_tsv("data/foundation_genes.tsv")$Hugo_Symbol,
               "MSK-I" = read_tsv("data/msk_impact_genes.tsv")$Hugo_Symbol)#,
               #"TSO-500" = read_tsv("data/tso_500_genes.tsv")$Hugo_Symbol)
models <- c("T", "Count", "OLM")

datasets <- list(nsclc_tables, nsclc_tables, nsclc_linear_tables)
names(datasets) <- models


non_ectmb_model_stats <- panels %>%
  purrr::map(function(panel) purrr::map(models, function(model) pred_refit_panel(model = model, pred_first = nsclc_pred_first_tmb, gene_lengths = ensembl_gene_lengths,
                                                                  genes = panel, training_data = datasets[[model]]$train,
                                                                  training_values = nsclc_tmb_values$train) %>%
                                          get_predictions(new_data = datasets[[model]]$val) %>%
                                          get_stats(biomarker_values = nsclc_tmb_values$val)) %>%
    bind_rows() %>%
    mutate(model = rep(models, each = 2))) %>%
    bind_rows() %>%
    mutate(panel = rep(names(panels), each = 6))

rm(datasets)

ectmb_preds <- list("TST-170" = read_tsv("data/results/val_pred_tst_170.tsv"),
                    "F1" = read_tsv("data/results/val_pred_f1.tsv"),
                    "MSK-I" = read_tsv("data/results/val_pred_msk.tsv"))#,
                    #"TSO-500" = read_tsv("data/results/val_pred_tso_500.tsv"))
ectmb_model_stats <- purrr::map(ectmb_preds, ~ list(predictions = select(column_to_rownames(., "Tumor_Sample_Barcode"), estimated_values),
                                                    panel_lengths = c(0))) %>%
  purrr::map(~ get_stats(., biomarker_values = nsclc_tmb_values$val, model = "ecTMB")) %>%
  bind_rows() %>%
  mutate(panel = rep(names(ectmb_preds), each = 2))

model_stats <- bind_rows(non_ectmb_model_stats, ectmb_model_stats) %>%
  group_by(panel) %>%
  mutate(panel_length = median(panel_length)) %>%
  ungroup() %>%
  mutate(metric = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(metric = factor(metric, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  mutate(panel = factor(panel, levels = c("Our Procedure", "TST-170", "F1", "MSK-I", "TSO-500"))) %>%
  mutate(model = factor(model, c("T", "ecTMB", "Count", "OLM"))) %>%
  mutate(panel_length = panel_length / 1000000)

fig7 <- bind_rows(refit_stats_tmb, first_stats_tmb) %>%
  mutate(metric = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(metric = factor(metric, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  mutate(panel = factor("Our Procedure", levels = c("Our Procedure", "TST-170", "F1", "MSK-I"))) %>% #, "TSO-500"))) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, colour = panel)) +
  geom_point(data = model_stats, aes(pch = model), size = 3, stroke = 1) +
  geom_line(aes(linetype = model), colour = "black", size = 1) +
  facet_wrap(~metric, labeller = label_parsed) +
  labs(x = "Panel Size (Mb)", y = "") +
  theme_minimal() + ylim(0.6,1) +
  scale_shape_discrete(solid = FALSE, name = "Models:", labels = list(TeX("$\\hat{T}$"), "ecTMB", "Count", "Linear")) +
  scale_linetype_discrete(name = TeX("Procedure:"), labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-fit $\\hat{T}$"))) +
  scale_colour_discrete(name = "Panels (left to right):") +
  theme(legend.position="bottom", legend.box = "horizontal") +
  guides(shape = guide_legend(label.position = "top"), colour = guide_legend(label.position = "top"), linetype = guide_legend(label.position = "top")) +
  scale_x_continuous(limits = c(0.2,1.5), breaks = c(0.4,0.8,1.2))

ggsave(filename = "figures/fig7.png", fig7, height = 5, width = 9)


### Figure 8
message("Creating Figure 8")



refit_predictions_tmb <- nsclc_pred_refit_tmb %>%
  get_predictions(new_data = nsclc_tables$test) %>%
  pred_intervals(pred_model = nsclc_pred_refit_tmb, biomarker_values = nsclc_tmb_values$test,
                 gen_model = nsclc_gen_model, training_matrix = nsclc_tables$train$matrix,
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 600000, biomarker = "TMB")

count_predictions_tmb <- nsclc_pred_count_tmb %>%
  get_predictions(new_data = nsclc_tables$test, max_panel_length = 600000) %>%
  {data.frame(estimated_value = .$predictions[nsclc_tmb_values$test[["Tumor_Sample_Barcode"]],],
              true_value = nsclc_tmb_values$test[["TMB"]],
              model = "Count", lower = NA, upper = NA)}

linear_predictions_tmb <- nsclc_pred_linear_tmb %>%
  get_predictions(new_data = nsclc_linear_tables$test, max_panel_length = 600000) %>%
  {data.frame(estimated_value = .$predictions[nsclc_tmb_values$test[["Tumor_Sample_Barcode"]],],
              true_value = nsclc_tmb_values$test[["TMB"]],
              model = "Linear", lower = NA, upper = NA)}

# write_tsv(data.frame(Hugo_Symbol = names(which(nsclc_pred_first_tmb$panel_genes[, max(which(nsclc_pred_first_tmb$panel_lengths <= 600000))]))),
#           file = "data/panel_0.6_genes.tsv")

# panel_0.6 <- "data/panel_0.6_bed.bed"
# sample_0.6_test = data.frame(SampleID = nsclc_tables$test$sample_list, BED = panel_0.6, stringsAsFactors = FALSE)
#
# testset_0.6 <- nsclc_maf_grch38 %>%
#   filter(Tumor_Sample_Barcode %in% nsclc_tables$test$sample_list) %>%
#   readData(exomef, covarf, mutContextf, ref, samplef = sample_0.6_test)
# write_rds(testset_0.6, "data/temporary_storage/testset_0.6")
# testset_0.6 <- read_rds("data/temporary_storage/testset_0.6")
#
# test_pred_0.6 <- pred_TMB(testset_0.6, WES = testset_WES, cores = 1, params = trainedModel, mut.nonsil = T,
#                              gid_nonsil_p = trainset$get_nonsil_passengers(0.95)) %>%
#   mutate(estimated_value = mean(nsclc_tmb_values$test$TMB) *ecTMB_panel_TMB/ mean(WES_TMB),
#          true_value = mean(nsclc_tmb_values$test$TMB)*WES_TMB / mean(WES_TMB),
#          Tumor_Sample_Barcode = sample) %>%
#   select(Tumor_Sample_Barcode, estimated_value, true_value)
# write_tsv(test_pred_0.6, "data/results/test_pred_0.6.tsv")

ectmb_predictions_tmb <- read_tsv("data/results/test_pred_0.6.tsv") %>%
  mutate(lower = NA, upper = NA, model = "ecTMB") %>%
  column_to_rownames("Tumor_Sample_Barcode")


fig8 <- bind_rows(refit_predictions_tmb$prediction_intervals, ectmb_predictions_tmb, count_predictions_tmb, linear_predictions_tmb) %>%
  mutate(model = factor(model, levels = c("Refitted T", "ecTMB", "Count", "Linear"))) %>%
  {ggplot() + geom_point(data = ., aes(x = true_value, y = estimated_value), size = 0.5) + facet_wrap(~model, nrow = 2) +
    geom_ribbon(data = refit_predictions_tmb$confidence_region, aes(x = x, ymin = y_lower, ymax = y_upper),
                alpha = 0.2, fill = "red") +
    geom_abline(colour = "blue", linetype = 2) +
    geom_hline(yintercept = 300, alpha = 0.5, linetype = 2) +
    geom_vline(xintercept = 300, alpha = 0.5, linetype = 2) +
    scale_x_log10() + scale_y_log10() + theme_minimal() + labs(x = "True TMB", y = "Predicted TMB")}

ggsave(fig8, filename = "figures/fig8.png", width = 8, height = 6)



### Section 3.2 stats
tst_170_length <- model_stats %>%
  filter(panel == "TST-170") %>%
  pull(panel_length) %>%
  unique()

tst_170_alt_length <- refit_stats_tmb %>%
  filter(panel_length <= 1000000 * tst_170_length) %>%
  pull(panel_length) %>%
  max()/1000000

tst_170_alt_r <- refit_stats_tmb %>%
  filter(panel_length == 1000000*tst_170_alt_length) %>%
  filter(metric == "R") %>%
  pull(stat)

tst_170_r <- model_stats %>%
  filter(panel == "TST-170") %>%
  filter(metric == "Regression ~ (R^2)") %>%
  filter(model == "OLM") %>%
  pull(stat)

n_tmb_h <- sum(nsclc_tmb_values$val$TMB >= 300)
n_tmb_l <- sum(nsclc_tmb_values$val$TMB < 300)

prop_tmb_h <- mean(nsclc_tmb_values$val$TMB >= 300)
prop_tmb_l <- mean(nsclc_tmb_values$val$TMB < 300)

t_0.6_r_tmb <- refit_predictions_tmb %>%
  {1 - sum((.$prediction_intervals$estimated_value - .$prediction_intervals$true_value)^2)/
    sum((.$prediction_intervals$true_value - mean(.$prediction_intervals$true_value))^2)}

count_0.6_r_tmb <- count_predictions_tmb %>%
  {1 - sum((.$estimated_value - .$true_value)^2)/
    sum((.$true_value - mean(.$true_value))^2)}

linear_0.6_r_tmb <- linear_predictions_tmb %>%
  {1 - sum((.$estimated_value - .$true_value)^2)/
    sum((.$true_value - mean(.$true_value))^2)}

ectmb_0.6_r_tmb <- ectmb_predictions_tmb %>%
  {1 - sum((.$estimated_value - .$true_value)^2)/
    sum((.$true_value - mean(.$true_value))^2)}

prop_confidence_tmb <- refit_predictions_tmb %>%
  {mean((.$prediction_intervals$upper >= .$prediction_intervals$estimated_value) &
        (.$prediction_intervals$lower <= .$prediction_intervals$estimated_value))}

s3.2.stats <- data.frame(tst_170_length = tst_170_length,
                         tst_170_alt_length = tst_170_alt_length,
                         tst_170_r = tst_170_r,
                         tst_170_alt_r = tst_170_alt_r,
                         n_tmb_h = n_tmb_h,
                         n_tmb_l = n_tmb_l,
                         prop_tmb_h = prop_tmb_h,
                         prop_tmb_l = prop_tmb_l,
                         t_0.6_r_tmb = t_0.6_r_tmb,
                         count_0.6_r_tmb = count_0.6_r_tmb,
                         linear_0.6_r_tmb = linear_0.6_r_tmb,
                         ectmb_0.6_r_tmb = ectmb_0.6_r_tmb,
                         prop_confidence_tmb = prop_confidence_tmb)

write_tsv(s3.2.stats, "data/results/s3.2.stats.tsv")



### Figure 9
message("Creating Figure 9")



first_stats_tib <- nsclc_pred_first_tib %>%
  get_predictions(new_data = nsclc_tables$val) %>%
  get_stats(biomarker_values = nsclc_tib_values$val, model = "First-fit T", threshold = 10)

refit_stats_tib <- nsclc_pred_refit_tib %>%
  get_predictions(new_data = nsclc_tables$val) %>%
  get_stats(biomarker_values = nsclc_tib_values$val, model = "Refitted T", threshold = 10)

fig9 <- bind_rows(refit_stats_tib, first_stats_tib) %>%
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line(size = 1) + ylim(0, 1) +
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") + theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))

ggsave(filename = "figures/fig9.png", plot = fig9, height = 4, width = 8)



### Figure 10
message("Creating Figure 10")

refit_predictions_tib <- nsclc_pred_refit_tib %>%
  get_predictions(new_data = nsclc_tables$test) %>%
  pred_intervals(pred_model = nsclc_pred_refit_tib, biomarker_values = nsclc_tib_values$test,
                 gen_model = nsclc_gen_model, training_matrix = nsclc_tables$train$matrix,
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 600000, biomarker = "TIB")

count_predictions_tib <- nsclc_pred_count_tib %>%
  get_predictions(new_data = nsclc_tables$test, max_panel_length = 600000) %>%
  {data.frame(estimated_value = .$predictions[nsclc_tib_values$test[["Tumor_Sample_Barcode"]],],
              true_value = nsclc_tib_values$test[["TIB"]],
              model = "Count", lower = NA, upper = NA)}

linear_predictions_tib <- nsclc_pred_linear_tib %>%
  get_predictions(new_data = nsclc_linear_tables$test, max_panel_length = 600000) %>%
  {data.frame(estimated_value = .$predictions[nsclc_tib_values$test[["Tumor_Sample_Barcode"]],],
              true_value = nsclc_tib_values$test[["TIB"]],
              model = "Linear", lower = NA, upper = NA)}

fig10 <- bind_rows(refit_predictions_tib$prediction_intervals, count_predictions_tib, linear_predictions_tib) %>%
  {ggplot() + geom_point(data = ., aes(x = true_value, y = estimated_value), size = 0.5) + facet_wrap(~model, nrow = 2) +
  geom_ribbon(data = refit_predictions_tib$confidence_region, aes(x = x, ymin = y_lower, ymax = y_upper),
              alpha = 0.2, fill = "red") +
  geom_abline(colour = "blue", linetype = 2) +
  geom_hline(yintercept = 10, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 10, alpha = 0.5, linetype = 2) +
  scale_x_log10() + scale_y_log10() + theme_minimal() + labs(x = "True TIB", y = "Predicted TIB")}

ggsave(filename = "figures/fig10.png", plot = fig10, height = 6, width = 8)



### Section 3.3 stats
n_tib_h <- sum(nsclc_tib_values$val$TIB >= 10)
n_tib_l <- sum(nsclc_tib_values$val$TIB < 10)

prop_tib_h <- mean(nsclc_tib_values$val$TIB >= 10)
prop_tib_l <- mean(nsclc_tib_values$val$TIB < 10)

t_0.6_r_tib <- refit_predictions_tib %>%
  {1 - sum((.$prediction_intervals$estimated_value - .$prediction_intervals$true_value)^2)/
      sum((.$prediction_intervals$true_value - mean(.$prediction_intervals$true_value))^2)}

count_0.6_r_tib <- count_predictions_tib %>%
  {1 - sum((.$estimated_value - .$true_value)^2)/
    sum((.$true_value - mean(.$true_value))^2)}

linear_0.6_r_tib <- linear_predictions_tib %>%
  {1 - sum((.$estimated_value - .$true_value)^2)/
    sum((.$true_value - mean(.$true_value))^2)}

prop_confidence_tib <- refit_predictions_tib %>%
  {mean((.$prediction_intervals$upper >= .$prediction_intervals$estimated_value) &
          (.$prediction_intervals$lower <= .$prediction_intervals$estimated_value))}

s3.3.stats <- data.frame(n_tib_h = n_tib_h,
                         n_tib_l = n_tib_l,
                         prop_tib_h = prop_tib_h,
                         prop_tib_l = prop_tib_l,
                         t_0.6_r_tib = t_0.6_r_tib,
                         count_0.6_r_tib = count_0.6_r_tib,
                         linear_0.6_r_tib = linear_0.6_r_tib,
                         prop_confidence_tib = prop_confidence_tib)

write_tsv(s3.3.stats, "data/results/s3.3.stats.tsv")



### Section 3.4 Stats


tst_170_genes <- read_tsv("data/tst_170_genes.tsv")$Hugo_Symbol
n_genes_tst_170 <- length(tst_170_genes)

# nsclc_pred_first_tmb_aug <- pred_first_fit(gen_model = nsclc_gen_model, lambda = exp(seq(-18, -26, length.out = 100)),
#                                         gene_lengths = ensembl_gene_lengths, training_matrix = nsclc_tables$train$matrix, 
#                                         free_genes = tst_170_genes)
# write_rds(x = nsclc_pred_first_tmb_aug, file = "data/results/nsclc_pred_first_tmb_aug")

nsclc_pred_first_tmb_aug <- read_rds("")


