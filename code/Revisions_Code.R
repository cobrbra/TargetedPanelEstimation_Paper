library(tidyverse)
library(devtools)
library(latex2exp)
library(cowplot)
library(survival)
library(survminer)
load_all("../../ICBioMark/")

### Melanoma analysis
skcm_maf <- read_tsv("data/skcm_tcga_pub_2015/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
skcm_tables <- get_mutation_tables(maf = skcm_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 250, val = 0, test = 96))

message("Getting generative model")
skcm_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = skcm_tables$train,
                                 progress = TRUE)

skcm_pred_first_tmb <- pred_first_fit(gen_model = skcm_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = skcm_tables$train$matrix,
                                      marker_mut_types = c("NS"))

write_rds(x = skcm_pred_first_tmb, "skcm_pred_first_tmb")
write_rds(x = skcm_gen_model, "skcm_gen_model")

skcm_pred_refit_tmb <- pred_refit_range(pred_first = skcm_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

skcm_tmb_values <- get_biomarker_tables(skcm_maf, biomarker = "TMB", split =  c(train = 250, val = 0, test = 96))



first_stats_tmb <- skcm_pred_first_tmb %>%
  get_predictions(new_data = skcm_tables$test) %>%
  get_stats(biomarker_values = skcm_tmb_values$test, model = "First-fit T", threshold = 300)

refit_stats_tmb <- skcm_pred_refit_tmb %>%
  get_predictions(new_data = skcm_tables$test) %>%
  get_stats(biomarker_values = skcm_tmb_values$test, model = "Refitted T", threshold = 300)

skcm_version_fig6 <- bind_rows(refit_stats_tmb, first_stats_tmb) %>%
  mutate(model = factor(model, levels = c("Refitted T", "First-fit T"))) %>% 
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  filter(panel_length <= 2000000) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line(size = 1) + ylim(0, 1) +
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") +
  theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))
ggsave(skcm_version_fig6, filename = "results/figures/skcm_version_fig8.png", width = 8, height = 4)


skcm_refit_predictions_tmb <- skcm_pred_refit_tmb %>%
  get_predictions(new_data = skcm_tables$test) %>%
  pred_intervals(pred_model = skcm_pred_refit_tmb, biomarker_values = skcm_tmb_values$test,
                 gen_model = skcm_gen_model, training_matrix = skcm_tables$train$matrix, marker_mut_types = c("NS"),
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 600000, biomarker = "TMB")

skcm_r2 <- 1 - sum((skcm_refit_predictions_tmb$prediction_intervals$estimated_value - skcm_refit_predictions_tmb$prediction_intervals$true_value)^2)/
  sum((skcm_refit_predictions_tmb$prediction_intervals$true_value - mean(skcm_refit_predictions_tmb$prediction_intervals$true_value))^2)

skcm_version_fig8 <- skcm_refit_predictions_tmb$prediction_intervals %>%
  {ggplot() + geom_point(data = ., aes(x = true_value, y = estimated_value), size = 0.5) + 
      geom_ribbon(data = skcm_refit_predictions_tmb$confidence_region %>% mutate(model = factor(model, levels = c("Refitted T", "ecTMB", "Count", "Linear"))), 
                  aes(x = x, ymin = y_lower, ymax = y_upper),
                  alpha = 0.2, fill = "red")  +
      geom_abline(colour = "blue", linetype = 2) +
      geom_hline(yintercept = 300, alpha = 0.5, linetype = 2) +
      geom_vline(xintercept = 300, alpha = 0.5, linetype = 2) +
      scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3))) +
      scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3)), limit = c(0,NA)) +
      theme_minimal() + labs(x = "True TMB", y = "Predicted TMB")}
ggsave(skcm_version_fig8, filename = "results/figures/skcm_version_fig8.png", width = 8, height = 6)


### Colorectal analysis

coadread_maf <- read_tsv("data/coadread_dfci_2016/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
coadread_tables <- get_mutation_tables(maf = coadread_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 500.1, val = 0, test = 119)) #don't ask about the .1, this needs fixing

message("Getting generative model")
coadread_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = coadread_tables$train,
                                progress = TRUE)

coadread_pred_first_tmb <- pred_first_fit(gen_model = coadread_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = coadread_tables$train$matrix,
                                      marker_mut_types = c("NS"))

coadread_pred_refit_tmb <- pred_refit_range(pred_first = coadread_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

coadread_tmb_values <- get_biomarker_tables(coadread_maf, biomarker = "TMB", split =  c(train = 500.1, val = 0, test = 119))



first_stats_tmb <- coadread_pred_first_tmb %>%
  get_predictions(new_data = coadread_tables$test) %>%
  get_stats(biomarker_values = coadread_tmb_values$test, model = "First-fit T", threshold = 300)

refit_stats_tmb <- coadread_pred_refit_tmb %>%
  get_predictions(new_data = coadread_tables$test) %>%
  get_stats(biomarker_values = coadread_tmb_values$test, model = "Refitted T", threshold = 300)

coadread_version_fig6 <- bind_rows(refit_stats_tmb, first_stats_tmb) %>%
  mutate(model = factor(model, levels = c("Refitted T", "First-fit T"))) %>% 
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  filter(panel_length <= 2000000) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line(size = 1) + ylim(0, 1) +
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") +
  theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))
ggsave(coadread_version_fig6, filename = "results/figures/coadread_version_fig8.png", width = 8, height = 4)


coadread_refit_predictions_tmb <- coadread_pred_refit_tmb %>%
  get_predictions(new_data = coadread_tables$test) %>%
  pred_intervals(pred_model = coadread_pred_refit_tmb, biomarker_values = coadread_tmb_values$test,
                 gen_model = coadread_gen_model, training_matrix = coadread_tables$train$matrix, marker_mut_types = c("NS"),
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 600000, biomarker = "TMB")

coadread_r2 <- 1 - sum((coadread_refit_predictions_tmb$prediction_intervals$estimated_value - coadread_refit_predictions_tmb$prediction_intervals$true_value)^2)/
  sum((coadread_refit_predictions_tmb$prediction_intervals$true_value - mean(coadread_refit_predictions_tmb$prediction_intervals$true_value))^2)

coadread_version_fig8 <- coadread_refit_predictions_tmb$prediction_intervals %>%
  {ggplot() + geom_point(data = ., aes(x = true_value, y = estimated_value), size = 0.5) + 
      geom_ribbon(data = coadread_refit_predictions_tmb$confidence_region %>% mutate(model = factor(model, levels = c("Refitted T", "ecTMB", "Count", "Linear"))), 
                  aes(x = x, ymin = y_lower, ymax = y_upper),
                  alpha = 0.2, fill = "red")  +
      geom_abline(colour = "blue", linetype = 2) +
      geom_hline(yintercept = 300, alpha = 0.5, linetype = 2) +
      geom_vline(xintercept = 300, alpha = 0.5, linetype = 2) +
      scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3))) +
      scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3)), limit = c(0,NA)) +
      theme_minimal() + labs(x = "True TMB", y = "Predicted TMB")}

ggsave(coadread_version_fig8, filename = "results/figures/coadread_version_fig8.png", width = 8, height = 6)

### Bladder analsyis

blca_maf <- read_tsv("data/blca_tcga_pan_can_atlas_2018/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
blca_tables <- get_mutation_tables(maf = blca_maf, 
                                       include_synonymous = FALSE,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       split = c(train = 300.1, val = 0, test = 109)) #don't ask about the .1, this needs fixing

message("Getting generative model")
blca_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = blca_tables$train,
                                    progress = TRUE)

blca_pred_first_tmb <- pred_first_fit(gen_model = blca_gen_model, 
                                          lambda = exp(seq(-18, -26, length.out = 100)),
                                          gene_lengths = ensembl_gene_lengths, 
                                          training_matrix = blca_tables$train$matrix,
                                          marker_mut_types = c("NS"))

blca_pred_refit_tmb <- pred_refit_range(pred_first = blca_pred_first_tmb, 
                                            gene_lengths = ensembl_gene_lengths,
                                            marker_mut_types = c("NS"))

blca_tmb_values <- get_biomarker_tables(blca_maf, biomarker = "TMB", split =  c(train = 500.1, val = 0, test = 119))



first_stats_tmb <- blca_pred_first_tmb %>%
  get_predictions(new_data = blca_tables$test) %>%
  get_stats(biomarker_values = blca_tmb_values$test, model = "First-fit T", threshold = 300)

refit_stats_tmb <- blca_pred_refit_tmb %>%
  get_predictions(new_data = blca_tables$test) %>%
  get_stats(biomarker_values = blca_tmb_values$test, model = "Refitted T", threshold = 300)

blca_version_fig6 <- bind_rows(refit_stats_tmb, first_stats_tmb) %>%
  mutate(model = factor(model, levels = c("Refitted T", "First-fit T"))) %>% 
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  filter(panel_length <= 2000000) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line(size = 1) + ylim(0, 1) +
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") +
  theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))
ggsave(blca_version_fig6, filename = "results/figures/blca_version_fig8.png", width = 8, height = 4)


blca_refit_predictions_tmb <- blca_pred_refit_tmb %>%
  get_predictions(new_data = blca_tables$test) %>%
  pred_intervals(pred_model = blca_pred_refit_tmb, biomarker_values = blca_tmb_values$test,
                 gen_model = blca_gen_model, training_matrix = blca_tables$train$matrix, marker_mut_types = c("NS"),
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 600000, biomarker = "TMB")

blca_r2 <- 1 - sum((blca_refit_predictions_tmb$prediction_intervals$estimated_value - blca_refit_predictions_tmb$prediction_intervals$true_value)^2)/
  sum((blca_refit_predictions_tmb$prediction_intervals$true_value - mean(blca_refit_predictions_tmb$prediction_intervals$true_value))^2)

blca_version_fig8 <- blca_refit_predictions_tmb$prediction_intervals %>%
  {ggplot() + geom_point(data = ., aes(x = true_value, y = estimated_value), size = 0.5) + 
      geom_ribbon(data = blca_refit_predictions_tmb$confidence_region %>% mutate(model = factor(model, levels = c("Refitted T", "ecTMB", "Count", "Linear"))), 
                  aes(x = x, ymin = y_lower, ymax = y_upper),
                  alpha = 0.2, fill = "red")  +
      geom_abline(colour = "blue", linetype = 2) +
      geom_hline(yintercept = 300, alpha = 0.5, linetype = 2) +
      geom_vline(xintercept = 300, alpha = 0.5, linetype = 2) +
      scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3))) +
      scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3)), limit = c(0,NA)) +
      theme_minimal() + labs(x = "True TMB", y = "Predicted TMB")}

ggsave(blca_version_fig8, filename = "results/figures/blca_version_fig8.png", width = 8, height = 6)


#### Joint Plots
joint_version_fig6 <- plot_grid(skcm_version_fig6 + ggtitle("Melanoma (SKCM)"),
                                coadread_version_fig6 + ggtitle("Colorectal Cancer (COADREAD)"),
                                blca_version_fig6 + ggtitle("Bladder Cancer (BLCA)"),
                                nrow = 3,
                                labels = "AUTO")
ggsave(joint_version_fig6, filename = "results/figures/joint_version_fig6.png", height = 10, width = 5)


joint_version_fig8 <- plot_grid(skcm_version_fig8 + ggtitle("Melanoma (SKCM)") + geom_label(data = data.frame(x = 1000, y = 10, label = round(sqrt(skcm_r2), 2)), aes(x =x, y= y, label = paste("R =", label))),
                                coadread_version_fig8 + ggtitle("Colorectal Cancer (COADREAD)") + geom_label(data = data.frame(x = 750, y = 10, label = round(sqrt(coadread_r2), 2)), aes(x =x, y= y, label = paste("R =", label))),
                                blca_version_fig8 + ggtitle("Bladder Cancer (BLCA)") + geom_label(data = data.frame(x = 470, y = 10, label = round(sqrt(blca_r2), 2)), aes(x =x, y= y, label = paste("R =", label))),
                                nrow = 3,
                                labels = "AUTO")
ggsave(joint_version_fig8, filename = "results/figures/joint_version_fig8.png", height = 9, width = 4)


### MSK-IMPACT External Validation
tmb_val_maf <- read_tsv("data/tmb_mskcc_2018/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
tmb_val_patient <- read_tsv("data/tmb_mskcc_2018/data_clinical_patient.txt", comment = "#")
tmb_val_sample <- read_tsv("data/tmb_mskcc_2018/data_clinical_sample.txt", comment = "#")
tmb_val_clinical <- inner_join(tmb_val_sample, tmb_val_patient, by = "PATIENT_ID") %>% 
  select(SAMPLE_ID, CANCER_TYPE, TMB_SCORE, OS_MONTHS, OS_STATUS, DRUG_TYPE)
msk_impact_genes <- read_tsv("data/msk_impact_genes.tsv")

skcm_val_tables <- get_mutation_tables(tmb_val_maf, 
                                   split = c(train = 0, val = 0, test = 320), 
                                   sample_list = filter(tmb_val_clinical, CANCER_TYPE == "Melanoma")$SAMPLE_ID, 
                                   gene_list = skcm_tables$train$gene_list,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol, 
                                   for_biomarker = "TMB",
                                   include_synonymous = FALSE)

skcm_msk_predictions <- pred_refit_panel(pred_first = skcm_pred_first_tmb, 
                 gene_lengths = ensembl_gene_lengths, 
                 genes = msk_impact_genes$Hugo_Symbol,
                 marker_mut_types = c("NS")) %>% 
  get_predictions(new_data = skcm_val_tables$test) %>% 
  {as.data.frame(.$predictions)}

colnames(skcm_msk_predictions) <- c("Estimated_TMB")
skcm_msk_predictions$Tumor_Sample_Barcode <- rownames(skcm_msk_predictions) 
skcm_msk_predictions <- skcm_msk_predictions %>% 
  inner_join(tmb_val_clinical, by = c("Tumor_Sample_Barcode"="SAMPLE_ID")) 



skcm_surv_25 <- skcm_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.25), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.25), "High", "Low"))
skcm_surv_50 <- skcm_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.5), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.5), "High", "Low"))
skcm_surv_75 <- skcm_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.75), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.75), "High", "Low"))

forest_25_est <- coxph(surv ~ `Estimated TMB`, data = skcm_surv_25) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "25% Quantile")
forest_25_true <- coxph(surv ~ `Original TMB`, data = skcm_surv_25) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_50_est <- coxph(surv ~ `Estimated TMB`, data = skcm_surv_50) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "50% Quantile")
forest_50_true <- coxph(surv ~ `Original TMB`, data = skcm_surv_50) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_75_est <- coxph(surv ~ `Estimated TMB`, data = skcm_surv_75) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "75% Quantile")
forest_75_true <- coxph(surv ~ `Original TMB`, data = skcm_surv_75) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")

skcm_forests <- plot_grid(forest_25_est, forest_25_true, forest_50_est, forest_50_true, forest_75_est, forest_75_true, ncol = 1)
ggsave(skcm_forests, filename = "results/figures/skcm_forests.png", width = 7, height = 12)

coadread_val_tables <- get_mutation_tables(tmb_val_maf, 
                                       split = c(train = 0, val = 0, test = 110), 
                                       sample_list = filter(tmb_val_clinical, CANCER_TYPE == "Colorectal Cancer")$SAMPLE_ID, 
                                       gene_list = coadread_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol, 
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)

coadread_msk_predictions <- pred_refit_panel(pred_first = coadread_pred_first_tmb, 
                                         gene_lengths = ensembl_gene_lengths, 
                                         genes = msk_impact_genes$Hugo_Symbol,
                                         marker_mut_types = c("NS")) %>% 
  get_predictions(new_data = coadread_val_tables$test) %>% 
  {as.data.frame(.$predictions)}

colnames(coadread_msk_predictions) <- c("Estimated_TMB")
coadread_msk_predictions$Tumor_Sample_Barcode <- rownames(coadread_msk_predictions) 
coadread_msk_predictions <- coadread_msk_predictions %>% 
  inner_join(tmb_val_clinical, by = c("Tumor_Sample_Barcode"="SAMPLE_ID")) 



coadread_surv_25 <- coadread_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.25), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.25), "High", "Low"))
coadread_surv_50 <- coadread_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.5), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.5), "High", "Low"))
coadread_surv_75 <- coadread_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.75), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.75), "High", "Low"))

forest_25_est <- coxph(surv ~ `Estimated TMB`, data = coadread_surv_25) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "25% Quantile")
forest_25_true <- coxph(surv ~ `Original TMB`, data = coadread_surv_25) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_50_est <- coxph(surv ~ `Estimated TMB`, data = coadread_surv_50) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "50% Quantile")
forest_50_true <- coxph(surv ~ `Original TMB`, data = coadread_surv_50) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_75_est <- coxph(surv ~ `Estimated TMB`, data = coadread_surv_75) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "75% Quantile")
forest_75_true <- coxph(surv ~ `Original TMB`, data = coadread_surv_75) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")

coadread_forests <- plot_grid(forest_25_est, forest_25_true, forest_50_est, forest_50_true, forest_75_est, forest_75_true, ncol = 1)
ggsave(coadread_forests, filename = "results/figures/coadread_forests.png", width = 7, height = 12)

blca_val_tables <- get_mutation_tables(tmb_val_maf, 
                                           split = c(train = 0, val = 0, test = 110), 
                                           sample_list = filter(tmb_val_clinical, CANCER_TYPE == "Colorectal Cancer")$SAMPLE_ID, 
                                           gene_list = blca_tables$train$gene_list,
                                           acceptable_genes = ensembl_gene_lengths$Hugo_Symbol, 
                                           for_biomarker = "TMB",
                                           include_synonymous = FALSE)

blca_msk_predictions <- pred_refit_panel(pred_first = blca_pred_first_tmb, 
                                             gene_lengths = ensembl_gene_lengths, 
                                             genes = msk_impact_genes$Hugo_Symbol,
                                             marker_mut_types = c("NS")) %>% 
  get_predictions(new_data = blca_val_tables$test) %>% 
  {as.data.frame(.$predictions)}

colnames(blca_msk_predictions) <- c("Estimated_TMB")
blca_msk_predictions$Tumor_Sample_Barcode <- rownames(blca_msk_predictions) 
blca_msk_predictions <- blca_msk_predictions %>% 
  inner_join(tmb_val_clinical, by = c("Tumor_Sample_Barcode"="SAMPLE_ID")) 



blca_surv_25 <- blca_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.25), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.25), "High", "Low"))
blca_surv_50 <- blca_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.5), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.5), "High", "Low"))
blca_surv_75 <- blca_msk_predictions %>% 
  {mutate(., surv = Surv(time = .$OS_MONTHS, event = (.$OS_STATUS == "1:DECEASED")))} %>% 
  mutate("Estimated TMB" = if_else(Estimated_TMB >= quantile(Estimated_TMB, 0.75), "High", "Low"),
         "Original TMB" = if_else(TMB_SCORE >= quantile(TMB_SCORE, 0.75), "High", "Low"))

forest_25_est <- coxph(surv ~ `Estimated TMB`, data = blca_surv_25) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "25% Quantile")
forest_25_true <- coxph(surv ~ `Original TMB`, data = blca_surv_25) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_50_est <- coxph(surv ~ `Estimated TMB`, data = blca_surv_50) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "50% Quantile")
forest_50_true <- coxph(surv ~ `Original TMB`, data = blca_surv_50) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")
forest_75_est <- coxph(surv ~ `Estimated TMB`, data = blca_surv_75) %>% 
  ggforest(cpositions = c(0.0, 0.15,0.38), main = "75% Quantile")
forest_75_true <- coxph(surv ~ `Original TMB`, data = blca_surv_75) %>% 
  ggforest(cpositions = c(0.0,0.15,0.38), main = "")

blca_forests <- plot_grid(forest_25_est, forest_25_true, forest_50_est, forest_50_true, forest_75_est, forest_75_true, ncol = 1)
ggsave(blca_forests, filename = "results/figures/blca_forests.png", width = 7, height = 12)


