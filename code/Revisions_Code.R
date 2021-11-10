library(tidyverse)
library(devtools)
library(latex2exp)
library(cowplot)
library(survival)
library(survminer)
load_all("../../ICBioMark/")

### Melanoma fits
skcm_maf <- read_tsv("data/skcm_tcga_pub_2015/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
skcm_tables <- get_mutation_tables(maf = skcm_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 250, val = 0, test = 96))

message("Getting generative model")
# skcm_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = skcm_tables$train,
#                                  progress = TRUE)
# 
# write_rds(x = skcm_gen_model, "data/pre_loaded/skcm_gen_model")
skcm_gen_model <- read_rds("data/pre_loaded/skcm_gen_model")
skcm_pred_first_tmb <- read_rds("data/pre_loaded/skcm_pred_first_tmb")

skcm_pred_first_tmb <- pred_first_fit(gen_model = skcm_gen_model,
                                      lambda = exp(seq(-18, -30, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths,
                                      training_matrix = skcm_tables$train$matrix,
                                      marker_mut_types = c("NS"))

skcm_pred_refit_tmb <- pred_refit_range(pred_first = skcm_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

skcm_tmb_values <- get_biomarker_tables(skcm_maf, biomarker = "TMB", split =  c(train = 250, val = 0, test = 96))


### Colorectal fits

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
write_rds(coadread_gen_model, "data/pre_loaded/coadread_gen_model")
coadread_gen_model <- read_rds("data/pre_loaded/coadread_gen_model")

coadread_pred_first_tmb <- pred_first_fit(gen_model = coadread_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = coadread_tables$train$matrix,
                                      marker_mut_types = c("NS"))

coadread_pred_refit_tmb <- pred_refit_range(pred_first = coadread_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

coadread_tmb_values <- get_biomarker_tables(coadread_maf, biomarker = "TMB", split =  c(train = 500.1, val = 0, test = 119))



### Bladder fits

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
write_rds(blca_gen_model, "data/pre_loaded/blca_gen_model")
blca_gen_model <- read_rds("data/pre_loaded/blca_gen_model")
blca_pred_first_tmb <- pred_first_fit(gen_model = blca_gen_model, 
                                          lambda = exp(seq(-18, -26, length.out = 100)),
                                          gene_lengths = ensembl_gene_lengths, 
                                          training_matrix = blca_tables$train$matrix,
                                          marker_mut_types = c("NS"))

blca_pred_refit_tmb <- pred_refit_range(pred_first = blca_pred_first_tmb, 
                                            gene_lengths = ensembl_gene_lengths,
                                            marker_mut_types = c("NS"))

blca_tmb_values <- get_biomarker_tables(blca_maf, biomarker = "TMB", split =  c(train = 300.1, val = 0, test = 109)) 


### Renal fits
kirc_maf <- read_tsv("data/kirc_tcga/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
kirc_tables <- get_mutation_tables(maf = kirc_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 350.1, val = 0, test = 101)) #don't ask about the .1, this needs fixing

message("Getting generative model")
kirc_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = kirc_tables$train,
                                progress = TRUE)
write_rds(kirc_gen_model, "data/pre_loaded/kirc_gen_model")
kirc_gen_model <- read_rds("data/pre_loaded/kirc_gen_model")

kirc_pred_first_tmb <- pred_first_fit(gen_model = kirc_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 200)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = kirc_tables$train$matrix,
                                      marker_mut_types = c("NS"))

kirc_pred_refit_tmb <- pred_refit_range(pred_first = kirc_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

kirc_tmb_values <- get_biomarker_tables(kirc_maf, biomarker = "TMB", split =  c(train = 350.1, val = 0, test = 101)) 

### Prostate fits

prad_maf <- read_tsv("data/prad_p1000/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
prad_tables <- get_mutation_tables(maf = prad_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 700.1, val = 0, test = 312)) #don't ask about the .1, this needs fixing

message("Getting generative model")
# prad_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = prad_tables$train,
#                                 progress = TRUE)
# write_rds(prad_gen_model, "data/pre_loaded/prad_gen_model")
prad_gen_model <- read_rds("data/pre_loaded/prad_gen_model")

prad_pred_first_tmb <- pred_first_fit(gen_model = prad_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = prad_tables$train$matrix,
                                      marker_mut_types = c("NS"))

prad_pred_refit_tmb <- pred_refit_range(pred_first = prad_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

prad_tmb_values <- get_biomarker_tables(prad_maf, biomarker = "TMB", split =  c(train = 700.1, val = 0, test = 312)) 

### Breast fits

brca_maf <- read_tsv("data/brca_tcga_pan_can_atlas_2018/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

message("Getting tables")
brca_tables <- get_mutation_tables(maf = brca_maf, 
                                   include_synonymous = FALSE,
                                   acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                   for_biomarker = "TMB",
                                   split = c(train = 700.1, val = 0, test = 309)) #don't ask about the .1, this needs fixing

message("Getting generative model")
brca_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = brca_tables$train,
                                progress = TRUE)
write_rds(brca_gen_model, "data/pre_loaded/brca_gen_model")
brca_gen_model <- read_rds("data/pre_loaded/brca_gen_model")

brca_pred_first_tmb <- pred_first_fit(gen_model = brca_gen_model, 
                                      lambda = exp(seq(-18, -26, length.out = 100)),
                                      gene_lengths = ensembl_gene_lengths, 
                                      training_matrix = brca_tables$train$matrix,
                                      marker_mut_types = c("NS"))

brca_pred_refit_tmb <- pred_refit_range(pred_first = brca_pred_first_tmb, 
                                        gene_lengths = ensembl_gene_lengths,
                                        marker_mut_types = c("NS"))

brca_tmb_values <- get_biomarker_tables(brca_maf, biomarker = "TMB", split =  c(train = 700.1, val = 0, test = 309)) 



### SKCM WES external validation
skcm_val_maf <- read_tsv("data/skcm_yale/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
skcm_val_tables <- get_mutation_tables(skcm_val_maf, 
  split = c(train = 0, val = 0, test = 91),
  gene_list = skcm_tables$train$gene_list,
  acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
  for_biomarker = "TMB",
  include_synonymous = FALSE)
skcm_val_tmb_values <- get_biomarker_tables(skcm_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 91))

skcm_refit_stats <- skcm_pred_refit_tmb %>%
  get_predictions(new_data = skcm_tables$test) %>%
  get_stats(biomarker_values = skcm_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Melanoma")

skcm_refit_stats_val <- skcm_pred_refit_tmb %>%
  get_predictions(new_data = skcm_val_tables$test) %>%
  get_stats(biomarker_values = skcm_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Melanoma")



### COADREAD WES external validation
coadread_val_maf <- read_tsv("data/coadread_genentech/data_mutations_extended.txt", comment = "#") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position) %>% 
  mutate(Tumor_Sample_Barcode = paste0("T", Tumor_Sample_Barcode))
coadread_val_tables <- get_mutation_tables(coadread_val_maf, 
                                       split = c(train = 0, val = 0, test = 72),
                                       gene_list = coadread_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)
coadread_val_tmb_values <- get_biomarker_tables(coadread_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 72))


coadread_refit_stats <- coadread_pred_refit_tmb %>%
  get_predictions(new_data = coadread_tables$test) %>%
  get_stats(biomarker_values = coadread_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Colorectal")

coadread_refit_stats_val <- coadread_pred_refit_tmb %>%
  get_predictions(new_data = coadread_val_tables$test) %>%
  get_stats(biomarker_values = coadread_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Colorectal")



### BLCA WES external validation
blca_val_maf <- read_tsv("data/blca_bgi/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
blca_val_tables <- get_mutation_tables(blca_val_maf, 
                                       split = c(train = 0, val = 0, test = 99),
                                       gene_list = blca_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)
blca_val_tmb_values <- get_biomarker_tables(blca_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 99))


blca_refit_stats <- blca_pred_refit_tmb %>%
  get_predictions(new_data = blca_tables$test) %>%
  get_stats(biomarker_values = blca_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Bladder")

blca_refit_stats_val <- blca_pred_refit_tmb %>%
  get_predictions(new_data = blca_val_tables$test) %>%
  get_stats(biomarker_values = blca_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Bladder")


### KIRC WES external validation
kirc_val_maf <- read_tsv("data/kirc_bgi/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
kirc_val_tables <- get_mutation_tables(kirc_val_maf, 
                                       split = c(train = 0, val = 0, test = 91),
                                       gene_list = kirc_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)
kirc_val_tmb_values <- get_biomarker_tables(kirc_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 91))

kirc_refit_stats <- kirc_pred_refit_tmb %>%
  get_predictions(new_data = kirc_tables$test) %>%
  get_stats(biomarker_values = kirc_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Renal")

kirc_refit_stats_val <- kirc_pred_refit_tmb %>%
  get_predictions(new_data = kirc_val_tables$test) %>%
  get_stats(biomarker_values = kirc_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Renal")

### PRAD WES external validation
prad_val_maf <- read_tsv("data/prad_fhcrc/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
prad_val_tables <- get_mutation_tables(prad_val_maf, 
                                       split = c(train = 0, val = 0, test = 141),
                                       gene_list = prad_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)
prad_val_tmb_values <- get_biomarker_tables(prad_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 141))

prad_refit_stats <- prad_pred_refit_tmb %>%
  get_predictions(new_data = prad_tables$test) %>%
  get_stats(biomarker_values = prad_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Prostate")

prad_refit_stats_val <- prad_pred_refit_tmb %>%
  get_predictions(new_data = prad_val_tables$test) %>%
  get_stats(biomarker_values = prad_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Prostate")


### BRCA WES external validation
brca_val_maf <- read_tsv("data/brca_smc_2018/data_mutations_extended.txt") %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)
brca_val_tables <- get_mutation_tables(brca_val_maf, 
                                       split = c(train = 0, val = 0, test = 185),
                                       gene_list = brca_tables$train$gene_list,
                                       acceptable_genes = ensembl_gene_lengths$Hugo_Symbol,
                                       for_biomarker = "TMB",
                                       include_synonymous = FALSE)
brca_val_tmb_values <- get_biomarker_tables(brca_val_maf, biomarker = "TMB", split =  c(train = 0, val = 0, test = 185))

brca_refit_stats <- brca_pred_refit_tmb %>%
  get_predictions(new_data = brca_tables$test) %>%
  get_stats(biomarker_values = brca_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "Internal Validation", cancer_type = "Breast")

brca_refit_stats_val <- brca_pred_refit_tmb %>%
  get_predictions(new_data = brca_val_tables$test) %>%
  get_stats(biomarker_values = brca_val_tmb_values$test, metrics = c("R"), model = "Refitted T", threshold = 300) %>% 
  mutate(Dataset = "External Test", cancer_type = "Breast")




external_validation_fig <- bind_rows(skcm_refit_stats, skcm_refit_stats_val, 
                                     coadread_refit_stats, coadread_refit_stats_val, 
                                     blca_refit_stats, blca_refit_stats_val,
                                     kirc_refit_stats, kirc_refit_stats_val,
                                     brca_refit_stats, brca_refit_stats_val,
                                     prad_refit_stats, prad_refit_stats_val) %>%
  mutate(Dataset = factor(Dataset, levels = c("Internal Validation", "External Test"))) %>% 
  mutate(model = factor(model, levels = c("Refitted T", "First-fit T"))) %>% 
  mutate(type = if_else(metric == "R", "Regression ~ (R^2)", "Classification ~ (AUPRC)")) %>%
  mutate(type = factor(type, levels = c("Regression ~ (R^2)", "Classification ~ (AUPRC)"))) %>%
  filter(panel_length <= 2000000) %>%
  mutate(panel_length = panel_length / 1000000) %>%
  ggplot(aes(x = panel_length, y = stat, colour = Dataset)) + geom_line(size = 1) + ylim(0, 1) + xlim(0.2, 1.5) +
  theme_minimal() + facet_wrap(~cancer_type, labeller = label_parsed, strip.position = "top") +
  theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = TeX("$R^2$")) +
  scale_color_manual(name = "Dataset:", values = c("black", "blue"), labels = list("Internal Validation", "External Test"))

ggsave(filename = "results/figures/external_validation_fig.png", external_validation_fig, width = 10, height = 6)


# Exploration of prostrate cancer TMB values between internal and external datasets
prad_val_internal_predictions <- prad_pred_refit_tmb %>%
  get_predictions(new_data = prad_tables$test) %>%
  pred_intervals(pred_model = prad_pred_refit_tmb, biomarker_values = prad_tmb_values$test,
                 gen_model = prad_gen_model, training_matrix = prad_tables$train$matrix,
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 1000000, marker_mut_types = c("NS")) 
prad_val_internal_predictions$prediction_intervals <- mutate(prad_val_internal_predictions$prediction_intervals, Dataset = "Internal (Armenia et al, 2018)")

prad_val_external_predictions <- prad_pred_refit_tmb %>%
    get_predictions(new_data = prad_val_tables$test) %>%
  pred_intervals(pred_model = prad_pred_refit_tmb, biomarker_values = prad_val_tmb_values$test,
                 gen_model = prad_gen_model, training_matrix = prad_tables$train$matrix,
                 gene_lengths = ensembl_gene_lengths, max_panel_length = 1000000, marker_mut_types = c("NS")) 
prad_val_external_predictions$prediction_intervals <- mutate(prad_val_external_predictions$prediction_intervals, Dataset = "External (Kumar et al, 2016)")

library(ggExtra)
p <- bind_rows(prad_val_internal_predictions$prediction_intervals,
               prad_val_external_predictions$prediction_intervals) %>%
  {ggplot(., aes(x = true_value, y = estimated_value, colour = Dataset)) + geom_point()  + theme_minimal() +
      scale_x_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3))) +
      scale_y_continuous(trans = scales::pseudo_log_trans(), breaks = c(0,10**(1:3)), limit = c(0,NA)) +
      geom_abline(colour = "grey", linetype = 2) + labs(x = "True Value", y = "Estimated Value")}

prad_int_ext_val_fig <- ggMarginal(p, margins = "x", groupFill = TRUE)


ggsave(filename = "results/figures/prad_int_ext_val.png", plot = prad_int_ext_val_fig, width = 8, height = 4)

prad_val_internal_predictions$prediction_intervals %>% 
  {1 - sum((.$true_value - .$estimated_value)^2)/sum((.$true_value - mean(.$true_value))^2)}

prad_val_external_predictions$prediction_intervals %>% 
  {1 - sum((.$true_value - .$estimated_value)^2)/sum((.$true_value - mean(.$true_value))^2)}



                

