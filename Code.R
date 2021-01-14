## devtools Package
#install.packages("devtools")
library(devtools)

## ICBioMark Package
#devtools::install_github("cobrbra/ICBioMark", upgrade = "ask")
library(ICBioMark)

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

## Figures path
fig_path <- "figures/"



### Main workflow
nsclc_tables <- get_mutation_tables(maf = nsclc_maf, include_synonymous = FALSE,
                                    acceptable_genes = ensembl_gene_lengths$Hugo_Symbol)
# nsclc_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = nsclc_tables$train,
#                                  progress = TRUE)
# write_rds(x = nsclc_gen_model, file = "data/results/nsclc_gen_model")

nsclc_gen_model <- read_rds("data/results/nsclc_gen_model")

# nsclc_pred_first_tmb <- pred_first_fit(gen_model = nsclc_gen_model, lambda = exp(seq(-18, -26, length.out = 100)), 
#                                         gene_lengths = ensembl_gene_lengths, training_matrix = nsclc_tables$train$matrix)
# write_rds(x = nsclc_pred_first_tmb, file = "data/results/nsclc_pred_first_tmb")

nsclc_pred_first_tmb <- read_rds("data/results/nsclc_pred_first_tmb")
nsclc_pred_refit_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb, gene_lengths = ensembl_gene_lengths)

nsclc_tmb_values <- get_biomarker_tables(nsclc_maf, biomarker = "TMB")


# nsclc_pred_refit_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb,
                                         # gene_lengths = ensembl_gene_lengths)


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
  ggplot(aes(x = panel_length, y = stat, linetype = model)) + geom_line() + ylim(0, 1) + 
  theme_minimal() + facet_wrap(~type, labeller = label_parsed, strip.position = "top") + theme(legend.position = "bottom") + labs(x = "Panel Size (Mb)", y = "") +
  scale_linetype(name = "Procedure:", labels = list(TeX("Refitted $\\hat{T}$"), TeX("First-Fit $\\hat{T}$")))

ggsave(filename = "figures/fig6.png", plot = fig6, height = 4, width = 8)
