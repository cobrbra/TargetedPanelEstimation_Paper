## devtools Package
#install.packages("devtools")
library(devtools)

## ICBioMark Package
#devtools::install_github("cobrbra/ICBioMark", upgrade = "ask")
library(ICBioMark)

## Other R packages needed
#install.packages("cowplot")
library(cowplot)




## Figures path
fig_path <- "InProgress/TargetedPanelDesign/figures/"



### Main workflow
nsclc_tables <- get_mutation_tables(maf = nsclc_maf, include_synonymous = FALSE)
nsclc_gen_model <- fit_gen_model(gene_lengths = ensembl_gene_lengths, table = nsclc_tables$train,
                                 progress = TRUE)
nsclc_pred_first_tmb <- pred_first_fit(gen_model = nsclc_gen_model, gene_lengths = ensembl_gene_lengths,
                                       training_matrix = nsclc_tables$train$matrix)
nsclc_pred_refit_tmb <- pred_refit_range(pred_first = nsclc_pred_first_tmb,
                                         gene_lengths = ensembl_gene_lengths)


### Figure 1



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
plot(fig1p1)

fig1p2 <- nsclc_survival %>% 
  group_by(SMOKING_HISTORY, STAGE) %>% 
  count() %>% 
  ungroup() %>% 
  mutate(Stage = stage_groups[STAGE]) %>% 
  filter(!is.na(Stage)) %>% 
  mutate(Stage = factor(Stage, levels = rev(c("I","II","III","IV")))) %>% 
  mutate(Smoking_History = smoking_groups[SMOKING_HISTORY])  %>%  
  filter(!is.na(Smoking_History)) %>% 
  select(Stage, Smoking_History, n) %>% 
  group_by(Stage, Smoking_History) %>% 
  mutate(n = sum(n)) %>% 
  distinct() %>% 
  ggplot(aes(x = Smoking_History, y = n, fill = Stage, order = Stage)) + 
  geom_col(position = position_stack(reverse = TRUE), alpha = 0.8, colour = "black") + scale_fill_brewer(palette = "Greys") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  labs(y = "Frequency") + theme(axis.title.x=element_blank())

plot(fig1p2)
fig1 <- plot_grid(fig1p1, fig1p2, labels = "AUTO", align = "h")
print(fig1)
ggsave(paste0(fig_path, "fig1.png"), fig1, width = 12, height = 5)



### Figure 2



tmb_tib_train <- inner_join(get_biomarker_tables(maf = nsclc_maf, biomarker = "TMB")$train, get_biomarker_tables(maf = nsclc_maf, biomarker = "TIB")$train)

## Subfigure 1
fig2p1 <- tmb_tib_train %>%
  gather(key = "Biomarker", value = "Count", -Tumor_Sample_Barcode) %>%
  mutate(Biomarker = factor(Biomarker, levels = c("TMB", "TIB"))) %>%
  ggplot(aes(x = Biomarker, y = Count)) +
  geom_violin(fill = "black", alpha = 0.5, colour = "black") +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Mutation Count") + theme(axis.title.x = element_blank())

print(fig2p1)

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

print(fig2p2)

fig2 <- plot_grid(fig2p1, fig2p2, labels = "AUTO", align = "h")
print(fig2)
ggsave(paste0(fig_path, "fig2.png"), fig2, width = 12, height = 5)




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
print(paste("TMB Mean:", tmb_mean))
print(paste("TIB Mean:", tib_mean))



### Figure 3