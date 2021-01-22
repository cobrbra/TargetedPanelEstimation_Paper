#' Foundation Panel Gene List
#' 
#' Gene list for the FoundationOne targeted gene panel.
#' 
#' @format A dataframe with one column:
#'   \describe{
#'     \item{Hugo_Symbol}{The genes constituting the panel}
#'   }
#' @source \url{https://assets.ctfassets.net/w98cd481qyp0/YqqKHaqQmFeqc5ueQk48w/0a34fcdaa3a71dbe460cdcb01cebe8ad/F1CDx_Technical_Specifications_072020.pdf}
"foundation_genes"

#' TruSight Tumor 170 Panel Gene List
#' 
#' Gene list for the TST-170 targeted gene panel.
#' 
#' @format A dataframe with one column:
#'   \describe{
#'     \item{Hugo_Symbol}{The genes constituting the panel}
#'   }
#' @source \url{https://emea.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/trusight-tumor-170-data-sheet-1170-2016-017.pdf}
"tst_170_genes"

#' TruSight Oncology 500 Panel Gene List
#' 
#' Gene list for the TSO-500 targeted gene panel.
#' 
#' @format A dataframe with one column:
#'   \describe{
#'     \item{Hugo_Symbol}{The genes constituting the panel}
#'   }
#' @source \url{http://albiogen.ru/upload/iblock/8e6/8e66914c1b78f53509fe9186d80f8d75.pdf}
"tso_500_genes"

#' Memorial-Sloan Kettering IMPACT Panel Gene List
#' 
#' Gene list for the MSK-IMPACT targeted gene panel.
#' 
#' @format A dataframe with one column:
#'   \describe{
#'     \item{Hugo_Symbol}{The genes constituting the panel}
#'   }
#' @source \url{https://www.oncokb.org/cancerGenes}
"tso_500_genes"


#' Foundation Panel BED
#' 
#' BED file (exon locations) for the FoundationOne targeted gene panel.
#' 
#' @format A dataframe with four columns:
#'   \describe{
#'     \item{Chromosome}{The chromosomal location of the exon.}
#'     \item{Start}{The start nucleotide base of the exon}
#'     \item{End}{The end nucleotide base of the exon}
#'     \item{name}{Exon name}
#'   }
#' @source Produced from foundation_genes dataframe using the UCSC table 
#' browser at \url{https://genome.ucsc.edu/cgi-bin/hgTables}
"foundation_bed"

#' TruSight Tumor 170 Panel BED
#' 
#' BED file (exon locations) for the TST-170 targeted gene panel.
#' 
#' @format A dataframe with four columns:
#'   \describe{
#'     \item{Chromosome}{The chromosomal location of the exon.}
#'     \item{Start}{The start nucleotide base of the exon}
#'     \item{End}{The end nucleotide base of the exon}
#'     \item{name}{Exon name}
#'   }
#' @source Produced from tso_170_genes dataframe using the UCSC table 
#' browser at \url{https://genome.ucsc.edu/cgi-bin/hgTables}
"tst_170_bed"

#' TruSight Oncology Panel BED
#' 
#' BED file (exon locations) for the TSO-500 targeted gene panel.
#' 
#' @format A dataframe with four columns:
#'   \describe{
#'     \item{Chromosome}{The chromosomal location of the exon.}
#'     \item{Start}{The start nucleotide base of the exon}
#'     \item{End}{The end nucleotide base of the exon}
#'     \item{name}{Exon name}
#'   }
#' @source Produced from tso_500_genes dataframe using the UCSC table 
#' browser at \url{https://genome.ucsc.edu/cgi-bin/hgTables}
"tso_500_bed"

#' MSK-IMPACT Panel BED
#' 
#' BED file (exon locations) for the MSK-IMPACT targeted gene panel.
#' 
#' @format A dataframe with four columns:
#'   \describe{
#'     \item{Chromosome}{The chromosomal location of the exon.}
#'     \item{Start}{The start nucleotide base of the exon}
#'     \item{End}{The end nucleotide base of the exon}
#'     \item{name}{Exon name}
#'   }
#' @source Produced from msk_impact_genes dataframe using the UCSC table 
#' browser at \url{https://genome.ucsc.edu/cgi-bin/hgTables}
"msk_impact_bed"

#' GRCh38 Build Non-Small Cell Lung Cancer Data
#' 
#' The same data as is loaded with the ICBioMark package,
#' with some extra columns and lifted to the GRCh38 genome 
#' build for use with ecTMB.
#' 
#' Produced with the following code:
#' #install.pacakages("TCGAretriever")
#' library(TCGAretriever)
#' library(ICBioMark)
#' library(dplyr)
#' library(magrittr)
#' library(readr)
#' library(stringr)
#' 
#' fetch_all_tcgadata(case_id = "nsclc_tcga_broad_2016_sequenced", gprofile_id = "nsclc_tcga_broad_2016_mutations", glist = ensembl_gene_lengths$Hugo_Symbol, mutations = TRUE) %>%
#'   select(Hugo_Symbol = gene_symbol, Entrez_Gene_Id = entrez_gene_id, Gene = entrez_gene_id, Chromosome = chr,  Start_Position = start_position, End_Position = end_position,
#'          Consequence = mutation_type, Variant_Classification = mutation_type, Reference_Allele = reference_allele, Tumor_Seq_Allele1 = variant_allele, Tumor_Seq_Allele2 = variant_allele, 
#'          Tumor_Sample_Barcode = case_id) %>%
#'   mutate(Protein_position = NA, Strand = NA, Consequence = tolower(Consequence)) %>%
#'   mutate(Variant_Type = if_else(Variant_Classification %in% c("Frame_Shift_Del", "In_Frame_Del"), "DEL", "SNP")) %>%
#'   mutate(Variant_Type = if_else(Variant_Classification %in% c("Frame_Shift_Ins", "In_Frame_Ins"), "INS", Variant_Type)) %>% 
#'   mutate(Variant_Type = if_else(Variant_Type == "SNP" & str_length(Reference_Allele) == 2, "DNP", Variant_Type)) %>% 
#'   mutate(Variant_Type = if_else(Variant_Type == "SNP" & str_length(Reference_Allele) == 3, "TNP", Variant_Type)) %>% 
#'   mutate(Chromosome = if_else(Chromosome == "23", "X", Chromosome), Chromosome = if_else(Chromosome == 24, "Y", Chromosome)) %>% 
#'   filter(Chromosome != "NA") %>% 
#'   write_tsv("data/nsclc_maf_full.tsv")
#' 
#' system("python3 data/maf_conversion.py")

#' 