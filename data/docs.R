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
#' Genomic data from the same study as is loaded with the ICBioMark package,
#' with some extra columns, including synonymous variants (which are not available 
#' from the automated download pipeline used in the ICBioMark package) and lifted 
#' to the GRCh38 genome build for use with ecTMB. Switch to GRCh38 was done with 
#' the python package LiftOver. Some very minor alterations had to be made to be usable with the
#' ecTMB package, namely excluding mutations with Start_Position labels larger than
#' their End_Position labels (affecting fewer than 100 mutations).
#' 
#' @source \url{https://www.cbioportal.org/study/summary?id=nsclc_tcga_broad_2016}