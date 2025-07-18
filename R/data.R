#' Multiplexed single cell RNA sequencing of healthy controls
#'
#' A subset of 15751 T cells from 8 human PBMCs samples. 2 of the samples belong to the same individual.
#'
#'
#' @format ## `dummysce`
#' A `SingleCellExperiment` object with 3000 rows (genes) and 15751 columns(cells), with 19 technical and biological variables available in the colData slot.
#' \describe{
#'   \item{nCount_originalexp}{}
#'   \item{nFeature_originalexp}{}
#'   \item{batch_cov}{Batch identifier}
#'   \item{ind_cov}{Subject identifier, from which the sample was taken}
#'   \item{Processing_Cohort}{}
#'   \item{louvain}{}
#'   \item{cg_cov}{}
#'   \item{ct_cov}{}
#'   \item{L3}{}
#'   \item{ind_cov_batch_cov}{Unique sample identifier}
#'   \item{Age}{Age}
#'   \item{Sex}{Sex}
#'   \item{pop_cov}{Ethnicity}
#'   \item{Status}{}
#'   \item{SLE_status}{}
#'   \item{mitoPercent}{Percentage of }
#'   \item{lab}{Source of the sample identifier, }
#'   \item{cell_id}{Cell type identifier}
#'   \item{sample_cell}{Sample-celltype identifier}
#' }
#' @source The subset was taken from the dataset "Multiplexed RNA-sequencing of 1M immune cells reveals the cellular, molecular, and genetic correlates of systemic lupus erythematosus", available at <https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2>
"dummysce"
