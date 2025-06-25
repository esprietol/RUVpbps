setGeneric("gen_PBPS", function(counts, ...) standardGeneric("gen_PBPS"))

.gen_PBPS <- function(counts,
                      PBC,
                      metadata,
                      ctype = "cg_cov",
                      BioVar,
                      NVar,
                      id_pb = "sample_cell",
                      id_sub = "ind_cov",
                      cell_id = "cell_id",
                      n = 1,
                      Seed = NULL) {
  # Identify sample-level covariates
  unique_pb <- metadata |>
    dplyr::group_by(!!rlang::sym(id_pb)) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::n_distinct), .groups = "drop") |>
    dplyr::select(-!!rlang::sym(id_pb)) |>
    colSums()

  metacovs <- c(id_pb, names(unique_pb)[unique_pb == dplyr::n_distinct(metadata[[id_pb]])])

  # Grouping metadata
  meta_samples <- unique(metadata[, metacovs])
  groups_info <- meta_samples |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(BioVar, NVar)))) |>
    dplyr::summarise(nsample = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      bsg = paste('bsg',!!!rlang::syms(BioVar), sep = "_"),
      group = paste('g', !!!rlang::syms(BioVar), !!!rlang::syms(NVar), sep = "_")
    )

  bsg.keep <- names(table(groups_info$bsg)[table(groups_info$bsg) > 1])

  metadata <- metadata |>
    dplyr::mutate(
      bsg = paste('bsg',!!!rlang::syms(BioVar), sep = "_"),
      group = paste('g', !!!rlang::syms(BioVar), !!!rlang::syms(NVar), sep = "_")
    ) |>
    dplyr::filter(bsg %in% bsg.keep)

  # Sample cells per group
  cells_to_pbps <- cellspbps(subsample = metadata, seed = Seed, n = n, pb_id = id_pb, cell_id = cell_id)

  # Aggregate pseudobulk pseudosamples

  pscounts <- lapply(cells_to_pbps, function (x) counts[, x, drop = FALSE])
  pbpscounts <- sapply(pscounts, function(x) DelayedArray::rowSums(x,drop=FALSE))

  # Metadata reconstruction
  pbps_info <- tibble::tibble(!!rlang::sym(id_pb) := colnames(pbpscounts)) |>
    dplyr::mutate(
      pseudosample = sub("\\..*$", "", .data[[id_pb]]),
      group = sub("^.*?\\.", "", .data[[id_pb]])
    ) |>
    dplyr::left_join(groups_info, by = "group") |>
    dplyr::mutate(!!rlang::sym(id_sub) := paste0("pbps_", bsg)) |>
    dplyr::select(-nsample)

  # Add missing metadata columns
  miss_var <- setdiff(colnames(meta_samples), colnames(pbps_info))
  for (v in miss_var) pbps_info[[v]] <- NA

  pbps_info <- dplyr::mutate(pbps_info, pbps = 1)

  meta_samples <- meta_samples |>
    dplyr::mutate(pbps = 0, pseudosample = "orig") |>
    dplyr::bind_rows(pbps_info)

  # Combine count matrices

  full_counts <- cbind(PBC, pbpscounts[rownames(PBC), ,drop = FALSE])
  rownames(meta_samples) <- meta_samples[,id_pb]
  meta_samples <- meta_samples[colnames(full_counts),]

  return(list(
    counts = full_counts,
    metadata = meta_samples
  ))
}

#' Generate a pseudobulk dataset with pseudobulk pseudosamples to be used as negative control samples in the removal of unwanted variation
#'
#' @description This function orchestrates the process of creating pseudobulk pseudosamples (PBPS) from a `SingleCellExperiment` object and associated metadata.
#' It combines sampling, aggregation, and metadata reconstruction, and returns a list, were each element is a  `SummarizedExperiment` object with the pseudobulk counts of a particular cell type.
#' @rdname gen_PBPS
#' @param counts A `SingleCellExperiment` object with the single-cell counts and cell level metadata stored in the `colData` slot.
#' @param ctype A character string specifying the column in the `colData` `DataFrame` used to define the cell type.
#' @param BioVar A character vector specifying column names in the `colData` `DataFrame` with the biological variables (e.g., treatment).
#' @param NVar A character vector specifying column names in the `colData` `DataFrame` with technical/nuisance variables (e.g., batch ID).
#' @param id_pb A character string specifying the column in the `colData` `DataFrame` used as sample ID.
#' @param id_sub A character string specifying the column in the `colData` `DataFrame` used as "subject-level" ID (e.g., patient ID).
#' @param cell_id A character string indicating the column in the `colData` `DataFrame` used as cell ID.
#' @param n Number of pseudo-replicates to generate per group. Default is 1.
#' @param Seed Integer used to set the random seed for reproducibility.
#'
#' @return A list of `SummarizedExperiment` objects.
#'
#' @importFrom dplyr group_by summarise mutate select left_join filter bind_rows n_distinct everything all_of across
#' @importFrom tibble tibble
#' @importFrom rlang sym syms :=
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom Matrix Matrix
#' @importFrom DelayedArray rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom BiocGenerics counts
#' @examples
#'
#' data(dummysce)
#'
#' pbpsc <- gen_PBPS(dummysce, ctype='cg_cov', BioVar='cg_cov', NVar='Processing_Cohort', id_pb = 'sample_cell', id_sub = 'ind_cov', cell_id = 'cell_id', n = 1 )
#' @export

setMethod(f = "gen_PBPS",
          signature = c(counts = "SingleCellExperiment"),
          definition = function(counts,
                                ctype='cg_cov',
                                BioVar,
                                NVar,
                                id_pb = 'sample_cell',
                                id_sub = 'ind_cov',
                                cell_id = 'cell_id',
                                n = 2,
                                Seed = NULL){

            metadata <- as.data.frame(counts@colData)
            ct <- levels(metadata[,ctype])

            PBC <- scuttle::aggregateAcrossCells(counts,metadata[,id_pb])
            PBC <- Matrix::Matrix(BiocGenerics::counts(PBC), sparse = TRUE)
            counts <- BiocGenerics::counts(counts)

            full_pbc <- .gen_PBPS(counts = counts,
                                 PBC = PBC,
                                 metadata = metadata,
                                 ctype = ctype,
                                 BioVar = BioVar,
                                 NVar = NVar,
                                 id_pb = id_pb,
                                 id_sub = id_sub,
                                 cell_id = cell_id,
                                 n = n,
                                 Seed = Seed)

            counts_pbc <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=full_pbc$counts),colData=full_pbc$metadata)

            counts_pbc_list <- lapply(ct,function(x) counts_pbc[,full_pbc$metadata[,ctype]==x])

            return(counts_pbc_list)
          }
)

