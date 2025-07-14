
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
                      seed = NULL) {
  # Identify sample-level covariates
  unique_pb <- metadata |>
    dplyr::group_by(!!rlang::sym(id_pb)) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::n_distinct), .groups = "drop") |>
    dplyr::select(-!!rlang::sym(id_pb)) |>
    colSums()

  metacovs <- c(id_pb, names(unique_pb)[unique_pb == dplyr::n_distinct(metadata[[id_pb]])])
  c_vars <- c(ctype,BioVar,NVar,id_pb,id_sub)
  problematic_vars <- c_vars[!c_vars %in% metacovs]

  if( length(problematic_vars)>0){
    stop("Consider providing unique sample identifiers. The following variables have multiple categories/levels associated with the same sample ID: ", paste(problematic_vars, collapse = ", "))
  }

  # Grouping metadata

  meta_samples <- unique(metadata[, metacovs])
  groups_info <- meta_samples |> # TODO: should meta_samples be replaced by metadata?
    dplyr::group_by(dplyr::across(dplyr::all_of(c(BioVar, NVar)))) |>
    dplyr::summarise(nsample = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      bsg = paste('bsg',!!!rlang::syms(BioVar), sep = "_"),
      group = paste('g', !!!rlang::syms(BioVar), !!!rlang::syms(NVar), sep = "_")
    )

  bsg.keep <- names(table(groups_info$bsg)[table(groups_info$bsg) > 1])

  if (any(table(groups_info$bsg) < 2)) {
    warning("Some biological groups have fewer than 2 samples and will not have pseudoreplicates.")
  }


  metadata <- metadata |>
    dplyr::mutate(
      bsg = paste('bsg',!!!rlang::syms(BioVar), sep = "_"),
      group = paste('g', !!!rlang::syms(BioVar), !!!rlang::syms(NVar), sep = "_")
     )  |>
    dplyr::filter(bsg %in% bsg.keep)

  # Sample cells per group
  cells_to_pbps <- cellspbps(subsample = metadata, seed = seed, n = n, pb_id = id_pb, cell_id = cell_id)

  # Aggregate pseudobulk pseudosamples
  pscounts <- lapply(cells_to_pbps, function (x) counts[, x, drop = FALSE])

  pbpscounts <- sapply(pscounts, function(x) Matrix::rowSums(x))


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

  # TODO: this errors for me.
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

#' @title Generate a pseudobulk dataset with pseudobulk pseudosamples to be used as negative control samples in the removal of unwanted variation
#'
#' @description This function orchestrates the process of creating pseudobulk pseudosamples (PBPS) from a `SingleCellExperiment` object and associated cell-level metadata.
#' It combines sampling, aggregation, and metadata reconstruction, and returns a list, were each element is a  `SummarizedExperiment` object containing the pseudobulk counts of a particular cell type.
#' @rdname gen_PBPS
#' @param counts A `SingleCellExperiment` object with the single-cell counts and cell-level metadata stored in the `colData` slot.
#' @param ctype A character string specifying the variable name in the `colData` `DataFrame` used to define the cell type.
#' @param BioVar A character vector specifying the biological variable names in the `colData` `DataFrame` (e.g., treatment).
#' @param NVar A character vector specifying the technical/nuisance variable names in the `colData` `DataFrame` (e.g., batch ID).
#' @param id_pb A character string specifying the variable name in the `colData` `DataFrame` used to pseudobulk the data, usually a combination of the sample ID and cell type.
#' @param id_sub A character string specifying the variable name in the `colData` `DataFrame` used as "subject-level" ID (e.g., patient ID).
#' If only one sample per subject is measured, then `id_sub` equals the sample ID.
#' @param cell_id A character string indicating the variable name in the `colData` `DataFrame` used as cell ID. Defaults to `"cell_id"`
#' @param n Number of pseudo-replicates to generate per group. Default is 1.
#' @param seed Integer used to set the random seed for reproducibility.
#'
#' @return A list of `SummarizedExperiment` objects.
#'
#' @importFrom dplyr group_by summarise mutate select left_join filter bind_rows n_distinct everything all_of across
#' @importFrom tibble tibble
#' @importFrom rlang sym syms :=
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom Matrix Matrix rowSums
#' @importFrom SummarizedExperiment SummarizedExperiment colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @importFrom BiocGenerics counts
#' @importFrom methods setGeneric setMethod
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
                                id_pb,
                                id_sub,
                                cell_id = 'cell_id',
                                n = 1,
                                seed = NULL){

            if (!inherits(counts, "SingleCellExperiment")) {
              stop("'counts' must be a SingleCellExperiment object.")
            }

            if(!ctype %in% BioVar){
              BioVar <- c(ctype,BioVar)

              message("Adding cell type to biological covariates")
            }

            metadata <- as.data.frame(SummarizedExperiment::colData(counts))

            required_cols <- c(id_pb, id_sub, cell_id, BioVar, NVar)
            missing_cols <- setdiff(required_cols, colnames(metadata))

            if (length(missing_cols) > 0) {
              stop(paste("The following columns are missing from metadata: ", paste(missing_cols, collapse = ", ")))
            }

            metadata <- metadata |>
              dplyr::mutate(dplyr::across(all_of(required_cols), as.character))

            ct <- unique(metadata[,ctype])

            PBC <- scuttle::aggregateAcrossCells(counts,metadata[,id_pb])
            PBC <- Matrix::Matrix(BiocGenerics::counts(PBC), sparse = TRUE)
            counts <- BiocGenerics::counts(counts)

            if(!all(metadata[,cell_id] %in% colnames(counts))){
              stop("The provided cell_id's do not match with the column names of the counts.")
            }

            unique_pb <- metadata |>
              dplyr::group_by(!!rlang::sym(id_pb)) |>
              dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::n_distinct), .groups = "drop") |>
              dplyr::select(-!!rlang::sym(id_pb)) |>
              colSums()

            metacovs <- c(id_pb, names(unique_pb)[unique_pb == dplyr::n_distinct(metadata[[id_pb]])])
            c_vars <- c(ctype,BioVar,NVar,id_pb,id_sub)
            problematic_vars <- c_vars[!c_vars %in% metacovs]

            if( length(problematic_vars)>0){
              msg <- paste("Consider providing unique sample identifiers. The following variables have multiple categories/levels associated with the same sample ID: ", paste(problematic_vars, collapse = ", "))
              stop(msg)
            }


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
                                 seed = seed)

            counts_pbc <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=full_pbc$counts),colData=full_pbc$metadata)

            counts_pbc_list <- lapply(ct,function(x) counts_pbc[,full_pbc$metadata[,ctype]==x])

            return(counts_pbc_list)
          }
)

