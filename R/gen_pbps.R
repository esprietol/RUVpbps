
#' Generate a pseudobulk dataset with pseudobulk pseudosamples to be used as negative control samples in the removal of unwanted variation
#'
#' @description This function orchestrates the process of creating pseudobulk pseudosamples (PBPS) from a Seurat object and associated metadata.
#' It combines sampling, aggregation, and metadata reconstruction, and returns an edgeR-compatible `DGEList` object.
#'
#' @param ds A Seurat object containing the single-cell data.
#' @param ctype A character string specifying the column in `ds@meta.data` used to define the cell type.
#' @param BioVar A character vector of column names in `ds@meta.data` representing biological variables (e.g., treatment).
#' @param NVar A character vector of column names in `ds@meta.data` representing technical/nuisance variables (e.g., batch ID).
#' @param id_pb Name of the column to use as sample ID (e.g., "sample_id").
#' @param id_sub Name of the column to store the "subject-level" ID (e.g., patient).
#' @param cell_id A character string indicating the column in `ds@meta.data` with the cell ID.
#' @param n Number of pseudo-replicates to generate per group. Default is 2.
#' @param Seed Integer used to set the random seed for reproducibility.
#'
#' @return A `DGEList` object (from edgeR) containing counts and metadata.
#'
#' @importFrom dplyr group_by summarise mutate select left_join filter bind_rows n_distinct everything all_of across
#' @importFrom tibble tibble
#' @importFrom Seurat ScaleData AggregateExpression
#' @importFrom edgeR DGEList
#' @importFrom rlang sym syms :=
#' @export

gen_PBPS <- function(ds, ctype='cg_cov', BioVar, NVar, id_pb = 'sample_cell', id_sub = 'ind_cov', cell_id = 'cell_id', n = 2, Seed = 2) {

  ds <- Seurat::ScaleData(ds)

  # Aggregate to pseudobulk level
  PBC <- Seurat::AggregateExpression(ds, group.by = id_pb, return.seurat = FALSE)

  # Identify sample-level covariates
  unique_pb <- ds@meta.data |>
    dplyr::group_by(!!rlang::sym(id_pb)) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::n_distinct), .groups = "drop") |>
    dplyr::select(-!!rlang::sym(id_pb)) |>
    colSums()

  metacovs <- c(id_pb, names(unique_pb)[unique_pb == dplyr::n_distinct(ds@meta.data[[id_pb]])])

  # Create DGEList
  PBC <- edgeR::DGEList(counts = PBC[[1]], samples = unique(ds@meta.data[, metacovs]))

  # Grouping metadata
  groups_info <- PBC$samples |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(BioVar, NVar)))) |>
    dplyr::summarise(nsample = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      bsg = paste(!!!rlang::syms(BioVar), sep = "_"),
      group = paste(bsg, !!!rlang::syms(NVar), sep = "_")
    )

  bsg.keep <- names(table(groups_info$bsg)[table(groups_info$bsg) > 1])

  PBC$samples <- PBC$samples |>
    dplyr::mutate(
      bsg = paste(!!!rlang::syms(BioVar), sep = "_"),
      group = paste(bsg, !!!rlang::syms(NVar), sep = "_")
    )

  sc.ref <- ds@meta.data |>
    dplyr::mutate(
      bsg = paste(!!!rlang::syms(BioVar), sep = "_"),
      group = paste(bsg, !!!rlang::syms(NVar), sep = "_")
    ) |>
    dplyr::filter(bsg %in% bsg.keep)

  # Cell sampling
  cells_to_pbps <- cellspbps(subsample = sc.ref, seed = Seed, n = n)

  # Count aggregation
  pbpscounts <- sapply(names(cells_to_pbps), function(x) pbps(cells_to_pbps, Gr = x, ds = ds, cell_id = cell_id))

  # Metadata reconstruction
  pbps_info <- tibble::tibble(!!rlang::sym(id_pb) := colnames(pbpscounts)) |>
    dplyr::mutate(
      pseudosample = sub("\\..*$", "", .data[[id_pb]]),
      group = sub("^.*?\\.", "", .data[[id_pb]])
    ) |>
    dplyr::left_join(groups_info, by = "group") |>
    dplyr::mutate(!!rlang::sym(id_sub) := paste0("pbps_", bsg)) |>
    dplyr::select(-nsample)

  miss_var <- setdiff(colnames(PBC$samples), colnames(pbps_info))
  for (v in miss_var) pbps_info[[v]] <- NA

  pbps_info <- dplyr::mutate(pbps_info, pbps = 1)

  sifull <- PBC$samples |>
    dplyr::mutate(pbps = 0, pseudosample = "orig") |>
    dplyr::bind_rows(pbps_info)

  fullcount <- cbind(PBC$counts, pbpscounts[rownames(PBC$counts), , drop = FALSE])

  PBPSC <- edgeR::DGEList(
    counts = fullcount,
    samples = dplyr::select(sifull, -lib.size, -norm.factors)
  )

  return(PBPSC)
}

