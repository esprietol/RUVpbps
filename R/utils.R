

#' select the cells to create pseudobulk pseudosamples.
#'
#' @description `cellspbps()` samples with replacement sets of cells from groups that share the same biological and technical conditions.
#' Preserving the average number of cells per biological subgroup.
#'
#' @param subsample A `data.frame` or `tibble` containing single-cell data with columns for group, bsg (biological subgroup), sample identifiers and cell identifiers.
#' @param pb_id A character string specifying the column name for biological replicate/sample identifier. Default is `"sample_cell"`.
#' @param cell_id A character string specifying the column name for cell identifiers. Default is `"cell_id"`.
#' @param seed A numeric value to set the random seed for reproducibility. Default is 1.
#' @param n An integer specifying how many pseudobulk pseudosamples to generate from each group. Default is 1.
#'
#' @return A list of character vectors containing sampled cell IDs for each pseudobulk pseudosample.
#'
#' @importFrom dplyr group_by summarise select left_join n
#' @importFrom tidyr all_of
#' @keywords internal

cellspbps <- function(subsample, pb_id = 'sample_cell', cell_id = 'cell_id', seed = 1, n = 1) {

  # compute the average number of cells per biological group
  total_cell <- subsample |>
    dplyr::group_by(dplyr::across(dplyr::all_of(pb_id))) |>
    dplyr::summarise(ncell = dplyr::n(), .groups = "drop")

  subsample <- dplyr::left_join(subsample, total_cell, by = pb_id)

  avg_cell <- subsample |>
    dplyr::group_by(bsg) |>
    dplyr::summarise(avg_ncell = mean(ncell), .groups = "drop")

  avg_cell <- subsample |>
    dplyr::select(group, bsg) |>
    dplyr::distinct() |>
    dplyr::left_join(avg_cell, by = "bsg")

  rownames(avg_cell) <- avg_cell$group

  # select the cells from each group
  groups <- split(subsample[[cell_id]], subsample$group)
  avg_cell <- avg_cell[names(groups), ]

  set.seed(seed)

  # sampling with replacement
  pseudo_samp <- replicate(n, list(
    mapply(function(x, y) base::sample(x = x, size = round(y), replace = TRUE),
           groups, avg_cell$avg_ncell, SIMPLIFY = FALSE)
  ))

  names(pseudo_samp) <- paste0("pa", 1:n)
  pseudo_samp <- unlist(pseudo_samp, recursive = FALSE)

  return(pseudo_samp)
}
