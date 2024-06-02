#' Spatial prioritization of the response to a biological perturbation and identification of spatially relevant features
#'
#' Vespucci first prioritize spatial barcodes involved in a complex biological process by training (1) a machine-learning model to predict sample labels (e.g., disease vs. control, treated vs. untreated, or time post-stimulus), and evaluate the performance of the model in cross-validation framework, and (2) a meta-learning model that trains on predictions from the prior machine learning models and distance metrics obtained from the spatial transcriptomics. Finally, Vespucci identifies spatially relevant features by fitting the expression and AUC values to a negative binomial generalized linear mixed model or a generalized linear mixed model. \cr Essentially \code{run_vespucci} runs the following functions in this particular order: \code{countsplit_matrix}, \code{calculate_spatial_auc} and \code{find_spatial_de}.
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param meta a data frame containing metadata about the \code{input}
#'   gene-by-cell matrix, containing the minimum of the coordinates columns specified (default \code{c('x', 'y')}), \code{barcode_col} column, \code{replicate_col} column and \code{label_col} column
#' @param seed seed for pseudorandom sampling; defaults to \code{42}
#'
#' @return a list containing the following elements:
#' \itemize{
#'   \item \code{spatial_auc_result}: list output from \code{calculate_spatial_auc}; most relevant element is \code{spatial_auc_results$aucs}, which is a barcode-to-auc data frame.
#'   \item \code{de_feature_result}: data frame output from \code{find_spatial_de}
#' }
#'
#' @importFrom dplyr do sample_n group_by ungroup tibble mutate select bind_rows pull rename n_distinct arrange desc filter summarise row_number left_join n all_of
#' @importFrom tibble repair_names rownames_to_column
#' @importFrom magrittr %>% %<>% extract extract2 set_rownames set_colnames
#' @importFrom stats setNames predict sd
#' @importFrom methods is
#' @importFrom sparseMatrixStats colVars
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel mclapply
#' @importFrom Augur select_variance calculate_auc
#' @importFrom RANN nn2
#' @importFrom tidyr gather
#' @importFrom nebula nebula
#' @importFrom lme4 lmerControl
#' @importFrom lmerTest lmer
#' @import Matrix
#'
#' @export
run_vespucci = function(
    input,
    meta,
    seed = 42,
    ...
) {

    # perform count splitting first
    message('Splitting matrix')
    split_mat_list = countsplit_matrix(input, seed=seed)

    spatial_auc_res = calculate_spatial_auc(
        input = split_mat_list$run_auc,
        meta = meta,
        seed = seed,
        ...
    )

    de_output_res = find_spatial_de(
        input = split_mat_list$run_de,
        meta = meta,
        auc_df = spatial_auc_res$aucs,
        ...
    )

    return (
        list(
            spatial_auc_result = spatial_auc_res,
            de_feature_result = de_output_res
        )
    )

}
