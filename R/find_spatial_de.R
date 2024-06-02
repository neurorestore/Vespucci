#' Spatially relevant differential gene/GO term analysis
#'
#' Find differentially expressed genes/GO terms that corresponds to AUC derived from \code{navigate_space} function. This is calculated by fitting either a negative binomial generalized mixed model or a generalized linear mixed model  of the expression values to the AUC values while accounting for biological replicates.
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param meta a data frame containing metadata about the \code{input}
#'   gene-by-cell matrix, containing the minimum of the \code{barcode_col} column and \code{replicate_col} column
#' @param auc_df a data frame containing \code{barcode} (matching meta \code{barcode_col} column) and \code{auc} (AUC values) columns
#' @param de_method method to run differential expression (defaults to \code{nbgmm}).
#'   \describe{
#'     \item{\code{nbgmm}}{Negative binomial generalized mixed model}
#'     \item{\code{glm}}{Generalized linear model}
#'   }
#' @param barcode_col the column of the \code{meta} data frame that
#'   contains barcode; defaults to \code{barcode}
#' @param replicate_col the column of the \code{meta} data frame that
#'   contains the biological replicate id; defaults to \code{replicate}
#' @param go_list A list of vectors of features for expression programs; each entry should be a vector of feature names; used as input into \code{Seurat::AddModuleScore}
#'
#' @return a data frame containing the following columns:
#' \itemize{
#'   \item \code{feature}: gene/GO term
#'   \item \code{effect_size}: effect size based on AUC values
#'   \item \code{p_val}: p-value
#'   \item \code{p_val_adj}: BH adjusted p-value
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
#' @importFrom nebula nebula group_cell
#' @importFrom RANN nn2
#' @importFrom lme4 lmerControl
#' @importFrom lmerTest lmer
#' @importFrom tidyr gather
#' @import Matrix
#'
#' @export
find_spatial_de = function(
    input,
    meta,
    auc_df,
    de_method = 'nbgmm',
    barcode_col = 'barcode',
    replicate_col = 'replicate',
    go_list = NULL,
    ...
) {
    barcodes = colnames(input)
    if (!is.null(go_list)) {
        go_start_time = Sys.time()
        message('Gene Ontology specified. Running AddModuleScore')

        input = get_GO_matrix(input, go_list)

        go_end_time = Sys.time()
        go_time = difftime(go_end_time, go_start_time, units='mins')
        message(paste0('Gene Ontology matrix obtained. Time taken: ', go_time, ' minutes'))

        de_method ='glm'
        message('Defaulting to glm since GO matrix is non-integer.')
    } else {
        go_time = -1
    }

    meta %<>%
        left_join(auc_df, by=barcode_col) %>%
        mutate(label = auc) %>%
        filter(!is.na(label))
    meta$replicate = meta[,replicate_col]

    if (de_method == 'nbgmm') {
        df = model.matrix(~label, data=meta)
        data_g = group_cell(count=input,id=meta$replicate,pred=df)
        if (!is.null(data_g)){
            expr0 = data_g$count
            id = data_g$id
            pred = data_g$pred
        } else {
            id = meta0$replicate
            pred = df
        }
        res = nebula(
            input,
            id,
            pred=pred,
            model = 'NBGMM',
            cpc = 0
        )
        de_output = res$summary %>% as_tibble() %>%
            dplyr::rename(
                feature = gene,
                p_val = p_label,
                effect_size = logFC_label
            ) %>%
            select(feature, effect_size, p_val)
        de_output$p_val_adj = p.adjust(de_output$p_val, method = 'BH')
        de_output %<>%
            arrange(p_val_adj, p_val)
    } else if (de_method == 'glm') {
        ctrl = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
        genes = rownames(input)
        de_output = data.frame(do.call(rbind, lapply(genes, function(x){
            message(x)
            p_val = NA
            effect_size = NA
            tryCatch({
                dat0 = data.frame(
                    GENE = input[x,],
                    label = meta$label,
                    replicate = factor(meta$replicate)
                )
                fit = lmerTest::lmer(GENE ~ label + (1 | replicate), dat0, REML = TRUE, control = ctrl)
                coefs = coef(summary(fit))
                p_val = coefs[2, 'Pr(>|t|)']
                effect_size = coefs[2, 'Estimate']
            },
            error = function(e){message(e)})
            return(
                data.frame(
                    'feature' = x,
                    'effect_size' = effect_size,
                    'p_val' = p_val
                )
            )
        })))
        de_output$p_val_adj = p.adjust(de_output$p_val, method = 'BH')
        de_output %<>%
            arrange(p_val_adj, p_val)
    }
    return (de_output)
}
