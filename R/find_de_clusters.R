#' Spatially relevant differential gene/GO term analysis
#'
#' Find clusters based on gene expression and DE genes
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param meta a data frame containing metadata about the \code{input}
#'   gene-by-cell matrix, containing the minimum of the \code{barcode_col} column and \code{replicate_col} column
#' @param ves_res results from run_vespucci
#' @param de_genes vector of DE genes chosen
#' @param barcode_col the column of the \code{meta} data frame that
#'   contains barcode; defaults to \code{barcode}
#' @param label_col the column of the \code{meta} data frame that
#'   contains label; defaults to \code{label}
#' @param norm if NormalizeData should be applied, defaults to true
#' @param label_to_run dataset to subset the label into; if null, use the full dataset
#' @param k k-means to run; if null, finds K with a naive heuristic of ELBO calculation from k-means
#'
#' @return a data frame containing the following columns:
#' \itemize{
#'   \item \code{gene}: gene
#'   \item \code{cluster}: which cluster it belongs to
#' }
#'
#' @importFrom dplyr do sample_n group_by ungroup tibble mutate select bind_rows pull rename n_distinct arrange desc filter summarise row_number left_join n all_of
#' @importFrom tibble repair_names rownames_to_column
#' @importFrom magrittr %>% %<>% extract extract2 set_rownames set_colnames
#' @importFrom stats setNames predict sd kmeans
#' @importFrom methods is
#' @importFrom tidyr gather
#' @import matrix
#'
#' @export
find_de_clusters = function(
    input,
    meta,
    ves_res=NULL,
    de_genes=NULL,
    barcode_col = 'barcode',
    label_col = 'label',
    norm = T,
    label_to_run = NULL,
    k = NULL,
    ...
) {
    if (is.null(ves_res) & is.null(de_genes)) {
         stop('Either use vespucci results as ves_res argument or use de_genes as de_genes arguments.')
    }
    if (norm == T) {
        input %<>% NormalizeData()
    }
    
    if (!is.null(ves_res)) {
        de_genes = ves_res$de_feature_result %>% filter(p_val_adj < 0.05) %>% pull(feature)
    }
    
    if (!is.null(label_to_run)) {
        barcodes = meta[,barcode_col][meta[,label_col] == label_to_run]
    } else{
        barcodes = meta[,barcode_col]
    }
    
    input0 = input[de_genes,barcodes]
    
    if (is.null(k)) {
        # automatically find clusters
        k_values = 3:10
        wss_values = sapply(k_values, function(k){kmeans(input0, centers = k, nstart = 10)$tot.withinss})
        # Compute second differences
        second_diff = diff(diff(wss_values))
        k = which.max(-second_diff) + 1 
        
    }
    clusts = kmeans(input0, centers = k)
    genes_cluster_df = clusts$cluster %>% data.frame() %>% set_colnames('cluster') %>% rownames_to_column('gene')
    return (genes_cluster_df)
}
