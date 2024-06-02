#' Derive nearest neighbors within spatial transcriptomics
#'
#' Find nearest neighbors of each barcode in a spatial transcriptomic object
#'
#' @param meta a data frame containing metadata of the spatial transcriptomics object,
#'   containing the minimum of the coordinates columns specified (default
#'   \code{c('x', 'y')}), \code{barcode_col} column and \code{label_col} column
#' @param nn_radius determines the radius of the nearest neighbors to search for; defaults to \code{1000}.
#' @param coord_cols the coordinate columns that are found in metadata; defaults to
#'  \code{('x', 'y')}
#' @param label_col the column of the \code{meta} data frame that
#'   contains condition labels (e.g., disease, timepoint); defaults to \code{label}
#' @param barcode_col the column of the \code{meta} data frame that
#'   contains barcode; defaults to \code{barcode}
#'
#' @return a data frame containing the following columns:
#' \itemize{
#'   \item \code{source_barcode}: the source barcode
#'   \item \code{source_label}: the label that the source barcode belongs to
#'   \item \code{target_barcode}: the target barcode
#'   \item \code{target_group}: the label that the target barcode belongs to
#'   \item \code{rank}: the rank based on the nearest neighbors within each label,
#'      the lower the number, the closer the target barcode is to the source barcode
#'   \item \code{neighbor_idx}: the rank based on the nearest neighbors,
#'      the lower the number, the closer the target barcode is to the source barcode
#' }
#'
#' @importFrom dplyr mutate row_number select rename
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %<>%
#' @importFrom RANN nn2
#'
#' @export
get_nn = function(
    meta,
    nn_radius = 1000, # 50 * 20
    coord_cols = c('x', 'y'),
    label_col = 'label',
    barcode_col = 'barcode'
  ) {
    coords = meta %>%
        mutate(
            idx = row_number(),
            label = get(label_col),
            barcode = get(barcode_col)
        )

    # get the nearest neighbors based on spatial coordinates
    nn = coords %>%
        dplyr::select(all_of(coord_cols)) %>%
        nn2(k = nn_radius) %>%
        extract("nn.idx") %>%
        as.data.frame() %>%
        mutate(idx = row_number()) %>%
        dplyr::select(-1) %>%
        gather(neighbor, row_index, -idx) %>%
        left_join(coords, by = 'idx') %>%
        # add in the neighbors and re-code everything
        dplyr::select(-idx, -all_of(coord_cols)) %>%
        dplyr::rename(source_barcode = barcode, source_label = label, idx = row_index) %>%
        left_join(coords %>% dplyr::select(idx, barcode, label),
                  by = 'idx') %>%
        dplyr::rename(target_barcode = barcode, target_label = label) %>%
        dplyr::select(-idx) %>%
        group_by(source_barcode, target_label) %>%
        mutate(
            rank = seq_len(n()),
            neighbor_idx=as.numeric(gsub('nn\\.idx\\.', '', neighbor))
        ) %>%
        ungroup() %>%
        dplyr::select(source_barcode, source_label, target_barcode, target_label, rank, neighbor_idx)
    return(nn)
}

#' Distance metrics within spatial transcriptomics
#'
#' Get the distance metrics of each barcode to its surrounding nearest neighbors
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param meta a data frame containing metadata about the \code{input}
#'   gene-by-cell matrix, containing the minimum of the coordinates columns specified
#'   (default \code{c('x', 'y')}), \code{barcode_col} column and \code{label_col} column
#' @param nn the nearest neighbors data frame derived from \link{get_nn}
#' @param k the nearest neighbors to calculate distance metrics for; defaults to
#' \code{50}.
#' @param label_col the column of the \code{meta} data frame that
#'   contains condition labels (e.g., disease, timepoint); defaults to \code{label}
#' @param barcode_col the column of the \code{meta} data frame that
#'   contains barcode; defaults to \code{barcode}
#' @param dismay_methods distance metrics to calculate from \href{https://github.com/skinnider/dismay}{dismay} package
#' @param select_var boolean to run \code{Augur}::select_variance, defaults to \code{True}
#' @param ncores number of cores to use
#'
#' @return a data frame containing the following columns:
#' \itemize{
#'   \item \code{barcode}: the barcode of each cell
#'   \item \code{dismay summary statistics}: for each dismay method in \code{dismay_methods}, the following summary
#'                  statistics will be calculated
#'                  \itemize{
#'                      \item \code{mean}
#'                      \item \code{median}
#'                      \item \code{sd}: standard deviation
#'                      \item \code{q1}: 25th quantile
#'                      \item \code{q3}: 75th quantile
#'                  }
#' }
#'
#'
#' @importFrom dplyr mutate row_number select rename
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %<>%
#' @importFrom dismay dismay
#' @importFrom Augur select_variance
#' @importFrom purrr flatten
#' @importFrom furrr future_map_dfr
#' @importFrom future multicore
#'
#' @export
get_distance_metrics = function(
    input,
    meta,
    nn,
    k = 50,
    label_col = 'label',
    barcode_col = 'barcode',
    dismay_methods = c('pearson', 'spearman', 'binomial'),
    select_var = T,
    ncores = 8
) {
    output_df = data.frame()

    barcodes = colnames(input)

    summary_stats = c('mean', 'median', 'q1', 'q3', 'sd')

    features = unlist(map(summary_stats, function(x) {
        paste0(dismay_methods, '_', x)
    }))
    # do the variance selection on the entire matrix

    if (select_var){
        input %<>% Augur::select_variance()
    }

    plan(multicore,workers=ncores)
    output_df = future_map_dfr(barcodes, function(barcode){
        barcodes_keep = nn %>%
            dplyr::filter(source_barcode == barcode) %>%
            group_by(target_label) %>%
            arrange(neighbor_idx) %>%
            mutate(rank = seq_len(n())) %>%
            filter(rank %in% seq_len(k)) %>%
            pull(target_barcode)

        input0 = input %>% extract(, barcodes_keep)
        meta0 = meta %>% extract(barcodes_keep,) %>%
            mutate(
                label = paste0('label', as.integer(factor(label)))
            )

        # calculate distance metrics
        dist_metrics = future_map_dfr(dismay_methods, ~ {
            metric = .
            # message(metric)
            label1_idx = meta0$label == 'label1'
            label2_idx = meta0$label == 'label2'
            dist = dismay(as.matrix(input0), metric = metric)
            dist = dist[label2_idx, label1_idx]
            # define output vector
            out = data.frame(
                metric = metric,
                mean = mean(dist),
                median = median(dist),
                q1 = quantile(dist, 0.25),
                q3 = quantile(dist, 0.75),
                sd = sd(dist)
            )
            out
        }
        )
        # flatten dist metrics into vector
        cell_metrics = dist_metrics %>%
            select(-metric) %>%
            flatten() %>%
            unlist()
        cell_metrics[is.nan(cell_metrics)] = 0
        cell_metrics[is.infinite(cell_metrics)] = 0
        cell_metrics = data.frame(t(cell_metrics))
        colnames(cell_metrics) = features
        cell_metrics
    })
    output_df$barcode = barcodes
    return(output_df)
}

#' Get Gene Ontology (GO) level score matrix
#'
#' Get the GO score matrix from expression matrix
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param go_list A list of vectors of features for expression programs; each entry should be a vector of feature names; used as input into \code{Seurat::AddModuleScore}
#' @param ncores number of cores to use
#'
#' @return a matrix containing GO level scores
#'
#' @importFrom dplyr mutate row_number select rename
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %<>%
#' @importFrom Seurat AddModuleScore
#' @importFrom stats xtabs
#' @importFrom furrr future_map_dfr
#'
#' @export
get_GO_matrix = function(
    input,
    go_list,
    ncores = 8
) {
    full_df = data.frame()
    sc0 = CreateSeuratObject(input)
    sc0 %<>% NormalizeData()

    plan(multicore,workers=ncores)
    full_df = future_map_dfr(names(go_list), function(go){
        go_features = go_list[[go]]
        dat_check = F
        dat = tryCatch({
            sc1 = AddModuleScore(sc0, features = list(go=go_features))
            dat = data.frame(go = go, barcode = colnames(sc1), module_score = sc1$Cluster1)
            dat_check = T
            dat
        },
        error = function(e){
        })

        if (!dat_check) {
            dat = data.frame(go = go, barcode = colnames(sc0), module_score = NA)
        }
        dat
    })
    out_mat = xtabs(module_score ~ go + barcode, full_df)
    return(out_mat)
}

#' Count split matrix
#'
#' split count matrix with Poisson spliting
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns)
#' @param seed seed for pseudorandom sampling; defaults to \code{42}
#'
#' @return a list containing two matrices (one for running \code{calculate_spatial_auc}, the other for running \code{find_spatial_de})
#'
#' @importFrom dplyr mutate row_number select rename
#' @importFrom tidyr gather
#' @importFrom magrittr %>% %<>%
#' @importFrom countsplit countsplit
#'
#' @export
countsplit_matrix = function(
        input,
        seed = 42
) {
    set.seed(seed)
    split = countsplit(input, epsilon=0.5)
    auc = split$train
    auc@x = as.numeric(auc@x)
    de = split$test
    out_list = list(
        'run_auc' = auc,
        'run_de' = de
    )
    return(out_list)
}
