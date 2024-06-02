#' Spatial prioritization of the response to a biological perturbation
#'
#' Prioritize spatial barcodes involved in a complex biological process by training (1) a machine-learning model to predict sample labels (e.g., disease vs. control, treated vs. untreated, or time post-stimulus), and evaluate the performance of the model in cross-validation framework, and (2) a meta-learning model that trains on predictions from the prior machine learning models and distance metrics obtained from the spatial transcriptomics.
#'
#' @param input a matrix containing gene expression values
#'   (genes in rows, cells in columns); columns matching the barcodes in meta
#' @param meta a data frame containing metadata about the \code{input}
#'   gene-by-cell matrix, containing the minimum of the coordinates columns specified (default \code{c('x', 'y')}), \code{barcode_col} column and \code{label_col} column
#' @param label_col the column of the \code{meta} data frame that
#'   contains condition labels (e.g., disease, timepoint); defaults to \code{label}
#' @param barcode_col the column of the \code{meta} data frame that
#'   contains barcode; defaults to \code{barcode}
#' @param coord_cols the coordinate columns that are found in metadata; defaults to
#'  \code{('x', 'y')}
#' @param nn the nearest neighbors data frame derived from \link{get_nn}
#' @param distance_metrics_df distance metrics data frame derived from \link{get_distance_metrics}
#' @param dismay_methods distance metrics tures to calculate from \href{https://github.com/skinnider/dismay}{dismay} package
#' @param k (1) the nearest neighbors to calculate distance metrics for and (2) number of neighbours to use for Augur; defaults to \code{50}.
#' @param n_subsamples the number of random subsamples of fixed size to draw from the data; defaults to \code{10}.
#' @param max_barcodes hyperparamter for convergence; the maximum number of barcodes to run Augur on before stopping; defaults to \code{1000}.
#' @param barcodes_per_step hyperparamter for convergence; the number of barcodes to run Augur on per convergence step, analogous to learning rate; defaults to \code{10}.
#' @param epsilon hyperparamter for convergence; \eqn{\epsilon}, maximum \eqn{\Delta} between correlation between convergence steps; defaults to \code{5e-4}.
#' @param steps_tracking hyperparamter for convergence; number of convergence steps to track to assume convergence --  converged if number of steps with \eqn{\Delta} < \eqn{\epsilon}; defaults to \code{5e-4}.
#' @param cor_method correlation method to use for convergence; defaults to \code{pearson}
#' @param min_cor minimum correlation before assumption of convergence; defaults to \code{0.8}
#' @param seed seed for pseudorandom sampling; defaults to \code{42}
#' @param ncores number of cores to use, mainly for \code{Augur}; defaults to \code{8}
#' @param test_sim set to true only if testing Vespucci; load pre-loaded calculated_auc_df for spatial_sim data; defaults to false
#'
#' @return a list containing the following elements:
#' \itemize{
#'   \item \code{calculated_auc_df}: returns a data frame containing \code{barcode} and corresponding \code{auc} from \code{Augur::calculate_auc} that are used as training data
#'   \item \code{distance_metrics_df}: if argument \code{distance_metrics_df} is \code{NULL}, returns a data frame from \link{get_distance_metrics}; else if \code{distance_metrics_df} is not \code{NULL}, returns the same data frame as given in the argument
#'   \item \code{aucs}: a data frame with a \code{barcode} column and \code{auc} column after convergence
#'   \item \code{auc_tracking}: a data frame that tracks the aucs over the convergence step; contains a \code{barcode} column and their corresponding \code{auc} columns calculated for each convergence step
#'   \item \code{cor_tracking}: a data frame that tracks the correlation over convergence steps; contains \code{step}, \code{cor_method}: correlation method used and \code{cor}: correlation between step \code{i} and step \code{i-1}
#'
#'   \item \code{time_tracking}: a list that contains
#'     \itemize{
#'          \item \code{global}: time taken for \link{get_nn} and \link{get_distance_metrics}
#'          \item \code{aucs}: time taken for \link{calculate_auc} for each barcode
#'          \item \code{cor_convergence}: time taken for building model and predicting AUCs for each convergence step
#'      }
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
#' @importFrom randomForest randomForest
#' @import Matrix
#'
#' @export
calculate_spatial_auc = function(
	input,
	meta,
	label_col = 'label',
	barcode_col = 'barcode',
	coord_cols = c('x', 'y'),
	nn = NULL,
	distance_metrics_df = NULL,
	dismay_methods = c('pearson', 'spearman', 'binomial'),
	k = 50,
	n_subsamples = 10,
	max_barcodes = 1000,
	barcodes_per_step = 10,
	epsilon = 5e-4,
	steps_tracking = 5,
	cor_method = 'pearson',
	min_cor = 0.8,
	seed = 42,
	ncores = 8,
	test_sim = F,
	...
) {

	# set seed
	set.seed(seed)

	barcodes = colnames(input)

	# check if nn is set
	if (is.null(nn)) {
		nn_start_time = Sys.time()
		message('Getting nearest neighbours based on spatial coordinates')
		nn = get_nn(
			meta,
			label_col = label_col,
			barcode_col = barcode_col,
			coord_cols = coord_cols,
			...
		)
		nn_end_time = Sys.time()
		nn_time = difftime(nn_end_time, nn_start_time, units='mins')
		message(paste0('Nearest neighbours obtained. Time taken: ', nn_time, ' minutes'))
	} else {
		message('Nearest neighbours given')
		nn_time = -1
	}

	# get feature names
	features = unlist(map(dismay_methods, function(x) {
		paste0(x, '_', c('mean', 'median', 'q1', 'q3', 'sd'))
	}))

	# check if distance metrics are calculated
	if (is.null(distance_metrics_df)) {
		message('Calculating distance metrics')
		distance_metrics_start_time = Sys.time()
		distance_metrics_df = get_distance_metrics(
			input,
			meta,
			nn,
			k = k,
			label_col = label_col,
			barcode_col = barcode_col,
			dismay_methods = dismay_methods,
			select_var = F,
			ncores = ncores
		)
		distance_metrics_end_time = Sys.time()
		distance_metrics_time = difftime(distance_metrics_end_time, distance_metrics_start_time, units='mins')
		message(paste0('Distance metrics calculated. Time taken: ', distance_metrics_time, ' minutes'))
	} else {
		message('Distance metrics given')
		distance_metrics_time = -1
	}

	message('Running Augur::select_variance')
	input %<>% Augur::select_variance()

	curr_barcode_size = barcodes_per_step
	step_count = 1
	converged = FALSE
	random_sampled_barcodes = c()

	# make all data types are set correct
	distance_metrics_df %<>% type_convert()

	if (test_sim == TRUE) {
		message('Pre-calculated auc data frame is loaded. Only for spatial_sim dataset')
		data('spatial_sim_calculated_auc_df')
		barcodes_to_sample = spatial_sim_calculated_auc_df$barcode
		max_barcodes = nrow(spatial_sim_calculated_auc_df)
	} else {
		barcodes_to_sample = barcodes
	}

	calculated_auc_df = data.frame()

	while (curr_barcode_size < (max_barcodes + 1) & !converged) {
		
		curr_sampled_barcodes = sample(
			barcodes_to_sample[!barcodes_to_sample %in% random_sampled_barcodes], barcodes_per_step)
		random_sampled_barcodes = c(random_sampled_barcodes, curr_sampled_barcodes)

		if (test_sim == FALSE) {
			for (curr_barcode in curr_sampled_barcodes) {
				barcodes_nn = nn %>%
					dplyr::filter(source_barcode == curr_barcode) %>%
					group_by(target_label) %>%
					arrange(neighbor_idx) %>%
					mutate(rank = seq_len(n())) %>%
					filter(rank %in% seq_len(k))

				barcodes_keep = barcodes_nn$target_barcode
				target_groups_n = as.numeric(table(barcodes_nn$target_label))

				# check if min cell counts for each group is > 20
				check = length(target_groups_n) == 2 & all(target_groups_n > 20)

				if (check) {
					input0 = input %>% extract(,barcodes_keep)
					meta0 = meta %>% extract(barcodes_keep, ) %>% set_rownames(barcodes_keep)
					meta0$cell_type = 'vespucci'

					start_time = Sys.time()
					augur_res = Augur::calculate_auc(input = input0,
										  meta = meta0,
										  n_subsamples = n_subsamples,
										  var_quantile = 1,
										  n_threads = ncores,
										  show_progress = T)
					end_time = Sys.time()
					time_taken = difftime(end_time, start_time, units='secs')
					out_row = data.frame('barcode'=curr_barcode, 'auc' = as.numeric(augur_res$AUC$auc), 'time' = as.numeric(time_taken))
				} else {
					out_row = data.frame('barcode'=curr_barcode, 'auc' = NA, 'time' = NA)
				}
				calculated_auc_df %<>% rbind(out_row)
			}   
		} else {
			calculated_auc_df %<>% rbind(
				spatial_sim_calculated_auc_df %>% filter(barcode %in% curr_sampled_barcodes)
			)
		}

		training_data = distance_metrics_df %>%
			as.data.frame() %>%
			left_join(calculated_auc_df %>% dplyr::select(barcode, auc)) %>%
			filter(
				barcode %in% random_sampled_barcodes,
				!is.na(auc)
				) %>%
			rename(y=auc)
		training_data = training_data[,c(features,'y')]

		message(paste0('Running step ', step_count))
		# build rf
		rf_start_time = Sys.time()
		rf_fit = randomForest(y ~ ., data=training_data)
		rf_end_time = Sys.time()
		rf_time = difftime(rf_end_time, rf_start_time, units='mins')

		# predict time
		predict_start_time = Sys.time()
		full_data = as.data.frame(distance_metrics_df[,features])
		rownames(full_data) = distance_metrics_df$barcode
		auc_vals = as.numeric(predict(rf_fit, full_data))
		predict_end_time = Sys.time()
		predict_time = difftime(predict_end_time, predict_start_time, units='mins')

		# first run
		curr_step = paste0('step_', step_count)
		if (step_count == 1) {
			auc_tracking = data.frame(
				distance_metrics_df$barcode,
				auc_vals
			) %>% set_colnames(c('barcode', curr_step))
			cor_tracking = data.frame(
				'step' = curr_step,
				'model_time' = rf_time,
				'predict_time' = predict_time,
				'cor_method' = cor_method,
				'cor' = 0
			)
			prev_auc_vals = auc_vals
			prev_cor = NULL
			message(paste0('Step ', step_count, '; no correlation'))
		} else {
			auc_tracking %<>%
				left_join(
					data.frame(
						distance_metrics_df$barcode,
						auc_vals
					) %>% set_colnames(c('barcode', curr_step)),
					by = 'barcode'
				)
			cor = cor(auc_vals, prev_auc_vals, method = cor_method, use = "complete.obs")
			message(paste0('Step ', step_count, '; correlation: ', round(cor, 3)))
			cor_tracking %<>%
				rbind(
					data.frame(
						'step' = curr_step,
						'model_time' = rf_time,
						'predict_time' = predict_time,
						'cor_method' = cor_method,
						'cor' = cor
					)
				)
			if (is.null(prev_cor)) {
				prev_cor = cor
				steps_tracking_count = 0
			} else {
				cor_delta = cor - prev_cor
				if (cor_delta <= epsilon) {
					steps_tracking_count = steps_tracking_count + 1
				} else {
					steps_tracking_count = 0
				}
				if (steps_tracking_count > steps_tracking & cor >= min_cor) {
					converged = TRUE
				}
			}
			prev_auc_vals = auc_vals
			prev_cor = cor
		}
		curr_barcode_size = curr_barcode_size + barcodes_per_step
		step_count = step_count + 1
	}

	if (!converged) {
		message(paste0('Failed to converge over ', step_count, ' steps'))
	}

	time_tracking = list(
		'global' = data.frame(
			'step' = c('nn', 'distance_metrics'),
			'time' = c(nn_time, distance_metrics_time)
		),
		'auc' = calculated_auc_df %>% dplyr::select(barcode, time),
		'cor_convergence' = cor_tracking %>% dplyr::select(step, model_time, predict_time)
	)

	return (list(
		'calculated_auc_df' = calculated_auc_df,
		'distance_metrics_df' = distance_metrics_df,
		'aucs' = auc_tracking[,c('barcode', curr_step)] %>%
			set_colnames(c('barcode', 'auc')),
		'auc_tracking' = auc_tracking,
		'cor_tracking' = cor_tracking %>% dplyr::select(step, cor_method, cor),
		'time_tracking' = time_tracking
	))
}
