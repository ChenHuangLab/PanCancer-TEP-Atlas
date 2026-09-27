#!/usr/bin/env Rscript
# Platelet bulk RNA-seq benchmark of five MHC-I genes against a geometric-mean score.
# Adapted from Platelet_bulk_diagnosis_AUC.R, 101_machine_learning_main.R,
# 9_ML_method.R and 21_ML_COM_method_new_version.R. See README.md for input data.
# Run with: Rscript Platelet_bulk_diagnosis_ML_benchmark.R

features <- c("HLA_A", "HLA_B", "HLA_C", "HLA_E", "HLA_F")
input_suffix <- "_5_mhc_i_expression_for_training_df.rds"
train_id <- "GSE183635"
seed <- 1314L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript so relative paths resolve from its directory.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L != 0L ||
    (length(args) > 0L && any(!args[seq(1L, length(args), by = 2L)] %in%
                              c("--data-dir", "--output-dir")))) {
  stop("Usage: Rscript Platelet_bulk_diagnosis_ML_benchmark.R [--data-dir PATH] [--output-dir PATH]")
}
options <- list("--data-dir" = file.path(script_dir, "data"),
                "--output-dir" = file.path(script_dir, "results"))
if (length(args)) {
  for (i in seq(1L, length(args), by = 2L)) options[[args[i]]] <- args[i + 1L]
}
data_dir <- options[["--data-dir"]]
output_dir <- options[["--output-dir"]]

packages <- c("pROC", "glmnet", "randomForest", "e1071", "xgboost", "gbm",
              "lightgbm", "superpc", "pls", "ggplot2", "writexl")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the missing R packages before running: ",
       paste(missing_packages, collapse = ", "))
}

if (!dir.exists(data_dir)) stop("Input directory does not exist: ", data_dir)
files <- list.files(data_dir, full.names = TRUE)
files <- files[endsWith(basename(files), input_suffix)]
if (!length(files)) stop("No matching RDS inputs in: ", data_dir)
ids <- substr(basename(files), 1L, nchar(basename(files)) - nchar(input_suffix))
if (anyDuplicated(ids)) stop("Duplicate dataset IDs in input filenames.")
names(files) <- ids
if (!train_id %in% ids) stop("Training dataset ", train_id, " is missing from ", data_dir)
validation_ids <- setdiff(ids, train_id)
if (!length(validation_ids)) stop("At least one other dataset is needed for the benchmark.")

prepare_dataset <- function(file, id) {
  raw <- readRDS(file)
  if (!is.data.frame(raw)) stop(id, ": RDS must contain a data.frame.")
  names(raw) <- gsub("[.-]", "_", names(raw))

  # The source analysis corrected a duplicated HLA-E in this particular RDS.
  # Retain that correction only when HLA-B is absent; verify it against source data.
  if (identical(id, "GSE232027") && !"HLA_B" %in% names(raw) &&
      sum(names(raw) == "HLA_E") == 2L) {
    names(raw)[which(names(raw) == "HLA_E")[2L]] <- "HLA_B"
    warning(id, ": second HLA_E column treated as HLA_B, following the source script; ",
            "verify the gene annotation in the input data.")
  }
  if (anyDuplicated(names(raw))) stop(id, ": duplicated column names; correct the RDS first.")
  required <- c(features, "Disease")
  if (!all(required %in% names(raw))) {
    stop(id, ": missing columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
  }

  x <- raw[, features, drop = FALSE]
  for (gene in features) {
    if (!is.numeric(x[[gene]])) stop(id, ": ", gene, " must be numeric.")
  }
  if (any(!is.finite(as.matrix(x))) || any(as.matrix(x) < 0)) {
    stop(id, ": gene expression must be finite and nonnegative.")
  }
  disease <- trimws(as.character(raw$Disease))
  label <- rep(NA_integer_, length(disease))
  label[grepl("Control|healthy|nonMalignant", disease, ignore.case = TRUE)] <- 0L
  if (id == "GSE232027") {
    label[tolower(disease) == "low"] <- 0L
    label[tolower(disease) == "high"] <- 1L
  } else if (id == "GSE126446") {
    label[tolower(disease) == "mature"] <- 0L
    label[tolower(disease) == "immature"] <- 1L
  } else {
    label[is.na(label) & !is.na(disease) & nzchar(disease)] <- 1L
  }
  if (anyNA(label)) stop(id, ": missing or unmapped Disease labels.")
  if (!identical(sort(unique(label)), c(0L, 1L))) {
    stop(id, ": both control (0) and case (1) samples are required.")
  }
  data.frame(x, class = factor(label, levels = c(0L, 1L)), check.names = FALSE)
}

datasets <- Map(prepare_dataset, files, names(files))
names(datasets) <- names(files)
train <- datasets[[train_id]]
validation <- datasets[validation_ids]
if (min(table(train$class)) < 3L) {
  stop("Training data need at least three samples per class for cross-validation.")
}
message("Training: ", train_id, " (", nrow(train), " samples); external datasets: ",
        paste(validation_ids, collapse = ", "))

# Each sample's half-minimum positive count is the pseudocount used in the
# source analysis. The score is 0 when all five expression values are zero.
geometric_mean <- function(values) {
  positive <- values[values > 0]
  if (!length(positive)) return(0)
  exp(mean(log(values + min(positive) / 2)))
}

auc_for_scores <- function(data, scores) {
  scores <- as.numeric(scores)
  if (length(scores) != nrow(data) || any(!is.finite(scores))) {
    stop("Prediction must contain one finite score per sample.")
  }
  as.numeric(pROC::auc(pROC::roc(response = data$class, predictor = scores,
                                levels = c("0", "1"), direction = "<", quiet = TRUE)))
}

results <- list()
add_result <- function(model_name, predict_one) {
  message("Evaluating ", model_name)
  values <- vapply(validation, function(data) auc_for_scores(data, predict_one(data)), numeric(1))
  results[[model_name]] <<- values
}

# Geometric mean is a fixed score and is never fitted on the training cohort.
add_result("Geometric_Mean", function(data) {
  apply(as.matrix(data[, features, drop = FALSE]), 1L, geometric_mean)
})

set.seed(seed)
x_train <- as.matrix(train[, features, drop = FALSE])
y_train <- as.integer(as.character(train$class))
foldid <- integer(nrow(train))
for (cls in c(0L, 1L)) {
  positions <- which(y_train == cls)
  foldid[positions] <- sample(rep(seq_len(10L), length.out = length(positions)))
}
if (any(tabulate(foldid, nbins = 10L) == 0L)) {
  stop("Training cohort is too small for 10-fold glmnet cross-validation.")
}

# Ridge, nine elastic-net mixing parameters, and LASSO: 11 models.
for (alpha in seq(0, 1, by = 0.1)) {
  set.seed(seed)
  fit <- glmnet::cv.glmnet(x_train, y_train, family = "binomial", alpha = alpha,
                           foldid = foldid, nfolds = 10L)
  model_name <- if (alpha == 0) "Ridge" else if (alpha == 1) "LASSO" else
    paste0("Enet", format(alpha, nsmall = 1L))
  add_result(model_name, function(data) {
    as.numeric(predict(fit, newx = as.matrix(data[, features, drop = FALSE]),
                       s = "lambda.min", type = "response"))
  })
}

set.seed(seed)
forest <- randomForest::randomForest(class ~ ., data = train, ntree = 3000L)
best_trees <- which.min(forest$err.rate[, "OOB"])
set.seed(seed)
forest <- randomForest::randomForest(class ~ ., data = train, ntree = best_trees)
add_result("RF", function(data) {
  as.numeric(predict(forest, newdata = data[, features, drop = FALSE],
                     type = "prob")[, "1"])
})

set.seed(seed)
svm_tune <- e1071::tune.svm(class ~ ., data = train, kernel = "linear",
                            cost = c(0.001, 0.01, 0.1, 1, 5, 10), probability = TRUE)
svm_fit <- svm_tune$best.model
add_result("SVM", function(data) {
  pred <- predict(svm_fit, newdata = data[, features, drop = FALSE], probability = TRUE)
  as.numeric(attr(pred, "probabilities")[, "1"])
})

set.seed(seed)
xgb_fit <- xgboost::xgboost(data = x_train, label = y_train, nrounds = 150L,
                            objective = "binary:logistic", eta = 0.05,
                            max_depth = 4L, gamma = 0.1, colsample_bytree = 0.8,
                            min_child_weight = 1, subsample = 0.8,
                            eval_metric = "auc", nthread = 1L, verbose = 0)
add_result("XGB", function(data) {
  predict(xgb_fit, newdata = as.matrix(data[, features, drop = FALSE]))
})

gbm_train <- train
gbm_train$class <- y_train
set.seed(seed)
gbm_fit <- gbm::gbm(class ~ ., data = gbm_train, distribution = "bernoulli",
                    n.trees = 1000L, shrinkage = 0.1, interaction.depth = 3L,
                    n.minobsinnode = 10L, cv.folds = 10L, n.cores = 1L,
                    train.fraction = 0.7, verbose = FALSE)
gbm_trees <- gbm::gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
add_result("GBM", function(data) {
  as.numeric(predict(gbm_fit, newdata = data[, features, drop = FALSE],
                     n.trees = gbm_trees, type = "response"))
})

set.seed(seed)
lgb_fit <- lightgbm::lgb.train(
  params = list(objective = "binary", metric = "binary_error", num_leaves = 31L,
                learning_rate = 0.05, max_depth = -1L, boosting_type = "gbdt",
                num_threads = 1L),
  data = lightgbm::lgb.Dataset(data = x_train, label = y_train), nrounds = 100L)
add_result("lgbm", function(data) {
  predict(lgb_fit, as.matrix(data[, features, drop = FALSE]))
})

full_glm <- stats::glm(class ~ ., family = stats::binomial(), data = train)
null_glm <- stats::glm(class ~ 1, family = stats::binomial(), data = train)
step_fits <- list(
  step_both = stats::step(full_glm, direction = "both", trace = 0),
  step_forward = stats::step(null_glm, scope = list(lower = ~1, upper = class ~ .),
                             direction = "forward", trace = 0),
  step_backward = stats::step(full_glm, direction = "backward", trace = 0)
)
for (model_name in names(step_fits)) {
  fit <- step_fits[[model_name]]
  add_result(model_name, function(data) {
    as.numeric(predict(fit, newdata = data[, features, drop = FALSE], type = "response"))
  })
}

super_train <- list(x = t(x_train), y = y_train, featurenames = features)
set.seed(seed)
super_fit <- superpc::superpc.train(super_train, type = "regression")
super_cv <- superpc::superpc.cv(super_fit, super_train, n.threshold = length(features),
                                n.fold = 10L, n.components = 3L, min.features = 1L,
                                max.features = length(features), compute.fullcv = TRUE,
                                compute.preval = TRUE)
scores <- as.numeric(super_cv$scor[1L, ])
if (!length(scores) || all(!is.finite(scores))) stop("superpc CV returned no valid threshold.")
scores[!is.finite(scores)] <- -Inf
super_threshold <- super_cv$thresholds[which.max(scores)]
add_result("super", function(data) {
  test <- list(x = t(as.matrix(data[, features, drop = FALSE])),
               y = as.integer(as.character(data$class)), featurenames = features)
  as.numeric(superpc::superpc.predict(super_fit, super_train, test,
                                      threshold = super_threshold, n.components = 1L)$v.pred)
})

pls_train <- train
pls_train$class <- y_train
set.seed(seed)
pls_fit <- pls::plsr(class ~ ., data = pls_train, validation = "CV")
pls_components <- pls::selectNcomp(pls_fit, method = "randomization", plot = FALSE)
pls_components <- max(1L, min(length(features), as.integer(pls_components)))
add_result("plsr", function(data) {
  as.numeric(predict(pls_fit, newdata = data[, features, drop = FALSE],
                     ncomp = pls_components, type = "response"))
})

auc_matrix <- do.call(rbind, results)
auc_matrix <- auc_matrix[order(rowMeans(auc_matrix), decreasing = TRUE), , drop = FALSE]
auc_table <- data.frame(Model_Combination = rownames(auc_matrix),
                        as.data.frame(auc_matrix, check.names = FALSE),
                        Mean_External_AUC = rowMeans(auc_matrix),
                        check.names = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(auc_table, file.path(output_dir, "all_model_auc_results.csv"), row.names = FALSE)
writexl::write_xlsx(auc_table, file.path(output_dir, "all_model_auc_results.xlsx"))

heatmap_data <- data.frame(Model = rep(rownames(auc_matrix), times = ncol(auc_matrix)),
                           Dataset = rep(colnames(auc_matrix), each = nrow(auc_matrix)),
                           AUC = as.vector(auc_matrix))
heatmap_data$Model <- factor(heatmap_data$Model, levels = rev(rownames(auc_matrix)))
plot <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Dataset, y = Model, fill = AUC)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", AUC)), size = 2.5) +
  ggplot2::scale_fill_gradient2(low = "#0084A7", mid = "#F5FACD", high = "#E05D00",
                                midpoint = 0.75, limits = c(0, 1)) +
  ggplot2::labs(title = "Five MHC-I genes: external-cohort AUC benchmark",
                x = "External dataset", y = "Method") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(output_dir, "AUC_Heatmap.pdf"), plot = plot,
                width = 11, height = 10, units = "in")
message("Wrote results to: ", normalizePath(output_dir))
