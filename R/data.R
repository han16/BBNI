#' Yeast Cell-Cycle Gene Expression Data
#'
#' A binary empirical dataset with the gene expression states of 14 yeast
#' cell-cycle genes across various experimental conditions. The dataset is used
#' in the original Han et al. (2014) paper to demonstrate the BBNI method for
#' inferring Boolean gene regulatory networks.
#'
#' @format A numeric matrix with 14 rows (representing individual genes/nodes)
#' and 385 columns (representing independent samples or sequential time points).
#'
#' @source Han, S., Wong, R. K. W., Lee, T. C. M., Shen, L., Li, S.-Y. R., & Fan, X. (2014).
#' A Full Bayesian Approach for Boolean Genetic Network Inference. *PLOS ONE*, 9(12), e115806.
#' \url{https://doi.org/10.1371/journal.pone.0115806}
#'
#' @usage data(yeast_data)
"yeast_data"
