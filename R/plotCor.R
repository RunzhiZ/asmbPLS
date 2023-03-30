#' Graphical output for the asmbPLS-DA framework
#'
#' Function to visualize correlations between PLS components from different blocks 
#' using the model fitted by the function \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' 
#' @param fit.results The output of \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' @param ncomp Which component to plot from each block. Should not be larger 
#' than the number of PLS components used (\code{PLS.comp}) in the function
#' \code{\link[asmbPLS]{asmbPLSDA.fit}}. The default is 1.
#' @param block.name A vector containing the named character for each block. It
#' must be ordered and match each block. 
#' @param group.name A vector containing the named character for each sample 
#' group. For \code{binary} outcome, first group name matches \code{Y.matrix} = 0, 
#' second group name matches \code{Y.matrix} = 1. For \code{multiclass} outcome,
#'  \code{i}th group name matches \code{i}th column of \code{Y.matrix} = 1. 
#' @param legend A logical value indicating whether the legend should be added.
#' The default is TRUE.
#' 
#' @return 
#' none
#' 
#' @details 
#' The function returns a plot to show correlations between PLS components from 
#' different blocks. The lower triangular panel indicates Pearson's 
#' correlation coefficient, and the upper triangular panel the scatter plot.
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' Y.matrix.binary = asmbPLSDA.example$Y.matrix.binary
#' Y.matrix.multiclass = asmbPLSDA.example$Y.matrix.morethan2levels
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb = asmbPLSDA.example$quantile.comb
#'  
#' ## asmbPLSDA fit for binary outcome
#' asmbPLSDA.fit.binary <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                       Y.matrix = Y.matrix.binary, 
#'                                       PLS.comp = PLS.comp, 
#'                                       X.dim = X.dim, 
#'                                       quantile.comb = quantile.comb,
#'                                       outcome.type = "binary")
#' 
#' ## asmbPLSDA fit for categorical outcome with more than 2 levels
#' asmbPLSDA.fit.multiclass <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                           Y.matrix = Y.matrix.multiclass,
#'                                           PLS.comp = PLS.comp, 
#'                                           X.dim = X.dim, 
#'                                           quantile.comb = quantile.comb,
#'                                           outcome.type = "multiclass")
#' 
#' ## visualization with default block.name and group.name using the first PLS component 
#' plotCor(asmbPLSDA.fit.binary, 1)
#' plotCor(asmbPLSDA.fit.multiclass, 1)
#' ## custom block.name and group.name
#' plotCor(asmbPLSDA.fit.binary, 
#'         ncomp = 1, 
#'         block.name = c("mRNA", "protein"), 
#'         group.name = c("control", "case"))
#' plotCor(asmbPLSDA.fit.multiclass, 
#'         ncomp = 1, 
#'         block.name = c("mRNA", "protein"), 
#'         group.name = c("healthy", "mild", "severe"))
#' 
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @import ggplot2 ggpubr
#' @importFrom stats cor

plotCor <- function(fit.results, ncomp = 1, block.name = NULL, group.name = NULL, legend = T) {
  if(is.null(block.name)) {
    n.block <- length(fit.results$X_dim)
    block.name <- paste0("block.", 1:n.block)
  }
  ## extract the Y information
  group_info <- fit.results$Y_group 
  
  if(is.null(group.name)) {
    if(ncol(group_info) == 1) {n.group = 2} else {n.group = ncol(group_info)}
    group.name <- paste0("group.", 1:n.group)
  }
  
  Group <- matrix(NA, nrow = nrow(group_info), ncol = 1)
  ## binary
  if(ncol(group_info) == 1) { 
    Group[which(group_info == 0)] <- group.name[1]
    Group[which(group_info == 1)] <- group.name[2]
  } else { 
    for(i in 1:ncol(group_info)) {
      Group[which(group_info[, i] == 1)] <- group.name[i]
    }
  }
  n_length <- length(block.name)
  
  ## create the layout matrix
  layout_matrix <- matrix(c(1:n_length^2), nrow = n_length, byrow = T)
  layout_diag <- diag(layout_matrix)
  layout_diag_index <- 1:n_length
  layout_lower <- sort(layout_matrix[lower.tri(layout_matrix)])
  layout_upper <- sort(layout_matrix[upper.tri(layout_matrix)])
  data_lower <- as.data.frame(matrix(NA, nrow = length(layout_lower), ncol = 3))
  data_upper <- as.data.frame(matrix(NA, nrow = length(layout_upper), ncol = 3))
  colnames(data_lower) <- colnames(data_upper) <- c("index", "X", "Y")
  data_lower$index <- layout_lower
  n = 1
  n_lower = (n_length*n_length - n_length)/2
  X <- 2
  Y <- 1
  data_lower[n, 2:3] <- c(X, Y)
  while(n < n_lower) {
    n = n + 1
    if(X <= Y + 1) {
      X = X + 1
      Y = 1
    } else {
      Y = Y + 1
    }
    data_lower[n, 2:3] <- c(X, Y)
  }
  data_upper$index <- layout_upper
  data_upper$X <- data_lower$X
  data_upper$Y <- data_lower$Y
  
  ## check if the corresponding index in diag, upper or lower part
  for(i in 1:n_length^2) {
    ## diag part: block.name
    if(i %in% layout_diag) {
      eval(parse(text = paste0("p", i, " <- ggplot(data.frame(x = 1, y = 1, text = \"", block.name[which(layout_diag == i)], "\"), aes(x, y)) + geom_text(aes(label = text), size = 10) + theme_bw() + 
  theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())")))
    }
    ## upper part: pair-wise point plot + ellipse
    if(i %in% layout_upper) {
      X_temp <- data_upper$X[which(data_upper$index == i)]
      Y_temp <- data_upper$Y[which(data_upper$index == i)]
      df <- as.data.frame(cbind(fit.results$X_score[[X_temp]][, ncomp], fit.results$X_score[[Y_temp]][, ncomp], Group))
      colnames(df) <- c("X", "Y", "Group")
      df$X <- as.numeric(df$X)
      df$Y <- as.numeric(df$Y)
      df$Group <- as.factor(df$Group)
      eval(parse(text = paste0("p", i, " <- ggplot(data = df, aes(X, Y, color = Group)) + geom_point(show.legend = ", legend, ") + stat_ellipse(show.legend = FALSE) + xlab(\"\") + ylab(\"\") + theme_bw()")))
    }
    ## lower part: pair-wise correlation
    if(i %in% layout_lower) {
      X_temp <- data_lower$X[which(data_lower$index == i)]
      Y_temp <- data_lower$Y[which(data_lower$index == i)]
      cor_value <- cor(fit.results$X_score[[X_temp]][, ncomp], fit.results$X_score[[Y_temp]][, ncomp])
      eval(parse(text = paste0("p", i, " <- ggplot(data.frame(x = 1, y = 1, text = round(cor_value,2)), aes(x, y)) + geom_text(aes(label = text), size = ", 25*sqrt(cor_value), ") + theme_bw() + 
  theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())")))
    }
  }
  p_output <- NULL
  eval(parse(text = paste0("p_output <- ggarrange(", paste(paste("p", 1:n_length^2, sep = ""), collapse = ", "), ", ncol = ", n_length, ", nrow = ", n_length, ", align = \"hv\", 
          widths = c(1, 1), heights = c(1, 1),
          common.legend = TRUE)")))
  return(p_output)
}
