#' PLS plot for asmbPLS-DA
#'
#' Function to visualize cluster of samples using super score of different PLS
#' components.
#' 
#' @param fit.results The output of \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' @param comp.X A integer indicating which PLS component to be used for the 
#' X.axis. The default is 1.
#' @param comp.Y A integer indicating which PLS component to be used for the 
#' Y.axis. The default is 2.
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
#' The function returns a plot to show cluster of samples using super score of 
#' different PLS components. 
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
#' ## visualization to show the cluster of samples using the first and the second super score
#' plotPLS(asmbPLSDA.fit.binary, comp.X = 1, comp.Y = 2)
#' plotPLS(asmbPLSDA.fit.multiclass, comp.X = 1, comp.Y = 2)
#' ## custom group.name
#' plotPLS(asmbPLSDA.fit.binary, 
#'         comp.X = 1, 
#'         comp.Y = 2, 
#'         group.name = c("control", "case"))
#' plotPLS(asmbPLSDA.fit.multiclass, 
#'         comp.X = 1, 
#'         comp.Y = 2, 
#'         group.name = c("healthy", "mild", "severe"))
#' 
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @import ggplot2 ggpubr

plotPLS <- function(fit.results, comp.X = 1, comp.Y = 2, group.name = NULL, legend = TRUE) {
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
  
  df <- as.data.frame(cbind(fit.results$X_super_score[, comp.X], fit.results$X_super_score[, comp.Y], Group))
  colnames(df) <- c("X", "Y", "Group")
  df$X <- as.numeric(df$X)
  df$Y <- as.numeric(df$Y)
  df$Group <- as.factor(df$Group)
  p_output <- NULL
  eval(parse(text = paste0("p_output <- ggplot(data = df, aes(X, Y, color = Group)) + geom_point(show.legend = ", legend, ") + 
                           stat_ellipse(show.legend = FALSE) + xlab(\"PLS ", comp.X, "\") + ylab(\"PLS ", comp.Y, "\") + theme_bw()")))
  return(p_output)
}