#' Relevance plot for asmbPLS-DA
#'
#' Function to visualize the most relevant features (relevant to outcome) in 
#' each block.
#' 
#' @param fit.results The output of \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' @param n.top A integer indicating the number of the most relevant features to 
#' be displayed for each block. The default is 10. If the number of selected 
#' features in block is smaller than \code{n.top}, all the selected features in
#' that block will be displayed.
#' @param ncomp Which component to plot from each block. Should not be larger 
#' than the number of PLS components used (\code{PLS.comp}) in the function
#' \code{\link[asmbPLS]{asmbPLSDA.fit}}. The default is 1.
#' @param block.name A vector containing the named character for each block. It
#' must be ordered and match each block.
#' 
#' @return 
#' none
#' 
#' @details 
#' The function returns a plot to show the most relevant features for each 
#' block.
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
#' ## visualization to show the most relevant features in each block
#' plotRelevance(asmbPLSDA.fit.binary)
#' plotRelevance(asmbPLSDA.fit.multiclass)
#' ## custom n.top and block.name
#' plotRelevance(asmbPLSDA.fit.binary, 
#'               n.top = 5,
#'               block.name = c("mRNA", "protein"))
#' plotRelevance(asmbPLSDA.fit.multiclass, 
#'               n.top = 7,
#'               block.name = c("miRNA", "protein"))
#' 
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @import ggplot2 ggpubr

plotRelevance <- function(fit.results, n.top = 10, ncomp = 1, block.name = NULL) {
  if(is.null(block.name)) {
    n.block <- length(fit.results$X_dim)
    block.name <- paste0("block.", 1:n.block)
  } else{n.block <- length(block.name)}
  X_weight <- fit.results$X_weight
  X_super_weight <- fit.results$X_super_weight
  
  for(i in 1:n.block) {
    non_zero_index <- which(X_weight[[i]][, ncomp] != 0)
    df <- as.data.frame(cbind(names(non_zero_index), X_weight[[i]][non_zero_index, ncomp]))
    colnames(df) <- c("feature", "value")
    df$value <- as.numeric(df$value)
    df$weight <- ifelse(df$value > 0, "positive", "negative")
    df$weight <- factor(df$weight, levels = c("positive", "negative"))
    df$feature <- factor(df$feature, levels = df$feature[order(abs(df$value))])
    df$block <- block.name[i]
    df <- df[order(df$value, decreasing = T), ]
    ## Take Top n
    if(n.top <= nrow(df)) {
      df <- df[1:n.top, ]
    }
    eval(parse(text = paste0("p", i, "<- ggplot(df, aes(feature, value, fill = weight)) + 
                           geom_bar(stat = \"identity\", show.legend = FALSE) + coord_flip() + xlab(\"\") + ylab(\"Feature weight\") + 
                           ggtitle(\"Top features in \'", block.name[i], "\'\\n Block weight: ", 
                             round(X_super_weight[i, ncomp], 2), "\") + theme(plot.title = element_text(hjust = 0.5, size = 20),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = \"black\"),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = \"black\"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = \"black\"),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = \"black\", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()
        )")))
    
  }
  p_output <- NULL
  eval(parse(text = paste0("p_output <- ggarrange(", paste(paste("p", 1:n.block, sep = ""), collapse = ", "), ", common.legend = TRUE)")))
  return(p_output)
}