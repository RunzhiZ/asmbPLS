result_compare <- function (prediction, true) {
  result_abs = sum(abs(prediction-true))/length(true)
  result_sq = sum((prediction-true)^2)/length(true)
  return(list(result_sq = result_sq,
              result_abs = result_abs))
}