#' filter Neutral Loss
#' @param df data.frame (needs a column Neutral.LossType)
#' @return a data.frame that contains only rows with Neutral.LossType == "None"
.filterNeutralLoss <- function(df) {
  return(df[which(df$Neutral.LossType == "None"), ])
}
