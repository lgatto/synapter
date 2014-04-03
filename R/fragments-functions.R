#' filter Neutral Loss
#' @param df data.frame (needs a column Neutral.LossType)
#' @return a data.frame that contains only rows with Neutral.LossType == "None"
.filterNeutralLoss <- function(df) {
  return(df[which(df$Neutral.LossType == "None"), ])
}

#' read final_fragments.csv
#' @param file file path
#' @param removeNeutralLoss remove rows with neutral loss != "none"?
#' @return a data.frame
.readFragements <- function(file, removeNeutralLoss=FALSE) {
  df <- read.csv(file, stringsAsFactors=FALSE)

  if (removeNeutralLoss) {
    return(.filterNeutralLoss(df))
  } else {
    return(df)
  }
}

