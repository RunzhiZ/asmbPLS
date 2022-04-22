KM.estimator <- function(survival_data) {
  colnames(survival_data)[1:2] <- c("Survival_time", "Event_indicator")
  ## event indicator = 1 with no repeat
  survival_unique_order <- unique(survival_data[which(survival_data[,"Event_indicator"] == 1), 1:2])[,"Survival_time"]
  ## KM estimator of failure times
  KM_table <- matrix(NA, nrow = length(survival_unique_order), ncol = 5)
  colnames(KM_table) <- c("Time", "Number_of_failure", "Number_at_risk", "Survival_at_this_point", "Survival_total")
  for (i in 1:length(survival_unique_order)) {
    if(i == 1){
      previous_survival <- 1
    } else {
      previous_survival <- KM_table[i-1, "Survival_total"]
    }
    KM_table[i, "Time"] <- survival_unique_order[i]
    KM_table[i, "Number_of_failure"] <- length(which((survival_data[, "Survival_time"] == KM_table[i, "Time"]) & 
                                                       (survival_data[, "Event_indicator"] == 1)))
    KM_table[i, "Number_at_risk"] <- length(which(survival_data[, "Survival_time"] >= KM_table[i, "Time"]))
    KM_table[i, "Survival_at_this_point"] <- 1 - KM_table[i, "Number_of_failure"]/KM_table[i, "Number_at_risk"]
    KM_table[i, "Survival_total"] <- KM_table[i, "Survival_at_this_point"] * previous_survival
  }
  return(KM_table)
}