library(dplyr)

cDNA_synth <- function(data, mass) {
  data <- data %>%
    mutate(RNA_uL = mass/Concentration) %>%
    mutate(iScript_5X_uL = 4) %>%
    mutate(RTase_uL = 1) %>%
    mutate(MQ_H2O_uL = 20 - (RNA_uL + iScript_5X_uL + RTase_uL)) %>%
    mutate(total_uL = 20)
  return(data)
}

generate_qPCR <- function(num_sample, trials, genotypes, genes) {
  num_sample_error <- (num_sample*trials*genotypes) +1
  total_vol <- num_sample_error * 10
  cocktail <- tibble(matrix(ncol = 5, nrow = 0))
  i <- length(genes)
  for(x in 1:i) {
    new_row <- c(genes[x], total_vol*0.5, total_vol*0.1, total_vol*0.2, total_vol*0.2)
    cocktail <- rbind(cocktail, new_row)
  }
  colnames(cocktail) <- c("gene", "tb_green", "primer", "mq_h2o", "template")
  return(cocktail)
}

