analysis_differences <- function(data, variables_of_interest, factor) {

    # Select relevant columns
  filtered_data <- data[,c(factor,variables_of_interest)] # escogemos las columnas de interes

  # Remove rows with NA values
    
  # Apply linear models and perform ANOVA
  results <- lapply(variables_of_interest, function(var) {


    model <- lm(as.formula(paste(factor, "~", var)), data = filtered_data)    
    anova_result <- anova(model)

    
    return(list( p.value = anova_result[1,"Pr(>F)"]))
  })
  
  names(results) <- variables_of_interest
  results <- unlist(results)
  return(results)
}

matrix.pvalues <- function(data, variables_of_interest){


factors <- colnames(data)[grep("Factor",colnames(data))]

result <- lapply(factors,function(factor){

res <- analysis_differences(data = meta,
variables_of_interest = variables_of_interest, 
factor = factor)

}
)

result <- do.call(rbind,result)
rownames(result) <- factors

return(result)
}
    