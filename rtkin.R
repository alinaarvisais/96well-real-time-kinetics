## This script is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This script is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this script.  If not, see <https://www.gnu.org/licenses/>.


## 0.0 Load Required Packages ==================================================


source("formatr.R") # not technically a package, but still the ideal spot

library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(dplyr)

## 1.0 Make the Concentration Curve ============================================


## Returns the constructed linear regression model for the glucose standard curve 
##     data for future use by the script.  Also exports .csv file of the
##     parameters if export = TRUE (default value), but does not if False.  Also
##     exports a plot of the standard curve is plot = True (the default value), 
##     but does not if False.
##
##    std_table must be standard curve data formatted as by make_standard_table
##    export must be one of True or False or left blank
##    plot must be one of True or False or left blank

make_concentration_curve <- function(std_table, export = TRUE, plot = TRUE){
  
  ## Prepares the replicate absorbances to be plotted
  std_tbl <- std_table
  
  melted_table <- std_tbl[,1]
  
  for (i in seq(2,length(colnames(std_tbl)))){
    row <- as.vector(std_tbl[,i])
    melted_table <- rbind(melted_table, row)
  }
  
  melted_table <- cbind(as.numeric(colnames(std_tbl)),melted_table)
  
  row.names(melted_table) <- c()
  
  colnames(melted_table) <- c("glucose_conc", 1, 2, 3)
  
  melted_table <- as.data.frame(melted_table)
  
  melt <- melt(melted_table, id.vars = "glucose_conc")
  
  ## Stores the table of standard curve data
  standard_curve_data <- make_curve_data(std_table)
  
  ## Stores the coefficients of the linear model (mean_abs = m * g_conc + b)
  standard_curve_params <- build_linear_model(standard_curve_data)
  
  if (export == TRUE) {
    ## Exports the std data if the user wishes to analyze by hand.
    
    write.csv(standard_curve_params,
              "standard_curve_params.csv",
              row.names = TRUE)
  }
  
  if (plot == TRUE){
    ## Produces the standard curve_plot
    
    std_curve_plot <- ggplot(data = standard_curve_data, 
                             mapping = aes(x = glucose_conc, y = mean_abs)) + 
      geom_point(size = 2) +
      geom_smooth(method = lm, se = FALSE, size = 1) +
      geom_point(data = melt, mapping = aes(x = glucose_conc, y = value), size = 2) +
      labs(title = "Standard Curve for Glucose Absorbance at 450nm", #check wavelength
           x = "Concentration of Glucose (mM)",
           y = "Absorbance at 450 nm") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Saves the plot
    
    ggsave(file="standard_curve.png", plot = std_curve_plot)
  }
  standard_curve_params
}


## Performs simple linear regression on the standard curve data.  Returns the 
##    coefficients describing the linear model fitting the data
##
##    std_crv_dt must be a data frame containing two columns ('mean_abs' and 
##    'glucose_conc'), one representing the means and the second containing  
##    their respective glucose concentration

build_linear_model <- function(std_crv_dt){
  
  model <- lm(mean_abs ~ glucose_conc, data = std_crv_dt)
  
  coefficients <- as.vector(model[[1]])
  
  intercept <- coefficients[1]
  
  slope <- coefficients[2]
  
  r_sqr <- summary(model)$r.squared
  
  model_coeffs <- data.frame(slope, intercept, r_sqr)
  
  colnames(model_coeffs) <- c("m", "b", "r_sqr")
  
  model_coeffs
}


## Creates a data frame with two columns, the first containing the mean of each
##      replicate for a glucose standards and the second containing the glucose
##      concentration for that standard.  Returns this data frame
##
##      std_tbl must be a data frame containing the absorbance data for all
##      replicates of a standard organized by col = concentration, row = rep#

make_curve_data <- function(std_tbl) {
  
  s1_mean <- mean(std_tbl$`2`) #change the column names if they are diff
  
  s2_mean <- mean(std_tbl$`1`)
  
  s3_mean <- mean(std_tbl$`0.5`)
  
  s4_mean <- mean(std_tbl$`0.25`)
  
  vector_mean <- c(s1_mean, s2_mean, s3_mean, s4_mean)
  
  vector_conc <- c(2, 1, 0.5, 0.25)
  
  std_curve_data <- data.frame(vector_mean, vector_conc)
  
  colnames(std_curve_data) <- c("mean_abs", "glucose_conc")
  
  std_curve_data
}


## 2.0 Analyze Kinetics ========================================================


## Analyzes the test data for a list of tests (combination of enzyme and subs 
##    conc.), to give the R_sqr, slope, and normalized glucose vs. time 
##    information.
##
##    test_ids must be a list of test_identifiers, i.e. a list of
##            list(test_id, subs conc, list(wellIds))
##    std_curve must be linear model parameters as prepared by make_concentration_curve
##    plot = TRUE will export the plot as a pdf and if FALSE will not
##    export = TRUE will export the data as a pdf and if FALSE will not

calc_test_kinetics <- function(test_ids, std_curve, export = TRUE, plot = TRUE){
  
  kinetics_analyzed <- calc_kinetics(test_ids, std_curve)
  
  ##kinetics_deviation <- kinetics_analyzed[[1]]
  ##kinetics_normalized <- kinetics_analyzed[[2]]
  
  if (export == TRUE) {
    ## Exports the analyzed data if the user wishes to view.
    
    write.csv(kinetics_analyzed[[1]],
              "analyzed_test.csv",
              row.names = FALSE)
    
    write.csv(kinetics_analyzed[[2]],
              "test_params.csv",
              row.names = TRUE)
  }
  
  if (plot == TRUE) {
    ## Exports the real-time kinetics plot for the enzyme tests.
    
    kinetics_plot <- plot_kinetics(kinetics_analyzed)
    
    ggsave(file="kinetics_plot.png", plot = kinetics_plot)
  }
  kinetics_analyzed
}


## Plots the enzyme test kinetics data (glucose concentration over time)
##
##    kin_data must be a valid series of [glucose_conc] over time data

plot_kinetics <- function(kin_data){
  kin_time_series <- as.data.table(kin_data[[3]])
  
  df <- melt(data = kin_time_series, id.vars = "Time")
  
  df.summary <- df %>% group_by(Time, variable) %>% summarise(sd = sd(value), value = mean(value))
  
enzyme_plot <-  ggplot(data = df, aes(x = Time, y = value, color = variable)) +
    geom_point(size = 0.75) + 
    geom_smooth(method = "lm", se = FALSE, size = 1) +
    #geom_errorbar(aes(ymin = value-sd, ymax = value+sd), data = df.summary, width = 2)+
    labs(title = "Progress Curve Plot for Enzyme Tests", #check title name
         x = "Time (min)",
         y = "Product, [Glucose] (mM)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


## Analyzes the test data for a list of tests (combination of enzyme and subs 
##    conc.), to give the R_sqr, slope, and normalized glucose vs. time 
##    information.  Returns a list where the first item is a data frame of 
##    averaged and normalized glucose vs. time data, and the second item 
##    is the linear regression coeffecients
##
##    test_ids must be a list of test_identifiers, i.e. a list of
##    list(test_id, subs conc, list(wellIds))
##    std_curve must be linear model parameters as prepared by make_concentration_curve

calc_kinetics <-function(test_ids, std_curve){
  
  test_kin_params <- function(ids){
    test_kinetics(ids, std_curve)
  }
  
  list_kin <- lapply(test_ids, test_kin_params)
  
  gluc_merged <- ((list_kin[[1]])[[2]])[,1]
  
  lm_params <- data.frame()
  
  row_names <- c()
  
  col_names <- c("Time")
  
  plot_names <-c("Time")
  
  plot_merged <- ((list_kin[[1]])[[2]])[,1]
  
  for (i in 1:length(list_kin))
  {
    test_i <- list_kin[[i]]
    test_i_name <- (test_ids[[i]])[[1]]
    test_i_conc <- (test_ids[[i]])[[2]]
    gluc_i <- test_i[[2]][,2:6]
    plot_i <- test_i[[2]]
    param_i <- test_i[[3]]
    
    gluc_merged <- cbind(gluc_merged, gluc_i)
    col_names <- c(col_names, colnames(gluc_i))
    
    plot_merged <- cbind(plot_merged, plot_i[,2:4])
    series_name <- paste(test_i_name, test_i_conc, sep = "_")
    plot_names <- c(plot_names, series_name, series_name, series_name)
    
    lm_params <- rbind(lm_params, param_i)
    row_names = c(row_names, paste(test_i_name, test_i_conc, sep = "_"))
  }
  
  colnames(gluc_merged) <- col_names
  
  colnames(plot_merged) <- plot_names
  
  rownames(lm_params) <- row_names
  
  list(gluc_merged, lm_params, plot_merged)
}


## Analyzes the test data for one test (combination of enzyme and subs conc.),
##    to give the R_sqr, slope, and normalized glucose vs. time information
##    Returns a list(test_name, test_table, lm_parameters)
##
##    test_id must be a test_identifier, i.e. list(test_id, subs conc, list(wellIds))
##    std_curve must be linear model parameters as prepared by make_concentration_curve

test_kinetics <- function(test_id, std_curve){
  
  data_table <- make_test_table(test_id)
  
  test_str <- paste(test_id[[1]], test_id[[2]], sep = "_")
  
  if (length(test_id) == 4){
    
    glc_table <- prep_glucose_data(data_table, std_curve, test_id[[4]])
    
  } else {
    
    glc_table <- prep_glucose_data(data_table, std_curve, 0)
  }
  
  kin_table <- cbind(as.numeric(row.names(glc_table)), glc_table)
  
  colnames(kin_table) <- c("Time", prepare_col_glc(test_str, glc_table))
  
  lm <- build_linear_model_test(as.data.frame(kin_table))
  
  list(test_id[[1]], kin_table, lm)
}


## Prepares column names such that each column is appended with the test_name
##    appropriate for that set of columns.  Returns the vector of column names.
##
##    test_id_str must be a string representing the test_id
##    glc_table may be any data frame containing at least one column

prepare_col_glc <- function(test_id_str, glc_table){
  
  col_names <- colnames(glc_table)
  
  new_names <- c(col_names[1:3])
  
  for (i in seq(4, 5)){
    new_names <- c(new_names, paste(col_names[[i]], test_id_str, sep = " "))
  }
  new_names
}


## Performs simple linear regression on the glucose kinetics data for a single
##    enzyme test (enzyme + concentration parameters).  Returns the coefficients
##    describing the linear model
##
##    norm_test_data must be a data frame where the first column is time and the
##    fifth column contains the normalized glucose concentration

build_linear_model_test <- function(norm_test_data){
  
  df <- norm_test_data
  
  model <- lm(df[,5] ~ df[,1], data = df)
  
  coefficients <- as.vector(model[[1]])
  
  intercept <- coefficients[1]
  
  slope <- coefficients[2]
  
  r_sqr <- summary(model)$r.squared
  
  model_coeffs <- data.frame(slope, intercept, r_sqr)
  
  colnames(model_coeffs) <- c("m", "b", "r_sqr")
  
  model_coeffs
}


## Prepares the kinetic data in terms of [glucose] and time for an enzyme test.
##    Returns the two column dataframe containing the time in one column and
##    the normalized test_data.
##
##    test_table must be a 3 col table with at least one row, where each col 
##    represents different replicates for each timepoint
##    std_curve must be linear model parameters as prepared by make_concentration_curve


prep_glucose_data <- function(test_table, std_curve, offset){
  
  glucose_table <- convert_to_glucose(test_table, std_curve, offset)
  
  normalized_table <- convert_to_normalized(glucose_table)
  
  averaged_table <- average_abs(normalized_table)
  
  dev_table <- calc_std_dev(averaged_table)
  
  dev_table
}


## Calculcates the standard deviation corresponding to the [glucose] (mM)
##    for each measured timepoint.
##    Returns these values in a new column at the right of the table.
##
##    avg_table must be a data frame with the first three column containing the 
##    normalized data and the fourth column contains the averaged concentration.

calc_std_dev <- function(avg_table){
  
  std_dev <- c()
  
  for(i in 1:num_timepoints) 
  {
    std_dev <- c(std_dev, sd(avg_table[i, 1:3]))
  }
  
  dev_table <- cbind(avg_table, std_dev)
  
  colnames(dev_table) <- c(colnames(avg_table), "Std. Dev.")
  
  dev_table
} 


## Averages the absorbance values for all three replicates of a test a each
##    measured timepoint.  Returns the data frame with these values in a new col 
##    at the end of the table
##    
##    test_table must be a data frame containing at least one row and one column

average_abs <- function(test_table){
  
  averages <- c()
  
  for(i in 1:num_timepoints) 
  {
    averages <- c(averages, mean(test_table[i,]))
  }
  
  averages
  
  averaged_table <- cbind(test_table, averages)
  
  colnames(averaged_table) <- c(colnames(test_table), "Mean [Glucose] (mM)")
  
  averaged_table
}


## Calculcates the normalized [glucose] (mM) corresponding to the [glucose] (mM)
##    for each measured timepoint. Normalized to the measurement at time 0. 
##    Returns these values in place of the current values for each well column.
##
##    glc_table must be a data frame

convert_to_normalized <- function(glc_table){
  
  for (c in seq(length(colnames(glc_table))))
  {
    norm <- glc_table[1, c]
    
    for (r in seq(length(row.names(glc_table))))
    {
      glc_conc <- glc_table[r, c]
      norm_conc <- (glc_conc - norm)
      glc_table[r, c] <- norm_conc
    }
  }
  glc_table
} 


## Calculcates the [glucose] (mM) corresponding to the  absorbance values
##    for each measured timepoint according to the linear regression with std 
##    values.  Returns the data frame with these values replacing the current
##    absorbance values
##
##    test_table must be a data frame containing at least one row and one column
##    std_curve must be linear model parameters as prepared by make_concentration_curve
##    offset must be an integer value denoting the blank offset to be applied

convert_to_glucose <- function(test_table, std_curve_params, offset){
  
  m <- std_curve_params[1,1]
  
  b <- std_curve_params[1,2]
  
  for (c in seq(length(colnames(test_table))))
  {
    for (r in seq(length(row.names(test_table))))
    {
      if (!(is.na(test_table[r, c]))){
        abs <- as.numeric(test_table[r, c]) - offset
        glucose_conc <- ((abs - b) / m)
        test_table[r, c] <- glucose_conc
      }
    }
  }
  test_table
}


## 3.0 Apply Effect of Blanks ==================================================


##  Apply_blanks applies the effect of the blanks onto the test identifiers.  If
##    the blank is greater than the threshold, it applies the blank values to 
##    its respective test identifiers (if an enzyme blank) or removes its
##    test_identifiers from the list (if a substrate blank).  If it is lower
##    than the threshold, no action is take to the list of test identifiers.
##
##    test_identifiers must be a list where each item is a test_identifier
##    blanks must be a list where each item is a blank_identifier
##    threshold must be a number greater than 0 (the default is 0.1)

apply_blanks <- function(test_identifiers, blanks, threshold = 0.1){
  
  for (i in blanks){
    
    blank_table <- set_blank(i)
    
    blank_offset <- check_blank(blank_table, threshold)
    
    if (typeof(blank_offset) == "logical"){
      
      print(paste("Effects from the ", i[[1]], "blank, ", i[[2]], ",are being considered."))
      
      if (i[[1]] == "substrate"){
        
        test_identifiers <- remove_tests(test_identifiers, i[[2]])
        
        print(paste("The", i[[1]], "blank, ", i[[2]], ", is above the threshold.  It's associated tests have been removed."))
        
      } else {
        print(paste("Offset from the ", i[[1]], "blank, ", i[[2]], ",has been applied."))
        
        test_identifiers <- add_blank_offset(test_identifiers, i[[2]], blank_offset)
      }
    } else {
      print(paste("No action from", i[[1]], "blank", i[[2]], "are taken."))
    }
  }
  test_identifiers
}


## Set blank returns a data table of the absrobance values of the blank wells in 
##    individual columns with the first column containing the time points.
##
##    blank_id should be a blank_description list(type, name, wellIDs) where
##    type is one of "substrate" or "enzyme"

set_blank <- function(blank_id){
  
  make_test_table(blank_id)
}


## Check blanks returns a data table of averaged absorbance values for the
##    replicates of a blank if the difference between the absorbance at the 
##    first and last timepoint is higher than 'threshold' and False otherwise.
##    
##    threshold must be a number greater than 0
##    blank_data must be a data table containing the absorbance values of the 
##    blank wells

check_blank  <- function(blank_data, threshold){
  
  averaged <- average_abs(blank_data)
  
  last_timepoint <- length(row.names(averaged))
  
  abs_difference <- as.numeric(averaged[last_timepoint, 4]) - as.numeric(averaged[1, 4])
  
  if (abs_difference >= threshold){
    
    return(averaged)
    
  } else {
    return(FALSE)
  }
}


## remove_tests returns a the list of test_identifiers 'test_ids' with all tests 
##    containing the substrate name 'subs' removed.
##
##    test_identifiers must be a list where each item is a test_identifier
##    subs must be a string containing a substrate name

remove_tests <- function(test_identifiers, subs){
  
  retained <- list()
  
  for (i in test_identifiers){
    if (!(grepl(subs, i[[1]], fixed = TRUE))){
      retained <- append(retained, i)
    }
  }
  
  retained
  
}


## add_blank_offset returns a list of test_identifiers 'test_ids' with all tests
##    containing the enzyme name 'enzyme' having the blank_data table attached.
##
##    test_identifiers must be a list where each item is a test_identifier
##    enzyme must be a string containing an enzyme name
##    offset must be a data table containing the average abs value of the blank

add_blank_offset <- function(test_identifiers, enzyme, offset){
  
  modified <- list()
  
  for (i in test_identifiers){
    
    if (!(grepl(enzyme, i[[1]], fixed = TRUE))){
      
      modified <- append(modified, append(i, offset))
      
    } else {
      
      modified <- append(modified, i)
    }
  }
  modified
}
