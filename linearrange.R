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

source("rtkin.R") # not technically a package, but still the ideal spot

library(zoo)


## 1.0 Determining a Linear Range ==============================================


## Formats user-input into a vector of timepoints that do not represent the 
##    desired linear range.
##    range_description is list(lin_min, lin_max, outliers)
##    dataset must be a dataframe wherein the first column contains the time at
##    which observations are taken, and other columns contain one test of
##    observations

format_non_lin_range <- function(range_description, data){
  min <- range_description[[1]]
  max <- range_description[[2]]
  out <- range_description[[3]]
  lin_range <- format_lin_range(min, max, out, data)
  
  timepoints <- as.vector(data[,1])
  
  non_lin <- c()
  
  for (i in timepoints){
    if (!(i %in% lin_range)){
      non_lin <- c(non_lin, ((i/2) + 1))
    }
  }
  non_lin
}


## Formats user-input into a vector of timepoints to be fit to a linear model 
##    and used to determine the velocity.
##    
##    lin_min must be an integer, representing a time-point
##    lin_max must be an integer, representing a time-point
##    outliers must be a vector of integers, representing a time-point
##    dataset must be a dataframe wherein the first column contains the time at
##    which observations are taken, and other columns contain one test of
##    observations

format_lin_range <- function(lin_min, lin_max, outliers, dataset){
  
  if(validate_time(lin_min, dataset) & validate_time(lin_max, dataset)) {
    
    if(length(outliers) == 0) {
      
      return(seq(lin_min, lin_max, 2))
      
    } else {
      
      times <- seq(lin_min, lin_max, 2)
      
      range <- c()
      
      for (i in times) {
        if(!(i %in% outliers)){
          range <- c(range, i)
        }
      }
      return(range)
    }
  }
}


## Suggests a range of timepoints for the linear increasing portion of the
##    progress curve based on clustering according to slope between data points
##
##    dataset must be a dataframe wherein the first column contains the time at
##    which observations are taken, and other columns contain one test of
##    observations
##    test must be the column number of the test for which a linear 
##    range is to be suggested, it must be a valid test_name


estimate_lin_range <- function(dataset, test){
  
  test_name <- paste("Mean [Glucose] (mM)", test) 
  
  if (validate_col(test_name, dataset)){
    
    data_table <- data.frame(x=as.vector(dataset[,1]), y=as.vector(dataset[[test_name]]))
    
    # This sets up the model to be used to determine the slope
    f <- function(table){
      
      model <- lm(y~x, as.data.frame(table))
      
      coef(model)[2]
    }
    
    co <- rollapply(data_table, 3, f, by.column = F)
    co.cl <- kmeans(co, 2)
    b.points <- which(co.cl$cluster == match(max(co.cl$centers), co.cl$centers))+1
    
    return(b.points)
  }
}


## 2.0 Identifying the Range and Removing Data Not Included ====================


## Replaces the non-linear timepoints of data in the progress curve set with N/A
##
##    lin_ranges must be a list(list(test_name, range_description))
##    where a range_description is list(lin_min, lin_max, outliers)
##    dataset must be a dataframe wherein the first column contains the time at
##    which observations are taken, and other columns contain one test of
##    observations

remove_non_lin <- function(lin_ranges, data){
  
  data <- as.data.frame(data)
  
  for (i in lin_ranges){
    test <- i[[1]]
    non_lin_time <- format_non_lin_range(i[[2]], data)
    
    for (i in non_lin_time){
      
      test_name <-paste("Mean [Glucose] (mM)", test)
      
      data[[test_name]][i] <- NA
    }
  }
  data
}


## 3.0 Verifying Whether the Input is Valid ====================================


## Validates the linear range input given for the assay data.  Returns True if
##    the input is valid and False otherwise.
##
##    lin_ranges must be a list(list(test_name, range_description)) where a 
##    range_description is list(lin_min, lin_max, outliers)

validate_lin_range <- function(lin_ranges, data_set){
  
  for (i in seq(length(lin_ranges))){
    
    test <- lin_ranges[[i]][[1]]
    
    test_name <-paste("Mean [Glucose] (mM)", test)
    
    if (!(validate_col(test_name, data_set))){
      
      print(paste("The inputted test_name,", lin_ranges[[i]][[1]],", is not a valid test_name."))
      
      return(FALSE)
      
    } else {
      
      time_points <- append(lin_ranges[[i]][[2]][[1]],lin_ranges[[i]][[2]][[2]])
      
      time_points <- append(lin_ranges[[i]][[2]][[3]], time_points)
      
      for (t in time_points){
        
        if (!(validate_time(t, data_set))){
          print(paste("The inputted time_point ", t," is not a valid time_point."))
          
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}


## Validates whether a test_name is a valid test_name for this assay data.  
##    Returns True if it is valid, and False otherwise.
##    
##    cnum can be any
##    dataset must be a dataframe whereing the first column contains the time
##    at which the measurements were taken

validate_col <- function(cnum, dataset){
  
  if (cnum %in% colnames(dataset)){
    TRUE
  } else {
    FALSE
  }
}


## Validates whether a timepoint is a valid timepoint for this assay data.  
##    Returns True if it is valid, and False otherwise.
##    
##    time can be any
##    dataset must be a dataframe whereing the first column contains the time
##    at which the measurements were taken

validate_time <- function(time, dataset){
  
  valid_timepoints <- as.vector(dataset[,1])
  
  if (time %in% valid_timepoints){
    TRUE
  } else {
    FALSE
  }
}


## 4.0 Knitting All Steps Together =============================================


## Calculates the parameters of modified kinetics (contains NA for items in the
##    non-linear ranges).
##
##    kin_data is normalized and averaged kinetics data wherein the rows are
##    different timepoints and the columns different test series
##    lin_ranges must be a list(list(test_name, range_description)) where a 
##    range_description is list(lin_min, lin_max, outliers)

calc_mod_kinetics <- function(analyzed_data, lin_ranges, export = TRUE, plot = TRUE){
  
  kinetics_analyzed_linear <- calc_kinetics_lin(analyzed_data, lin_ranges)
  
  if (export == TRUE) {
    ## Exports the analyzed data if the user wishes to view.
    
    write.csv(kinetics_analyzed_linear[[1]],
              "analyzed_test_linear.csv",
              row.names = FALSE)
    
    write.csv(kinetics_analyzed_linear[[2]],
              "test_params_linear.csv",
              row.names = TRUE)
  }
  
  if (plot == TRUE) {
    ## Exports the real-time kinetics plot for the enzyme tests.
    
    kinetics_plot_linear <- plot_kinetics_lin(kinetics_analyzed_linear[[1]])
    
    ggsave(file="kinetics_plot_linear.png", plot = kinetics_plot_linear)
  }
  
  kinetics_analyzed_linear
}


## Plots the enzyme test kinetics data (glucose concentration over time)
##
##    kin_data must be a valid series of [glucose_conc] over time data
##    as formatted by calc_kinetics_lin

plot_kinetics_lin <- function(kin_data){
  
  plot_data <- gather_data_lin_plot(kin_data)
  
  indiv_data <- as.data.frame(plot_data[[1]])
  
  avg_data <- as.data.frame(plot_data[[2]])
  
  df_indiv <- melt(data = indiv_data, id.vars = "Time")
  
  df_avg <- melt(data = avg_data, id.vars = "Time")
  
  enzyme_plot <-  ggplot(data = df_indiv, aes(x = Time, y = value, color = variable)) +
    geom_point(size = 0.75) +
    geom_smooth(data = df_avg, method = "lm", se = FALSE, size = 1) +
    labs(title = "Progress Curve Plot for Enzyme Tests", #check title name
         x = "Time (min)",
         y = "Product, [Glucose] (mM)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}


## Formats the data given by calc_kinetics_lin into two different data.tables
##    to be plotted by plot_kinetics_lin
##
##    kin_data must be a valid series of [glucose_conc] over time data
##    as formatted by calc_kinetics_lin

gather_data_lin_plot <- function(kin_data){
  
  avg_merged <- kin_data[,1]
  
  indiv_merged <- kin_data[,1]
  
  avg_names <-c("Time")
  
  indiv_names <- c("Time")
  
  for (i in seq(2, length(colnames(kin_data)), 5)){
    
    avg_merged <- cbind(avg_merged, kin_data[,(i+3)])
    avg_names <- c(avg_names, colnames(kin_data)[(i + 3)])
    
    indiv_merged <- cbind(indiv_merged, kin_data[,i:(i+2)])
    indiv_names <- c(indiv_names, colnames(kin_data)[i:(i + 2)])
  }
  
  colnames(avg_merged) <- avg_names
  
  colnames(indiv_merged) <- indiv_names
  
  list(indiv_merged, avg_merged)
}

## Replaces the non-linear data points for each test with NA and recalculates
##    the linear regression parameters based on this.
##
##    kin_data is a list with the first item being normalized and averaged 
##    kinetics data wherein the rows are different timepoints and the columns 
##    different test series, and the second item are linear regression params
##    lin_ranges must be a list(list(test_name, range_description)) where a 
##    range_description is list(lin_min, lin_max, outliers)

calc_kinetics_lin <- function(analyzed_kin, lin_ranges){
  
  kin_data <- as.data.frame(analyzed_kin[[1]])
  
  ## Removes the non-linear ranges
  
  if (validate_lin_range(lin_ranges, kin_data)){
    
    linear_data <- remove_non_lin(lin_ranges, kin_data)
    
    lm_params <- data.frame()
    
    row_names <- c()
    
    ## Calculates regression parameters
    
    for (i in seq(5, length(colnames(linear_data)), 5)){
      
      test_data_i <- as.data.frame(cbind(linear_data[,1], linear_data[,i]))
      
      test_i_name <- colnames(kin_data)[i]
      
      params_i <- build_linear_model_linear(test_data_i)
      
      lm_params <- rbind(lm_params, params_i)
      
      row_names <- c(row_names, test_i_name)
      
      rownames(lm_params) <- row_names
    }
  } else{
    return(NULL)
  }
  list(linear_data, lm_params)
}

## Performs simple linear regression on the glucose kinetics data for a single
##    linearized enzyme test (enzyme + concentration parameters).  Returns the 
##    coefficients describing the linear model
##
##    lin_test_data must be a data frame where the first column is time and the
##    second column contains the normalized glucose concentration

build_linear_model_linear <- function(lin_test_data){
  
  df <- lin_test_data
  
  model <- lm(df[,2] ~ df[,1], data = df)
  
  coefficients <- as.vector(model[[1]])
  
  intercept <- coefficients[1]
  
  slope <- coefficients[2]
  
  r_sqr <- summary(model)$r.squared
  
  model_coeffs <- data.frame(slope, intercept, r_sqr)
  
  colnames(model_coeffs) <- c("m", "b", "r_sqr")
  
  model_coeffs
}
