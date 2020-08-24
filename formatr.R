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


library(purrr)


## 1.0 Organize the Standard Curve Data ========================================


## Returns the table for the glucose standard curve data for future analysis by 
##    the script.  Also exports .csv file of the data for future analysis if 
##    export = True, but does not if False (the default value)
##
##    s1, s2, s3, and s4 must be lists in the format list("WellID", "WellID", "WellID")
##    data must be tidied data such as that returned by tidy_plate_data
##    export must be one of True or False

make_standard_table <- function(s1, s2, s3, s4, data, export = FALSE){
  
  ## Stores the table for the glucose standard curve data.
  
  if (validate_ids(s1) == FALSE){
    print("Well IDs for standard 1 are invalid.")
    return()
    
  } else if (validate_ids(s2) == FALSE){
    print("Well IDs for standard 2 are invalid.")
    return()
    
  } else if (validate_ids(s3) == FALSE){
    print("Well IDs for standard 3 are invalid.")
    return()
    
  } else if (validate_ids(s4) == FALSE){
    print("Well IDs for standard 4 are invalid.")
    return()}
  
  standard_table <- create_std_table(s1, 
                                     s2, 
                                     s3, 
                                     s4,
                                     data)
  
  if (export == TRUE) {
    
    ## Exports the std data if the user wishes to analyze by hand.
    
    write.csv(standard_table,
              "standards.csv", 
              row.names = TRUE)
  }
  standard_table
}


## Makes a data frame (table) for the concentration curve data.  The columns
##    represent the concentrations of glucose and the rows represent different
##    replicates.  Returns the data frame.
##
##    s1, s2, s3, and s4 must be lists in the format list("WellID", "WellID", "WellID")
##    data must be tidied data such as that returned by tidy_plate_data

create_std_table <- function(s1, s2, s3, s4, data){
  
  data_s1 <- get_standard_data(s1, data)
  
  data_s2 <- get_standard_data(s2, data)
  
  data_s3 <- get_standard_data(s3, data)
  
  data_s4 <- get_standard_data(s4, data)
  
  std_table <- data.frame(data_s1, data_s2, data_s3, data_s4)
  
  colnames(std_table) <- c(2,1,0.5,0.25) # This is based on the concentrations for a PGO assay
  
  row.names(std_table) <- c("Replicate 1", "Replicate 2", "Replicate 3")
  
  std_table
}


## Fetches the absorbance data for the wells corresponding to the replicates of
##    one glucose standard at time = 0 and stores it as a vector of absorbance 
##    data.  Returns the vector of absorbances.
##
##    std_list in the format list("WellID", "WellID", "WellID")
##    data must be tidied data such as that returned by tidy_plate_data

get_standard_data <- function(std_list, data) {
  
  col = c()
  
  for(i in seq(length(std_list))) 
  {
    if (is.na(std_list[[i]])){
      
      col <- c(col, NA)
      
    } else {
      
      well_id <- std_list[[i]]
      
      col <- c(col, 
               as.numeric(data[which(data$WellIds == well_id),
                               "0"])) 
    }
  }
  col
}


## 2.0 Organize the Enzyme Test Data ===========================================


## Returns the tables for the enzyme test data for future analysis by the script.
##    Also exports .csv file of the data for future analysis if 
##    export = True (the default value), but does not if False
##
##    test_ids must be a list of test_identifiers list("Name", subs_conc, list("WellID"))
##    export must be one of True or False or left blank

make_test_data_tables <- function(test_ids, export = TRUE){
  
  # Future work will validate the test identifiers first
  
  test_tables <- make_test_tables(test_ids)
  
  if (export == TRUE) {
    
    ## Exports the std data if the user wishes to analyze by hand.
    
    export_table <- make_export_table(test_tables)
    
    export_test_data(export_table)
  }
  
  test_tables
}


## Exports a data frame representing all the enzyme test data into a .csv file
##
##    export_table must be a data frame representing the data to be exported

export_test_data <- function(export_table){
  
  write.csv(export_table,
            "test_data.csv",
            row.names = TRUE)
}


## Merges a list of data frames each representing an individual test identifier
##    into one large data frame.  Returns the merged data frame
##
##    test_tables_lst must be a list of data frames (representing a test
##    identifier each)

make_export_table <- function(test_tables_lst){
  
  export_table <- NULL
  
  for(i in 1:length(test_tables_lst)) 
  {
    export_table <- cbind(export_table, test_tables_lst[[i]])
  }
  
  export_table
}


## Makes a data frame of test data for each test in a list of test identifiers.
##    Returns a list where each item is a data frame representing data for each
##    test identifier
##
##    test_descriptions must be a list of test identifiers, described previously

make_test_tables <- function(test_descriptions){
  
  test_tables <- lapply(test_descriptions, make_test_table)
  
  test_tables
}


## Makes a data frame containing all the absorbance data for all replicates of a 
##    set of enzyme + substrate conditions.  Returns the created data frame
##
##    test_description must be a test identifier (three-item list described 
##    previously)

make_test_table <- function(test_description){
  
  test_name <- paste(test_description[[1]], test_description[[2]], sep = "_")
  
  well_list <- test_description[[3]]
  
  df_rows <- rows_to_df(well_list)
  
  df <- assign_x_y(df_rows, test_name)
}


## Assigns column and row names to a data frame containing the absorbance data
##    for all three replicates of a set of enzyme and substrate concentration
##    conditions  The columns names are the timepoints and the row names are the
##    test name appended with the replicate number.  Returns the data frame
##
##    df_rows must be a data frame (containing the rows of absorbance values)
##    test_name must be a name for a set of enzyme + substrate conditions

assign_x_y <- function(df_rows, test_name) {
  
  rownames(df_rows) <- timepoint_vector(num_timepoints)
  
  colnames(df_rows) <- make_ids_rep(test_name)
  
  df_rows
}


## Creates a new list of a test ID (ex: NtSI_2) concatenated with _'replicate#'
##    Returns the list of test ID's concatenated with the replicate number.
##
##    test_id must be a character

make_ids_rep <- function(test_id) {
  
  ## Number of replicates done for each test (combination of enzyme and susbstrate
  ##    concentration)
  
  replicate_num <- 3
  
  assign("replicate_num", replicate_num, envir = .GlobalEnv)
  
  (paste(test_id, seq_len(replicate_num), sep = "_"))
  
}


## Creates a vector containing all the timepoints at which reads where taken in
##    starting a 0 and in minutes
##
##    num_timepoints must be a natural number

timepoint_vector <- function(num_timepoints) {
  
  times <- c()
  
  for(i in 1:num_timepoints) 
  {
    times <- c(times, ((i-1)*timepoint_dif))
  }
  times
}


## Creates a data frame where each col is all the absorbance data for a well ID
##    assorted by increasing time
##
##    list_ids must be a list containing valid well IDs for a 96-well plate A-H 
##    concatenated with 1-12, ex: A12

rows_to_df <- function(list_ids) {
  
  df <- NULL
  
  for(i in 1:length(list_ids)) 
  {
    df <- cbind(df, extract_row(list_ids[[i]]))
  }
  
  df
}


## Extracts all absorbance values for a well ID over all the timepoints
##    contained in the tidied data table as a 'row'(vector).
##
##    well_id must be a valid well ID for a 96-well plate A-H concatenated with 
##    1-12, ex: A12

extract_row <- function(well_id) {
  
  if (is.na(well_id)){
    
    extracted_row <- c()
    
    for (i in seq(num_timepoints)){
      extracted_row <- c(extracted_row, NA)
    }
  } else {
    extracted_row <- tidy_data[which(tidy_data$WellIds==well_id),
                               2:(num_timepoints + 1)]
  }
  
  row_values <- as.vector(as.numeric(extracted_row))
}


## 3.0 Validate Well ID input ==================================================


## Validates whether a list of Well IDs are valid Well IDs.  Returns True if
##    they are valid well IDs and False otherwise.
##    
##    well_ids must be a list of any

validate_ids <- function(list_ids){
  
  ## Sets the valid well IDs for the standards locations
  
  valid_ids <- gen_well_ids(8, 12) # set for 96-well plate
  
  assign("valid_ids", valid_ids, envir = .GlobalEnv)
  
  for (i in (1:length(list_ids))){
    
    id <- list_ids[[i]]
    
    if (is.na(id)){
      return(TRUE)
      
    } else if (!(typeof(id) == "character")){
      return(FALSE)
      
    } else if (!(validate_id(id))){
      return(FALSE)
    }
  }
  TRUE
}


## Validates whether a well ID is a valid Well ID.  Returns True if it is valid,
##    and False otherwise.
##    
##    id can be any

validate_id <- function(id){
  if (id %in% valid_ids){
    TRUE
  } else {
    FALSE
  }
}

## Generates a vector of well IDs for a well plate.  In the vector, the
##    the rows are identified alphabetically and the columns are identified
##    numerically.  The resulted is ordered alphabetically then numerically
##
##    nrow must be a natural number and < 26 (=8 for 96-well plate)
##    ncol must be a natural number (=12 for 96-well plate)

gen_well_ids <- function(nrow, ncol){
  
  well_col_names <- as.character(seq_len(ncol))
  
  well_row_names <- as.character(LETTERS[1:nrow])
  
  paste_ids <- function(row_id) {
    (paste(row_id, well_col_names, sep = ""))
  }
  
  ids <- c()
  for (i in 1:length(well_row_names)) {
    ids <- c(ids, paste_ids(well_row_names[i]))
  }
  
  ids
}
