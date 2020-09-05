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

library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)

## 1. Importing the file with the right formatting =============================

## This creates a .csv file with the proper formatting for analysis, if it 
##    takes too long to run, can later omit this and do it manually with Excel.

## USER-INPUT: location of original data file.
##    example path  : C:/Users/Alina/Desktop/example_data.txt

raw_path <- readline(prompt = "Enter original data file path: ")

raw_lines <- readLines(raw_path)

## The code converts all tabs in a document to space so that there is a
##    single delimiter across the whole document

lines_tab_delim <- gsub("\\t$", "", raw_lines)

raw_table <- read.delim(text = lines_tab_delim,
                        skip = 3, 
                        header = FALSE)

## This identifies and selects the columns in the raw data table that contains
##    absorbance data.

last_col <- 14 #last col with abs data 

first_col <- 3 #first col with abs data

abs_data_cols <- raw_table[, first_col:last_col]

abs_data_cols[is.na(abs_data_cols)] <- "" #removes NA from data

## USER-INPUT: location to export csv file
##    example path  : C:\\Users\\Alina\\Desktop\\comma_sep.csv
##    if such file already exists with this path it WILL BE OVER-WRITTEN

csv_path <- readline(prompt = "Enter path to save csv: ")

## Exports the data to a csv file that will be re-read in as a table.

write.csv(abs_data_cols,
          csv_path,
          row.names = FALSE)


## 2. Separate the plates and read plate-view data into lines ==================

data_file_csv <- csv_path

plate_size <- 96 # set for a 96-well plate

num_rows <- 8 #there are 8 rows for a 96-well plate

num_col <- 12 # there are 12 columns for a 96-well plate

## Reads in all plate reads into one formatted tidy data table.  Each plate read 
##    represents absorbances of the plate wells at a different timepoint.
##
##    data_file must be a .csv file containing only plate reads arranged in 
##              an 8x12 size rectangle, with different plate reads of the 
##              file separated by a blank line.  Also exists a single row at the
##              top of the file containing column names

read_plates <- function(data_file) {
  list_timepoint_plates <- make_list_of_timepoints(data_file)
  list_plates_collapsed <- collapse_plates(list_timepoint_plates)
  plates_data_frame <- combine_plates_list(list_plates_collapsed)
  class(plates_data_frame) <- c("tbl_df", "tbl", "data.frame")
  plates_data_frame
}

## Combines the list of data frames corresponding to every timepoint into one
##    data frame.
##
##    list_plates must be a list of plates where each plate is a two column data 
##                frame containing absorbance data.  The first column contains a
##                list of Well ID's and the second contains absorbance values.

combine_plates_list <- function(list_plates) {
  if (length(list_plates) == 1) {
    df_plates <- list_plates[[1]] #if only one timepoint return the input
  } else {
    #combine result into one data frame
    df_plates <- Reduce(function(x,y) merge(x, y, by = "WellIds", all = TRUE),
                        list_plates)
  }
  # only return rows which have value for more than the well ID
  keep <- rowSums(!is.na(df_plates)) > 1
  df_plates <- df_plates[keep, ]
}

## Applies plate_to_column to each plate in the list of plate reads.
##    In this manner collapses plate-view data to columnized data
##
##    list_plate_reads must be a list of plate_view data contained as a data
##                     frame.  Each item in the list is one plate worth of data 

collapse_plates <- function(list_plate_reads) {
  num_timepoints <- length(list_plate_reads)
  map2(list_plate_reads,
       1:num_timepoints,
       plate_to_column)
}

## Converts a single plate (timepoint read) to two-column data frame.  The first
##    column contains the Well ID's and the second the absorbance data.
##
##    plate_read is a group of 8 lines each containing 12 comma-separated values
##               representing each well in a 96-well plate.
##
##    timepoint_num is the number of plate_reads contained in the data set

timepoint_dif <- 2 # there is 2 minutes between each timepoint

plate_to_column <- function(plate_read, timepoint_num) {
  plate <- plate_text_to_data_frame(plate_read)
  column_name <- (timepoint_num - 1) * timepoint_dif
  #convert the plate to a vector
  vect_plate <- unlist(lapply(seq_len(num_rows),
                              function(i) unname(plate[i, ])))
  well_ids <- gen_well_ids(num_rows, num_col)
  df <- data.frame(well_ids, vect_plate, stringsAsFactors = FALSE)
  names(df) <- c("WellIds", column_name)
  return(df)
}

## Converts lines from a plate read to a data frame containing the absorbance
##    data from that plate read stored in an 8x12 data frame
##
##    plate_read is a group of 8 lines each containing 12 comma-separated values
##               representing each well in a 96-well plate.

plate_text_to_data_frame <- function(plate_read) {
  connection <- textConnection(plate_read)
  on.exit(close.connection(connection))
  utils::read.table(connection, sep = ",",
                    na.strings = "", stringsAsFactors = FALSE,
                    comment.char = "", colClasses = "character")
}

## Makes a list of all the plate reads representing each timepoint from a data
##    file
##
##    data_file must be a .csv file containing only plate reads arranged in 
##              an 8x12 size rectangle, with different plate reads of the 
##              file separated by a blank line.  Also exists a single row at the
##              top of the file containing column names

make_list_of_timepoints <- function(data_file) {
  # import data as lines
  data_lines <- readLines(data_file)
  #make a list of data frames (plate reads)
  num_timepoints <- calc_timepoint_num(data_lines)
  list_timepoint_plates <- lapply(1:num_timepoints, 
                                  FUN = function(plate) {
                                    frst_rw <- (plate - 1) * (num_rows + 1) + 2
                                    lst_rw <- frst_rw + (num_rows - 1)
                                    data_lines[frst_rw:lst_rw]
                                  }
  )
  list_timepoint_plates
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

## Calculates the number of timepoints/plate-reads contained inside a given data
##    input
##
##    data_lines is a vector of text lines containing comma separated values

calc_timepoint_num <- function(data_lines) {
  # determines if a number is an integer
  is_integer <- function(x) x %% 1 == 0
  # take one away due to the header line
  quotient <- (length(data_lines) - 2) / (num_rows + 1)
  if (is_integer(quotient)) {
    return(quotient)
  } else {
    quotient <- (length(data_lines) - 1) / (num_rows + 1)
    if (data_file[length(data_lines)] == "" || is_integer(quotient)) {
      return(quotient)
    } else {
      stop(paste0("File length is incorrect.  It must be a multiple of:",
                  "the number of rows in the plate + one blank + a header row,",
                  "the previous + an extra blank row at the end of the file."),
           call. = FALSE)
    }
  }
}

## Stores the number of timepoints/plate-reads in the data file.

num_timepoints <- calc_timepoint_num(readLines(data_file_csv))

## Stores the data_file as the tidied data organized in a data frame.

tidy_data <- read_plates(data_file_csv)


## 3. Organizes the Standard Curve Data ========================================

## Prompts user to input the well ID corresponding to each standard concentration

glucose_2mM <- list(readline(prompt = "Well ID for 2mM Glucose Std - Replicate 1: "),
                    readline(prompt = "Well ID for 2mM Glucose Std - Replicate 2: "),
                    readline(prompt = "Well ID for 2mM Glucose Std - Replicate 3: "))

glucose_1mM <- list(readline(prompt = "Well ID for 1mM Glucose Std - Replicate 1: "),
                    readline(prompt = "Well ID for 1mM Glucose Std - Replicate 2: "),
                    readline(prompt = "Well ID for 1mM Glucose Std - Replicate 3: "))

glucose_0.5mM <- list(readline(prompt = "Well ID for 0.5mM Glucose Std - Replicate 1: "),
                      readline(prompt = "Well ID for 0.5mM Glucose Std - Replicate 2: "),
                      readline(prompt = "Well ID for 0.5mM Glucose Std - Replicate 3: "))

glucose_0.25mM <- list(readline(prompt = "Well ID for 0.25mM Glucose Std - Replicate 1: "),
                        readline(prompt = "Well ID for 0.25mM Glucose Std - Replicate 2: "),
                        readline(prompt = "Well ID for 0.25mM Glucose Std - Replicate 3: "))

## Stores all glucose standard well IDs as a list of lists

standards <- list(glucose_2mM, glucose_1mM, glucose_0.5mM, glucose_0.25mM)

## Fetches the absorbance data for the wells corresponding to the replicates of
##    one glucose standard at time = 0 and stores it as a vector absorbance data
##
##    well_list is a list of well IDs stored as characters(ex: list("A1", "A2"))

get_standard_data <- function(well_list) {
  col = c()
  for(i in 1:length(well_list)) 
  {
    well_id <- well_list[[i]]
    col <- c(col, as.numeric(tidy_data[which(tidy_data$WellIds == well_id), "0"])) 
  }
  col
}

## Makes a data frame (table) for the concentration curve data.  The columns
##    represent the concentrations of glucose and the rows represent different
##    replicates.
##
##    glucose_2 must be a list of Well IDs corresponding to the location of the
##              2 mM glucose standards.  Same for glucose_1 (1 mM), glucose_0.5
##              (0.5 mM), and glucose_0.1 (0.125 mM)

make_standard_table <- function(glucose_2, glucose_1, glucose_0.5, glucose_0.2){
  g2 <- get_standard_data(glucose_2)
  g1 <- get_standard_data(glucose_1)
  g0.5 <- get_standard_data(glucose_0.5)
  g0.25 <- get_standard_data(glucose_0.2)
  std_table <- data.frame(g2, g1, g0.5, g0.25)
  colnames(std_table) <- c(2,1,0.5,0.25)
  row.names(std_table) <- c("Replicate 1", "Replicate 2", "Replicate 3")
  std_table
}

## Stores the table for the glucose standard curve data.

standard_table <- make_standard_table(glucose_2mM, 
                                      glucose_1mM, 
                                      glucose_0.5mM, 
                                      glucose_0.25mM)

## USER-INPUT: location to export csv file for glucose standard data.
##    example path  : C:\\Users\\Alina\\Desktop\\standards.csv
##    if such file already exists with this path it WILL BE OVER-WRITTEN

standard_path <- readline(prompt = "Enter path to save std csv: ")

## Exports the std data if the user wishes to analyze by hand.

write.csv(standard_table,
          standard_path,
          row.names = TRUE)


## 4. Organizes the Enzyme Test Data ===========================================

## Obtains names for all the different Enzymes being tested.  These can be
##    entered with or without double quotes around the names.  But once entered
##    this code must be re-run to change the names.  Any number of names > 1    
##    can be entered.

ask_test_names <- function() {
  print("Enter a short name for each enzyme you are testing on this plate: ")
  return(as.list(scan(what = character())))
}

## Calls ask_test_names() function to store the list of enzymes being tested.

test_names <- ask_test_names()

##  Nests the substrate concentration (in mM) at which an enzyme is being 
##    tested as a list of two-item lists (ex: list(list("NtSI", 2), 
##    list("NtSI", 1))).  Requires input from user to determine the substrate 
##    concentrations tested.
##
##    test_name must be a character describing an enzyme being assayed

nest_subs_conc <- function(test_name){
  print(paste("Enter Substrate Concentration in mM for test: ", 
              test_name, 
              sep=""))
  subs_concs <- as.list(scan(what = numeric()))
  test_ids <- list()
  for (i in 1:length(subs_concs)) {
    test_ids <- append(list(append(test_name, subs_concs[[i]])), test_ids)
  }
  return(test_ids)
}

## Nests the substrate concentration(s) at which each enzyme on the plate is 
##    being tested.
##    
##    list_tests must be a list of characters (representing enzyme names)

ask_subs_conc <- function(list_tests) {
  lst <- list()
  for (i in 1:length(list_tests)){
    test <- list_tests[[i]]
    lst <- append(nest_subs_conc(test), lst)
  }
  return(lst)
}

## Stores the list of test identifiers (list("enzyme_name", substrate_conc))

test_conc <- ask_subs_conc(test_names)

## Adds the list of well IDs where the replicates of each enzyme test is
##    located.  This is stored as a third item in the test identifier; a list
##    of well IDs as characters.  This returns a list of three-item lists.  Ex:
##    (list("enzyme_name", substrate_conc, list("A11", "A12", "A4"))).
##
##    Three well IDs should be entered for each test identifier.  Which are 
##    defined as list("enzyme_name", substrate_conc)
##    list_test_id must be a list of test identifiers 

ask_test_wells <- function(list_test_id) {
  lst <- list()
  for (i in 1:length(list_test_id)){
    test <- list_test_id[[i]]
    test_name <- test[[1]]
    subs_conc <- test[[2]]
    print(
      paste("Enter Well Numbers for test: ", 
            paste(test_name, subs_conc, sep="_"),
            sep=""))
    lst <- append(list(append(test, list(scan(what = character())))), lst)
  }
  return(lst)
}

## Stores the list of test identifiers with the well numbers for each test.

test_conc_wells <- ask_test_wells(test_conc)

## Makes a test ID where the enzyme_name is replaced by the enzyme name appended
##    with the substrate concentration.
##
##    test_descriptor must be a three-item list/test identifier.  Ex: 
##    list("enzyme_name", substrate_conc, list("A11", "A12", "A4"))

make_test_id <- function(test_descriptor){
  enzyme_name <- test_descriptor[[1]]
  subs_conc <- test_descriptor[[2]]
  wells <- test_descriptor[[3]]
  test_id <- paste(enzyme_name, subs_conc, sep = "_")
  test_id_descriptor <- list(test_id, subs_conc, wells)
  return(test_id_descriptor)
}

## Applies make_test_id to every item in a list of test identifiers.
##
##    test_list must be a list of test identifiers (three item list described
##    above)

make_test_ids <- function(test_list){
  lapply(test_list, make_test_id)
}

## Stores the list of test identifiers for the imported data set

test_identifiers <- make_test_ids(test_conc_wells)

## Number of replicates done for each test (combination of enzyme and susbstrate
##    concentration)
replicate_num <- 3

## Extracts all absorbance values for a well ID over all the timepoints
##    contained in the tidied data table as a 'row'(vector).
##
##    well_id must be a valid well ID for a 96-well plate A-H concatenated with 
##    1-12, ex: A12

extract_row <- function(well_id) {
  extracted_row <- tidy_data[which(tidy_data$WellIds==well_id),
                             2:(num_timepoints + 1)]
  row_values <- as.vector(as.numeric(extracted_row))
}

## Creates a new list of a test ID (ex: NtSI_2) concatenated with _'replicate#'
##
##    test_id must be a character

make_ids_rep <- function(test_id) {
  (paste(test_id, seq_len(replicate_num), sep = "_"))
}

read_increments <- 2 # there are 2 minutes in between each read

## Creates a vector containing all the timepoints at which reads where taken in
##    starting a 0 and in minutes
##
##    num_timepoints must be a natural number

timepoint_vector <- function(num_timepoints) {
  times <- c()
  for(i in 1:num_timepoints) 
  {
    times <- c(times, ((i-1)*read_increments))
  }
  times
}

## Creates a data frame where each row is all the absorbance data for a well ID
##    assorted by increasing time
##
##    list_ids must be a list containing valid well IDs for a 96-well plate A-H 
##    concatenated with 1-12, ex: A12

rows_to_df <- function(list_ids) {
  df <- NULL
  for(i in 1:length(list_ids)) 
  {
    df <- rbind(df, extract_row(list_ids[[i]]))
  }
  df
}

## Assigns column and row names to a data frame containing the absorbance data
##    for all three replicates of a set of enzyme and substrate concentration
##    conditions  The columns names are the timepoints and the row names are the
##    test name appended with the replicate number
##
##    df_rows must be a data frame (containing the rows of absorbance values)
##    test_name must be a name for a set of enzyme + substrate conditions

assign_x_y <- function(df_rows, test_name) {
  colnames(df_rows) <- timepoint_vector(num_timepoints)
  rownames(df_rows) <- make_ids_rep(test_name)
  df_rows
}

## Makes a data frame containing all the absorbance data for all replicates of a 
##    set of enzyme + substrate conditions
##
##    test_description must be a test identifier (three-item list described 
##    previously)

make_test_table <- function(test_description){
  test_name <- test_description[[1]]
  well_list <- test_description[[3]]
  df_rows <- rows_to_df(well_list)
  df <- assign_x_y(df_rows, test_name)
}

## Makes a data frame of test data for each test in a list of test identifiers.
##    It returns this within a list where each item is a data frame representing
##    on test identifier
##
##    test_descriptions must be a list of test identifiers, described previously

make_test_tables <- function(test_descriptions){
  test_tables <- lapply(test_descriptions, make_test_table)
  test_tables
}

## Stores a list of data frames representing the list of test identifiers input

test_tables_list <- make_test_tables(test_identifiers)

## Merges a list of data frames each representing an individual test identifier
##    into one large data frame
##
##    test_tables_lst must be a list of data frames (representing a test
##    identifier each)

make_export_table <- function(test_tables_lst){
  export_table <- NULL
  for(i in 1:length(test_tables_lst)) 
  {
    export_table <- rbind(export_table, test_tables_lst[[i]])
  }
  return(export_table)
}

## Makes a table containing the absorbance data for all the enzyme tests
##    contained on the assay plate

test_table <- make_export_table(test_tables_list)

## USER-INPUT: location to export csv file for the enzyme test data.
##    example path  : C:\\Users\\Alina\\Desktop\\test_data.csv
##    if such file already exists with this path it WILL BE OVER-WRITTEN

tests_path <- readline(prompt = "Enter path to save test data: ")

write.csv(test_table,
          tests_path,
          row.names = TRUE)


## 5. Make the Concentration Curve =============================================

## Creates a data frame with two columns, the first containing the mean of each
##      replicate for a glucose standards and the second containing the glucose
##      concentration for that standard
##
##      std_tbl must be a data frame containing the absorbance data for all
##      replicates of a standard organized by col = concentration, row = rep#

make_curve_data <- function(std_tbl) {
  g2_mean <- mean(std_tbl$`2`)
  g1_mean <- mean(std_tbl$`1`)
  g0.5_mean <- mean(std_tbl$`0.5`)
  g0.25mean <- mean(std_tbl$`0.25`)
  vector_means <- c(g2_mean, g1_mean, g0.5_mean, g0.25mean)
  vector_conc <- c(2, 1, 0.5, 0.25)
  std_curve_data <- data.frame(vector_means, vector_conc)
  colnames(std_curve_data) <- c("mean_abs", "glucose_conc")
  std_curve_data
}

## Stores the table of standard curve data

standard_curve_data <- make_curve_data(standard_table)

## Performs simple linear regression on the standard curve data
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

## Stores the coefficients of the linear model (mean_abs = m * glucose_conc + b)

standard_curve_params <- build_linear_model(standard_curve_data)

## Produces a scatter plot with linear regression line for the standard curve
##    data

standard_curve_plot <- ggplot(data = standard_curve_data, 
                              mapping = aes(x = glucose_conc, y = mean_abs)) + 
  geom_point(color = "orangered3", size = 2) +
  geom_smooth(method = lm, se = TRUE, color = "orangered3", size = 1) +
  labs(title = "Standard Curve for Glucose Absorbance at 450nm", #check wavelength
       x = "Concentration of Glucose (mM)",
       y = "Absorbance at 450 nm") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## Saves the standard_curve_plot to a filed named standard_curve.df
##    USER-INPUT: location to export pdf file for the standard curve plot.
##    example path  : C:\\Users\\Alina\\Desktop\\standard_curve.pdf
##    IF SUCH A FILE ALREADY EXISTS IT WILL BE OVERWRITTEN

std_plot_path <- readline(prompt = "Enter path to save standard curve plot: ")
pdf(std_plot_path)
standard_curve_plot
dev.off()

## 6. Analyze Kinetics =========================================================

## Analyzes the test data for a list of tests (combination of enzyme and subs 
##    conc.), to give the R_sqr, slope, and normalized glucose vs. time 
##    information
##
##    test_ids must be a list of test_identifiers, i.e. a list of
##            list(test_id, subs conc, list(wellIds))

calc_kinetics <-function(test_ids){
  list_kin <- lapply(test_ids, calc_test_kinetics)
  gluc_merged <- ((list_kin[[1]])[[2]])[,1]
  lm_params <- data.frame()
  row_names <- c()
  col_names <- c("Time")
  for (i in 1:length(list_kin))
    {
    test_i <- list_kin[[i]]
    test_i_name <- (test_ids[[i]])[[1]]
    gluc_i <- test_i[[2]]
    param_i <- test_i[[3]]
    
    gluc_merged <- cbind(gluc_merged, gluc_i[,2])
    col_names <- c(col_names, test_i_name)
    
    lm_params <- rbind(lm_params, param_i)
    row_names = c(row_names, test_i_name)
  }
  colnames(gluc_merged) <- col_names
  rownames(lm_params) <- row_names
  list(gluc_merged, lm_params)
}

## Analyzes the test data for one test (combination of enzyme and subs conc.),
##    to give the R_sqr, slope, and normalized glucose vs. time information
##
##    test_id must be a test_identifier, i.e. 
##            list(test_id, subs conc, list(wellIds))

calc_test_kinetics <- function(test_id){
  data_table <- make_test_table(test_id)
  test_id_str <- test_id[[1]]
  glc_table <- prep_glucose_data(data_table)
  kin_table <- data.frame(as.numeric(colnames(glc_table)), 
                          as.numeric(glc_table[6,]))
  colnames(kin_table) <- c("time", test_id_str)
  lm <-build_linear_model_test(kin_table)
  return(list(test_id[[1]], kin_table, lm))
}

## Averages the absorbance values for all three replicates of a test a each
##    measured timepoint.  Returns these values in a new row at the bottom of
##    the table
##
##    test_table must be a data frame containing at least one row and one column

average_abs <- function(test_table) {
  averages <- c()
  for(i in 1:num_timepoints) 
  {
    averages <- c(averages, mean(test_table[, i]))
  }
  averages
  averaged_table <- rbind(test_table, averages)
  row.names(averaged_table) <- c(row.names(test_table), "Mean Abs")
  averaged_table
}

## Calculcates the [glucose] (mM) corresponding to the average absorbance values
##    for each measured timepoint according to the linear regression with std 
##    values.  Returns these values in a new row at the bottom of the table
##
##    avg_table must be a data frame containing at four rows and at least one 
##               column the last row must contain the averaged absorbance data
##    lm_params  must be a data frame formatted as in build_linear_model

convert_to_glucose <- function(avg_table, lm_params) {
  m <- lm_params[1,1]
  b <- lm_params[1,2]
  glucose_conc <- c()
  for(i in 1:num_timepoints) 
  {
    avg_abs <- avg_table[4, i]
    glucose_conc <- c(glucose_conc, ((avg_abs - b) / m))
  }
  glucose_conc
  glucose_table <- rbind(avg_table, glucose_conc)
  row.names(glucose_table) <- c(row.names(avg_table), "[Glucose] (mM)")
  glucose_table
}

## Calculcates the normalized [glucose] (mM) corresponding to the [glucose] (mM)
##    for each measured timepoint. Normalized to the measurement at time 0. 
##    Returns these values in a new row at the bottom of the table
##
##    glc_table must be a data frame containing fives rows and at least on 
##              column. The last row must contain the [glucose] data

convert_to_normalized <- function(glc_table) {
  norm <- glc_table[5,1]
  norm_conc <- c()
  for(i in 1:num_timepoints) 
  {
    norm_glc <- (glc_table[5, i] - norm)
    norm_conc <- c(norm_conc, norm_glc)
  }
  norm_conc
  norm_table <- rbind(glc_table, norm_conc)
  row.names(norm_table) <- c(row.names(glc_table), "Normalized [Glucose] (mM)")
  norm_table
} 

## Prepares the kinetic data in terms of [glucose] and time for an enzyme test.
##
##    norm_table must be a 3 row table with at least one column, where each row 
##               represents different replicates for each timepoint

prep_glucose_data <- function(test_table) {
  averaged_table <- average_abs(test_table)
  glucose_table <- convert_to_glucose(averaged_table, standard_curve_params)
  normalized_table <- convert_to_normalized(glucose_table)
}

## Performs simple linear regression on the glucose kinetics data for a single
##    enzyme test (enzyme + concentration parameters)
##
##    glc_test_data must be a data frame prepared as by prep_glucose_data

build_linear_model_test <- function(glc_test_data){
  df <- glc_test_data
  model <- lm(df[,2] ~ df[,1], data = df)
  coefficients <- as.vector(model[[1]])
  intercept <- coefficients[1]
  slope <- coefficients[2]
  r_sqr <- summary(model)$r.squared
  model_coeffs <- data.frame(slope, intercept, r_sqr)
  colnames(model_coeffs) <- c("m", "b", "r_sqr")
  model_coeffs
}

## Stores all the analyzed kinetics data for the plate:

kinetics_data <- calc_kinetics(test_identifiers)

## USER-INPUT: location to export csv file for the enzyme test data normalized.
##    example path  : C:\\Users\\Alina\\Desktop\\test_data_analyzed.csv
##    if such file already exists with this path it WILL BE OVER-WRITTEN
##    the columns of the table include the normalized [glucose] (mM) at 
##    each timepoint

kinetics_path_1 <- readline(prompt = "Enter path to save test data: ")

write.csv(kinetics_data[[1]],
          kinetics_path_1,
          row.names = FALSE)

## USER-INPUT: location to export csv file for the enzyme test data normalized.
##    example path  : C:\\Users\\Alina\\Desktop\\test_data_params.csv
##    if such file already exists with this path it WILL BE OVER-WRITTEN
##    the columns of the table include the velocity and R_sqr for each test

kinetics_path_2 <- readline(prompt = "Enter path to save test data: ")

write.csv(kinetics_data[[2]],
          kinetics_path_2,
          row.names = TRUE)

## Plots the enzyme test kinetics data (glucose concentration over time)

kin_time_series <- as.data.table(kinetics_data[[1]])

enzyme_melted <- melt(data = kin_time_series, id.vars = "Time")

enzyme_plot <- 
  ggplot(data = enzyme_melted, aes(x = Time, y = value, color = variable)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  labs(title = "Kinetics/Velocity Plot for Enzyme Tests", #check title name
      x = "Time (min)",
     y = "Product, [Glucose] (mM)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
  
## Saves the enzyme_kinetics_plot to a filed named kinetics_plot.pdf
##    USER-INPUT: location to export pdf file for enzyme kinetics plot.
##    example path  : C:\\Users\\Alina\\Desktop\\kinetics_plot.pdf
##    IF SUCH A FILE ALREADY EXISTS IT WILL BE OVERWRITTEN

kinetics_plot_path <- readline(prompt = "Enter path to save standard curve plot: ")

pdf(kinetics_plot_path)
enzyme_plot
dev.off()
