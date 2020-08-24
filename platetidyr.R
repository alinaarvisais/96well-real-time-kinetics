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


## 1.0 Importing the file with the right formatting ============================


## Creates a .csv file with the proper formatting for analysis, if it takes
##   too long to run, can later omit this and do it manually with Excel.
##    Returns the path to the created .csv file path
##
##    raw_data is the file path for the raw .txt file exported from the
##    spectrophotometer

read_data <- function(raw_data){
  
  raw_lines <- readLines(raw_data)
  
  ## The code converts all tabs in a document to a space character so that there 
  ## is a single delimiter across the whole document
  
  lines_tab_delim <- gsub("\\t$", "", raw_lines)
  
  raw_table <- read.delim(text = lines_tab_delim,
                          skip = 3, 
                          header = FALSE)
  
  ## Identifies and selects the columns in the raw data table that contains
  ## absorbance data.
  
  last_col <- 14 #last col with abs data 
  
  first_col <- 3 #first col with abs data
  
  abs_data_cols <- raw_table[, first_col:last_col]
  
  abs_data_cols[is.na(abs_data_cols)] <- "" #removes NA from data
  
  ## Exports the data to a csv file that will be re-read in as a table.
  
  write.csv(abs_data_cols,
            "comma_sep.csv",
            row.names = FALSE)
  
  "comma_sep.csv"
}


## 2.0 Separate the plates and read plate-view data into data lines ============


## Wrapper for read_plates containing plate dimensions and data file path.  Runs
##    read_plates on the data contained in a comma-separated plate-view file and 
##    returns the outcome.
##
##    raw_data_csv must be a valide file pathname containing the data to be 
##    converted to a tidy-data format

tidy_plate_data <- function(raw_data_csv){
  
  num_rows <- 8 #there are 8 rows for a 96-well plate
  
  assign("num_rows", num_rows, envir = .GlobalEnv)
  
  num_col <- 12 # there are 12 columns for a 96-well plate
  
  assign("num_col", num_col, envir = .GlobalEnv)
  
  num_timepoints <- calc_timepoint_num(readLines(raw_data_csv))
  
  assign("num_timepoints", num_timepoints, envir = .GlobalEnv)
  
  read_plates(raw_data_csv)
}


## Reads in all plate reads into one formatted tidy data table.  Each plate read 
##    represents absorbances of the plate wells at a different timepoint.
##    Returns the tidied data table.
##
##    data_file must be a .csv file containing only plate reads arranged in 
##              an 8x12 size rectangle, with different plate reads of the 
##              file separated by a blank line.  Also exists a single row at the
##              top of the file containing column names

read_plates <- function(data_file){
  
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
##    In this manner collapses plate-view data to columnized data.  Returns a
##    list of the collapsed plates.
##
##    list_plate_reads must be a list of plate-view data contained as a data
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
##    timepoint_read is a group of 8 lines each containing 12 comma-separated 
##    values representing each well in a 96-well plate.
## 
##    timepoint_num is the number corresponding to which timepoint_read it
##    is (indexed from 1:number of abs reads)

plate_to_column <- function(timepoint_read, timepoint_num) {
  
  timepoint_dif <- 2 # there are 2 minutes between each timepoint
  
  assign("timepoint_dif", timepoint_dif, envir = .GlobalEnv) 
  
  plate <- plate_text_to_data_frame(timepoint_read)
  
  column_name <- (timepoint_num - 1) * timepoint_dif
  
  #convert the plate to a vector
  vect_plate <- unlist(lapply(seq_len(num_rows),
                              function(i) unname(plate[i, ])))
  
  well_ids <- gen_well_ids(num_rows, num_col)
  
  df <- data.frame(well_ids, vect_plate, stringsAsFactors = FALSE)
  
  names(df) <- c("WellIds", column_name)
  
  df
}


## Converts lines from a plate read to a data frame containing the absorbance
##    data from that plate read stored in an 8x12 data frame.  Returns the
##    8x12 data frame.
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
##    file.  Returns the list of timepoints.
##
##    data_file must be a .csv file containing only plate reads arranged in 
##              an 8x12 size rectangle, with different plate reads of the 
##              file separated by a blank line.  Also exists a single row at the
##              top of the file containing column names

make_list_of_timepoints <- function(data_file) {
  
  # import data as lines
  data_lines <- readLines(data_file)
  
  #make a list of data frames (plate reads)
  
  list_timepoint_plates <- lapply(1:num_timepoints, 
                                  FUN = function(plate) {
                                    frst_rw <- (plate - 1) * (num_rows + 1) + 2
                                    lst_rw <- frst_rw + (num_rows - 1)
                                    data_lines[frst_rw:lst_rw]
                                  })
  list_timepoint_plates
}


## Generates a vector of well IDs for a well plate.  In the vector, the
##    the rows are identified alphabetically and the columns are identified
##    numerically.  Returns this vector
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
                  "or the previous + an extra blank row at the end of the file."),
           call. = FALSE)
    }
  }
}
