## 0.0 Importing the Required Functions ========================================


## NOTE: these files must be in your current working directory, and every time
##    you change data within it, you will have to re-run the source line.

source("platetidyr.R")

source("formatr.R")

source("rtkin.R")


## 1.0 Importing and Tidying the Plate-Shaped Data =============================


raw_data <- "C:/Users/Alina/Desktop/example_data.txt" # USER-INPUT

imported_data <- read_data(raw_data)

## Convert the raw_data data to tidy data:

tidy_data <- tidy_plate_data(imported_data)


## 2.0 Organizing the standard curve data ======================================

## WellID's should be entered manually based on the well map in the accompanying
##    word document.  Replace the characters with valid well IDs.

std1 <- list("A1", "A2", "A3")

std2 <- list("A4", "A5", "A6")

std3 <- list("A7", "A8", "A9")

std4 <- list("A10", "A11", "A12")


## To have the standard table exported as a .csv use export = TRUE, if not, use
##    export = FALSE.  The default is export = TRUE.

standard_table <- make_standard_table(std1, std2, std3, std4, tidy_data, export = TRUE)


## 3.0 Organizing the enzyme assay data ========================================


## Build test identifiers for each different set of tests conducted on the plate.
##    A test identifier is a: list("Enzyme_subs", subs_conc_mM, list(well IDs)))
##    You should form as many or as few as are required.

t1 <- list("NtSI_malt", 1, list("B1", "B2", "B3"))

t2 <- list("NtSI_malt", 0.5, list("B4", "B5", "B6"))

t3 <- list("NtSI_malt", 0.25, list("B7", "B8", "B9"))

t4 <- list("NtSI_malt", 0.125, list("B10", "B11", "B12"))

test_identifiers <- list(t1, t2, t3, t4) # add as many or few as required

## To have the standard table exported as a .csv use export = TRUE, if not, use
##    export = FALSE.  The default is export = TRUE.

test_data <- make_test_data_tables(test_identifiers, export = TRUE)


## 4.0 Analyze the Test Data ===================================================


## To have the standard curve parameters exported as a .csv use export = TRUE, 
##    and to have the plot saved and viewed use plot = TRUE if not, use
##    export = FALSE/plot = FALSE for either or.  The default for both is TRUE.

std_curve <- make_concentration_curve(standard_table, export = TRUE, plot = TRUE)


## To have both the analyzed kinetics data and linear fit parameters as a .csv 
##    use export = TRUE,  and to have the plot saved and viewed use plot = TRUE 
##    if not, use export = FALSE/plot = FALSE for either or.  
##    The default for both is TRUE.

analyzed_kinetics <- calc_test_kinetics(test_identifiers, 
                                        std_curve, 
                                        export = TRUE, 
                                        plot = TRUE)


## 5.0 Modify to Remove Non-Linear Portions and Outliers =======================


## There exists the option to have the program suggest a linear range, however,
##    it is often more accurate to do this by eye.  To do so call the following
##    for each test_name you'd like to have a linear range estimated for:

test_name <- "NtSI_malt_1" #change as required

suggested_points <- estimate_lin_range(analyzed_kinetics[[1]], test_name)

print(suggested_points)  # will allow visualization of which points it suggests


## Next for each test for which a linear range wants to be assigned (that isn't
##    the whole test range), create a list called linear where each item is a 
##    list(test_name, range_description) where a range_description is 
##    list(lin_min, lin_max, outliers).  This DOES NOT have to be done if the 
##    whole timepoint range wants to be taken into account for the analysis.

linear <- list(list("NtSI_malt_1", list(0, 60, list())),
               list("NtSI_malt_0.5", list(10, 60, list())),
               list("NtSI_malt_0.25", list(0, 30, list())),
               list("NtSI_malt_0.125", list(0, 60, list(2, 4))))


## Apply these linear range changes to the data set to generate new progress
##    curves and new data and parameter .csv files.
##    To have both the modified kinetics data and linear fit parameters as a .csv 
##    use export = TRUE,  and to have the plot saved and viewed use plot = TRUE 
##    if not, use export = FALSE/plot = FALSE for either or.  
##    The default for both is TRUE.

modified_kinetics <- calc_mod_kinetics(analyzed_kinetics, 
                                       linear, 
                                       export = TRUE, 
                                       plot = TRUE)
