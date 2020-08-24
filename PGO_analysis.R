## 0.0 Importing the Required Functions ========================================


## NOTE: these files must be in your current working directory, and every time
##    you change data within it, you will have to re-run the source line.

source("platetidyr.R")

source("formatr.R")

source("rtkin.R")

source("linearrange.R")


## 1.0 Importing and Tidying the Plate-Shaped Data =============================


raw_data <- "path_here" # USER-INPUT

imported_data <- read_data(raw_data)

## Convert the raw_data data to tidy data:

tidy_data <- tidy_plate_data(imported_data) # must be named tidy_data


## 2.0 Organizing the standard curve data ======================================


## WellID's should be entered manually based on the well map in the accompanying
##    word document.  Replace the characters with valid well IDs.

std1 <- list("WellID Glucose 2mM Rep1", "Rep2", "Rep3") # 2mM Glucose

std2 <- list("WellID Glucose 1mM Rep1", "Rep2", "Rep3") # 1 mM Glucose

std3 <- list("WellID Glucose 0.5mM Rep1", "Rep2", "Rep3") # 0.5 mM Glucose

std4 <- list("WellID Glucose 0.25mM Rep1", "Rep2", "Rep3") # 0.25 mM Glucose


## To have the standard table exported as a .csv use export = TRUE, if not, use
##    export = FALSE.  The default is export = TRUE.

standard_table <- make_standard_table(std1, std2, std3, std4, tidy_data, export = TRUE)


## 3.0 Organizing the enzyme assay data ========================================


## Build test identifiers for each different set of tests conducted on the plate.
##    A test identifier is a: list("Enzyme_subs", subs_conc_mM, list(well IDs))
##    You should form as many or as few as are required.

t1 <- list("Enzyme1_subs1", 2, list("WellId", "WellId", "WellId"))

t2 <- list("Enzyme1_subs1", 1, list("WellId", "WellId", "WellId"))

t3 <- list("Enzyme1_subs2", 2, list("WellId", "WellId", "WellId"))

test_identifiers <- list(t1, t2, t3) # add as many or few as required

## To have the standard table exported as a .csv use export = TRUE, if not, use
##    export = FALSE.  The default is export = TRUE.

test_data <- make_test_data_tables(test_identifiers, export = TRUE)


## 4.0 Assign the Blank Wells (Enzyme + Substrate Blanks) ======================

## Build blank identifiers for each different blank conducted on the plate.
##    A blank identifier is a: list("blank_type", "name", list(well IDs))
##    You should form as many or as few as are required.
##    "blank_type" must be one of "substrate" or "enzyme"
##    name should be the substrate name or enzyme it relates to, match the
##    formatting used in the test_identifiers

b1 <- list("substrate", "subs1", list("WellId", "WellId", "WellId"))

b2 <- list("enzyme", "Enzyme1", list("WellId", "WellId", "WellId"))

blanks <- list(b1, b2) # add as many or as few as required

## Analyze the blanks to see if their absorbance value is greater than a given
##    threshold (the default value is 0.1).  If the blank values are lower than
##    threshold, no actions are taken.  If it is an enzyme blank and the value
##    is above the threshold, it applies a blank offset (avg abs of blanks at that
##    same time for each timepoint) to the test_identifiers that contain the 
##    blank enzyme name.  If it is a substrate blank and the value is above the
##    threshold, it removes the test_identifiers that contains its name
##    from the list.  It will print a message saying which actions have been
##    taken.

test_identifiers_blanked <- apply_blanks(test_identifiers, blanks, 0.1) # 0.1 = threshold


## 5.0 Analyze the Test Data ===================================================


## To have the standard curve parameters exported as a .csv use export = TRUE,
##    and to have the plot saved and viewed use plot = TRUE if not, use
##    export = FALSE/plot = FALSE for either or.  The default for both is TRUE.

std_curve <- make_concentration_curve(standard_table, export = TRUE, plot = TRUE)


## To have both the analyzed kinetics data and linear fit parameters as a .csv
##    use export = TRUE,  and to have the plot saved and viewed use plot = TRUE
##    if not, use export = FALSE/plot = FALSE for either or.
##    The default for both is TRUE.

analyzed_kinetics <- calc_test_kinetics(test_identifiers_blanked, # or test_identifiers if blanks are N/A
                                        std_curve,
                                        export = TRUE,
                                        plot = TRUE)


## 6.0 Modify to Remove Non-Linear Portions and Outliers =======================


## There exists the option to have the program suggest a linear range, however,
##    it is often more accurate to do this by eye.  To do so call the following
##    for each test_name you'd like to have a linear range estimated for:

test_name <- "Enzyme1_subs1_1" #change as required

suggested_points <- estimate_lin_range(analyzed_kinetics[[1]], test_name)

print(suggested_points)  # will allow visualization of which points it suggests

## Next for each test for which a linear range wants to be assigned (that isn't
##    the whole test range), create a list called linear where each item is a 
##    linear_identifiers where a linear_identifier is a
##    list("test_name", list(lin_min, lin_max, outliers)).  
##    This DOES NOT have to be done if the whole timepoint range wants to be 
##    taken into account for the analysis.

lin1 <- list("Enzyme1_subs1_2", list(0, 60, list()))

lin2 <- list("Enzyme1_subs1_1", list(10, 60, list()))

lin3 <- list("Enzyme1_subs2_2", list(0, 30, list()))

linear <- list(lin1, lin2, lin3)

## Apply these linear range changes to the data set to generate new progress
##    curves and new data and parameter .csv files.
##    To have both the modified kinetics data and linear fit parameters as a .csv 
##    use export = TRUE,  and to have the plot saved and viewed use plot = TRUE 
##    if not, use export = FALSE/plot = FALSE for either or.  
##    The default for both is TRUE.

linear_kinetics <- calc_mod_kinetics(analyzed_kinetics, 
                                     linear, 
                                     export = TRUE, 
                                     plot = TRUE)
