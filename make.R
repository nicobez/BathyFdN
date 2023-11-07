##################
# BATHYMETRY Fernando de Norona
# nicolas.bez@ird.fr
# May 2023
###################

rm(list=ls())

# Upload required packages
source("R/Libraries.R")
# Load specific functions
source("R/Adds_on_functions.R")
# Read the different polygons used for the analysis
source("R/Polygons.R")

# # Raw data have been cleaned from two problems
# #  1- isolated points (errors in the coordinates)
# #  2- truncation of some observation to various threshold 100, 200 or 300 m 
# # The following code allows traceability and allows to redo it if needed
# # Input : FAROFA123_logSa_1ping / Output : FAROFA123_logSa_1ping_cleaned
# source("R/Cleaning_Database.R") 

# # Geographical coordinates get not enough digits 
# # ==> high numbers of duplicates 
# # ==> regularisation over small pixels to solve the problem (and reduce the size of the data set)
# # Input : FAROFA123_logSa_1ping_cleaned / Output : AROFA123_depth_cleaned_regul
# source("R/pingDatabase.R") 

# Loading the dataframe and building the Db
source("R/Building_myDb.R")

# Computing order 1 and 2 as functions of distAcross and the standardized depth
source("R/Drift_and_residuals.R")

# Grid and mesh definition + unfolding grid coordinates
source("R/Grid_definition_and_meshing.R")

# Building the A matrix for projecting mesh points towards data points 
A = cs_toTL(ProjMatrix_create(myDb,myMesh)$getAproj())

### SPDE of the standardized depth
compute <- F # upload results ; must be turned 'T' to compute the likelihood 
source("R/mLogL_Stationary_depthStd.R")
source("R/SPDE_Stationary_depthStd.R")
source("R/mLogL_Non_Stationary_depthStd.R")
source("R/SPDE_Non_Stationary_depthStd.R")

source("R/Xvalid_Farofa3_depthStd.R")



### SPDE of th centered depth
myDb$setLocator("depthC",ELoc_Z(),0)
compute <- F # upload results ; must be turned 'T' to compute the likelihood 
source("R/mLogL.R")

source("R/SPDE_Stationary.R")
source("R/SPDE_Non_Stationary.R")

source("R/Basic_graphics.R")


