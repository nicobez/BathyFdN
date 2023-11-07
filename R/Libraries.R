##############
# Libraries
##############
# CRAN packages
CRAN_package_list <- c("Matrix",
                       "wesanderson",
                       'ggplot2',
                       'ggpubr')


not_in_list <- CRAN_package_list[!(CRAN_package_list %in% installed.packages())]
lapply(not_in_list,install.packages,dependencies=TRUE)

lapply(CRAN_package_list, require, character.only = TRUE)
       
# NON-CRAN packages
if("gstlearn" %in% installed.packages()) library(gstlearn)
if(!("gstlearn" %in% installed.packages())) cat("gstlearn must be loaded and installed from GitHub \n")

myPalette = wes_palette("Zissou1", 16, type = "continuous")

rm(CRAN_package_list,not_in_list)