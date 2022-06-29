####################################################################################################################
## Author: GREG CHISM
## Date: APR 2022
## email: gchism@arizona.edu
## Project: Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space
## Title: Site fidelity zone assignment, calculations, and distance functions
####################################################################################################################

# This code is to replicate the analyses and figures for the following:
# Zone assignment for site fidelity
# Fidelity zone calculations
# Fidelity zone distance functions (to nest entrance and brood center)

####################################################################################################################
# IMPORT ALL NECESSARY DATASETS 
# This code imports all necessary data sets for the script below
####################################################################################################################

# COLOR REFRENCE DATA SET
ColorRefFull <- read.csv("ColorRefFull.csv")

# NEST AREA REFERENCES
NestAreaFull <- read.csv("NestAreaFull.csv")

# SCALING REFERENCE DATA SET FOR CIRCLE NESTS (BASICALLY COLONY SIZE)
ScalingCircleSFZ <- read.csv("ScalingCircleSFZ.csv")

####################################################################################################################
# CORRECTING COLOR COORDINATES
# The following function corrects viable worker color identification
# This means that there is a missing color and only one possible combination
# The function uses the worker color combinations for each colony as a reference
# First the worker data set is set so that missing colors are "X"
# Then the colors are filled by looking at all possible combinations of one missing color mark for each colony
# The function then checks for any rows that didn't have a corrected combination, resulting most likely from the color being misinterpreted
####################################################################################################################

ColorCoords <- function(data.table){
  # Fixing missing marks to only by "X" 
  ColorSetting <- data.table %>%
    # Separate the AntID column into 4 columns, one for each mark 
    separate(AntID, c("Head", "Thorax", "Abd1", "Abd2"), sep = ",", remove = FALSE) %>% 
    # Changing "?" values to "X" in each column
    mutate(Head = ifelse(Head == "?", "X", Head),
           Thorax = ifelse(Thorax == "?", "X", Thorax),
           Abd1 = ifelse(Abd1 == "?", "X", Abd1),
           Abd2 = ifelse(Abd2 == "?", "X", Abd2)) %>%
    # Uniting the 4 columns to the original AntID column
    unite(AntID, c("Head", "Thorax", "Abd1", "Abd2"), sep = ",", remove = TRUE)
  # Creating a reference data set of all possible combinations for one missing color mark
  ColorRefFull <- ColorRefFull%>%
    # Unite the color reference columns
    unite(AntIDRef, c("Head", "Thorax", "Abd1", "Abd2"), sep = ",", remove = FALSE)
  # Creating individual data sets to create one column full of "X"
  X1 <- ColorRefFull
  X2 <- ColorRefFull
  X3 <- ColorRefFull
  X4 <- ColorRefFull
  # Filling one of the reference columns
  X1[3] <- "X"
  X2[4] <- "X"
  X3[5] <- "X"
  X4[6] <- "X"
  # Joining the reference data sets
  FullColorCoordRef <- full_join(X1, X2) %>%
    full_join(X3) %>%
    full_join(X4) %>%
    # Uniting the reference columns
    # This new column is the same name as in the real data set
    # So any column with the same missing value will now have the true color combination in the AntIDRef column
    unite(AntID, c("Head", "Thorax", "Abd1", "Abd2"), sep = ",", remove = TRUE)
  # Creating the final data set with the original combination and true color combination
  SFZDataFull <<- left_join(ColorSetting, FullColorCoordRef) %>%
    # Separating the true AntID column and removing rows with more than one "X"
    # This is done by creating a column that assigns a 1 when two columns are "X"
    separate(AntID, c("Head", "Thorax", "Abd1", "Abd2"), sep = ",", remove = FALSE) %>%
    mutate(XCount = ifelse(Head == "X" & Thorax == "X" | 
                             Abd1 == "X" & Abd2 == "X" | 
                             Head == "X" & Abd1 == "X" |
                             Head == "X" & Abd2 == "X" | 
                             Thorax == "X" & Abd1 == "X"|
                             Thorax == "X" & Abd2 == "X", 1, 0),
           # Changing NA values in AntIDRef to 0, which occurs when a AntID doesn't have just one missing color
           # The objective is to remove both rows with more than two missing colors and keep all fully marked workers 
           AntIDRef = ifelse(is.na(AntIDRef), 0, AntIDRef),
           AntID = ifelse(AntIDRef == "0", AntID, AntIDRef)) %>%
    # Filtering out rows with more than one "X", or missing color
    filter(XCount != 1) %>%
    # Selecting the desired columns 
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, AntID) %>%
    distinct()
}

# Run color correction function on the combined workers data set
ColorCoords(FullDataCoordWorkersRD1_RD2)

####################################################################################################################
# SITE FIDELITY ASSIGNMENT FUNCTION
# The code below bins x and y coordinate colony data into twenty-four even area nest sections
# Note that the code is specifically used below for the low density treatment colonies (11-20)
# This is because the high density treatment colonies (1-10) were done through excel but can be run through this script as well
# To do this, a reference data set of bin coordinates is used and coordinates are run through a series of conditional statements
# Where each conditional statement checks whether the coordinate is in one of twenty-four bins sequentially
# The reference data set for tube nest bins was used here as well, since the bins were divided by three to achieve the desired nest sections here
# The circle nest shape has a different set of reference coordinates
####################################################################################################################

# Creating the reference zones for the circle nest shape
BinCoord24Circ <- NestAreaFull %>%
  mutate(Radius = Diameter / 2) %>% # Radius of each circle
  left_join(ScalingCircleSFZ) %>% # Joining the scaling reference data set
  mutate(AddValue = Radius * Scaling) %>% # Value to add to the center of the circle nest for zones, where instead of manually assigning reference coordinates for every bin, a value is added from the nest center (3.75) to find it instead
  distinct() # Remove duplicates

# COLONY 1
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "1") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "1") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "1")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony1BinnedTubeSFZ <- data_table %>%
    filter(Colony == "1" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony1BinnedCircleSFZ <- data_table %>%
    filter(Colony == "1" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony1BinnedSFZ <<- full_join(Colony1BinnedTubeSFZ, Colony1BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 2
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "2") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "2") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "2")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony2BinnedTubeSFZ <- data_table %>%
    filter(Colony == "2" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony2BinnedCircleSFZ <- data_table %>%
    filter(Colony == "2" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony2BinnedSFZ <<- full_join(Colony2BinnedTubeSFZ, Colony2BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 3
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "3") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "3") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "3")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony3BinnedTubeSFZ <- data_table %>%
    filter(Colony == "3" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony3BinnedCircleSFZ <- data_table %>%
    filter(Colony == "3" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony3BinnedSFZ <<- full_join(Colony3BinnedTubeSFZ, Colony3BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 4
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "4") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "4") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "4")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony4BinnedTubeSFZ <- data_table %>%
    filter(Colony == "4" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony4BinnedCircleSFZ <- data_table %>%
    filter(Colony == "4" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony4BinnedSFZ <<- full_join(Colony4BinnedTubeSFZ, Colony4BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 5
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "5") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "5") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "5")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony5BinnedTubeSFZ <- data_table %>%
    filter(Colony == "5" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony5BinnedCircleSFZ <- data_table %>%
    filter(Colony == "5" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony5BinnedSFZ <<- full_join(Colony5BinnedTubeSFZ, Colony5BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 6
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "6") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "6") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "6")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony6BinnedTubeSFZ <- data_table %>%
    filter(Colony == "6" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony6BinnedCircleSFZ <- data_table %>%
    filter(Colony == "6" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony6BinnedSFZ <<- full_join(Colony6BinnedTubeSFZ, Colony6BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 7
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "7") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "7") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "7")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony7BinnedTubeSFZ <- data_table %>%
    filter(Colony == "7" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony7BinnedCircleSFZ <- data_table %>%
    filter(Colony == "7" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony7BinnedSFZ <<- full_join(Colony7BinnedTubeSFZ, Colony7BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 8
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "8") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "8") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "8")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony8BinnedTubeSFZ <- data_table %>%
    filter(Colony == "8" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony8BinnedCircleSFZ <- data_table %>%
    filter(Colony == "8" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony8BinnedSFZ <<- full_join(Colony8BinnedTubeSFZ, Colony8BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 9
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "9") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "9") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "9")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony9BinnedTubeSFZ <- data_table %>%
    filter(Colony == "9" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony9BinnedCircleSFZ <- data_table %>%
    filter(Colony == "9" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony9BinnedSFZ <<- full_join(Colony9BinnedTubeSFZ, Colony9BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 10
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "10") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "10") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "10")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony10BinnedTubeSFZ <- data_table %>%
    filter(Colony == "10" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony10BinnedCircleSFZ <- data_table %>%
    filter(Colony == "10" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony10BinnedSFZ <<- full_join(Colony10BinnedTubeSFZ, Colony10BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 11
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "11") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "11") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "11")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony11BinnedTubeSFZ <- data_table %>%
    filter(Colony == "11" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony11BinnedCircleSFZ <- data_table %>%
    filter(Colony == "11" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony11BinnedSFZ <<- full_join(Colony11BinnedTubeSFZ, Colony11BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 12
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "12") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "12") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "12")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony12BinnedTubeSFZ <- data_table %>%
    filter(Colony == "12" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony12BinnedCircleSFZ <- data_table %>%
    filter(Colony == "12" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony12BinnedSFZ <<- full_join(Colony12BinnedTubeSFZ, Colony12BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 13
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "13") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "13") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "13")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony13RD2BinnedTubeSFZ <- data_table %>%
    filter(Colony == "13" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                 ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                 ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                               # Zone 17
                                                                                                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                         ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                         ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                       # Zone 18
                                                                                                                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                 ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                 ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                               # Zone 19
                                                                                                                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                         ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                         ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                       # Zone 20
                                                                                                                                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                 ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                 ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                               # Zone 21
                                                                                                                                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                         ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                         ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                       # Zone 22
                                                                                                                                                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                 ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                 ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                               # Zone 23
                                                                                                                                                                                               if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                         ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                       # Zone 24
                                                                                                                                                                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                 ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                 ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                       )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony13RD2BinnedCircleSFZ <- data_table %>%
  filter(Colony == "13" & Nest == "Circle") %>%
  mutate(Zone =
           # Zone 1
           if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                     ScaledY <= BinCoordAssign$ScaledY[10], 1,
                   # Zone 2
                   if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                             ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                             ScaledY <= BinCoordAssign$ScaledY[10], 2,
                           # Zone 3
                           if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                     ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                   # Zone 4
                                   if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                             ScaledY >= BinCoordAssign$ScaledY[10] & 
                                             ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                           # Zone 5
                                           if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                     ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                     ScaledY >= BinCoordAssign$ScaledY[10] &
                                                     ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                   # Zone 6
                                                   if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                             ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                             ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                           # Zone 7
                                                           if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                     ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                     ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                   # Zone 8
                                                                   if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                             ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                             ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                             ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                           # Zone 9
                                                                           if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                     ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                     ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                   # Zone 10
                                                                                   if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                             ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                             ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                           # Zone 11
                                                                                           if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                     ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                     ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                     ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                   # Zone 12
                                                                                                   if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                             ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                             ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                           # Zone 13
                                                                                                           if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                     ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                   # Zone 14
                                                                                                                   if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                             ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                             ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                           # Zone 15
                                                                                                                           if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                     ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                   # Zone 16
                                                                                                                                   if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                             ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                           # Zone 17
                                                                                                                                           if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                     ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                     ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                   # Zone 18
                                                                                                                                                   if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                             ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                           # Zone 19
                                                                                                                                                           if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                     ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                   # Zone 20
                                                                                                                                                                   if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                             ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                             ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                           # Zone 21
                                                                                                                                                                           if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                     ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                   # Zone 22
                                                                                                                                                                                   if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                           # Zone 23
                                                                                                                                                                                           if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                     ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                     ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                   # Zone 24
                                                                                                                                                                                                   if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                             ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                           )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony13BinnedSFZ <<- full_join(Colony13RD2BinnedTubeSFZ, Colony13RD2BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 14
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "14") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "14") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "14")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony14BinnedTubeSFZ <- data_table %>%
    filter(Colony == "14" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony14BinnedCircleSFZ <- data_table %>%
    filter(Colony == "14" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony14BinnedSFZ <<- full_join(Colony14BinnedTubeSFZ, Colony14BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 15
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "15") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "15") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "15")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony15BinnedTubeSFZ <- data_table %>%
    filter(Colony == "15" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony15BinnedCircleSFZ <- data_table %>%
    filter(Colony == "15" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony15BinnedSFZ <<- full_join(Colony15BinnedTubeSFZ, Colony15BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 16
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "16") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "16") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "16")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony16BinnedTubeSFZ <- data_table %>%
    filter(Colony == "16" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony16BinnedCircleSFZ <- data_table %>%
    filter(Colony == "16" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony16BinnedSFZ <<- full_join(Colony16BinnedTubeSFZ, Colony16BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 17
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "17") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "17") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "17")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony17BinnedTubeSFZ <- data_table %>%
    filter(Colony == "17" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony17BinnedCircleSFZ <- data_table %>%
    filter(Colony == "17" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony17BinnedSFZ <<- full_join(Colony17BinnedTubeSFZ, Colony17BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 18
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "18") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "18") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "18")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony18BinnedTubeSFZ <- data_table %>%
    filter(Colony == "18" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony18BinnedCircleSFZ <- data_table %>%
    filter(Colony == "18" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony18BinnedSFZ <<- full_join(Colony18BinnedTubeSFZ, Colony18BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 19
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "19") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "19") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "19")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony19BinnedTubeSFZ <- data_table %>%
    filter(Colony == "19" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony19BinnedCircleSFZ <- data_table %>%
    filter(Colony == "19" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony19BinnedSFZ <<- full_join(Colony19BinnedTubeSFZ, Colony19BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# COLONY 20
CoordBinnedSFZ <- function(data_table){
  # Filtering the desired colony from the reference data sets
  SFZBinCoordAssign <- BinCoordFull %>%
    filter(Colony == "20") 
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "20") 
  BinCoord24Circ <- BinCoord24Circ %>%
    filter(Colony == "20")
  
  # Filtering the desired colony from the worker tube nest data sets
  # Note that the site fidelity Bins will be referred to as "Zones" from here on
  # TUBE NEST ZONES
  Colony20BinnedTubeSFZ <- data_table %>%
    filter(Colony == "20" & Nest == "Tube") %>%
    # The following columns create reference points from the eight even area bin reference data set
    # The goal is to produce a segment length that is along one of the eight bins, each encompassing three of the twenty-four zones
    # Since these bins were cut into thirds, it is more accurate to simply use these segments than to create a whole new set of reference coordinates
    mutate(Bin1.3 = SFZBinCoordAssign$ScaledY[2], #Length from Zone 1 to 3
           Bin4.6 = abs(SFZBinCoordAssign$ScaledX[2] - SFZBinCoordAssign$ScaledX[3]), #Length from Zone 4 to 6
           Bin7.9 = abs(SFZBinCoordAssign$ScaledY[4] - SFZBinCoordAssign$ScaledY[3]), #Length from Zone 7 to 9
           Bin10.12 = abs(SFZBinCoordAssign$ScaledX[5] - SFZBinCoordAssign$ScaledX[4]), #Length from Zone 10 to 12
           Bin13.15 = abs(SFZBinCoordAssign$ScaledX[6] - SFZBinCoordAssign$ScaledX[5]), #Length from Zone 13 to 15
           Bin16.18 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[6]), #Length from Zone 16 to 18
           Bin19.21 = abs(SFZBinCoordAssign$ScaledX[7] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 19 to 21
           Bin22.24 = abs(SFZBinCoordAssign$ScaledX[8] - SFZBinCoordAssign$ScaledX[9]), #Length from Zone 22 to 24
           # Binning the worker coordinates into a Zone
           Zone =
             # Zone 1
             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                       ScaledY <= (1/3) * (Bin1.3), 1, 
                     # Zone 2
                     if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                               ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                               ScaledY >= (1/3) * (Bin1.3) &
                               ScaledY <= (2/3) * (Bin1.3), 2,
                             # Zone 3
                             if_else(ScaledX >= SFZBinCoordAssign$ScaledX[2] & 
                                       ScaledX <= SFZBinCoordAssign$ScaledX[1] &
                                       ScaledY >= (2/3) * (Bin1.3) &
                                       ScaledY <= (Bin1.3), 3, 
                                     # Zone 4
                                     if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) & 
                                               ScaledX <= BinCoordAssign$ScaledX[3] + (Bin4.6) &
                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                               ScaledY <= BinCoordAssign$ScaledY[2], 4,
                                             # Zone 5
                                             if_else(ScaledX >= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) & 
                                                       ScaledX <= BinCoordAssign$ScaledX[3] + ((2/3) * (Bin4.6)) &
                                                       ScaledY >= BinCoordAssign$ScaledY[3] &
                                                       ScaledY <= BinCoordAssign$ScaledY[2], 5,
                                                     # Zone 6
                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[3] + ((1/3) * (Bin4.6)) &
                                                               ScaledX >= BinCoordAssign$ScaledX[3] &
                                                               ScaledY >= BinCoordAssign$ScaledY[3] &
                                                               ScaledY <= BinCoordAssign$ScaledY[2], 6,
                                                             # Zone 7
                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)), 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                               ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[3] + ((1/3) * (Bin7.9)) &
                                                                               ScaledY <= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)), 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                                                                       ScaledX <= BinCoordAssign$ScaledX[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[3] + ((2/3) * (Bin7.9)) & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[3] + (Bin7.9), 9,
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + ((1/3) * Bin10.12) & 
                                                                                                       ScaledX <= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) & 
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX <= BinCoordAssign$ScaledX[4] + Bin10.12 & 
                                                                                                               ScaledX >= BinCoordAssign$ScaledX[4] + ((2/3) * Bin10.12) &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 12,
                                                                                                             # Zone13
                                                                                                             if_else(ScaledX >= BinCoordAssign$ScaledX[4] + Bin10.12 &
                                                                                                                       ScaledX <= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[5], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= BinCoordAssign$ScaledX[5] + ((1/3) * Bin13.15) & 
                                                                                                                               ScaledX <= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) & 
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[5], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else( ScaledX >= BinCoordAssign$ScaledX[5] + ((2/3) * Bin13.15) &
                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 15,
                                                                                                                                      # Zone 16
                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[5] + Bin13.15 & 
                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 16,
                                                                                                                                              # Zone 17
                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((1/3) * Bin16.18) & 
                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) & 
                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[5], 17,
                                                                                                                                                      # Zone 18
                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[6] + ((2/3) * Bin16.18) &
                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[6] + Bin16.18 &
                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[4] &
                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[5], 18,
                                                                                                                                                              # Zone 19
                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] & 
                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21), 19,
                                                                                                                                                                      # Zone 20
                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[7] - (1/3) * (Bin19.21) & 
                                                                                                                                                                                ScaledY >= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21), 20,
                                                                                                                                                                              # Zone 21
                                                                                                                                                                              if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                                                                                                                                        ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[7] - (2/3) * (Bin19.21) & 
                                                                                                                                                                                        ScaledY >= BinCoordAssign$ScaledY[7] - (Bin19.21), 21,
                                                                                                                                                                                      # Zone 22
                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (Bin22.24) &
                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 22,
                                                                                                                                                                                              # Zone 23
                                                                                                                                                                                              if_else(ScaledX <= BinCoordAssign$ScaledX[9] + ((2/3) * (Bin22.24)) & 
                                                                                                                                                                                                        ScaledX >= BinCoordAssign$ScaledX[9] + ((1/3) * (Bin22.24)) &
                                                                                                                                                                                                        ScaledY <= BinCoordAssign$ScaledY[9], 23,
                                                                                                                                                                                                      # Zone 24
                                                                                                                                                                                                      if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                                                                                                                                                                ScaledX <= BinCoordAssign$ScaledX[9] + (1/3) * (Bin22.24) &
                                                                                                                                                                                                                ScaledY <= BinCoordAssign$ScaledY[9], 24, 0
                                                                                                                                                                                                      )))))))))))))))))))))))))
  # CIRCLE NEST ZONES
  Colony20BinnedCircleSFZ <- data_table %>%
    filter(Colony == "20" & Nest == "Circle") %>%
    mutate(Zone =
             # Zone 1
             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                       ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Zone 2
                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                               ScaledY <= BinCoordAssign$ScaledY[10], 2,
                             # Zone 3
                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                       ScaledY <= BinCoordAssign$ScaledY[10], 3,
                                     # Zone 4
                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                               ScaledY <= BinCoordAssign$ScaledY[11], 4,
                                             # Zone 5
                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                       ScaledY >= BinCoordAssign$ScaledY[10] &
                                                       ScaledY <= BinCoordAssign$ScaledY[11], 5,
                                                     # Zone 6
                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                               ScaledY >= BinCoordAssign$ScaledY[10] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[11], 6,
                                                             # Zone 7
                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 7,
                                                                     # Zone 8
                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                               ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                               ScaledY <= BinCoordAssign$ScaledY[12], 8,
                                                                             # Zone 9
                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                       ScaledY >= BinCoordAssign$ScaledY[11] & 
                                                                                       ScaledY <= BinCoordAssign$ScaledY[12], 9, 
                                                                                     # Zone 10
                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 10,
                                                                                             # Zone 11
                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                       ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                       ScaledY <= BinCoordAssign$ScaledY[13], 11,
                                                                                                     # Zone 12
                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                               ScaledY >= BinCoordAssign$ScaledY[12] & 
                                                                                                               ScaledY <= BinCoordAssign$ScaledY[13], 12,
                                                                                                             # Zone 13
                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 13,
                                                                                                                     # Zone 14
                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[14], 14,
                                                                                                                             # Zone 15
                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[4] &
                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[14], 15,
                                                                                                                                     # Zone 16
                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 16,
                                                                                                                                             # Zone 17
                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[15], 17,
                                                                                                                                                     # Zone 18
                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[3] &
                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[15], 18,
                                                                                                                                                             # Zone 19
                                                                                                                                                             if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 19,
                                                                                                                                                                     # Zone 20
                                                                                                                                                                     if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledX <= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                               ScaledY <= BinCoordAssign$ScaledY[16], 20,
                                                                                                                                                                             # Zone 21
                                                                                                                                                                             if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[2] &
                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                                                                                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 21,
                                                                                                                                                                                     # Zone 22
                                                                                                                                                                                     if_else(ScaledX <= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 22,
                                                                                                                                                                                             # Zone 23
                                                                                                                                                                                             if_else(ScaledX >= 3.75 - BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledX <= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                       ScaledY >= BinCoordAssign$ScaledY[16], 23,
                                                                                                                                                                                                     # Zone 24
                                                                                                                                                                                                     if_else(ScaledX >= 3.75 + BinCoord24Circ$AddValue[1] &
                                                                                                                                                                                                               ScaledY >= BinCoordAssign$ScaledY[16], 24, 0
                                                                                                                                                                                                     )))))))))))))))))))))))))
  
  # Join the site fidelity zone data sets
  Colony20BinnedSFZ <<- full_join(Colony20BinnedTubeSFZ, Colony20BinnedCircleSFZ) %>%
    select(Colony, Nest, Day, Bin, ScaledX, ScaledY, Zone, AntID)
}

# Run the zone assignment function on the usable color coordinates worker data set
CoordBinnedSFZ(SFZDataFull)

# Joining all of the separate site fidelity zone data sets
SFZDataFullZones <- 
  full_join(Colony1BinnedSFZ,Colony2BinnedSFZ) %>%
  full_join(Colony3BinnedSFZ) %>%
  full_join(Colony4BinnedSFZ) %>%
  full_join(Colony5BinnedSFZ) %>%
  full_join(Colony6BinnedSFZ) %>%
  full_join(Colony7BinnedSFZ) %>%
  full_join(Colony8BinnedSFZ) %>%
  full_join(Colony9BinnedSFZ) %>%
  full_join(Colony10BinnedSFZ) %>%
  full_join(Colony11BinnedSFZ) %>%
  full_join(Colony12BinnedSFZ) %>%
  full_join(Colony13BinnedSFZ) %>%
  full_join(Colony14BinnedSFZ) %>%
  full_join(Colony15BinnedSFZ) %>%
  full_join(Colony16BinnedSFZ) %>%
  full_join(Colony17BinnedSFZ) %>%
  full_join(Colony18BinnedSFZ) %>%
  full_join(Colony19BinnedSFZ) %>%
  full_join(Colony20BinnedSFZ) %>%
  distinct()

####################################################################################################################
# SPATIAL FIDELITY & OCCURRENCE ZONE CALCULATIONS
# The following functions are the site fidelity calculations
# The function uses a nest area reference and the color coordinate data set created above
# First, the reference data set is joined and a new column is created with a 1
# Then, we find the frequency of the color coordinate for the colony in each nest
# We require at least 3 observations of the color ID
# Next, we categorize spatial fidelity zones as zones having at least 15% of the total observations, where all obversations are included in the occurrence zone
# Last, each zone size is calculated by adding the number of zones in each category and multiplying that value by 1/24 the area of the nest 
####################################################################################################################
FidelityZones <- function(data_table){
  # Joining the area reference data set
  FidelityZonesDataRaw <- left_join(data_table, NestAreaFull) %>%
    # Creating a column filled with 1s
    mutate(AntIDNum = 1) %>%
    # Group by the Colony, Nest, Day, and AntID columns
    group_by(Colony, Nest, Day, AntID) %>%
    # Count the number of each group
    mutate(count = n()) %>%
    # There should only be one observation of each color ID
    # If this isn't true we don't want the color ID in case of error
    filter(count == 1) %>%
    # Group by the Colony, Nest, and AntID columns
    group_by(Colony, Nest, AntID) %>%
    # Finding the frequency of the color ID in each group
    mutate(Freq = sum(AntIDNum)) %>%
    # Removing all color IDs with fewer than three observations
    filter(Freq > 6) %>%
    # Group by the Colony, Nest, AntID, and Zone columns
    group_by(Colony, Nest, AntID, Zone) %>%
    # Count the number of observations in each group
    mutate(ZoneCount = n(),
           # Calculating the proportion of color ID observations in each observed zone 
           PropSFZ = (ZoneCount / Freq)) %>%
    select(-c(Day, ScaledX, ScaledY)) %>%
    # Remove any duplicates - i.e., days 2 and 10 can have two observations in zone 5, so there will be duplicate values of PropSFZ
    distinct() %>%
    # Group by the Colony, Nest, and AntID columns
    group_by(Colony, Nest, AntID) %>%
    # Determining whether the zone has at least 15% of total observations or not
    # If so, the zone will be included in fidelity zone
    mutate(FullZone = 1, 
           FidZone = ifelse(PropSFZ >= 0.15, 1, 0), 
           # Calculating the zone sizes
           Occur = sum(FullZone), # Scaled occurrence zones
           SFZ = sum(FidZone), # Scaled spatial fidelity zones
           Occur_Area = Area * (Occur / 24), # Unscaled occurrence zones
           SFZ_Area = Area * (SFZ / 24), # Unscaled spatial fidelity zones
           Density = ifelse(Colony < 11, "High", "Low")) %>%
    # Select the desired columns
    select(Colony, Nest, AntID, Density, SFZ, Occur, Occur_Area, SFZ_Area, Freq, Number.ants, FullZone, FidZone) %>%
    distinct()
  
  # Final data set
  FidelityZonesDataRD1_RD2 <<- left_join(data_table, FidelityZonesDataRaw) %>%
    distinct() %>%
    drop_na()
}

# Running the site fidelity calculation function on the SFZDataFullZones data set
FidelityZones(SFZDataFullZones)

####################################################################################################################
# SITE FIDELITY ZONE DISTANCE TO ENTRANCE
# The following function calculates the mean distance to the entrance for each worker to correspond with their site fidelity
# The function is very similar to the distance to the entrance functions used in the "DistanceFunctions.R" script. 
# In sum, the function uses reference shortest distances to the nest entrance for each of the eight nest sections
# The shortest distance to the nest entrance is then determined by finding the shortest pythagorean distance for each worker
# Which is then added to the appropriate reference distance
####################################################################################################################

DistanceCoordsFunctionSFZ<-function(data.table){
  Bin1 <- DistBinsFull %>% # Bin 1 has corners that are ignored unless there is a way to make the shortest distance to the corner first, this code produces the reference coordinates for the corners
    filter(Bin == 1) %>%
    group_by(Colony, Nest) %>%
    mutate(BinX1 = BinX, LeftCorner = BinX1 - 0.1, RightCorner = BinX1 + 0.1) %>%
    select(Colony, Nest, BinX1, LeftCorner, RightCorner)
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% # This first section of script makes "Bin4" a column that will make the shortest distance point the same as bin 3, which corrects the issue stated above
    group_by(Colony) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(Colony, BinY4 ,Distance4)
  Bin7 <- DistBinsFull %>% # Bin 7 has the same problem as above, one shortest distance that uses the x coordinate of bin 7, and another at Bin 4 where Scaled Y values are greater than the y coordinate of Bin 7
    filter(Nest == "Tube" & Bin == 4) %>%
    group_by(Colony) %>%
    mutate(BinX7 = BinX,
           Distance7 = Distance) %>% 
    select(Colony, BinX7, Distance7) # The distance that is added to the Pythagorean distance is also different
  Bin3 <- DistBinsFull %>% # Bin 3 has the same problem as above, one shortest distance that uses the x coordinate of bin 3, and another at Bin 2 where Scaled Y values are less than the y coordinate of Bin 3
    filter(Nest == "Tube" & Bin == 2) %>%
    group_by(Colony) %>%
    mutate(BinX3 = BinX,
           BinY3 = BinY,
           Distance3 = Distance) %>%
    select(Colony, BinY3, BinX3, Distance3) # The distance that is added to the Pythagorean distance is also different
  DistBinsFull <- left_join(DistBinsFull, Bin4) %>% #Full distance references, joining all alternative references from above
    left_join(Bin7) %>%
    left_join(Bin3) %>%
    left_join(Bin1)
  DistBins.tube <- DistBinsFull %>% # Subsetting tube distance references
    filter(Nest == "Tube") 
  DistBins.circle <- DistBinsFull %>% # Subsetting circle distance references
    filter(Nest == "Circle") 
  BinsTube <- data.table %>% # Subsetting tube nest data
    filter(Nest == "Tube")
  BinsCircle <- data.table %>% # Subsetting circle nest data
    filter(Nest == "Circle") 
  
  # TUBE NEST SHAPE
  DistBinsTube <- left_join(BinsTube, DistBins.tube) %>% # Joining tube nest data with associated reference coordinates
    group_by(Colony, Bin) %>% # Group by Colony and Bin columns
    filter(Bin != 1) %>% # Filter out Bin 1, since this will be handled in the code below
    mutate(DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference Bin4 coordinate
           DistanceX7 = ScaledX - BinX7, # Distance from each individual's x coordinate to the x reference Bin7 coordinate
           DistanceX3 = ScaledX - BinX3, # Distance from each individual's x coordinate to the x reference Bin2 coordinate
           DistanceY3 = ScaledY - BinY3, # Distance from each individual's x coordinate to the x reference Bin2 coordinate
           # Pythagorean theorm sqrt(X^2 + Y^2) = the hypotenuse (PythagDist)
           PythagDist = ifelse(Bin == 4 & Nest == "Tube" & ScaledX < BinX, 
                               # This uses an if else statement where if its coordinates from bin 4 in the tube nest, and the x coordinate is less than the reference for the bin, then it uses DistanceY4, else just DistanceY from above
                               sqrt((DistanceX^2) + (DistanceY4^2)), # Pythagorean distance from the alternative x and y axes distances for bin 4 calculated above
                               # This uses an if else statement where if its coordinates from bin 7 in the tube nest, and the y coordinate is greater than the reference for the, then it uses DistanceX7, else just DistanceX from above
                               ifelse(Bin == 7 & Nest == "Tube" & ScaledY > BinY, 
                                      sqrt((DistanceX7^2) + (DistanceY^2)), # Pythagorean distance from the alternative x and y axes distances for bin 7 calculated above
                                      # This uses an if else statement where if its coordinates from bin 3 in the tube nest, and the y coordinate is greater than the reference for the, then it uses DistanceX7, else just DistanceX from above
                                      ifelse(Bin == 3 & Nest == "Tube" & ScaledY < BinY4,
                                             sqrt((DistanceX3^2) + (DistanceY3^2)), # Pythagorean distance from the alternative x and y axes distances for bin 3 calculated above
                                             sqrt((DistanceX^2) + (DistanceY^2))))), # Pythagorean distance from the x and y axes distances calculated above 
           # Total distance, which includes the the closest bin to the nest entrance and reference distance for that bin to the nest entrance
           DistanceTotal = ifelse(Bin == 4 & Nest == "Tube" & ScaledX < BinX, # Nested ifelse statements that produce the shortest distance to the nest entrance for each coordinate
                                  PythagDist + Distance4, # Shortest distance for coordinates using the alternative bin 4 distance 
                                  ifelse(Bin == 7 & Nest == "Tube" & ScaledY > BinY, # Shortest distance for coordinates using the alternative bin 7 distance 
                                         PythagDist + Distance7, # Shortest distance for coordinates using the alternative bin 7 distance
                                         ifelse(Bin == 3 & Nest == "Tube" & ScaledY < BinY4,
                                                PythagDist + Distance3, # Shortest distance for coordinates using the alternative bin 3 distance
                                                PythagDist + Distance)))) # Shortest distance for coordinates without alternatives
  # Bin 1 for the tube nest, this considers corners at the entrance
  DistBinsTubeBin1 <- left_join(BinsTube, DistBins.tube) %>% # Joining tube nest data with associated reference coordinates
    group_by(Colony, Bin) %>% # Group by Colony and Bin columns
    filter(Bin == 1) %>% # Filter out all but Bin 1
    mutate(SegEntLeft = sqrt((BinX1 - LeftCorner)^2 + (0.2 - 0)^2), # Distance from the left corner to the entrance
           SegEntRight = sqrt((RightCorner - BinX1)^2 + (0.2 - 0)^2), # Distance from the right corner to the entrance
           AngleEnt1 = tan((0.2 - 0) / (BinX1 - LeftCorner)), # Angle between the distance from a corner to the entrance and the vertical line segment created from the entrance to the left corner
           AngleEnt2 = tan((0.2 - 0) / (RightCorner - BinX1)), # Angle between the distance from a corner to the entrance and the vertical line segment created from the entrance to the right corner
           DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           PythagDistEnt = ifelse(ScaledX < LeftCorner, # Pythagorean distance to the left or right corner, based on whether the x coordinate is greater than or less than the left corner x reference
                                  sqrt((LeftCorner - ScaledX)^2 + (ScaledY - 0.2)^2), # Pythagorean distance to the left corner
                                  sqrt((ScaledX - RightCorner)^2 + (ScaledY - 0.2)^2)), # Pythagorean distance to the right corner
           AngleEntExp = ifelse(ScaledX < LeftCorner, # ifelse statement that determines if the coordinate is on the left or right side, from the perspective of the corners 
                                cos(abs(LeftCorner - ScaledX) / PythagDistEnt), # Angle between the coordinate and the horizontal segment from the left corner to the left side nest wall in bin 1
                                cos(abs(ScaledX - RightCorner) / PythagDistEnt)), # Angle between the coordinate and the horizontal segment from the right corner to the right side nest wall in bin 1
           AngleEntExp = ifelse(ScaledX >= LeftCorner & ScaledX <= RightCorner, # ifelse statement that determines if coordinates are not in the field considered for the corners
                                180, AngleEntExp), # Assigns 180 to those coordinates, else it keeps the AngleEntExp from above
           # Checks is the total angle from AngleEnt1 or AngleEnt1 + AngleEntExp + 90 degrees is less than 180 degrees
           AngleCheckEnt = ifelse(((AngleEnt1 + AngleEntExp + 90) < 180) | ((AngleEnt2 + AngleEntExp + 90) < 180), "Yes", "No"),
           # Shortest distance to the entrance, where distances are either just the pythagorean distance from the coordinate to the nest entrance or first to a corner + the corner reference distance to the entrance
           DistanceTotal = ifelse(ScaledX <= BinX, # ifelse statement asking whether a coordinates x value is less than the x coordinate reference for bin 1
                                  # ifelse statement either calculating the alternative distance to the entrance for coordinates that would cut a corner
                                  ifelse(AngleCheckEnt == "Yes",  PythagDistEnt + SegEntLeft, # Left corner, if the condition is true
                                         sqrt((DistanceX^2) + (DistanceY^2))), # Uses a coordinates pythagorean distance to the nest entrance if the corner will not be cut ("No" from above)
                                  ifelse(AngleCheckEnt == "Yes",  PythagDistEnt + SegEntRight, # Right corner, if the condition is true
                                         sqrt((DistanceX^2) + (DistanceY^2))))) # Uses a coordinates pythagorean distance to the nest entrance if the corner will not be cut ("No" from above)
  # Joins together the tube nest shortest distances to the nest entrance
  DistBinsTube <- full_join(DistBinsTube, DistBinsTubeBin1) 
  
  # CIRCLE NEST SHAPE
  # Circle nest distances are the same as above except all coordinates run through the angle test 
  DistBinsCircle <- left_join(BinsCircle, DistBins.circle) %>%
    group_by(Colony) %>% # Group by the Colony column
    mutate(DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           SegEntLeft = sqrt((BinX1 - LeftCorner)^2 + (0.2 - 0)^2), # Distance from the left corner to the entrance
           SegEntRight = sqrt((RightCorner - BinX1)^2 + (0.2 - 0)^2), # Distance from the right corner to the entrance
           AngleEnt1 = tan((0.2 - 0) / (BinX1 - LeftCorner)), # Angle between the distance from a corner to the entrance and the vertical line segment created from the entrance to the left corner
           AngleEnt2 = tan((0.2 - 0) / (RightCorner - BinX1)), # Angle between the distance from a corner to the entrance and the vertical line segment created from the entrance to the right corner
           PythagDistEnt = ifelse(ScaledX < LeftCorner, # Pythagorean distance to the left or right corner, based on whether the x coordinate is greater than or less than the left corner x reference
                                  sqrt((LeftCorner - ScaledX)^2 + (ScaledY - 0.2)^2), # Pythagorean distance to the left corner
                                  sqrt((ScaledX - RightCorner)^2 + (ScaledY - 0.2)^2)), # Pythagorean distance to the right corner
           AngleEntExp = ifelse(ScaledX < LeftCorner, # ifelse statement that determines if the coordinate is on the left or right side, from the perspective of the corners 
                                cos(abs(LeftCorner - ScaledX) / PythagDistEnt), # Angle between the coordinate and the horizontal segment from the left corner to the left side nest wall in bin 1
                                cos(abs(ScaledX - RightCorner) / PythagDistEnt)), # Angle between the coordinate and the horizontal segment from the right corner to the right side nest wall in bin 1
           AngleEntExp = ifelse(ScaledX >= LeftCorner & ScaledX <= RightCorner, # ifelse statement that determines if coordinates are not in the field considered for the corners
                                180, AngleEntExp), # Assigns 180 to those coordinates, else it keeps the AngleEntExp from above
           # Checks is the total angle from AngleEnt1 or AngleEnt1 + AngleEntExp + 90 degrees is less than 180 degrees
           AngleCheckEnt = ifelse(((AngleEnt1 + AngleEntExp + 90) < 180) | ((AngleEnt2 + AngleEntExp + 90) < 180), "Yes", "No"),
           # Shortest distance to the entrance, where distances are either just the pythagorean distance from the coordinate to the nest entrance or first to a corner + the corner reference distance to the entrance
           DistanceTotal = ifelse(ScaledX <= BinX, # ifelse statement asking whether a coordinates x value is less than the x coordinate reference for bin 1
                                  # ifelse statement either calculating the alternative distance to the entrance for coordinates that would cut a corner
                                  ifelse(AngleCheckEnt == "Yes",  PythagDistEnt + SegEntLeft, # Left corner, if the condition is true
                                         sqrt((DistanceX^2) + (DistanceY^2))), # Uses a coordinates pythagorean distance to the nest entrance if the corner will not be cut ("No" from above)
                                  ifelse(AngleCheckEnt == "Yes",  PythagDistEnt + SegEntRight, # Right corner, if the condition is true
                                         sqrt((DistanceX^2) + (DistanceY^2))))) # Uses a coordinates pythagorean distance to the nest entrance if the corner will not be cut ("No" from above)
  
  # Has all distances and zones
  WorkerDistScaledRD1_RD2SFZFull <<- full_join(DistBinsTube, DistBinsCircle) %>%
    # Group by the columns Colony, Nest and AntID
    group_by(Colony, Nest, Day, AntID) %>%
    # Finding the mean distance of workers within the groupings
    mutate(ScaledDist = DistanceTotal / (MaxDist),
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist))%>%
    ungroup()%>%
    select(Colony, Nest, AntID, SFZ, Occur, ScaledDist, Density, Day, SFZ_Area, Occur_Area, Freq, Number.ants) %>%
    distinct()

  # Contains all mean distances
  WorkerDistScaledRD1_RD2SFZWorking <<- WorkerDistScaledRD1_RD2SFZFull %>%
    group_by(Colony, Nest, AntID) %>%
    mutate(MeanScaledDist = mean(ScaledDist)) %>%
    select(Colony, Nest, AntID, SFZ, Occur, MeanScaledDist, Density, SFZ_Area, Occur_Area, Freq, Number.ants)%>%
    distinct()
}

# Run the mean distance to the nest entrance function for the color marked worker data set FidelityZonesDataRD1_RD2
DistanceCoordsFunctionSFZ(FidelityZonesDataRD1_RD2)

####################################################################################################################
# CALCULATE THE BROOD CENTER 
# This code calculated the brood center for the nest section with the most brood in each observation
# Below first calculates the brood proportions in each nest section, then finds which nest section had the highest proportion
# The code then finds the centroid of the brood within the nest section with the highest proportion of brood in each observation
####################################################################################################################

# PROPORTIONS OF BROOD IN EACH NEST SECTION
# BROOD
# High density treatment
Prop_functionBrood <- function(data.table) {
  BroodProp <- data.table %>% # Creating the data set of brood proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of brood in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each brood in each bin in each observation
           PropBrood = (BinCount / count)) %>% # Calculate the proportion of brood in each bin 
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null brood proportions in each nest section
  BroodPropNull <- BroodProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(BroodPropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullBrood <<- full_join(NestArchNullBins, BroodProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(PropBrood = ifelse(is.na(PropBrood), 0, PropBrood),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropBrood)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% #Select only the desired columns
    left_join(CornerFull) %>%# Joins with a data set that assigned corner presence to each nest section
    distinct()
}

# Run the proportions of brood in nest sections function for the FullDataCoordBrood data set 
Prop_functionBrood(FullDataCoordBrood)

# BROOD
# Low density treatment
Prop_functionBrood <- function(data.table) {
  BroodProp <- data.table %>% # Creating the data set of brood proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of brood in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each brood in each bin in each observation
           PropBrood = (BinCount / count)) %>% # Calculate the proportion of brood in each bin 
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null brood proportions in each nest section
  BroodPropNull <- BroodProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(BroodPropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullBroodRD2 <<- full_join(NestArchNullBins, BroodProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(PropBrood = ifelse(is.na(PropBrood), 0, PropBrood),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropBrood)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% #Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}

# Run the proportions of brood in nest sections function for the FullDataCoordBroodRD2 data set 
Prop_functionBrood(FullDataCoordBroodRD2)

# Join brood proportions in nest sections data sets
AntPropFullBroodRD1_RD2 <- full_join(AntPropFullBrood, AntPropFullBroodRD2) 

# FINDING THE NEST SECTION WITH THE HIGHEST PROPORTION OF BROOD IN EACH OBSERVATION
MeanBroodCoordProps <- AntPropFullBroodRD1_RD2 %>%
  group_by(Colony, Nest, Day) %>% # Group by the Colony, Nest, and Day columns
  mutate(MaxBrood = max(PropBrood), # Determine the largest proportion of brood in each observation
         MaxBin = ifelse(PropBrood == MaxBrood, Bin, NA)) %>% # If the nest section holds the largest proportion the row copies the Bin column, else with NA 
  ungroup() %>% # Ungroup the data set
  drop_na() %>% # Drop NAs created from forming the MaxBin 
  select(Colony, Nest, Density, Day, Bin) # Select desired columns

# CALCULATE THE CENTROID OF THE NEST SECTION WITH THE HIGHEST PROPORTION OF BROOD IN EACH OBSERVATION
MeanBroodCoordFull <- left_join(MeanBroodCoordProps, FullDataCoordBroodRD1_RD2) %>% # Join brood observations that match the bins left over from the above function
  drop_na() %>% # Drop any NAs from the join
  group_by(Colony, Nest, Day) %>% # Group by the Colony Nest and Day columns
  mutate(BroodX = mean(ScaledX), # Calculate the x coordinate for the brood centroids
         BroodY = mean(ScaledY)) %>% # Calculate the y coordinate for the brood centroids
  rename(BroodBin = Bin) %>% # Rename the bin column to "BroodBin", so that the name is different than the mobile colony member data sets
  select(Colony, Nest, Day, BroodBin, BroodX, BroodY, Density) %>% # Select the desired columns
  distinct() # Remove duplicates

####################################################################################################################
# TUBE NEST SHAPE REFERENCE DISTANCES
# This code calculated all reference distances used in the remaining code
# Reference distances are the shortest distances between nest sections
####################################################################################################################


# First we make data sets of distance segments either between a bin and the entrance or between two sections
# This is done outside of the actual function below to take less computing time
TubeDistances <- DistBinsFull %>% # Filtering out the Tube reference distances
  filter(Nest == "Tube")

# Filtering out only distances for bin 1
DistBin1 <- TubeDistances %>% 
  filter(Bin == 1) %>%
  mutate(Dist1 = Distance) %>% # Set reference distance
  select(Colony, Nest, Dist1) # Select desired columns

# Filtering and creating reference distances for bins 2 - 8
# Bin 2
DistBin2 <- TubeDistances %>%
  filter(Bin == 2) %>%
  mutate(Dist2 = Distance) %>%
  select(Colony, Nest, Dist2)

# Bin 3
DistBin3 <- TubeDistances %>%
  filter(Bin == 3) %>%
  mutate(Dist3 = Distance) %>%
  select(Colony, Nest, Dist3)

# Bin 4
DistBin4 <- TubeDistances %>%
  filter(Bin == 4) %>%
  mutate(Dist4 = Distance) %>%
  select(Colony, Nest, Dist4)

# Bin 5
DistBin5 <- TubeDistances %>%
  filter(Bin == 5) %>%
  mutate(Dist5 = Distance) %>%
  select(Colony, Nest, Dist5)

# Bin 6
DistBin6 <- TubeDistances %>%
  filter(Bin == 6) %>%
  mutate(Dist6 = Distance) %>%
  select(Colony, Nest, Dist6)

# Bin 7
DistBin7 <- TubeDistances %>%
  filter(Bin == 7) %>%
  mutate(Dist7 = Distance) %>%
  select(Colony, Nest, Dist7)

# Bin 8
DistBin8 <- TubeDistances %>%
  filter(Bin == 8) %>%
  mutate(Dist8 = Distance) %>%
  select(Colony, Nest, Dist8)

# Create a combined data set of reference coordinates and distances for each bin
DistBin1_8Full <- full_join(DistBin1, DistBin2) %>%
  full_join(DistBin3) %>%
  full_join(DistBin4) %>%
  full_join(DistBin5) %>%
  full_join(DistBin6) %>%
  full_join(DistBin7) %>%
  full_join(DistBin8) %>%
  mutate(
    Distance8_6 = abs(Dist8 - Dist7), # Distance from bin 8 to 6
    Distance4_2 = Distance8_6, # Distance from bin 4 to 2
    Distance3_5 = Distance8_6, # Distance from bin 3 to 5
    Distance7_5 = Distance8_6, # Distance from bin 7 to 5 
    Distance8_1 = abs(Dist8 - Dist2), # Distance from bin 8 to 1
    Distance7_1 = abs(Dist7 - Dist2), # Distance from bin 7 to 1
    Distance4_1 = abs(Dist4 - Dist2), # Distance from bin 4 to 1
    Distance3_1 = abs(Dist3 - Dist2), # Distance from bin 3 to 1
    Distance8_2 = abs(Dist8 - Dist3), # Distance from bin 8 to 2
    Distance7_2 = abs(Dist7 - Dist3), # Distance from bin 7 to 2
    Distance7_3 = abs(Dist7 - Dist4), # Distance from bin 7 to 3
    Distance8_3 = abs(Dist8 - Dist4), # Distance from bin 8 to 3
    Distance5_1 = abs(Distance4_1 + Distance8_6), # Distance from bin 5 to 1 
    Distance6_1 = abs(Distance7_1 - Distance8_6), # Distance from bin 6 to 1
    Distance5_6 = abs(Distance6_1 - Distance5_1), # Distance from bin 5 to 6 
    Distance7_4 = abs(Distance7_3 - Distance8_6), # Distance from bin 7 to 4
    Distance2_5 = abs(Distance8_6 * 2), # Distance from bin 2 to 5
    Distance2_6 = abs(Distance2_5 + Distance5_6), # Distance from bin 2 to 6
    Distance3_6 = abs(Distance8_6 + Distance5_6), # Distance from bin 3 to 6 
    Distance8_4 = abs(Distance8_3 - Distance8_6), # Distance from bin 8 to 4
    Distance8_5 = abs(Distance8_4 - Distance5_6) # Distance from bin 8 to 5
    )

# Function that uses the data sets of individual or simulation distances to the entrance
# The function uses these distances and either adds or subtracts reference distances to obtain the shortest distance to each nest section

# Reduce the distance bins reference data set to desired columns for only tube nest references
DistBinsTubeRefFull <- DistBinsFull %>%
  filter(Nest == "Tube") %>% 
  select(-c(Distance, Xmax, Ymax, TubeRatio)) 

# Make reference coordinates for each bin
# Bin 1
DistBinsTubeRef1 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 1) %>%
  rename(BinXRef1 = BinX, BinYRef1 = BinY) %>%
  select(-c(Bin))

# Bin 2
DistBinsTubeRef2 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 2) %>%
  rename(BinXRef2 = BinX, BinYRef2 = BinY) %>%
  select(-c(Bin))

# Bin 3
DistBinsTubeRef3 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 3) %>%
  rename(BinXRef3 = BinX, BinYRef3 = BinY)%>%
  select(-c(Bin))

# Bins 4-6
DistBinsTubeRef4_6 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 4) %>%
  rename(BinXRef4_6 = BinX, BinYRef4_6 = BinY)%>%
  select(-c(Bin))

# Bin 7
DistBinsTubeRef7 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 7) %>%
  rename(BinXRef7 = BinX, BinYRef7 = BinY)%>%
  select(-c(Bin))

# Bin 8
DistBinsTubeRef8 <- DistBinsTubeRefFull %>%
  filter(Nest == "Tube" & Bin == 8) %>%
  rename(BinXRef8 = BinX, BinYRef8 = BinY)%>%
  select(-c(Bin))

# Combine all of the reference distance data sets
DistBinsTubeRef <- left_join(DistBinsTubeRefFull, DistBinsTubeRef1) %>%
  left_join(DistBinsTubeRef2) %>%
  left_join(DistBinsTubeRef3) %>%
  left_join(DistBinsTubeRef4_6) %>%
  left_join(DistBinsTubeRef7) %>%
  left_join(DistBinsTubeRef8)

# Filtering tube nest data from the brood centroid data set and joining reference coordinate and distance columns
MeanBroodCoordFullTube <- MeanBroodCoordFull %>% 
  filter(Nest == "Tube") %>%
  left_join(DistBinsTubeRef)

# Filtering circle nest data from the brood centroid data set and joining reference coordinate and distance columns
MeanBroodCoordFullCircle <- MeanBroodCoordFull %>% 
  filter(Nest == "Circle") %>%
  left_join(DistBinsFull)

####################################################################################################################
# DISTANCE TO THE BROOD CENTER
# The code finds each mobile color marked worker's distance to the brood center, then scales it
# The scale is such that the longest shortest distance in the tube nest shape (entrance-to-back) is 1
# The code first finds these distances for the tube nest, using the references calculated above
# The code then finds these distances for the circle nest, which is just the Pythagorean distance between the points
####################################################################################################################

DistanceToBroodSFZFunction <- function(data.table){
  # Create the alternative bin 4 reference y coordinate
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% # This first section of script makes "Bin4" a column that will make the shortest distance point the same as bin 3, which corrects the issue stated above
    group_by(Colony) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(Colony, BinY4, Distance4)
  
  Bin7 <- DistBinsFull %>%
    filter(Nest == "Tube" & Bin == 4) %>%
    group_by(Colony) %>%
    mutate(BinX7 = BinX) %>%
    select(Colony, BinX7)
  
  # TUBE NEST
  DistanceToBroodSFZTube <- data.table %>% # Full worker data set
    filter(Nest == "Tube") %>% # Filter out tube nest shape
    left_join(MeanBroodCoordFullTube) %>% # Join the tube nest brood centroid data set
    left_join(DistBin1_8Full) %>% # Join the tube nest reference distance data set
    left_join(Bin4) %>%
    left_join(Bin7) %>%
    drop_na() %>% # Drop any NAs
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    # Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), # Shortest distance from bins 3 to 4
           # Creating the x and y distances from each individual to all Bin x and y references
           # Distances from the individual to the occupied Bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference bin 4 coordinate
           DistanceX7 = ScaledX - BinX7, # Distance from each worker x coordinate to the second bin 6 and 7 x coordinate creating a direct path to bin 4 
           # Calculating the shortest distance from each individual to each nest section
           # This uses the Pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           # A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each individual
           # Shortest distances from each individual to bins towards the nest entrance
           # If the individual is in bin 4 but to the left of the x reference, use the first formula, else use the second
           # This is only done for nest sections towards the entrance, the rest are done using reference coordinates
           PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               sqrt((DistanceX^2) + (DistanceY^2))),
           PythagDist4 = sqrt((DistanceX7^2) + (DistanceY^2)), # Creates the shortest distance to bin 4
           # Calculating the distance of each individual to the brood center
           # If the individual and brood centroid can be connected with a straight line, then this is done with the Pythagorean theorem 
           # If they aren't, the shortest distance to the bin closest to each other is calculated for each, than a reference distance is added where needed
           BroodDist = 
             # Distances from workers to brood centroids in bin 1
             ifelse(BroodBin == 1,
                    # Bin 1
                    ifelse(Bin == 1, sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                           # Bin 2
                           ifelse(Bin == 2, 
                                  ifelse(BroodY > BinYRef2, 
                                         sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                         PythagDist + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                  ifelse(Bin == 3, 
                                         ifelse(ScaledY < BinYRef3,
                                                ifelse(BroodY > BinYRef2, 
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                                ifelse(BroodY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2)) + Distance3_1)),
                                         ifelse(Bin == 4,
                                                ifelse(ScaledX < BinXRef4_6,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance3_1),
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance3_1)),
                                                # Bin 5
                                                ifelse(Bin == 5,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance3_1),
                                                       # Bin 6
                                                       ifelse(Bin == 6,
                                                              ifelse(BroodY > BinYRef2,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance3_1),
                                                              # Bin 7
                                                              ifelse(Bin == 7,
                                                                     ifelse(ScaledY > BinYRef7,
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist4 + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist4 + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance3_1),
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance7_1 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)))),
                                                                     # Bin 8
                                                                     ifelse(Bin == 8,
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance8_2 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2))),
                                                                            NA)
                                                              ))))))),
                    # Distances from workers to brood centroids in bin 2
                    ifelse(BroodBin == 2,
                           # Bin 1
                           ifelse(Bin == 1, 
                                  ifelse(ScaledY > BinYRef2,
                                         sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                         sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                  # Bin 2
                                  ifelse(Bin == 2,
                                         sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                         # Bin 3
                                         ifelse(Bin == 3,
                                                ifelse(ScaledY < BinYRef3,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                                # Bin 4
                                                ifelse(Bin == 4,
                                                       ifelse(ScaledX < BinXRef4_6,
                                                              PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                       # Bin 5
                                                       ifelse(Bin == 5,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              # Bin 6
                                                              ifelse(Bin == 6,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     # Bin 7
                                                                     ifelse(Bin == 7,
                                                                            ifelse(ScaledY > BinYRef7,
                                                                                   PythagDist4 + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                            # Bin 8
                                                                            ifelse(Bin == 8,
                                                                                   PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   NA)))))))),
                           # Distances from workers to brood centroids in bin 3
                           ifelse(BroodBin == 3,
                                  # Bin 1
                                  ifelse(Bin == 1,
                                         ifelse(BroodY < BinYRef3,
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance3_1)),
                                         # Bin 2
                                         ifelse(Bin == 2,
                                                ifelse(BroodY < BinYRef3,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                                # Bin 3
                                                ifelse(Bin == 3,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       # Bin 4
                                                       ifelse(Bin == 4,
                                                              ifelse(ScaledX < BinXRef4_6,
                                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                     PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                              # Bin 5
                                                              ifelse(Bin == 5,
                                                                     PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                     # Bin 6
                                                                     ifelse(Bin == 6,
                                                                            PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                            # Bin 7
                                                                            ifelse(Bin == 7,
                                                                                   ifelse(ScaledY > BinYRef7,
                                                                                          PythagDist4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          PythagDist + Distance7_3 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                   ifelse(Bin == 8,
                                                                                          PythagDist + Distance8_3 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          NA)))))))),
                                  # Distances from workers to brood centroids in bin 4
                                  ifelse(BroodBin == 4,
                                         # Bin 1
                                         ifelse(Bin == 1,
                                                ifelse(BroodX < BinXRef4_6,
                                                       ifelse(ScaledY > BinYRef2, 
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance3_1),
                                                       ifelse(ScaledY > BinYRef2,
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1)),
                                                # Bin 2
                                                ifelse(Bin == 2,
                                                       ifelse(BroodX < BinXRef4_6,
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4),
                                                       # Bin 3
                                                       ifelse(Bin == 3,
                                                              ifelse(BroodX < BinXRef4_6,
                                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                     sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                              # Bin 4
                                                              ifelse(Bin == 4,
                                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                     # Bin 5
                                                                     ifelse(Bin == 5,
                                                                            sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                            # Bin 6
                                                                            ifelse(Bin == 6,
                                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                   # Bin 7
                                                                                   ifelse(Bin == 7,
                                                                                          ifelse(ScaledY > BinYRef7,
                                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                 PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                          # Bin 8
                                                                                          ifelse(Bin == 8,
                                                                                                 PythagDist + Distance8_6 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)), 
                                                                                                 NA)))))))),
                                         # Distances from workers to brood centroids in bin 5
                                         ifelse(BroodBin == 5,
                                                # Bin 1
                                                ifelse(Bin == 1,
                                                       ifelse(ScaledY > BinYRef2,
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                       # Bin 2
                                                       ifelse(Bin == 2,
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                              # Bin 3
                                                              ifelse(Bin == 3,
                                                                     sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                     # Bin 4
                                                                     ifelse(Bin == 4,
                                                                            sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                            # Bin 5
                                                                            ifelse(Bin == 5,
                                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                   # Bin 6
                                                                                   ifelse(Bin == 6,
                                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                          # Bin 7
                                                                                          ifelse(Bin == 7,
                                                                                                 ifelse(ScaledY > BinYRef7,
                                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                        PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                                 # Bin 8
                                                                                                 ifelse(Bin == 8,
                                                                                                        PythagDist + Distance8_6 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                        NA)))))))),
                                                # Distances from workers to brood centroids in bin 6
                                                ifelse(BroodBin == 6,
                                                       # Bin 1
                                                       ifelse(Bin == 1,
                                                              ifelse(ScaledY > BinYRef2,
                                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                     sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                              # Bin 2
                                                              ifelse(Bin == 2,
                                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                     # Bin 3
                                                                     ifelse(Bin == 3,
                                                                            sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                            # Bin 4
                                                                            ifelse(Bin == 4,
                                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                   # Bin 5
                                                                                   ifelse(Bin == 5,
                                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                          # Bin 6
                                                                                          ifelse(Bin == 6,
                                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                 # Bin 7
                                                                                                 ifelse(Bin == 7,
                                                                                                        ifelse(ScaledY > BinYRef7,
                                                                                                               sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                               PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                                        # Bin 8
                                                                                                        ifelse(Bin == 8,
                                                                                                               PythagDist + Distance8_6 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                               NA)))))))),
                                                       # Distances from workers to brood centroids in bin 8, note that there weren't any brood centroids in bin 7, so it is skipped here
                                                       ifelse(BroodBin == 8,
                                                              # Bin 1
                                                              ifelse(Bin == 1,
                                                                     ifelse(ScaledY > BinYRef2,
                                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                            sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_1),
                                                                     # Bin 2
                                                                     ifelse(Bin == 2,
                                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                            # Bin 3
                                                                            ifelse(Bin == 3,
                                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_4,
                                                                                   # Bin 4
                                                                                   ifelse(Bin == 4,
                                                                                          sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_6,
                                                                                          # Bin 5
                                                                                          ifelse(Bin == 5,
                                                                                                 sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_6,
                                                                                                 # Bin 6
                                                                                                 ifelse(Bin == 6,
                                                                                                        ifelse(ScaledX > BinXRef7,
                                                                                                               sqrt(((ScaledX - BinXRef8)^2) + ((ScaledY - BinYRef8)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)),
                                                                                                               sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_6),
                                                                                                        ifelse(Bin == 7,
                                                                                                               ifelse(ScaledY < BinYRef8,
                                                                                                                      sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                                      sqrt(((ScaledX - BinXRef8)^2) + ((ScaledY - BinYRef8)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2))),
                                                                                                               # Bin 8
                                                                                                               ifelse(Bin == 8,
                                                                                                                      sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                                      NA)
                                                                                                        )
                                                                                                 )
                                                                                          )
                                                                                   )
                                                                            )
                                                                     )
                                                              ),
                                                              NA
                                                       )))))))) %>%
    group_by(Colony) %>% # Group by the colony column
    mutate(ToBrood = BroodDist / MaxDist, # Scale the shortest distance from each individual to the brood center
           ToBrood = ifelse(ToBrood > 1, 1, ToBrood)) %>% # As a precaution, if any scaled distances are greater than 1, the value is converted to 1
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, ToBrood, AntID, SFZ, Occur, SFZ_Area, Occur_Area) %>% # Select desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs

# CIRCLE NEST
DistanceToBroodSFZCircle <- data.table %>% # Full worker data set 
  filter(Nest == "Circle") %>% # Filtering out circle nest coordinates
  left_join(MeanBroodCoordFullCircle) %>% # Joining the circle brood centroid data set
  group_by(Colony, Day) %>% # Group by the Colony and Day columns
  mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)), # Shortest distance from each individual to the brood center
         ToBrood = BroodDist / MaxDist) %>% # Select the desired columns 
  select(Colony, Nest, Day, ScaledX, ScaledY, Density, ToBrood, AntID, SFZ, Occur, SFZ_Area, Occur_Area) %>% # Select the desired columns
  distinct() %>% # Remove any duplicates (shouldn't exist)
  drop_na() # Remove NAs

BroodCentDistWorkersSFZFull <<- full_join(DistanceToBroodSFZTube, DistanceToBroodSFZCircle) %>%
  unite('AntIDColNest', c(Colony, AntID, Nest), remove = FALSE) %>%
  mutate(AntIDColNest = as.factor(AntIDColNest))

BroodCentDistWorkersSFZ <<- full_join(DistanceToBroodSFZTube, DistanceToBroodSFZCircle) %>%
  group_by(Colony, Nest, AntID) %>%
  mutate(MeanToBrood = mean(ToBrood)) %>%
  select(-c(Day, ScaledX, ScaledY, ToBrood)) %>%
  distinct() 
}

# Run the mean distance to the brood center for the color marked workers data set FidelityZonesDataRD1_RD2
DistanceToBroodSFZFunction(FidelityZonesDataRD1_RD2)

####################################################################################################################
# SUPPLEMENTARY FOR NUMBER OF OBSERVATIONS & SITE FIDELITY FIGURE 
# SPATIAL FIDELITY & OCCURRENCE ZONE CALCULATIONS
# The following functions are the site fidelity calculations
# The function uses a nest area reference and the color coordinate data set created above
# First, the reference data set is joined and a new column is created with a 1
# Then, we find the frequency of the color coordinate for the colony in each nest
# We require at least 3 observations of the color ID
# Next, we categorize spatial fidelity zones as zones having at least 15% of the total observations, where all obversations are included in the occurrence zone
# Last, each zone size is calculated by adding the number of zones in each category and multiplying that value by 1/24 the area of the nest 
####################################################################################################################
FidelityZonesSupp <- function(data_table){
  # Joining the area reference data set
  FidelityZonesDataRaw <- left_join(data_table, NestAreaFull) %>%
    # Creating a column filled with 1s
    mutate(AntIDNum = 1) %>%
    # Group by the Colony, Nest, Day, and AntID columns
    group_by(Colony, Nest, Day, AntID) %>%
    # Count the number of each group
    mutate(count = n()) %>%
    # There should only be one observation of each color ID
    # If this isn't true we don't want the color ID in case of error
    filter(count == 1) %>%
    # Group by the Colony, Nest, and AntID columns
    group_by(Colony, Nest, AntID) %>%
    # Finding the frequency of the color ID in each group
    mutate(Freq = sum(AntIDNum)) %>%
    # Removing all color IDs with fewer than three observations
    filter(Freq > 2) %>%
    # Group by the Colony, Nest, AntID, and Zone columns
    group_by(Colony, Nest, AntID, Zone) %>%
    # Count the number of observations in each group
    mutate(ZoneCount = n(),
           # Calculating the proportion of color ID observations in each observed zone 
           PropSFZ = (ZoneCount / Freq)) %>%
    select(-c(Day, ScaledX, ScaledY)) %>%
    # Remove any duplicates - i.e., days 2 and 10 can have two observations in zone 5, so there will be duplicate values of PropSFZ
    distinct() %>%
    # Group by the Colony, Nest, and AntID columns
    group_by(Colony, Nest, AntID) %>%
    # Determining whether the zone has at least 15% of total observations or not
    # If so, the zone will be included in fidelity zone
    mutate(FullZone = 1, 
           FidZone = ifelse(PropSFZ >= 0.15, 1, 0), 
           # Calculating the zone sizes
           Occur = sum(FullZone), # Scaled occurrence zones
           SFZ = sum(FidZone), # Scaled spatial fidelity zones
           Occur_Area = Area * (Occur / 24), # Unscaled occurrence zones
           SFZ_Area = Area * (SFZ / 24), # Unscaled spatial fidelity zones
           Density = ifelse(Colony < 11, "High", "Low")) %>%
    # Select the desired columns
    select(Colony, Nest, AntID, Density, SFZ, Occur, Occur_Area, SFZ_Area, Freq, Number.ants, FullZone, FidZone) %>%
    distinct()
  
  # Final data set
  FidelityZonesDataRD1_RD2Supp <<- left_join(data_table, FidelityZonesDataRaw) %>%
    select(Colony, Nest, AntID, Density, SFZ, Occur, Freq) %>%
    distinct() %>%
    drop_na()
}

# Running the site fidelity calculation function on the SFZDataFullZones data set
FidelityZonesSupp(SFZDataFullZones)
