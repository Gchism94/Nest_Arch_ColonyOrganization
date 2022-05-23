####################################################################################################################
## Author: GREG CHISM
## Date: APR 2022
## email: gchism@arizona.edu
## Project: Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space
## Title: Nest section bin functions & nest density calculations
####################################################################################################################

# BIN ASSIGNMENT FUNCTION & NEST SECTION PROPORTIONS
# This code does the following:  
# (1) Processes empirical and netlogo data for nest section binning (eight even area bins) 
# (2) Bins the data into the eight even area nest sections
# (3) Calculates the proportions of individuals in each nest sections

install.packages("pacman") # Download package with function to load multiple packaged at once
####################################################################################################################
# Loading required packages for code below  
# p_load() will download packages that aren't in system library
####################################################################################################################

pacman::p_load(assertthat,
               forcats,
               formattable,
               ggdist,
               ggpointdensity,
               ggpubr,
               lme4,
               lmerTest,
               magrittr,
               MuMIn,
               RColorBrewer,
               scales,
               tidyverse,
               twosamples,
               viridis)

####################################################################################################################
# IMPORT ALL NECESSARY DATASETS 
# This code imports all necessary data sets for the script below
####################################################################################################################

# WORKERS
# High nest density
FullDataCoordWorkers <- read.csv("FullDataCoordWorkers.csv")

# Low nest density
FullDataCoordWorkersRD2 <- read.csv("FullDataCoordWorkersRD2.csv")

# BROOD
# High nest density
FullDataCoordBrood <- read.csv("FullDataCoordBrood.csv")

# Low nest density
FullDataCoordBroodRD2 <- read.csv("FullDataCoordBroodRD2.csv")

# QUEENS
# High nest density
FullDataCoordQueen <- read.csv("FullDataCoordQueen.csv")

# Low nest density
FullDataCoordQueenRD2 <- read.csv("FullDataCoordQueenRD2.csv")

# ALATES
# Low nest density
FullDataCoordAlate <- read.csv("FullDataCoordAlate.csv")

# NETLOGO 
ArchitectureMoveModelFull <- read.csv("ArchitectureMoveModelFull.csv")

# NULL BINS REFERENCE (EMPIRICAL)
BinsNullFull <- read.csv("BinsNullFull.csv")

# NULL BINS REFERENCE (NETLOGO SIMULATIONS)
BinsNullNetlogo <- read.csv("BinsNullNetlogo.csv")

# BIN REFERENCE COORDINATES (EMPIRICAL)
BinCoordFull <- read.csv("BinCoordFull.csv")

# BIN REFERENCE COORDINATES (NETLOGO SIMULATIONS)
BinCoordNetlogo <- read.csv("BinCoordNetlogo.csv")

# REFERENCE COORDINATES FOR CORNERS (EMPIRICAL)
CornerFull <- read.csv("CornerFull.csv")

# REFERENCE COORDINATES FOR CORNERS (NETLOGO SIMULATIONS)
CornerFullSim <- read.csv("CornerFullSim.csv")

####################################################################################################################
# EMPIRICAL DATA PROCESSING 
# This code processes empirical data for binning 
####################################################################################################################

# Joining all databases, to be subset later into individual colony members
# WORKERS
FullDataCoordWorkersRD1_RD2 <- full_join(FullDataCoordWorkers, FullDataCoordWorkersRD2) %>%
  mutate(ColonyMember = "Workers")

# BROOD
FullDataCoordBroodRD1_RD2 <- full_join(FullDataCoordBrood, FullDataCoordBroodRD2) %>%
  mutate(ColonyMember = "Brood")

# QUEENS
FullDataCoordQueenRD1_RD2 <- full_join(FullDataCoordQueen, FullDataCoordQueenRD2) %>%
  mutate(ColonyMember = "Queens")

# ALATES
FullDataCoordAlates <- FullDataCoordAlate %>%
  mutate(ColonyMember = "Alates") %>%
  separate(SexID, c("Alate", "Sex", "SexNumber", "TotalNumber"), sep = ",", remove = TRUE) %>%
  select(-c(Alate))

# Alate sex ratios
# Note, these ratios include "?", where the sex of the alate was uncertain.
# We chose males because they were always present when alates were in the nest.
FullDataCoordAlatesMales <- FullDataCoordAlates %>%
  filter(TotalNumber != "?" & SexNumber != "?") %>%
  group_by(Colony, Nest, Day, Sex) %>% #Group by the Colony, Nest, and Day columns, allowing us to find the below values for each photo.
  mutate(MaxTotal = max(as.numeric(TotalNumber)), #The number of alates in that observation
         MaxSex = max(as.numeric(SexNumber)), #The number of males in that observation
         Ratio = MaxSex / MaxTotal) %>% #The ratio of males / total alates
  filter(Sex == "M") %>%
  ungroup() %>%
  select(Colony, Nest, Day, Ratio) %>% #Select the desired columns
  distinct() #Remove duplicate rows

# Join the alate and male data sets, this is the data set used below
FullDataCoordAlatesRatio <- left_join(FullDataCoordAlates, FullDataCoordAlatesMales)

# ALL COLONY MEMBERS
FullDataCoordAll <- full_join(FullDataCoordWorkersRD1_RD2, FullDataCoordBroodRD1_RD2) %>%
  full_join(FullDataCoordQueenRD1_RD2) %>%
  full_join(FullDataCoordAlatesRatio)

####################################################################################################################
# NETLOGO DATA PROCESSING
# This code processes the netlogo random-walk simulated data set
####################################################################################################################

# Removing brackets from xcor and ycor list 
ArchitectureMoveModelFullCorrect <- ArchitectureMoveModelFull %>%
  mutate(xcor = gsub("\\[|\\]", "", xcor),
         ycor = gsub("\\[|\\]", "", ycor))

# Keep desired columns
ArchitectureMoveModelFullCorrectRed <- ArchitectureMoveModelFullCorrect %>%
  select(c(RunNumber, NestSize, Nest, TimeStep))

# Keeping only X and Y coordinates
# X coordinate
ArchitectureMoveModelFullCorrectRedXcor <- ArchitectureMoveModelFullCorrect %>%
  select(-c(ycor))

# Y coordinate
ArchitectureMoveModelFullCorrectRedYcor <- ArchitectureMoveModelFullCorrect %>%
  select(-c(xcor))

# Splitting the xcor and ycor lists
XCoordCorrect <- strsplit(ArchitectureMoveModelFullCorrectRedXcor$xcor, split = " ") # xcor
YCoordCorrect <- strsplit(ArchitectureMoveModelFullCorrectRedYcor$ycor, split = " ") # ycor

# Create a column that assigns unique IDs to each xcor and ycor row within each simulation run
# X coordinate
NetlogoTestX <- data.frame(RunNumber = rep(ArchitectureMoveModelFullCorrectRedXcor$RunNumber, sapply(XCoordCorrect, length)), xcor = unlist(XCoordCorrect)) %>%
  rowid_to_column(var = "id")

# Y coordinate
NetlogoTestY <- data.frame(RunNumber = rep(ArchitectureMoveModelFullCorrectRedYcor$RunNumber, sapply(YCoordCorrect, length)), ycor = unlist(YCoordCorrect)) %>%
  rowid_to_column(var = "id")

# Creating the working Netlogo data sets
NetlogoTestFull <- left_join(NetlogoTestX, NetlogoTestY) %>%
  mutate(xcor = as.numeric(xcor),
         ycor = as.numeric(ycor)) %>%
  select(-c(id)) %>%
  left_join(ArchitectureMoveModelFullCorrectRed) %>%
  mutate(xcor = xcor * 0.1,
         ycor = ycor * 0.1) %>% 
  rename(ScaledX = xcor, ScaledY = ycor)

####################################################################################################################
# EMPIRICAL DATA BINNING FUNCTION
# The code below bins x and y coordinate colony data into eight even area nest sections
# The code for Netlogo simulated x and y coordinate results is separate and below
# To do this, a reference data set of bin coordinates is used and coordinates are run through a series of conditional statements
# Where each conditional statement checks whether the coordinate is in one of eight bins sequentially

# The code is set up such that you can run each colony sequentially, which was done to avoid errors in loops that can result in losing an entire set of data 
####################################################################################################################

# COLONY 1
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "1")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony1BinnedTube <- data_table %>%
    filter(Colony == "1" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                                       )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony1BinnedCircle <- data_table %>%
    filter(Colony == "1" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                             )))))))))
  Colony1Binned <<- full_join(Colony1BinnedTube, Colony1BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 2
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "2")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony2BinnedTube <- data_table %>%
    filter(Colony == "2" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony2BinnedCircle <- data_table %>%
    filter(Colony == "2" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony2Binned <<- full_join(Colony2BinnedTube, Colony2BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 3
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "3")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony3BinnedTube <- data_table %>%
    filter(Colony == "3" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony3BinnedCircle <- data_table %>%
    filter(Colony == "3" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony3Binned <<- full_join(Colony3BinnedTube, Colony3BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 4
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "4")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony4BinnedTube <- data_table %>%
    filter(Colony == "4" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony4BinnedCircle <- data_table %>%
    filter(Colony == "4" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony4Binned <<- full_join(Colony4BinnedTube, Colony4BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 5
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "5")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony5BinnedTube <- data_table %>%
    filter(Colony == "5" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony5BinnedCircle <- data_table %>%
    filter(Colony == "5" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony5Binned <<- full_join(Colony5BinnedTube, Colony5BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 6
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "6")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony6BinnedTube <- data_table %>%
    filter(Colony == "6" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony6BinnedCircle <- data_table %>%
    filter(Colony == "6" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony6Binned <<- full_join(Colony6BinnedTube, Colony6BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 7
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "7")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony7BinnedTube <- data_table %>%
    filter(Colony == "7" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony7BinnedCircle <- data_table %>%
    filter(Colony == "7" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony7Binned <<- full_join(Colony7BinnedTube, Colony7BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 8
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "8")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony8BinnedTube <- data_table %>%
    filter(Colony == "8" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony8BinnedCircle <- data_table %>%
    filter(Colony == "8" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony8Binned <<- full_join(Colony8BinnedTube, Colony8BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 9
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "9")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony9BinnedTube <- data_table %>%
    filter(Colony == "9" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony9BinnedCircle <- data_table %>%
    filter(Colony == "9" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony9Binned <<- full_join(Colony9BinnedTube, Colony9BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 10
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "10")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony10BinnedTube <- data_table %>%
    filter(Colony == "10" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony10BinnedCircle <- data_table %>%
    filter(Colony == "10" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony10Binned <<- full_join(Colony10BinnedTube, Colony10BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 11
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "11")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony11BinnedTube <- data_table %>%
    filter(Colony == "11" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony11BinnedCircle <- data_table %>%
    filter(Colony == "11" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony11Binned <<- full_join(Colony11BinnedTube, Colony11BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 12
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "12")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony12BinnedTube <- data_table %>%
    filter(Colony == "12" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony12BinnedCircle <- data_table %>%
    filter(Colony == "12" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony12Binned <<- full_join(Colony12BinnedTube, Colony12BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 13
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "13")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony13BinnedTube <- data_table %>%
    filter(Colony == "13" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony13BinnedCircle <- data_table %>%
    filter(Colony == "13" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony13Binned <<- full_join(Colony13BinnedTube, Colony13BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 14
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "14")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony14BinnedTube <- data_table %>%
    filter(Colony == "14" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony14BinnedCircle <- data_table %>%
    filter(Colony == "14" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony14Binned <<- full_join(Colony14BinnedTube, Colony14BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 15
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "15")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony15BinnedTube <- data_table %>%
    filter(Colony == "15" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony15BinnedCircle <- data_table %>%
    filter(Colony == "15" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony15Binned <<- full_join(Colony15BinnedTube, Colony15BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 16
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "16")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony16BinnedTube <- data_table %>%
    filter(Colony == "16" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony16BinnedCircle <- data_table %>%
    filter(Colony == "16" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony16Binned <<- full_join(Colony16BinnedTube, Colony16BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 17
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "17")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony17BinnedTube <- data_table %>%
    filter(Colony == "17" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony17BinnedCircle <- data_table %>%
    filter(Colony == "17" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony17Binned <<- full_join(Colony17BinnedTube, Colony17BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll) 

# COLONY 18
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "18")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony18BinnedTube <- data_table %>%
    filter(Colony == "18" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony18BinnedCircle <- data_table %>%
    filter(Colony == "18" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony18Binned <<- full_join(Colony18BinnedTube, Colony18BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll)

# COLONY 19
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "19")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony19BinnedTube <- data_table %>%
    filter(Colony == "19" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony19BinnedCircle <- data_table %>%
    filter(Colony == "19" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony19Binned <<- full_join(Colony19BinnedTube, Colony19BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll)

# COLONY 20
CoordBinned <- function(data_table){
  # Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "20")
  
  # Binning the tube colony member coordinates
  # Filtering out the tube data for the selected colony
  Colony20BinnedTube <- data_table %>%
    filter(Colony == "20" & Nest == "Tube") %>%
    # Ifelse conditional statements for each bin
    mutate(Bin =
             # The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             # The order is always x coordinates first, then y
             # Bin 1
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               # Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       # Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               # Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       # Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               # Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       # Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               # Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                               )))))))))
  # Binning circle nest coordinates 
  # Filtering out the circle data for the selected colony
  # The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony20BinnedCircle <- data_table %>%
    filter(Colony == "20" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                     )))))))))
  Colony20Binned <<- full_join(Colony20BinnedTube, Colony20BinnedCircle)
}

# Run the eight nest section binning function for the combined colony member data set
CoordBinned(FullDataCoordAll)

# The following can be used to check how many coordinates have zero bin values, meaning they were not placed in a bin
# These coordinates need to be removed as they are most likely an error (a coordinate where no ant exists)
# This also checks whether the correct colony assignment was used throughout
# Each issue was examined in Excel 

# e.g.Colony1
Colony1Binned %>%
  filter(Bin == "0") %>%
  group_by(Bin) %>%
  summarise(n = n())

# Combine all data sets into the final one for workers, colonies 1-10
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordWorkers <- Colony1Binned %>%
  full_join(Colony2Binned) %>%
  full_join(Colony3Binned) %>%
  full_join(Colony4Binned) %>%
  full_join(Colony5Binned) %>% 
  full_join(Colony6Binned) %>% 
  full_join(Colony7Binned) %>%
  full_join(Colony8Binned) %>% 
  full_join(Colony9Binned) %>%
  full_join(Colony10Binned) %>%
  # Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Workers" & Colony < 11) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ColorID, Density)

# Combine all data sets into the final one for workers, colonies 11-20
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordWorkersRD2 <- Colony11Binned %>%
  full_join(Colony12Binned) %>%
  full_join(Colony13Binned) %>%
  full_join(Colony14Binned) %>%
  full_join(Colony15Binned) %>%
  full_join(Colony16Binned) %>%
  full_join(Colony17Binned) %>%
  full_join(Colony18Binned) %>%
  full_join(Colony19Binned) %>%
  full_join(Colony20Binned) %>%
  # Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Workers" & Colony > 10) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ColorID, Density)

# Join all data
FullDataCoordWorkersRD1_RD2 <- full_join(FullDataCoordWorkers, FullDataCoordWorkersRD2)

# Combine all data sets into the final one for brood, colonies 1-10
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordBrood <- Colony1Binned %>%
  full_join(Colony2Binned) %>%
  full_join(Colony3Binned) %>%
  full_join(Colony4Binned) %>%
  full_join(Colony5Binned) %>% 
  full_join(Colony6Binned) %>% 
  full_join(Colony7Binned) %>%
  full_join(Colony8Binned) %>% 
  full_join(Colony9Binned) %>%
  full_join(Colony10Binned) %>%
  # Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Brood" & Colony < 11) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Density)

# Combine all data sets into the final one for brood, colonies 11-20
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordBroodRD2 <- Colony11Binned %>%
  full_join(Colony12Binned) %>%
  full_join(Colony13Binned) %>%
  full_join(Colony14Binned) %>%
  full_join(Colony15Binned) %>%
  full_join(Colony16Binned) %>%
  full_join(Colony17Binned) %>%
  full_join(Colony18Binned) %>%
  full_join(Colony19Binned) %>%
  full_join(Colony20Binned) %>%
  #Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Brood" & Colony > 10) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Density)

# Join all data
FullDataCoordBroodRD1_RD2 <- full_join(FullDataCoordBrood, FullDataCoordBroodRD2)

# Combine all data sets into the final one for queens, colonies 1-10
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordQueen <- Colony1Binned %>%
  full_join(Colony2Binned) %>%
  full_join(Colony3Binned) %>%
  full_join(Colony4Binned) %>%
  full_join(Colony5Binned) %>% 
  full_join(Colony6Binned) %>% 
  full_join(Colony7Binned) %>%
  full_join(Colony8Binned) %>% 
  full_join(Colony9Binned) %>%
  full_join(Colony10Binned) %>%
  #Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Queens" & Colony < 11) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Density)

# Combine all data sets into the final one for queens, colonies 11-20
# This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordQueenRD2 <- Colony11Binned %>%
  full_join(Colony12Binned) %>%
  full_join(Colony13Binned) %>%
  full_join(Colony14Binned) %>%
  full_join(Colony15Binned) %>%
  full_join(Colony16Binned) %>%
  full_join(Colony17Binned) %>%
  full_join(Colony18Binned) %>%
  full_join(Colony19Binned) %>%
  full_join(Colony20Binned) %>%
  # Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Queens" & Colony > 10) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Density)

# Join all data
FullDataCoordQueenRD1_RD2 <- full_join(FullDataCoordQueen, FullDataCoordQueenRD2)

# Combine all data sets into the final one for alates
# Note that alates are only found in the low density treatment so there is only one script to combine alates
FullDataCoordAlates <- Colony11Binned %>%
  full_join(Colony12Binned) %>%
  full_join(Colony13Binned) %>%
  full_join(Colony14Binned) %>%
  full_join(Colony15Binned) %>%
  full_join(Colony18Binned) %>%
  full_join(Colony19Binned) %>%
  # Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(ColonyMember == "Alates" & Colony > 10) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Sex, Ratio)

####################################################################################################################
# NETLOGO SIMULATION RESULTS BIN FUNCTION
# The below code is the same as for the experimental coordinates
# There is also a bin coordinate reference data set used here 
# Except for a renaming line that changes "x" and "y" columns to "ScaledX" and "ScaledY" in order to work with the function
# This also makes later joining real data and simulated results easier 
# There is a separate set of code for each combination of nest shape (Tube, Circle) and size (Small, Large)
####################################################################################################################

CoordBinnedNetlogo <- function(data_table){
  
  # Separating out the high nest density simulations
  BinCoordAssignSmall <- BinCoordNetlogo %>%
    filter(NestSize == "Small")
  
  # Separating out the low nest density simulations
  BinCoordAssignLarge <- BinCoordNetlogo %>%
    filter(NestSize == "Large")
  
  # Binning small tube nest coordinates 
  NetlogoBinnedTubeSmall <<- data_table %>%
    filter(NestSize == "Small" & Nest == "Tube") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[2] & 
                       ScaledX <= BinCoordAssignSmall$ScaledX[1] &
                       ScaledY <= BinCoordAssignSmall$ScaledY[2], 1,
                     # Bin 2
                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[3] & 
                               ScaledX <= BinCoordAssignSmall$ScaledX[2] &
                               ScaledY >= BinCoordAssignSmall$ScaledY[3] &
                               ScaledY <= BinCoordAssignSmall$ScaledY[2], 2,
                             # Bin 3
                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[4] & 
                                       ScaledX <= BinCoordAssignSmall$ScaledX[3] &
                                       ScaledY >= BinCoordAssignSmall$ScaledY[3] &
                                       ScaledY <= BinCoordAssignSmall$ScaledY[4], 3,
                                     # Bin 4
                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[4] & 
                                               ScaledX <= BinCoordAssignSmall$ScaledX[5] &
                                               ScaledY >= BinCoordAssignSmall$ScaledY[4] &
                                               ScaledY <= BinCoordAssignSmall$ScaledY[5], 4,
                                             # Bin 5
                                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[5] & 
                                                       ScaledX <= BinCoordAssignSmall$ScaledX[6] &
                                                       ScaledY >= BinCoordAssignSmall$ScaledY[4] &
                                                       ScaledY <= BinCoordAssignSmall$ScaledY[5], 5,
                                                     # Bin 6
                                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[6] & 
                                                               ScaledX <= BinCoordAssignSmall$ScaledX[7] &
                                                               ScaledY >= BinCoordAssignSmall$ScaledY[7] &
                                                               ScaledY <= BinCoordAssignSmall$ScaledY[6], 6,
                                                             # Bin 7
                                                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[8] & 
                                                                       ScaledX <= BinCoordAssignSmall$ScaledX[7] &
                                                                       ScaledY <= BinCoordAssignSmall$ScaledY[7], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[9] & 
                                                                               ScaledX <= BinCoordAssignSmall$ScaledX[8] &
                                                                               ScaledY <= BinCoordAssignSmall$ScaledY[9], 8, 0
                                                                             )))))))))
  # Binning large tube nest coordinates 
  NetlogoBinnedTubeLarge <<- data_table %>%
    filter(NestSize == "Large" & Nest == "Tube") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[2] & 
                       ScaledX <= BinCoordAssignLarge$ScaledX[1] &
                       ScaledY <= BinCoordAssignLarge$ScaledY[2], 1,
                     # Bin 2
                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[3] & 
                               ScaledX <= BinCoordAssignLarge$ScaledX[2] &
                               ScaledY >= BinCoordAssignLarge$ScaledY[3] &
                               ScaledY <= BinCoordAssignLarge$ScaledY[2], 2,
                             # Bin 3
                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[4] & 
                                       ScaledX <= BinCoordAssignLarge$ScaledX[3] &
                                       ScaledY >= BinCoordAssignLarge$ScaledY[3] &
                                       ScaledY <= BinCoordAssignLarge$ScaledY[4], 3,
                                     # Bin 4
                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[4] & 
                                               ScaledX <= BinCoordAssignLarge$ScaledX[5] &
                                               ScaledY >= BinCoordAssignLarge$ScaledY[4] &
                                               ScaledY <= BinCoordAssignLarge$ScaledY[5], 4,
                                             # Bin 5
                                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[5] & 
                                                       ScaledX <= BinCoordAssignLarge$ScaledX[6] &
                                                       ScaledY >= BinCoordAssignLarge$ScaledY[4] &
                                                       ScaledY <= BinCoordAssignLarge$ScaledY[5], 5,
                                                     # Bin 6
                                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[6] & 
                                                               ScaledX <= BinCoordAssignLarge$ScaledX[7] &
                                                               ScaledY >= BinCoordAssignLarge$ScaledY[7] &
                                                               ScaledY <= BinCoordAssignLarge$ScaledY[6], 6,
                                                             # Bin 7
                                                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[8] & 
                                                                       ScaledX <= BinCoordAssignLarge$ScaledX[7] &
                                                                       ScaledY <= BinCoordAssignLarge$ScaledY[7], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[9] & 
                                                                               ScaledX <= BinCoordAssignLarge$ScaledX[8] &
                                                                               ScaledY <= BinCoordAssignLarge$ScaledY[9], 8, 0
                                                                     )))))))))
  # Binning small circle nest coordinates 
  NetlogoBinnedCircleSmall <<- data_table%>%
    filter(NestSize == "Small" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssignSmall$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[10] & 
                               ScaledY <= BinCoordAssignSmall$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[11] & 
                                       ScaledY <= BinCoordAssignSmall$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[12] & 
                                               ScaledY <= BinCoordAssignSmall$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssignSmall$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssignSmall$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssignSmall$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[16], 8, 0
                                                                             )))))))))
  # Binning large circle nest coordinates 
  NetlogoBinnedCircleLarge <<- data_table %>%
    filter(NestSize == "Large" & Nest == "Circle") %>%
    mutate(Bin =
             # Bin 1
             if_else(ScaledY <= BinCoordAssignLarge$ScaledY[10], 1,
                     # Bin 2
                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[10] & 
                               ScaledY <= BinCoordAssignLarge$ScaledY[11], 2,
                             # Bin 3
                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[11] & 
                                       ScaledY <= BinCoordAssignLarge$ScaledY[12], 3,
                                     # Bin 4
                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[12] & 
                                               ScaledY <= BinCoordAssignLarge$ScaledY[13], 4,
                                             # Bin 5
                                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssignLarge$ScaledY[14], 5,
                                                     # Bin 6
                                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssignLarge$ScaledY[15], 6,
                                                             # Bin 7
                                                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssignLarge$ScaledY[16], 7,
                                                                     # Bin 8
                                                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[16], 8, 0
                                                                     )))))))))
  
  NetlogoBinnedFull <<- full_join(NetlogoBinnedTubeSmall, NetlogoBinnedTubeLarge) %>%
    full_join(NetlogoBinnedCircleSmall) %>% 
    full_join(NetlogoBinnedCircleLarge)
}

# Run the eight nest section binning function for the Netlogo simulation data set
CoordBinnedNetlogo(NetlogoTestFull) #NetlogoDataset

####################################################################################################################
# NEST SECTION DENSITY CALCULATIONS
# The following scripts calculate the proportion of colony members in each nest section (Column Bin: 1-8), calculated in the file Bins_Working.R 
# The script uses a null data set that's just Colony, Nest, and Bin called BinsNullFull, which just shows bin 1-8 for each Colony and nest combination
# The following also uses a reference data set CornerFull, which assignes the presence of corners (Y/N) to nest sections
####################################################################################################################

# First we create the proportions data set, then a null one 
# This approach allows zeros to be present in the proportions data sets, since no workers in the nest section is relevant data

# WORKERS
# High density treatment
Prop_functionWorker <- function(data.table) {
  AntProp <- data.table %>% # Creating the data set of worker proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of workers in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each worker in each bin in each observation
           PropWorker = (BinCount / count)) %>% # Calculate the proportion of workers in each bin 
    select(Colony, Day, Nest, Bin, PropWorker, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null worker proportions in each nest section
  AntPropNull <- AntProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(AntPropNull, BinsNullFull) %>% # Join the two null data sets
    drop_na()
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullWorkers <<- full_join(NestArchNullBins, AntProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(WorkerType = "Obsv", # Create a column filled with "Obsv" which identifies this data as the real one vs. the Netlogo simunlations 
           PropWorker = ifelse(is.na(PropWorker), 0, PropWorker),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropWorker)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropWorker, Density, WorkerType) %>% # Select only the desired columns
    left_join(CornerFull)  # Joins with a data set that assigned corner presence to each nest section
}

# Run the proportions of workers in nest sections function for the FullDataCoordWorkers data set 
Prop_functionWorker(FullDataCoordWorkers)

# WORKERS
# Low density treatment
Prop_functionWorker <- function(data.table) {
  AntProp <- data.table %>% # Creating the data set of worker proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of workers in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each worker in each bin in each observation
           PropWorker = (BinCount / count)) %>% # Calculate the proportion of workers in each bin 
    select(Colony, Day, Nest, Bin, PropWorker, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null worker proportions in each nest section
  AntPropNull <- AntProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(AntPropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullWorkersRD2 <<- full_join(NestArchNullBins, AntProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(WorkerType = "Obsv", # Create a column filled with "Obsv" which identifies this data as the real one vs. the Netlogo simunlations 
           PropWorker = ifelse(is.na(PropWorker), 0, PropWorker),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropWorker)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropWorker, Density, WorkerType) %>% # Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}

# Run the proportions of workers in nest sections function for the FullDataCoordWorkersRD2 data set 
Prop_functionWorker(FullDataCoordWorkersRD2)

#Join worker proportions in nest sections data sets
AntPropFullWorkersRD1_RD2 <- full_join(AntPropFullWorkers, AntPropFullWorkersRD2)

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
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% # Select only the desired columns
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
    select(Colony, Day, Nest, Bin, PropBrood, Density) %>% # Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}

# Run the proportions of brood in nest sections function for the FullDataCoordBroodRD2 data set 
Prop_functionBrood(FullDataCoordBroodRD2)

# Join brood proportions in nest sections data sets
AntPropFullBroodRD1_RD2 <- full_join(AntPropFullBrood, AntPropFullBroodRD2) 

# QUEENS
# High density treatment
Prop_functionQueen <- function(data.table) {
  QueenProp <- data.table %>% # Creating the data set of queen proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of queen in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each queen in each bin in each observation
           PropQueen = (BinCount / count)) %>% # Calculate the proportion of queen in each bin 
    select(Colony, Day, Nest, Bin, PropQueen, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null queen proportions in each nest section
  QueenPropNull <- QueenProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(QueenPropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullQueen <<- full_join(NestArchNullBins, QueenProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(PropQueen = ifelse(is.na(PropQueen), 0, PropQueen),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropQueen)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropQueen, Density) %>% # Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}
# Run the proportions of queens in nest sections function for the FullDataCoordQueen data set 
Prop_functionQueen(FullDataCoordQueen)

# QUEENS
# Low density treatment
Prop_functionQueen <- function(data.table) {
  QueenProp <- data.table %>% # Creating the data set of queen proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of queen in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each queen in each bin in each observation
           PropQueen = (BinCount / count)) %>% # Calculate the proportion of queen in each bin 
    select(Colony, Day, Nest, Bin, PropQueen, Density) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null queen proportions in each nest section
  QueenPropNull <- QueenProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day, Density)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(QueenPropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullQueenRD2 <<- full_join(NestArchNullBins, QueenProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(PropQueen = ifelse(is.na(PropQueen), 0, PropQueen),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropQueen)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropQueen, Density) %>% #Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}

# Run the proportions of queens in nest sections function for the FullDataCoordQueenRD2 data set 
Prop_functionQueen(FullDataCoordQueenRD2)

# Join queen proportions in nest sections data sets
AntPropFullQueenRD1_RD2 <- full_join(AntPropFullQueen, AntPropFullQueenRD2)

# ALATES
# Only the high density treatment
Prop_functionAlate <- function(data.table) {
  AlateProp <- data.table %>% # Creating the data set of alate proportions in each nest section
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, Day
    mutate(count = n()) %>% # Count total number of alates in each observation
    group_by(Colony, Nest, Day, Bin) %>% # Group by columns Colony, Nest, Day, Bin
    mutate(BinCount = n(), # Counting the number of each alate in each bin in each observation
           PropAlate = (BinCount / count)) %>% # Calculate the proportion of alates in each bin 
    select(Colony, Day, Nest, Bin, PropAlate) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null alate proportions in each nest section
  AlatePropNull <- AlateProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(Colony, Nest, Day)) %>% # Select the desired columns
    drop_na() %>% # Remove any NAs
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(AlatePropNull, BinsNullFull) # Join the two null data sets
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullAlateRD2 <<- full_join(NestArchNullBins, AlateProp) %>%  
    group_by(Colony, Nest, Day) %>% # Group by columns Colony, Nest, and Day
    mutate(PropAlate = ifelse(is.na(PropAlate), 0, PropAlate),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropAlate)) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(Colony, Day, Nest, Bin, PropAlate) %>% # Select only the desired columns
    left_join(CornerFull) # Joins with a data set that assigned corner presence to each nest section
}

# Run the Prop_functionAlate function for the FullDataCoordAlates data set 
Prop_functionAlate(FullDataCoordAlates)

####################################################################################################################
# NETLOGO SIMULATION: NEST SECTION DENSITY CALCULATIONS
# The following scripts calculate the proportion of Netlogo simulated agents (NetlogoBinnedFull) in each nest section (Column Bin: 1-8), calculated in the file Bins_Working.R 
# The script uses a null data set that's just Colony, Nest, and Bin called BinsNullNetlogo, which just shows bin 1-8 for each Colony and nest combination
# The following also uses a reference data set CornerFull, which assignes the presence of corners (Y/N) to nest sections
####################################################################################################################

# NETLOGO SIMULATIONS
# All data at once, not split into density treatments
Prop_functionWorkerSim <- function(data.table) {
  NetlogoProp <- data.table %>% # Creating the data set of simulates result proportions in each nest section
    group_by(RunNumber, Nest, NestSize) %>% # Group by columns RunNumber, Nest, NestSize
    mutate(count = n()) %>% # Count total number of simulates result in each observation
    group_by(RunNumber, Nest, NestSize, Bin) %>% # Group by columns RunNumber, Nest, NestSize, Bin
    mutate(BinCount = n(), # Counting the number of each simulates result in each bin in each observation
           PropWorker = (BinCount / count)) %>% # Calculate the proportion of simulated results in each bin 
    select(NestSize, Nest, Bin, RunNumber, PropWorker, TimeStep) %>% # Select only the desired columns
    distinct() # Remove duplicate rows
  # Creating the data set of null simulates result proportions in each nest section
  NetlogoPropNull <- NetlogoProp %>%  
    ungroup() %>% # Ungroup the data set
    select(c(RunNumber, Nest, NestSize)) %>% # Select the desired columns
    distinct() # Remove duplicate rows
  NestArchNullBins <- full_join(NetlogoPropNull, BinsNullNetlogo) # Join the two null data sets
  
  # Joining the working data set to the null one, which keeps the zeros in the final data set
  AntPropFullSim <<- full_join(NestArchNullBins, NetlogoProp) %>%  
    group_by(RunNumber, Nest, NestSize) %>% # Group by columns Colony, Nest, and Day
    mutate(PropWorker = ifelse(is.na(PropWorker), 0, PropWorker),# NAs are produced in the join above, this makes them zeros
           Binsum = sum(PropWorker),
           WorkerType = "RandSim",
           Density = ifelse(NestSize == "Small", "High", "Low")) %>% # Create a column that sums the proportions
    filter(Binsum != 0) %>% # Removes any rows with zeros from the Binsum column. This is only a precaution 
    select(RunNumber, NestSize, Nest, Bin, PropWorker, WorkerType, TimeStep, Density) %>% # Select only the desired columns
    left_join(CornerFullSim) %>% # Joins with a data set that assigned corner presence to each nest section
    drop_na() # Remove any NAs, also a precaution only
}

# Run the Prop_functionWorkerSim function for the NetlogoBinnedFull data set 
Prop_functionWorkerSim(NetlogoBinnedFull)
