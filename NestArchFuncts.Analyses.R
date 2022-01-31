####################################################################################################################
## Author: GREG CHISM
## Date: Dec 2021
## email: gchism@email.arizona.edu
## Project: Nest shape influences colony organization in ants
## Title: Calculating colony member proportions in nest sections, importing and preparing Netlogo simulation data, all paper analyses and figures 
####################################################################################################################

### IMPORTANT ### 
## BEFORE RUNNING THIS CODE, PLEASE RUN ALL CODE FROM THE R SCRIPTS "Bins_Working.R", "DistanceFunctions.R", "DistToBrood.R", AND FINALLY THIS R SCRIPT 
# This code is to replicate the analyses and figures for the following in my first chapter:
# Density of colony members in nest sections
# Netlogo simulation & worker comparisons
# Density figures

install.packages("pacman") # Download package with function to load multiple packaged at once
####################################################################################################################
# Loading required packages for code below  
# p_load() will download packages that aren't in system library
####################################################################################################################

pacman::p_load(extrafonts, 
               forcats,
               ggpubr,
               Kendall,
               kuiper.2samp,
               lme4,
               lmerTest,
               magick,
               magrittr,
               MuMIn,
               RColorBrewer,
               tidyverse,
               wesanderson,
               assertthat,
               twosamples,
               RColorBrewer,
               ggpointdensity,
               readxl,
               scales,
               viridis)

# Importing fonts for plots
font_import()

# Yes, import all fonts (takes a few minutes)
y

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
FullDataCoordAlates <- read.csv("FullDataCoordAlates.csv")

# NETLOGO 
ArchitectureMoveModelFull <- read.csv("ArchitectureMoveModelFull.csv")

# NULL BINS REFERENCE (EMPIRICAL)
BinsNullFull <- read.csv("BinsNullFull.csv")

# NULL BINS REFERENCE (NETLOGO SIMULATIONS)
BinsNullNetlogo <- read.csv("BinsNullNetlogo.csv")

# REFERENCE COORDINATES FOR CORNERS (EMPIRICAL)
CornerFull <- read.csv("CornerFull.csv")

# REFERENCE COORDINATES FOR CORNERS (NETLOGO SIMULATIONS)
CornerFullSim <- read.csv("CornerFullSim.csv")

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
  select(c(RunNumber, NestSize, Nest, MovementRule, TimeStep))

# Keeping only X and Y coordinates
# X coordinate
ArchitectureMoveModelFullCorrectRedXcor <- ArchitectureMoveModelFullCorrectRed %>%
  select(-c(ycor))

# Y coordinate
ArchitectureMoveModelFullCorrectRedYcor <- ArchitectureMoveModelFullCorrectRed %>%
  select(-c(xcor))


# Splitting the xcor and ycor lists
XCoordCorrect <- strsplit(ArchitectureMoveModelFullCorrectRedXcor$xcor, split = " ") #xcor
YCoordCorrect <- strsplit(ArchitectureMoveModelFullCorrectRedYcor$ycor, split = " ") #ycor

# Create a column that assigns unique IDs to each xcor and ycor row within each simulation run
# X coordinate
NetlogoTestX <- data.frame(RunNumber = rep(ArchitectureMoveModelFullCorrectRedXcor$RunNumber, sapply(XCoordCorrect, length)), xcor = unlist(XCoordCorrect)) %>%
  rowid_to_column(var = "id")

# Y coordinate
NetlogoTestY <- data.frame(RunNumber = rep(ArchitectureMoveModelFullCorrectRedYcor$RunNumber, sapply(YCoordCorrect, length)), ycor = unlist(YCoordCorrect)) %>%
  rowid_to_column(var = "id")


# Creating the working Netlogo data sets
NetlogoTestFull <- left_join(NetlogoTestX,NetlogoTestY) %>%
  mutate(xcor = as.numeric(xcor),
         ycor = as.numeric(ycor)) %>%
  select(-c(id)) %>%
  left_join(ArchitectureMoveModelFull_GC_30_Sept_2021_TestExp_TableRed) %>%
  mutate(MovementRule = "Random",
         xcor = xcor * 0.1,
         ycor = ycor * 0.1) %>% 
rename(ScaledX = xcor, ScaledY = ycor)

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
    select(c(Colony, Nest, Day)) %>% # Select the desired columns
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
    select(c(Colony, Nest, Day)) %>% # Select the desired columns
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

####################################################################################################################
# SUBSETTING / COMPARING EMPIRICAL & NETLOGO SIMULATION DISTRIBUTIONS
# The below script separates the empirical worker and Netlogo simulated data sets so that each combination of nest, density treatment can be pairwise compared
# Distributions are compared using Cramer von mises tests, which compare points across entire distributions
# The script then determines which nest section holds the highest proportions of both empirical and Netlogo simulated workers 
####################################################################################################################

# SEPARATING NEST SHAPES AND DENSITIES
# NETLOGO SIMULATIONS
# High nest density
WorkerSim1PropAnalysisSmall <- AntPropFullSim %>% filter(NestSize == "Small")
# Tube nest shape
WorkerSim1TubePropAnalysis <- WorkerSim1PropAnalysisSmall%>%filter(Nest=="Tube")
# Circle nest shape
WorkerSim1CirclePropAnalysis<-WorkerSim1PropAnalysisSmall%>%filter(Nest=="Circle")

# Low nest density
WorkerSim2PropAnalysisLarge <- AntPropFullSim %>% filter(NestSize == "Large")
# Tube nest shape
WorkerSim2TubePropAnalysisRD2<-WorkerSim2PropAnalysisLarge%>%filter(Nest=="Tube")
# Circle nest shape
WorkerSim2CirclePropAnalysisRD2<-WorkerSim2PropAnalysisLarge%>%filter(Nest=="Circle")

# EMPIRICAL WORKERS (already separated into high and low nest densities - RD2 are low density)
# High nest density
# Tube nest shape
WorkerTubePropAnalysis <- AntPropFullWorkers %>% filter(Nest == "Tube")
# Circle nest shape
WorkerCirclePropAnalysis <- AntPropFullWorkers %>% filter(Nest == "Circle")

# Low nest density
# Tube nest shape
WorkerTubePropAnalysisRD2 <- AntPropFullWorkersRD2 %>% filter(Nest == "Tube")
# Circle nest shape
WorkerCirclePropAnalysisRD2 <- AntPropFullWorkersRD2 %>% filter(Nest == "Circle")

# SIMULATION & EXPERIMENTAL DISTRIBUTION COMPARISIONS - CRAMER VON MISES TESTS
# High density tube nest comparison
cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)
# Low density tube nest comparison
cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)
# High density circle nest comparison
cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)
# Low density circle nest comparison
cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)
# High density tube nest empirical workers v. low density tube nest simulations
cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)
# Low density tube nest empirical workers v. high density tube nest simulations
cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)
# High density circle nest empirical workers v. low density circle nest simulations
cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)
# Low density circle nest empirical workers v. high density circle nest simulations
cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)

# COMPARING SIMULATIONS
# High density vs low density tube nest
# Note that because this statistic is calculated through bootstrap resamplings there is a slightlu different p-value each time, but it is always > 0.95
cvm_test(WorkerSim1TubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)
# High density vs low density circle nest
cvm_test(WorkerSim1CirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)

# DETERMINING MAX WORKER AND NETLOGO SIMULATION PROPORTIONS IN NEST SECTIONS
# WORKERS 
MaxWorkerProp <- AntPropFullWorkersRD1_RD2 %>% #Full worker proportions data set
  group_by(Colony, Nest, Day) %>% # Group by Colony, Nest, and Day columns
  # Create a column with the max worker proportion and another confirming showing what the associated bin is
  mutate(MaxPropWorker = max(PropWorker), MaxBin = ifelse(MaxPropWorker == PropWorker, 1, 0)) %>%  
  # Remove unwanted columns
  select(-c(PropWorker, MaxPropWorker)) 

# NETLOGO SIMULATIONS
# Same as for empirical workers but with the Netlogo proportions data set
MaxWorkerPropSim <- AntPropFullSim %>% #Full Netlogo simulations proportions data set
  group_by(RunNumber, NestSize, Nest) %>%
  mutate(MaxPropWorker = max(PropWorker), MaxBin = ifelse(MaxPropWorker == PropWorker, 1, 0)) %>%
  select(-c(PropWorker, MaxPropWorker)) %>%
  mutate(Density = ifelse(NestSize == "Small", "High", "Low"))

# Joining the max proportions data sets 
MaxPropWorkersFull <- full_join(MaxWorkerProp, MaxWorkerPropSim) 

####################################################################################################################
# PLOTS AND ANALYSES: Colony member densities through the nest
# The scripts below are to analyze and visualize: 
# Worker, brood, queen, and alate densities through the nest (including comparing workers and Netlogo simulations)
####################################################################################################################

# WORKER AND NETLOGO SIMULATION DENSITIES IN NEST SECTIONS

# Creating the working data set for analysis and visualization
WorkerSimPropFullRD1 <- full_join(AntPropFullWorkers, WorkerSim1PropAnalysisSmall)
WorkerSimPropFullRD2 <- full_join(AntPropFullWorkersRD2, WorkerSim2PropAnalysisLarge)

# All proportions
AntPropFullWorkersSimRD1_RD2 <- full_join(WorkerSimPropFullRD1, WorkerSimPropFullRD2)# Join all worker proportions data sets

# BOXPLOTS 
# TUBE NEST
# High nest density
WorkerPlotRD1Tube <- ggplot(data = WorkerSimPropFullRD1 %>% filter(Nest == "Tube"),
                          aes(x = as.factor(Bin), y = PropWorker, 
                              fill = WorkerType)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      position = position_dodge(1),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab("Nest section") + 
  ylab("Proportions of workers") +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.87, vjust = -10),
        legend.position = "none") +
  scale_fill_manual(values = c("red", "white"), 
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.8)

# Low nest density
WorkerPlotRD2Tube <- ggplot(data = WorkerSimPropFullRD2 %>% filter(Nest == "Tube"),
                          aes(x = as.factor(Bin), y = PropWorker, 
                              fill = WorkerType)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      position = position_dodge(1),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab("Nest section") + 
  ylab("Proportions of workers") +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Worker type", family = "Arial", color = "black")) +
  scale_fill_manual(values = c("red", "white"),
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.8)

# Compiling the empirical worker density in nest sections box plots 
WorkerPropPlot <- ggarrange(WorkerPlotRD1Tube, WorkerPlotRD2Tube,
                            labels = c("(a)", "(b)"),
                            font.label = list(size = 18, family = "Arial", face = "plain"),
                            label.x = 0.9,
                            label.y = 1,
                            ncol = 2, nrow = 1,
                            common.legend = FALSE)

# Annotating the compiled tube nest plot to include a title
WorkerPropPlotAnnot <-annotate_figure(WorkerPropPlot,
                                      top = text_grob("Tube nest", color = "black",
                                                      size = 18, x = 0.08, y = -0.8, family = "Arial"),
                                      bottom = NULL,
                                      left = NULL,
                                      right = NULL, 
                                      fig.lab.size = 18
)

# CIRCLE NEST
# High density treatment
WorkerPlotRD1Circle <- ggplot(data = WorkerSimPropFullRD1 %>% filter(Nest == "Circle"),
                            aes(x = as.factor(Bin), y = PropWorker, 
                                fill = WorkerType)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      position = position_dodge(1),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab("Nest section") + 
  ylab("Proportions of workers") +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "white", hjust = 0.875, vjust = -10),
        legend.position = "none") +
  scale_fill_manual(values = c("blue", "white"),
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.5)

#Low density treatment
WorkerPlotRD2Circle <- ggplot(data = WorkerSimPropFullRD2 %>% filter(Nest == "Circle"),
                            aes(x = as.factor(Bin), y = PropWorker, 
                                fill = WorkerType)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      position = position_dodge(1),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab("Nest section") + 
  ylab("Proportions of workers") +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, color = "white", family = "Arial"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "white", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.75),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Worker type", family = "Arial", color = "black")) +
  scale_fill_manual(values = c("blue", "white"),
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.5)

# Compiling the two above circle nest box plots 
WorkerPropPlot2 <- ggarrange(WorkerPlotRD1Circle, WorkerPlotRD2Circle,
                           labels = c("(c)", "(d)"),
                           font.label = list(size = 18, family = "Arial", color = "black", face = "plain"),
                           label.x = 0.9,
                           label.y = 1,
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled circle nest plot to include a title
WorkerPropPlotAnnot2 <- annotate_figure(WorkerPropPlot2,
                                        top = text_grob("Circle nest", color = "black",
                                                        size = 18, x = 0.08, y = -0.8, family = "Arial"),
                                        bottom = NULL,
                                        left = NULL,
                                        right = NULL
)

# Compiling the empirical and Netlogo simulated worker densities in nest sections box plots 
WorkerPropPlotAnnotFull <- ggarrange(WorkerPropPlotAnnot, WorkerPropPlotAnnot2,
                             ncol = 1, nrow = 2,
                             common.legend = FALSE)

# Annotating the compiled plot to include shared axes lables
annotate_figure(WorkerPropPlotAnnotFull,
                         top = NULL,
                         bottom = text_grob("Nest section", color = "black",
                                            size = 18, x = 0.5, family = "Arial"),
                         left = text_grob("Proportions of workers", color = "black",
                                          size = 18, rot = 90, family = "Arial"),
                         right = NULL)

# Linear regression: observed & simulated workers
# RESPONSE VARIABLE
# PropWorker - Proportion of either empirical or Netlogo simulated workers 
# EFFECTS
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# WorkerType - Empirical or Netlogo random walk simulated worker (Worker / RandSim)
# Corner - Presence of a corner in the nest section (Y / N)
summary(lm(PropWorker ~ poly(Bin, degree = 2, raw = TRUE) * Nest * WorkerType * Density + Corner, data = AntPropFullWorkersSimRD1_RD2))

# LINEAR MIXED EFFECTS MODEL: Observed worker density through the nest
# RESPONSE VARIABLE
# PropWorker - Proportion of total brood found in each nest section
# FIXED EFFECTS 
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(PropWorker ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullWorkersRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(PropWorker ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullWorkersRD1_RD2))

# NEST SECTION THAT HOLDS THE MAXIMUM EMPIRICAL AND NETLOGO SIMULATED WORKER PROPORTION
# GLM LINE PLOTS 
# High nest density
MaxLog1 <- MaxPropWorkersFull %>% 
  filter(WorkerType == "Obsv") %>%
  arrange(Nest) %>%
  ggplot(aes(Bin, MaxBin)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(color = Nest, linetype = Density), 
              se = FALSE) +
  theme_pubclean() + 
  ggtitle("Observed") +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = -0.1, vjust = -10),
        legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle","Tube")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  ylim(0, 1)

# Low nest density
MaxLog2 <- MaxPropWorkersFull %>% 
  filter(WorkerType == "RandSim") %>%
  arrange(Nest) %>%
  ggplot(aes(Bin, MaxBin)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              aes(color = Nest, linetype = Density), 
              se = FALSE) +
  theme_pubclean() +
  ggtitle("Random walk") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = -0.1, vjust = -10),
        legend.key = element_blank(),
        legend.position = c(0.975, 0.95),
        legend.direction = "horizontal",
        legend.justification = c(1, 1),
        legend.text = element_text(size = 18, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.key.size = unit(1, 'cm')) +
  xlab(NULL) +
  ylab(NULL) +
  guides(color = guide_legend(title = "Nest", family = "Arial",
                              order = 1),
         linetype = guide_legend(order = 2)) +
  scale_color_manual(breaks = c("Tube", "Circle"), 
                     name = "Nest",
                     values = c("red", "blue")) +
  guides(lty = guide_legend(override.aes = list(col = 'black'))) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
  ylim(0, 1)

# Compiling the above GLM line plots
MaxLogPlot <- ggarrange(MaxLog1, MaxLog2,
                        labels = c("(a)", "(b)"),
                        label.x = 0.9,
                        font.label = list(size = 18, family = "Arial", face = "plain"),
                        ncol = 2, nrow = 1,
                        common.legend = FALSE)

# Annotating the compiled GLM line plots to include shared axes labels
annotate_figure(MaxLogPlot,
                                  top = NULL,
                                  bottom = text_grob("Nest section", color = "black",
                                                     size = 18, family = "Arial"),
                                  left = text_grob("Max worker proportion", color = "black",
                                                   size = 18, rot = 90, family = "Arial"),
                                  right = NULL
)

# GENERALIZED LINEAR REGRESSION (BINOMIAL FIT): Max worker and Netlogo simulation proportions in nest sections
# RESPONSE VARIABLE
# MaxBin - Successes / failures dictated by whether the nest section held the highest proportion of empirical or Netlogo simulated workers (0/1)
# EFFECTS
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# WorkerType - Empirical or Netlogo random walk simulated worker (Worker / RandSim)
summary(glm(MaxBin ~ poly(Bin, degree = 2, raw = TRUE) * WorkerType * Nest * Density, family = binomial, data = MaxPropWorkersFull))

# BROOD DENSITIES IN NEST SECTIONS
# BOXPLOTS
# High density treatment
BroodProp1 <- ggplot(data = AntPropFullBrood %>% arrange(Nest), 
                     aes(x = as.factor(Bin), y = PropBrood)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.87, vjust = -10),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

# Low density treatment
BroodProp2 <- ggplot(data = AntPropFullBroodRD2 %>% arrange(Nest), 
                     aes(x = as.factor(Bin), y = PropBrood)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest),
                      color = "grey25", 
                      alpha = 0.65) + 
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)


# Compiling the two above brood densities in nest sections boxplots 
BroodPropPlot <- ggarrange(BroodProp1, BroodProp2,
                           labels = c("(a)", "(b)"),
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           label.x = 0.9,
                           label.y = 1,
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled boxplots to include shared axes lables
BroodPropPlotFull <- annotate_figure(BroodPropPlot,
                                   top = text_grob("Brood", color = "black", size = 18, x = 0.055, y = -0.6, family = "Arial"),
                                   bottom = NULL,
                                   left = text_grob("Proportions of brood", color = "black",
                                                    size = 18, rot = 90, family = "Arial"),
                                   right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Brood density through the nest
# RESPONSE VARIABLE
# PropBrood - Proportion of total brood found in each nest section
# FIXED EFFECTS 
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(PropBrood ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullBroodRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(PropBrood ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullBroodRD1_RD2))

# QUEEN DENSITIES IN NEST SECTIONS
# BOXPLOTS
# High density treatment
QueenProp1 <- ggplot(data = AntPropFullQueen %>% arrange(Nest), 
                   aes(x = as.factor(Bin), y = PropQueen)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

#Low density treatment
QueenProp2 <- ggplot(data = AntPropFullQueenRD2 %>% arrange(Nest), 
                   aes(x = as.factor(Bin), y = PropQueen)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest),
                      color = "grey25", 
                      alpha = 0.65) +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        axis.text.y = element_text(size=18, color = "white", family="Arial"),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

# Compiling the two above queen densities in nest sections boxplots 
QueenPropPlot <- ggarrange(QueenProp1, QueenProp2,
                           labels = c("(c)", "(d)"),
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           label.x = 0.9,
                           label.y = 1.045,
                           ncol = 2, nrow = 1)

# Annotating the compiled boxplots to include common axes labels
QueenPropPlotFull <- annotate_figure(QueenPropPlot,
                top = text_grob("Queens", color = "black", size = 18, x = 0.0625, y = 0, family = "Arial"),
                bottom = text_grob("Nest section", color = "black",
                                   size = 18, x = 0.5, family = "Arial"),
                left = text_grob("Proportions of queens", color = "black",
                                 size = 18, rot = 90, family = "Arial"),
                right = NULL
)

# Compiling the brood and queen densities in nest sections plots
ggarrange(BroodPropPlotFull, QueenPropPlotFull,
          ncol = 1, nrow = 2,
          common.legend = TRUE)

# LINEAR MIXED EFFECTS MODEL: Queen density through the nest
# RESPONSE VARIABLE
# PropWorker - Proportion of total brood found in each nest section
# FIXED EFFECTS 
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(PropQueen ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullQueenRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(PropQueen ~ poly(Bin, degree = 2, raw = TRUE) * Nest * Density + Day + Corner + (1|Colony), data = AntPropFullQueenRD1_RD2))

# ALATE DENSITIES IN NEST SECTIONS
# BOXPLOTS 
# Note that alates were only in the low nest density
AlatePropFig <- ggplot(data = AntPropFullAlateRD2 %>% arrange(Nest), 
       aes(x = as.factor(Bin), y = PropAlate)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +
  xlab("Nest section") + 
  ylab("Proportions of alates") +
  theme_pubclean() +
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "black"),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

# Arranging the above alate densities in nest sections boxplot 
AlatePropPlot <- ggarrange(AlatePropFig,
                           labels = c("(e)"),
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           label.x = 0.9,
                           label.y = 1.0325,
                           ncol = 1, nrow = 1)

# Annotating the alate proportions boxplot to include a title
annotate_figure(AlatePropPlot,
                                   top = text_grob("Alates", color = "black", size = 18, x = 0.15, y = 0, family = "Arial"),
                                   bottom = NULL,
                                   left = NULL,
                                   right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Observed worker density through the nest
# RESPONSE VARIABLE
# PropWorker - Proportion of total brood found in each nest section
# FIXED EFFECTS 
# Bin - Nest section, transformed to raw polynomial term because of a priori assumptions 
# Nest - Nest shape (Tube / Circle)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification  
summary(lmer(PropAlate ~  poly(Bin, degree = 2, raw = TRUE) * Nest + Day + Corner + (1|Colony), data = AntPropFullAlateRD2))

ggplot(BroodCentDistWorkersSFZ %>% filter(Nest == "Tube"), aes(MeanToBrood, Occur)) +
  geom_smooth(method = "lm")

# Marginal and conditional R-squared values for the linear mixed effects model above
r.squaredGLMM(lmer(PropAlate ~ poly(Bin, degree = 2, raw = TRUE) * Nest + Day + Corner + (1|Colony), data = AntPropFullAlateRD2))

####################################################################################################################
# PLOTS AND ANALYSES: Colony member distances from the nest entrance
# The scripts below are to analyze and visualize: 
# Worker, brood, queen, and alate distances from the nest entrance (including comparing workers and Netlogo simulations)
####################################################################################################################

# WORKER AND NETLOGO SIMULATION SCALED DISTANCES FROM THE NEST ENTRANCE
# WORKER SCALED DISTANCES FROM THE NEST ENTRANCE
# HISTOGRAMS
# High density treatment
WorkerDist1 <- ggplot(WorkerDistScaled %>% arrange(Nest), 
                      aes(ScaledDist, 
                          fill = Nest)) + 
  geom_histogram(position = "identity",
                 alpha = 0.7,
                 binwidth = 0.0416666) +
  ggtitle("High density") +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 2500)

# Low density treatment
WorkerDist2 <- ggplot(WorkerDistScaledRD2 %>% arrange(Nest), 
                      aes(ScaledDist,
                          fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("Low density") + 
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 2500)

# Compiling the above empirical worker distances to the nest entrance histograms
WorkerDistPlot <- ggarrange(WorkerDist1, WorkerDist2,
                          labels = c("(a)", "(b)"),
                          label.x = 0.9,
                          font.label = list(size = 18, family = "Arial", face = "plain"),
                          ncol = 2, nrow = 1,
                          common.legend = FALSE)

# Annotating the compiled plot to add a title
WorkerDistPlotFull <- annotate_figure(WorkerDistPlot,
                top = text_grob("Workers", color = "black",
                                size = 18, x = 0.055, y = -0.6, family = "Arial"),
                bottom = NULL,
                left =  NULL,
                right = NULL
)

# NETLOGO SIMULATION SCALED DISTANCES FROM THE NEST ENTRANCE
# High density treatment
SimDist1 <- ggplot(SimDistScaled %>% filter(NestSize == "Small") %>% arrange(Nest), 
                      aes(ScaledDist, 
                          fill = Nest)) + 
  geom_histogram(position = "identity",
                 alpha = 0.7,
                 binwidth = 0.0416666) +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlim(0, 1) +
  ylim(0, 30000)

# Low density treatment
SimDist2 <- ggplot(SimDistScaled %>% filter(NestSize == "Large") %>% arrange(Nest), 
                   aes(ScaledDist, 
                       fill = Nest)) + 
  geom_histogram(position = "identity",
                 alpha = 0.7,
                 binwidth = 0.0416666) +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlim(0, 1) +
  ylim(0, 30000)

# Compiling the above Netlogo simulated worker distances to the nest entrance histograms
WorkerSimPlot <- ggarrange(SimDist1, SimDist2,
                          labels = c("(c)", "(d)"),
                          label.x = 0.9,
                          label.y = 0.99,
                          font.label = list(size = 18, family = "Arial", face = "plain"),
                          ncol = 2, nrow = 1)

# Annotating the compiled plot to add a title
WorkerSimPlotFull <- annotate_figure(WorkerSimPlot,
                top = text_grob("Random walk", color = "black",
                                size = 18, x = 0.089, y = -1, family = "Arial"),
                bottom = NULL,
                left =  NULL,
                right = NULL
)

# Compiling the worker and Netlogo simulation histograms
FullDistPlot<-ggarrange(WorkerDistPlotFull, WorkerSimPlotFull,
                        ncol = 1, nrow = 2,
                        common.legend = FALSE)

# Annotating the compiled histograms to include common axes labels 
annotate_figure(FullDistPlot,
                top = NULL,
                bottom = text_grob("Scaled distance to nest entrance", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = text_grob("Obsv / sim worker count", color = "black",
                                 size = 18,  rot = 90, family = "Arial"),
                right = NULL
)

# Combining the empirical and Netlogo simulated worker distances from the nest entrance data sets
AllDistScaledRD1_RD2 <- SimDistScaled %>%
  mutate(Density = ifelse(NestSize == "Small", "High", "Low")) %>%
  full_join(WorkerDistScaledRD1_RD2)

# LINEAR REGRESSION: Empirical and Netlogo simulated worker scaled distances from the nest entrance
# RESPONSE VARIABLE 
# ScaledDist - Worker scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# WorkerType - Empirical or Netlogo random walk simulated worker (Worker / RandSim)
# Density - Nest density (High / Low)
# Corner - Presence of a corner in the nest section (Y / N)
summary(lm(ScaledDist ~ Nest * WorkerType + Density + Corner, AllDistScaledRD1_RD2))

# LINEAR MIXED EFFECTS MODEL: Empirical worker distances from the nest entrance
# RESPONSE VARIABLE
# ScaledDist - Worker scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1 | Colony), data = WorkerDistScaledRD1_RD2))

#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1 | Colony), data = WorkerDistScaledRD1_RD2))

# BROOD SCALED DISTANCES FROM THE NEST ENTRANCE
# HISTOGRAMS
# High density treatment
BroodDist1 <- ggplot(BroodDistScaled %>% arrange(Nest),
                   aes(ScaledDist, 
                       fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +  
  ggtitle("High density") +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 8000) +
  xlim(0, 1)

# Low density treatment
BroodDist2 <- ggplot(BroodDistScaledRD2 %>% arrange(Nest), 
                     aes(ScaledDist, 
                         fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("Low density") +
  labs(color="Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.685),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 8000) +
  xlim(0, 1)

# Compiling the brood distance to the nest entrance histograms
BroodDistPlot <- ggarrange(BroodDist1, BroodDist2,
                           labels = c("(a)", "(b)"),
                           label.x = 0.9,
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled histograms to include a common x axis labels and  title
BroodFullDist <- annotate_figure(BroodDistPlot,
                                 top = text_grob("Brood", color = "black",
                                                 size = 18, x = 0.055, y = -0.6, family = "Arial"),
                                 bottom = NULL,
                                 left =  text_grob("Brood count", color = "black",
                                                   size = 18, x = 0.525, rot = 90, family = "Arial"),
                                 right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Brood scaled distances from the nest entrance
# RESPONSE VARIABLE
# ScaledDist - Brood scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodDistScaledRD1_RD2))

#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodDistScaledRD1_RD2))

# QUEEN SCALED DISTANCES FROM THE NEST ENTRANCE
# HISTOGRAMS
# High density treatment
QueenDist1 <- ggplot(QueenDistScaled %>% arrange(Nest),
                     aes(ScaledDist, 
                         fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 150) +
  xlim(0, 1)

# Low density treatment
QueenDist2 <- ggplot(QueenDistScaledRD2 %>% arrange(Nest),
                     aes(ScaledDist, 
                         fill= Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 150) +
  xlim(0, 1)

# Compiling the queens distance to the nest entrance histograms
QueenDistPlot <- ggarrange(QueenDist1, QueenDist2,
                           labels = c("(c)", "(d)"),
                           label.x = 0.9,
                           label.y = 0.99,
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           ncol = 2, nrow = 1)

# Annotating the compiled histograms to include a common x axis labels and  title
QueenFullDist <- annotate_figure(QueenDistPlot,
                               top = text_grob("Queens", color = "black",
                                               size = 18, x = 0.06, y = -1, family = "Arial"),
                               bottom = NULL,
                               left = text_grob("Queen count", color = "black",
                                                         size = 18, x = 0.525, rot = 90, family = "Arial"),
                               right = NULL
)

# Compiling the brood and queen scaled distance to the nest entrance plots
BroodQueenDist <- ggarrange(BroodFullDist, QueenFullDist,
          ncol = 1, nrow = 2)


annotate_figure(BroodQueenDist,
                top = NULL,
                bottom = text_grob("Scaled distance to nest entrance", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Queen scaled distances from the nest entrance
# RESPONSE VARIABLE
# ScaledDist - Queen scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Corner - Presence of a corner in the nest section (Y / N)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ScaledDist ~ Nest * Density + Corner + Day + (1 | Colony), data = QueenDistScaledRD1_RD2))

#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1 | Colony), data = QueenDistScaledRD1_RD2))

# ALATE SCALED DISTANCES TO THE NEST ENTRANCE
# Removing individuals with unknown alate sex
# In these plots and analyses, we care about the individuals with identifiable sex
# 10 individuals are removing the unknown sex individuals, representing < 1% of the data set  
AlateDistScaledRD2Plot <- AlateDistScaledRD2 %>%
  filter(Sex != "?")

# HISTOGRAM
AlateDist1 <- ggplot(AlateDistScaledRD2Plot %>% arrange(Nest), 
                     aes(ScaledDist, 
                         fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  xlab("Scaled distance to nest entrance") +
  ylab("Alate count") +
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 110)

# BOXPLOT SHOWING THE RELATIONSHIP BETWEEN ALATE SEX AND SCALED DISTANCE TO THE NEST ENTRANCE 
AlateDist2 <- ggplot(AlateDistScaledRD2Plot %>% arrange(Nest),
                     aes(x = Sex, y = ScaledDist, 
                         fill = Nest)) +
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) + 
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlab("Alate sex") +
  ylab("Scaled distance to nest entrance")

# Compiling the alate distance to the entrance plots 
AlateDistPlot <- ggarrange(AlateDist1, AlateDist2,
                           labels = c("(e)", "(f)"),
                           label.x = 0.9,
                           label.y = 0.965, 
                           font.label = list(size = 18, face = "plain", family = "Arial"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled plots to include a title
annotate_figure(AlateDistPlot,
                top = text_grob("Alates", color = "black",
                                size = 18, x = 0.06, y = -1.30, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Alate scaled distances from the nest entrance
# RESPONSE VARIABLE
# ScaledDist - Alate scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Sex - Alate sex (M / F)
# Corner - Presence of a corner in the nest section (Y / N)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Ratio - Ratio of male alates over total alates in the observation (0 - 1)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ScaledDist ~ Nest + Sex + Ratio + Day + Corner + (1 | Colony), data = AlateDistScaledRD2Plot))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest + Sex + Ratio + Day + Corner + (1 | Colony), data = AlateDistScaledRD2Plot))

####################################################################################################################
# PLOTS AND ANALYSES: Worker distances from the physical nest center
# The scripts below are to analyze and visualize: 
# Worker mean distance to all other nest sections but the one they occupy
####################################################################################################################

# HISTOGRAM
# High density treatment
WorkerMeanDist1 <- ggplot(WorkerDistScaledMeanDist %>% arrange(Nest) %>% filter (Density == "High"),
                          aes(ScaledDistMean, 
                              fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("High density") +
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlab("Mean scaled dist to nest sections") +
  ylab("Worker count") +
  xlim(0, 0.6) +
  ylim(0, 10000)

# Low density treatment
WorkerMeanDist2 <- ggplot(WorkerDistScaledMeanDist %>% arrange(Nest) %>% filter (Density == "Low"),
                          aes(ScaledDistMean, 
                              fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("Low density") +
  labs(color="Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.685),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlim(0, 0.6) +
  ylim(0, 10000)

# REGRESSION LINE PLOT SHOWING THE RELATIONSHIP BETWEEN NEST SHAPE AND WORKER MEAN SCALED DISTANCE TO NEST SECTIONS
WorkerMeanDist3 <- ggplot(WorkerDistScaledMeanDist %>% arrange(Nest),
                          aes(x = ScaledDist, y = ScaledDistMean,
                              color = Nest, 
                              linetype = Density)) +
  geom_smooth(method = lm, se = FALSE) +
  xlab("Scaled dist to the entrance") +
  ylab("Mean scaled dist to nest sections") +
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        plot.title = element_text(size = 18, face = "bold", family = "Arial"),
        legend.key = element_blank(),
        legend.justification = c(1, 0.5),
        legend.box = "verticle",
        legend.position = "right",
        legend.margin=margin(),
        legend.text = element_text(size = 18, family = "Arial"),
        legend.title = element_text(size = 18, family = "Arial"),
        legend.key.size = unit(1, 'cm')) +
  guides(color = guide_legend(title = "Nest", family = "Arial",
                              order = 1),
         linetype = guide_legend(order = 2)) +
  scale_color_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  guides(lty = guide_legend(override.aes = list(col = 'black'))) 

# Compiled worker distance to nest sections histograms
WorkerMeanDistPlot <- ggarrange(WorkerMeanDist1, WorkerMeanDist2,
                         labels = c("(a)", "(b)"),
                         label.x = 0.9,
                         label.y = 1,    
                         font.label = list(size = 18, face = "plain", family = "Arial"),
                         ncol = 2, nrow = 1,
                         common.legend = FALSE)

# Annotating the compiled plots to include common axes
annotate_figure(WorkerMeanDistPlot,
                top = NULL,
                bottom = text_grob("Mean scaled dist to nest sections", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = text_grob("Worker count", color = "black",
                                 size = 18, x = 0.525, rot = 90, family = "Arial"),
                right = NULL
)

# Compiling the worker distance to nest sections regression line plot 
ggarrange(WorkerMeanDist3,
                                labels = c("(c)"),
                                label.x = 0.725,
                                label.y = 0.97,    
                                font.label = list(size = 18, face = "plain", family = "Arial"),
                                ncol = 1, nrow = 1,
                                common.legend = FALSE)

# LINEAR MIXED EFFECTS MODEL: Worker mean scaled distances to nest sections
# RESPONSE VARIABLE
# ScaledDistMean - Worker mean scaled distances to nest sections (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# ScaledDist - Worker scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ScaledDistMean ~ Nest * Density * ScaledDist + Day + Corner + (1 | Colony), data = WorkerDistScaledMeanDist))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDistMean ~ Nest * Density * ScaledDist + Day + Corner + (1 | Colony), data = WorkerDistScaledMeanDist))

####################################################################################################################
# PLOTS AND ANALYSES: Mobile colony member distances from the brood center
# The scripts below are to analyze and visualize: 
# Worker, queen, and alate distances from the brood center
####################################################################################################################

# WORKER SCALED DISTANCE TO THE BROOD CENTER
# HISTOGRAMS
# High density treatment
WorkerBroodDist1 <- ggplot(BroodCentDistWorkersRD1 %>% arrange(Nest),
                           aes(ToBrood, 
                               fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("High density") +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) + 
  xlim(0, 1) + 
  ylim(0, 3000)

# Low density treatment
WorkerBroodDist2 <- ggplot(BroodCentDistWorkersRD2 %>% arrange(Nest), 
                           aes(ToBrood, 
                               fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) + 
  xlim(0, 1) + 
  ylim(0, 3000)

#Compiling worker scaled distances from brood center plots
WorkerBroodDistPlot <- ggarrange(WorkerBroodDist1, WorkerBroodDist2,
                                 labels = c("(a)", "(b)"),
                                 label.x = 0.9,
                                 font.label = list(size = 18, family = "Arial", face = "plain"),
                                 ncol = 2, nrow = 1,
                                 common.legend = FALSE)

#Annotating the compiled plot to include a common y axis and title
WorkerFullBroodDist <- annotate_figure(WorkerBroodDistPlot,
                                       top = text_grob("Workers", color = "black",
                                                       size = 18, x = 0.06, y = -0.675, family = "Arial"),
                                       bottom = NULL,
                                       left = text_grob("Worker count", color = "black",
                                                  size = 18, rot = 90, family = "Arial"),
                                       right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Worker scaled distances to the brood center
# RESPONSE VARIABLE
# ToBrood - Worker scaled distances from the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ToBrood ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodCentDistWorkersRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodCentDistWorkersRD1_RD2))

# QUEEN SCALED DISTANCE TO THE BROOD CENTER
# HISTOGRAMS
# High density treatment
QueenBroodDist1 <- ggplot(BroodCentDistQueensRD1 %>% arrange(Nest),
                        aes(ToBrood, 
                            fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("High density") +
  labs(color = "Nest") +
  xlab(NULL) + 
  ylab(NULL) +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 200)

#Low density treatment
QueenBroodDist2 <- ggplot(BroodCentDistQueensRD2 %>% arrange(Nest),
                          aes(ToBrood, fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 200)

# Compiling queens scaled distances to the brood center plots
QueenBroodDistPlot <- ggarrange(QueenBroodDist1, QueenBroodDist2,
                                labels = c("(c)", "(d)"),
                                label.x = 0.9,
                                font.label = list(size = 18, family = "Arial", face = "plain"),
                                ncol = 2, nrow = 1,
                                common.legend = FALSE)

# Annotating the compiled plots to include a common y axis and title
QueenBroodFullDist <- annotate_figure(QueenBroodDistPlot,
                                      top = text_grob("Queens", color = "black",
                                                      size = 18, x = 0.06, y = -0.75, family = "Arial"),
                                      bottom = NULL,
                                      left = text_grob("Queen count", color = "black",
                                                  size = 18, rot = 90, family = "Arial"),
                                      right = NULL
)

# Compiling the Worker and queen scaled distance to the brood center
WorkerQueenBroodDist<-ggarrange(WorkerFullBroodDist, QueenBroodFullDist,
                          ncol = 1, nrow = 2)

# Annotating the compiled plot to include a commpn x axis
annotate_figure(WorkerQueenBroodDist,
                top = NULL,
                bottom = text_grob("Scaled distance to brood center", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Queen scaled distances to the brood center
# RESPONSE VARIABLE
# ToBrood - Queen scaled distances from the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Corner - Presence of a corner in the nest section (Y / N)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ToBrood ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodCentDistQueensRD1_RD2))

#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest * Density + Day + Corner + (1 | Colony), data = BroodCentDistQueensRD1_RD2))

# ALATE SCALED DISTANCE TO THE BROOD CENTER
# Removing individuals with unknown alate sex
# In these plots and analyses, we care about the individuals with identifiable sex
# 10 individuals are removing the unknown sex individuals, representing < 1% of the data set  
BroodCentDistAlatesRD2Plot <- BroodCentDistAlatesRD2 %>% filter(Sex != "?")

# HISTOGRAM
AlateBroodPlot <- ggplot(BroodCentDistAlatesRD2Plot %>% arrange(Nest),
                         aes(ToBrood, 
                             fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  xlab("Scaled distance to brood center") +
  ylab("Alate count") +
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial"),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue"))

# BOXPLOT SHOWING THE RELATIONSHIP BETWEEN ALATE SEX AND SCALED DISTANCE TO THE BROOD CENTER
AlateBroodPlot2 <- ggplot(BroodCentDistAlatesRD2Plot %>% arrange(Nest),
                     aes(x = Sex, y = ToBrood, 
                         fill = Nest)) +
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) + 
  theme_pubclean() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlab("Alate sex") +
  ylab("Scaled distance to brood center")

# Compiling the alate distance to the the brood center plots 
AlateToBroodPlot <- ggarrange(AlateBroodPlot, AlateBroodPlot2,
                           labels = c("(e)", "(f)"),
                           label.x = 0.9,
                           label.y = 0.965, 
                           font.label = list(size = 18, face = "plain", family = "Arial"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled pltos to include a title
annotate_figure(AlateToBroodPlot,
                top = text_grob("Alates", color = "black",
                                size = 18, x = 0.06, y = -1.30, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Alate scaled distances from the nest entrance
# RESPONSE VARIABLE
# ToBrood - Alate scaled distances from the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Sex - Alate sex (M / F)
# Corner - Presence of a corner in the nest section (Y / N)
# Day - Experimental observation (From days 1-16, always days 3-14 in High density treatment)
# Ratio - Ratio of male alates over total alates in the observation (0 - 1)
# RANDOM EFFECT
# (1|Colony) - Colony identification 
summary(lmer(ToBrood ~ Nest + Sex + Day + Corner + Ratio + (1 | Colony), data = BroodCentDistAlatesRD2Plot))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest + Sex + Day + Corner + Ratio + (1 | Colony), data = BroodCentDistAlatesRD2Plot))


####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes & nest shape
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes, and consider how they relate to nest shape
####################################################################################################################

# SPATIAL FIDELITY ZONE SIZE AND NEST SHAPE
# BOXPLOTS
# High density treatment
FidZone.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>%
                      filter(Colony < 11) %>% arrange(Nest), 
                    aes(x = as.factor(Colony), y = SFZ),
                    position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +  
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 0.5)

# Low density treatment
FidZone.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>%
                      filter(Colony > 10) %>% arrange(Nest), 
                    aes(x = as.factor(Colony), y = SFZ),
                    position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +  
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.725),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 0.5)

#Compiling worker fidelity zone size plots
FidZonePlot <- ggarrange(FidZone.1, FidZone.2,
                         labels = c("(a)", "(b)"),
                         label.x = 0.9,
                         font.label = list(size = 18, family = "Arial", face = "plain"),
                         ncol = 2, nrow = 1,
                         common.legend = FALSE)

#Annotating the compiled plots to include a common y-axis
FidZonePlotFull <- annotate_figure(FidZonePlot,
                                   top = NULL,
                                   bottom = NULL,
                                   left = text_grob("Fidelity zone size", color = "black",
                                                    size = 18, rot = 90, family = "Arial"),
                                   right = NULL
)

# OCCURRENCE ZONE SIZE AND NEST SHAPE
# BOXPLOTS
# High density treatment
Occur.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZFull %>%
                    filter(Colony < 11) %>% arrange(Nest), 
                  aes(x = as.factor(Colony), y = Occur),
                  position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +  
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none",
        legend.key = element_blank()) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 0.5)

# Low density treatment
Occur.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZFull %>%
                    filter(Colony > 10) %>% arrange(Nest), 
                  aes(x = as.factor(Colony), y = Occur),
                  position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, color = "white", family = "Arial"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 0.5)

# Compiling worker occurrence zone size plots
OccurZonePlot <- ggarrange(Occur.1, Occur.2,
                           labels = c("(c)", "(d)"),
                           label.x = 0.9,
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotate the combined plots to include a common y-axis 
OccurZonePlotFull <- annotate_figure(OccurZonePlot,
                                     top = NULL,
                                     bottom = NULL,
                                     left = text_grob("Occurrence zone size", color = "black",
                                                      size = 18, rot = 90, family = "Arial"),
                                     right = NULL
)

# Compile the spatial fidelity and occurrence zone plots
FidOccurPlot<-ggarrange(FidZonePlotFull, OccurZonePlotFull,
                            ncol = 1, nrow = 2,
                            common.legend = TRUE)

# Annotate the compiled plots to include a common x-axis
annotate_figure(FidOccurPlot,
                top = NULL,
                bottom = text_grob("Colony identification", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(SFZ ~ Nest * Density + Colony + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))
  
#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ Nest * Density + Colony + (1 | ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# Occur - Worker occurrence zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(Occur ~ Nest * Density + Colony + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

#Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ Nest * Density + Colony + (1 | ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

# SUPPLEMENTAL SFZ AREA PLOTS AND ANALYSES
# SPATIAL FIDELITY ZONE SIZE (IN UNSCALED NEST AREA - cm^2) AND NEST SHAPE
# BOXPLOTS
# High density treatment
FidZoneArea.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>%
                          filter(Colony < 11) %>% arrange(Nest), 
                        aes(x = as.factor(Colony), y = SFZ_Area),
                        position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) + 
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 4)

# Low density treatment
FidZoneArea.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>%
                          filter(Colony > 10) %>% arrange(Nest), 
                        aes(x = as.factor(Colony), y = SFZ_Area),
                        position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +  
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.05, 0.725),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", family = "Arial", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 4)

# Compiling worker unscaled fidelity zone size plots
FidZoneAreaPlot <- ggarrange(FidZoneArea.1, FidZoneArea.2,
                         labels = c("(a)", "(b)"),
                         label.x = 0.9,
                         font.label = list(size = 18, family = "Arial", face = "plain"),
                         ncol = 2, nrow = 1,
                         common.legend = FALSE)

#Annotating the compiled plots to include a common y-axis
FidZoneAreaPlotFull <- annotate_figure(FidZoneAreaPlot,
                                   top = NULL,
                                   bottom = NULL,
                                   left = text_grob(expression(paste('Fidelity zone size ('*cm^2*')')), color = "black",
                                                    size = 18, rot = 90, family = "Arial"),
                                   right = NULL
)

# OCCURRENCE ZONE SIZE (IN UNSCALED NEST AREA - cm^2) AND NEST SHAPE
# BOXPLOTS
# High density treatment
OccurArea.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                        filter(Colony < 11) %>% arrange(Nest), 
                      aes(x = as.factor(Colony), y = Occur_Area),
                      position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +  
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none",
        legend.key = element_blank()) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 5)

#Low density treatment
OccurArea.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                        filter(Colony > 10) %>% arrange(Nest), 
                      aes(x = as.factor(Colony), y = Occur_Area),
                      position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65) +
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, color = "white", family = "Arial"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 5)

# Compiling worker occurrence zone size plots
OccurZoneAreaPlot <- ggarrange(OccurArea.1, OccurArea.2,
                               labels = c("(c)", "(d)"),
                               label.x = 0.9,
                               font.label = list(size = 18, family = "Arial", face = "plain"),
                               ncol = 2, nrow = 1,
                               common.legend = FALSE)

# Annotate the combined plots to include a common y-axis 
OccurZoneAreaPlotFull <- annotate_figure(OccurZoneAreaPlot,
                                         top = NULL,
                                         bottom = NULL,
                                         left = text_grob(expression(paste('Occurrence zone size ('*cm^2*')')), color = "black",
                                                          size = 18, rot = 90, family = "Arial"),
                                         right = NULL
)

# Compile the spatial fidelity and occurrence zone plots and include a common legend
FidOccurAreaPlot<-ggarrange(FidZoneAreaPlotFull, OccurZoneAreaPlotFull,
                            ncol = 1, nrow = 2,
                            common.legend = TRUE)

# Annotate the compiled plots to include a common x-axis
annotate_figure(FidOccurAreaPlot,
                top = NULL,
                bottom = text_grob("Colony identification", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (cm^2) and nest shape
# RESPONSE VARIABLE
# SFZ_Area - Worker spatial fidelity zone size (cm^2), zones have at least 3 observations 
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(SFZ_Area ~ Nest * Density + Colony + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ_Area ~ Nest * Density + Colony + (1 | ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (cm^2) and nest shape
# RESPONSE VARIABLE
# Occur_Area - Worker spatial fidelity zone size (cm^2), 
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(Occur_Area ~ Nest * Density + Colony + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur_Area ~ Nest * Density + Colony + (1 | ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

#Function to create large points in a geom_point legend
large_points <- function(data, params, size) {
  # Multiply by some number, it doesn't matter what value, but larger numbers = large sized points in the legend
  data$size <- data$size * 2.5
  draw_key_point(data = data, params = params, size = size)
}

####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes & distances from the nest entrance
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes, and consider how they relate to both nest shape and mean distance from the nest entrance
####################################################################################################################

# FIDELITY ZONE SIZE AND DISTANCE TO THE NEST ENTRANCE
# LINE PLOTS
# High density treatment
SFZDist1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>% filter(Density == "High") %>% arrange(Nest), 
                   aes(x = MeanScaledDist, y=SFZ, 
                       linetype = Nest,
                       color = Nest, 
                       shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.position = "none") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  ylim(0.05, 0.5) +
  xlim(0, 0.75)

# Low density treatment
SFZDist2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorkingFid %>% filter(Density == "Low") %>% arrange(Nest),
                 aes(x = MeanScaledDist, y = SFZ, 
                     linetype = Nest,
                     color = Nest,
                     shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  ylim(0.05, 0.5) +
  xlim(0, 0.75)

# Compiling worker fidelity zone size v. worker scaled distance to the nest entrance plots
SFZDistPlot <- ggarrange(SFZDist1, SFZDist2,
                         labels = c("(a)", "(b)"),
                         label.x = 0.9,
                         font.label = list(size = 18, family = "Arial", face = "plain"),
                         ncol = 2, nrow = 1,
                         common.legend = FALSE)

# Annotating the compiled plots to include a common y-axis
SFZFullDistPlot <- annotate_figure(SFZDistPlot,
                                   top = NULL,
                                   bottom = NULL,
                                   left = text_grob("Fidelity zone size", color = "black",
                                                    size = 18, rot = 90, family = "Arial"),
                                   right = NULL
)

# OCCURRENCE ZONE SIZE AND DISTANCE TO THE NEST ENTRANCE
# LINE PLOTS
# High density treatment
OccurDist1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "High") %>% arrange(Nest),
                   aes(x = MeanScaledDist, y = Occur,
                       linetype = Nest,
                       color = Nest,
                       shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  ylim(0.05, 0.5) +
  xlim(0, 0.75)


#Low density treatment
OccurDist2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "Low") %>% arrange(Nest),
                     aes(x = MeanScaledDist, y = Occur,
                         color = Nest, 
                         linetype = Nest,
                         shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  ylim(0.05, 0.5) +
  xlim(0, 0.75)

# Compiling worker occurrence zone size v. worker scaled distance to the nest entrance plots
OccurDistPlot <- ggarrange(OccurDist1, OccurDist2,
                           labels = c("(c)", "(d)"),
                           label.x = 0.9,
                           font.label = list(size = 18, family = "Arial", face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled plots to include a common y-axis
OccurFullDistPlot <- annotate_figure(OccurDistPlot,
                                       top = NULL,
                                       bottom = NULL,
                                       left = text_grob("Occurrence zone size", color = "black",
                                                        size = 18, rot = 90, family = "Arial"),
                                       right = NULL
)

# Compile the spatial fidelity and occurrence zone v. worker scaled distance to the nest entrance plots and include a common legend
FidOccurDistPlot<-ggarrange(SFZFullDistPlot, OccurFullDistPlot,
                                ncol = 1, nrow = 2,
                            common.legend = TRUE)

# Annotate the compiled plots to include a common x-axis
annotate_figure(FidOccurDistPlot,
                top = NULL,
                bottom = text_grob("Average scaled distance to nest entrance", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# MeanScaledDist - each worker's average distance to the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(SFZ ~ MeanScaledDist * Nest * Density + (1 | Colony) + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ MeanScaledDist * Nest * Density + (1 | Colony) + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorkingFid))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations 
# FIXED EFFECTS 
# MeanScaledDist - each worker's average distance to the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(Occur ~ MeanScaledDist * Nest * Density + (1 | Colony) + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ MeanScaledDist * Nest * Density + (1 | Colony) + (1|ColorID), data = WorkerDistScaledRD1_RD2SFZWorking))

####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes & distances from the nest entrance and brood center
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes, and consider how they relate to both nest shape and mean distance from the brood center
####################################################################################################################

# FIDELITY ZONE SIZE AND DISTANCE TO THE BROOD CENTER
# LINE PLOTS
# High density treatment
SFZBroodDist1 <- ggplot(data = BroodCentDistWorkersSFZFid %>% filter(Density == "High") %>% arrange(Nest), 
                        aes(x = MeanToBrood, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("High density") +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -20),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  ylim(0.05, 0.4) +
  xlim(0, 0.5)

# Low density treatment
SFZBroodDist2 <- ggplot(data = BroodCentDistWorkersSFZFid %>% filter(Density == "Low") %>% arrange(Nest),
                        aes(x = MeanToBrood, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 18, family = "Arial", color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", hjust = 0.875, vjust = -10),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 18, family = "Arial", color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  ylim(0.05, 0.4) +
  xlim(0, 0.5)

# Compiling worker fidelity zone size v. worker scaled distance to the brood center plots and include a common legend
SFZBroodDistPlot <- ggarrange(SFZBroodDist1, SFZBroodDist2,
                              labels = c("(a)", "(b)"),
                              label.x = 0.9,
                              font.label = list(size = 18, family = "Arial", face = "plain"),
                              ncol = 2, nrow = 1,
                              common.legend = TRUE)

# Annotating the compiled plots to include a common y-axis
SFZFullBroodDistPlot <- annotate_figure(SFZBroodDistPlot,
                                        top = NULL,
                                        bottom = NULL,
                                        left = text_grob("Fidelity zone size", color = "black",
                                                         size = 18, rot = 90, family = "Arial"),
                                        right = NULL
)

# OCCURRENCE ZONE SIZE AND DISTANCE TO THE BROOD CENTER
# LINE PLOTS
# High density treatment
OccurBroodDist1 <- ggplot(data = BroodCentDistWorkersSFZ %>% filter(Density == "High") %>% arrange(Nest),
                          aes(x = MeanToBrood, y = Occur,
                              color = Nest,
                              linetype = Nest,
                              shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  ylim(0.1, 0.5) +
  xlim(0, 0.5)

# Low density treatment
OccurBroodDist2 <- ggplot(data = BroodCentDistWorkersSFZ %>% filter(Density == "Low") %>% arrange(Nest),
                          aes(x = MeanToBrood, y = Occur, 
                              color = Nest, 
                              linetype = Nest,
                              shape = Nest)) +
  geom_point(key_glyph = large_points, size = 3, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 1.25, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, family = "Arial", color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", family = "Arial", color = "white", hjust = 0.75, vjust = -20),
        legend.position = "none") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  ylim(0.1, 0.5) +
  xlim(0, 0.5)

# Compiling worker occurrence zone size v. worker scaled distance to the nest entrance plots
OccurBroodDistPlot <- ggarrange(OccurBroodDist1, OccurBroodDist2,
                                labels = c("(c)", "(d)"),
                                label.x = 0.9,
                                font.label = list(size = 18, family = "Arial", face = "plain"),
                                ncol = 2, nrow = 1,
                                common.legend = FALSE)

# Annotating the compiled plots to include a common y-axis
OccurFullBroodDistPlot <- annotate_figure(OccurBroodDistPlot,
                                     top = NULL,
                                     bottom = NULL,
                                     left = text_grob("Occurrence zone size", color = "black",
                                                      size = 18, rot = 90, family = "Arial"),
                                     right = NULL
)

# Compile the spatial fidelity and occurrence zone v. worker scaled distance to the brood center plots and include a common legend
FidOccurBroodDistPlot<-ggarrange(SFZFullBroodDistPlot, OccurFullBroodDistPlot,
                            ncol = 1, nrow = 2,
                            common.legend = TRUE)

# Annotating the compiled plots to include a common x-axis
annotate_figure(FidOccurBroodDistPlot,
                top = NULL,
                bottom = text_grob("Average scaled distance to brood center", color = "black",
                                   size = 18, x = 0.525, family = "Arial"),
                left = NULL,
                right = NULL
)

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# BroodDist - each worker's average distance to the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(SFZ ~ MeanToBrood * Nest * Density + (1 | Colony) + (1|ColorID), data = BroodCentDistWorkersSFZFid))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ MeanToBrood * Nest * Density + (1 | Colony) + (1|ColorID), data = BroodCentDistWorkersSFZFid))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# Occur - Worker occurrence zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations 
# FIXED EFFECTS 
# BroodDist - each worker's average distance to the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
# (1|ColorID) - Worker color identification marks on the head, thorax, abdomen 1, abdomen 2 (e.g. W,G,W,G)
summary(lmer(Occur ~ MeanToBrood * Nest * Density + (1 | Colony) + (1|ColorID), data = BroodCentDistWorkersSFZ))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ MeanToBrood * Nest * Density + (1 | Colony) + (1|ColorID), data = BroodCentDistWorkersSFZ))

####################################################################################################################
# COLONY MEMBER SCATTERPLOTS
# The below script produces the colony member density plots
# Each scatter plot colors the points by K-nearest neighbors, which is defined by how many points are in proximity to one another
# First the scatter plots for the Netlogo simulation results are created, then each colony's worker, brood, queen, and alate (low nest density only) are created
# Each combined plot is assigned to an object that can be run to view 
####################################################################################################################

# SIMULATION SCATTERPLOTS
# High density treatment
# Note that the data table has the ScaledDist column which is used here to remove points at the very back of the nest for the Netlogo simulation scatter plot 
# Without these points being removed the back of the nest overlaps into the first nest section visually
# This correction is not used in any of the other plots and represents a very small portion of the data points
SmallRandWalkDensity <- ggplot(data = SimDistScaled %>% filter(NestSize == "Small" & ScaledDist < 0.995),
                               aes(ScaledX, ScaledY)) +
  geom_pointdensity() +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Random walk") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN") +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# Low density treatment
# Note that the data table has the ScaledDist column which is used here to remove points at the very back of the nest for the Netlogo simulation scatter plot 
# Without these points being removed the back of the nest overlaps into the first nest section visually
# This correction is not used in any of the other plots and represents a very small portion of the data points
LargeRandWalkDensity <- ggplot(data = SimDistScaled %>% filter(NestSize == "Large" & ScaledDist < 0.995),
                               aes(ScaledX, ScaledY)) +
  geom_pointdensity() +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Random walk") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN") +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5)

# EMPIRICAL SCATTERPLOTS
# Below are the colony member scatter plots for every colony

# MAIN MANUSCRIPT FIGURE
# This figure is the same a colony 11 in the supplementaty figures below
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 11),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 11),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 11),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 11),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 11 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity, 
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(FullDensityPlot,
                top = text_grob("Colony Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# SUPPLEMENTARY PLOTS
# HIGH NEST DENSITY
# COLONY 1  
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 1),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
# Queens are represented as red triangles instead of color scaled by K-nearest neighbor
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 1),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 1),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# Compiling the colony 1 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony1 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony1,
                top = text_grob("Colony 1 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 2
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 2),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 2),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 2),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# Compiling the colony 2 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony2 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony2,
                top = text_grob("Colony 2 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 3
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 3),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 3),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 3),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# Compiling the colony 3 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which can be used below or later
DensityPlotColony3 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony3,
                top = text_grob("Colony 3 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 4
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 4),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 4),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 4),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 4 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony4 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony4,
                top = text_grob("Colony 4 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 5
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 5),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 5),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 5),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# Compiling the colony 5 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony5 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony5,
                top = text_grob("Colony 5 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 6
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 6),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 6),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 6),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 6 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony6 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony6,
                top = text_grob("Colony 6 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 7
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 7),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 7),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 7),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 7 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony7 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony7,
                top = text_grob("Colony 7 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 8
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 8),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 8),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 8),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 8 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony8 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony8,
                top = text_grob("Colony 8 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 9
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 9),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 9),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 9),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 9 worker, brood & queen scatter plots with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony9 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony9,
                top = text_grob("Colony 9 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 10
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 10),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 10),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 10),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 10 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             SmallRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony10 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony10,
                top = text_grob("Colony 10 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# LOW NEST DENSITY
# COLONY 11
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 11),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 11),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 11),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 11),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 11 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity, 
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony11 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(FullDensityPlot,
                top = text_grob("Colony 11 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 12
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 12),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 12),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 12),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 12) %>% add_row(Nest = "Circle"),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest, drop = F) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 12 worker, brood & queen, and alate scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony12 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(FullDensityPlot,
                top = text_grob("Colony 12 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 13
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 13),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 13),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 13),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 13) %>% add_row(Nest = "Circle"),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_grid(~ Nest, drop = FALSE) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 13 worker, brood & queen, and alate scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony13 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony13,
                top = text_grob("Colony 13 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 14
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 14),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 14),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 14),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 14),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 14 worker, brood & queen, and alate scatter plots, with the high density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony14 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony14,
                top = text_grob("Colony 14 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 15
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 15),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 15),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 15),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 15) %>% add_row(Nest = "Circle"),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest, drop = F) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 15 worker, brood & queen, and alate scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony15 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony15,
                top = text_grob("Colony 15 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 16
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 16),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 16),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 16),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 16 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony16 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony16,
                top = text_grob("Colony 16 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 17
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 17),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 17),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 17),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 17 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony17 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony17,
                top = text_grob("Colony 17 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 18
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 18),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 18),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 18),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 18) %>% add_row(Nest = "Tube"),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 18 worker, brood & queen, and alate scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony18 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony18,
                top = text_grob("Colony 18 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 19
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 19),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 19),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 19),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# ALATE SCATTERPLOTS ONLY IN THE FOLLOWING COLONIES
#   |  Colony |
#   |  11     |
#   |  12     |
#   |  13     |
#   |  14     |
#   |  15     |
#   |  18     |
#   |  19     |

AlateDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 19) %>% add_row(Nest = "Tube"),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  facet_wrap(~ Nest) +
  labs(color = "KNN")+
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5)+
  ylim(-0.25, 5) 

# Compiling the colony 19 worker, brood & queen, and alate scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             AlateDensity,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony19 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(FullDensityPlot,
                top = text_grob("Colony 19 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)

# COLONY 20
# WORKERS
WorkerDensityColony <- ggplot(data = WorkerDistScaledRD1_RD2 %>% filter(Colony == 20),
                              aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Workers") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD & QUEENS
BroodQueenDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 20),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 20),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, family = "Arial", color = "black"),
        legend.title = element_text(size = 20, family = "Arial", color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, family = "Arial", color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN")+
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 


# Compiling the colony 20 worker, brood & queen scatter plots with the low density simulation scatter plot
FullDensityPlot <- ggarrange(WorkerDensityColony, BroodQueenDensityColony,
                             LargeRandWalkDensity,
                             font.label = list(size = 18, face = "plain", family = "Arial"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony20 <- FullDensityPlot

# Annotating the compiled plot and include a title
annotate_figure(DensityPlotColony20,
                top = text_grob("Colony 20 Member Densities", color = "black",
                                size = 18, x = 0.2, family = "Arial"),
                bottom = NULL,
                left = NULL,
                right = NULL
)
