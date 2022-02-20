####################################################################################################################
## Author: GREG CHISM
## Date: FEB 2022
## email: gchism@email.arizona.edu
## Project: Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space
## Title: Distance functions 
####################################################################################################################

# READ ME:
# YOU MUST RUN THE FOLLOWING SCRIPT BEFORE THIS ONE: Bins_Working.R

# These code is to replicate the data transformations for my first chapter:
# Distance to the nest entrance: search for DISTANCE TO THE NEST ENTRANCE
# Shortest distance to nest sections: search for SHORTEST DISTANCE TO NEST SECTIONS
# Distance to the brood center: starts at CALCULATE THE BROOD CENTER 

####################################################################################################################
# IMPORT ALL NECESSARY DATASETS 
# This code imports all necessary data sets for the script below
####################################################################################################################

# DISTANCE TO NEST ENTRANCE (EMPIRICAL) REFERENCES
DistBinsFull <- read.csv("DistBinsFull.csv")

# DISTANCE TO NEST ENTRANCE (NETLOGO SIMULATIONS) REFERENCES
DistBinsFullNetlogo <- read.csv("DistBinsFullNetlogo.csv")

# CIRCLE REFERENCES FOR WORKER SHORTEST DISTANCE TO NEST SECTIONS
NestBinCircle <- read.csv("NestBinCircle.csv")

# DISTANCE FROM THE BACK OF EACH CIRCLE NEST TO THE ENTRANCE 
MaxCircleDist <- read.csv("MaxCircleDist.csv")

# CORNER REFERENCE STATING IF A BIN HAS A CORNER PRESENT
CornerFull <- read.csv("CornerFull.csv")


####################################################################################################################
## DISTANCE TO THE NEST ENTRANCE
# The following scripts find each colony member's distance from the nest entrance 
# The code then scales these distances from 0 - 1, such that 0 is the entrance and 1 is the back of the tube nest shape. 
# This means that the farthest distance in the circle nest is ~ 0.24-0.32, which varies due to different methods in cutting out nests
# These reference distances have the shortest possible distance from the beginning of that section to the entrance, represented in the "Distance" column.
# The reference data set also has the maximum possible distance from each colony's tube nest to the entrance: The "MaxDist" column. 
####################################################################################################################

# WORKERS
# High density treatment
DistanceCoordsFunction <- function(data.table) {
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
  WorkerDistScaled <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist),
           WorkerType = "Obsv") %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, WorkerType) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not 
}

# Run the shortest distance from the bin function for the high density workers data set
DistanceCoordsFunction(FullDataCoordWorkers)


# Low density treatment
DistanceCoordsFunctionRD2 <- function(data.table) {
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
  WorkerDistScaledRD2 <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist),
           WorkerType = "Obsv") %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, WorkerType) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not.
}

# Run the shortest distance from the bin function for the low density workers data set
DistanceCoordsFunctionRD2(FullDataCoordWorkersRD2)

# Join the worker distance to the entrance data sets 
WorkerDistScaledRD1_RD2 <- full_join(WorkerDistScaled, WorkerDistScaledRD2) %>%
  mutate(WorkerType = "Worker")

# BROOD
# High density treatment
DistanceCoordsFunctionBrood <- function(data.table) {
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
  BroodDistScaled <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist)) %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not.
}

# Run the shortest distance from the bin function for the high density brood data set
DistanceCoordsFunctionBrood(FullDataCoordBrood)

# BROOD
# Low density treatment
DistanceCoordsFunctionBroodRD2 <- function(data.table) {
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
  BroodDistScaledRD2 <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist)) %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not
}

# Run the shortest distance from the bin function for the low density brood data set
DistanceCoordsFunctionBroodRD2(FullDataCoordBroodRD2)

# Join brood distance to the entrance data sets
BroodDistScaledRD1_RD2 <- full_join(BroodDistScaled, BroodDistScaledRD2) 

# QUEENS
# High density treatment
DistanceCoordsFunctionQueen <- function(data.table) {
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
  QueenDistScaled <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist)) %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not
}

# Run the shortest distance from the bin function for the high density queens data set
DistanceCoordsFunctionQueen(FullDataCoordQueen)

# QUEENS
# Low density treatment
DistanceCoordsFunctionQueenRD2 <- function(data.table) {
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
  QueenDistScaledRD2 <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist)) %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not
}

# Run the shortest distance from the bin function for the high density queens data set
DistanceCoordsFunctionQueenRD2(FullDataCoordQueenRD2)

# Join queen distance data sets
QueenDistScaledRD1_RD2 <- full_join(QueenDistScaled, QueenDistScaledRD2) 

# ALATE SEX RATIOS
# Note, these ratios incluse "?", where the sex of the alate was uncertain.
# We chose males because they were always present when alates were in the nest.
FullDataCoordAlatesMales <- FullDataCoordAlates %>%
  filter(Sex == "M") %>% # Filter only males
  group_by(Colony, Nest, Day) %>% # Group by the Colony, Nest, and Day columns, allowing us to find the below values for each photo.
  mutate(MaxTotal = max(TotalNumber), # The number of alates in that observation
         MaxSex = max(SexNumber), # The number of males in that observation
         Ratio = (MaxSex) / (MaxTotal)) %>% # The ratio of males / total alates
  select(Colony, Nest, Day, Ratio) %>% # Select the desired columns
  distinct() # Remove duplicate rows

# Join the alate and male data sets
FullDataCoordAlatesRatio <- left_join(FullDataCoordAlates, FullDataCoordAlatesMales)

# ALATES
# Alates are only in the Low density treatment
DistanceCoordsFunctionAlate <- function(data.table) {
  Bin1 <- DistBinsFull %>% # Bin 1 has corners that are ignored unless there is a way to make the shortest distance to the corner first, this code produces the reference coordinates for the corners
    filter(Bin == 1) %>%
    group_by(Colony, Nest) %>%
    mutate(BinX1 = BinX, LeftCorner = BinX1 - 0.1, RightCorner = BinX1 + 0.1) %>%
    select(Colony, Nest, BinX1, LeftCorner, RightCorner)
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% 
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
  DistBinsFull <- left_join(DistBinsFull, Bin4) %>% # Full distance references, joining all alternative references from above
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
  AlateDistScaledRD2 <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(Colony) %>% # Group by the Colony column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist)) %>% 
    ungroup() %>% # Ungroup the data set
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Sex, Ratio) %>% # Select the desired columns. 
    left_join(CornerFull) # Join with data set that associates nest sections with the presence of a corner or not
}

# Run the shortest distance from the bin function for the alate data set
DistanceCoordsFunctionAlate(FullDataCoordAlatesRatio)

# NETLOGO SIMULATIONS
DistanceCoordsFunctionNetlogo <- function(data.table) {
  Bin1 <- DistBinsFullNetlogo %>% # Bin 1 has corners that are ignored unless there is a way to make the shortest distance to the corner first, this code produces the reference coordinates for the corners
    filter(Bin == 1) %>%
    group_by(NestSize, Nest) %>%
    mutate(BinX1 = BinX, LeftCorner = 3.65, RightCorner = 3.85) %>%
    select(NestSize, Nest, BinX1, LeftCorner, RightCorner)
  Bin4 <- DistBinsFullNetlogo %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% 
    group_by(NestSize) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(NestSize, BinY4 ,Distance4)
  Bin7 <- DistBinsFullNetlogo %>% # Bin 7 has the same problem as above, one shortest distance that uses the x coordinate of bin 7, and another at Bin 4 where Scaled Y values are greater than the y coordinate of Bin 7
    filter(Nest == "Tube" & Bin == 4) %>%
    group_by(NestSize) %>%
    mutate(BinX7 = BinX,
           Distance7 = Distance) %>% 
    select(NestSize, BinX7, Distance7) # The distance that is added to the Pythagorean distance is also different
  Bin3 <- DistBinsFullNetlogo %>% # Bin 3 has the same problem as above, one shortest distance that uses the x coordinate of bin 3, and another at Bin 2 where Scaled Y values are less than the y coordinate of Bin 3
    filter(Nest == "Tube" & Bin == 2) %>%
    group_by(NestSize) %>%
    mutate(BinX3 = BinX,
           BinY3 = BinY,
           Distance3 = Distance) %>%
    select(NestSize, BinY3, BinX3, Distance3) # The distance that is added to the Pythagorean distance is also different
  DistBinsFullNetlogo <- left_join(DistBinsFullNetlogo, Bin4) %>% #Full distance references, joining all alternative references from above
    left_join(Bin7) %>%
    left_join(Bin3) %>%
    left_join(Bin1)
  DistBins.tube <- DistBinsFullNetlogo %>% # Subsetting tube distance references
    filter(Nest == "Tube") 
  DistBins.circle <- DistBinsFullNetlogo %>% # Subsetting circle distance references
    filter(Nest == "Circle") 
  BinsTube <- data.table %>% # Subsetting tube nest data
    filter(Nest == "Tube")
  BinsCircle <- data.table %>% # Subsetting circle nest data
    filter(Nest == "Circle") 
  
  # TUBE NEST SHAPE
  DistBinsTube <- left_join(BinsTube, DistBins.tube) %>% # Joining tube nest data with associated reference coordinates
    group_by(NestSize, Bin) %>% # Group by NestSize and Bin columns
    filter(Bin != 1) %>% # Filter out Bin 1, since this will be handled in the code below
    mutate(DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference bin 4 coordinate
           DistanceX7 = ScaledX - BinX7, # Distance from each individual's x coordinate to the x reference bin 7 coordinate
           DistanceX3 = ScaledX - BinX3, # Distance from each individual's x coordinate to the x reference bin 2 coordinate
           DistanceY3 = ScaledY - BinY3, # Distance from each individual's x coordinate to the x reference bin 2 coordinate
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
    group_by(NestSize, Bin) %>% # Group by NestSize and Bin columns
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
    group_by(NestSize) %>% # Group by the NestSize column
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
  SimDistScaled <<- full_join(DistBinsTube, DistBinsCircle) %>% # Full join the Tube and Circle nest data sets 
    group_by(NestSize) %>% # Group by the NestSize column
    mutate(ScaledDist = DistanceTotal / (MaxDist), # Divide the distance of each individual to the entrance by the maximum possible distance in the tube nest shape. 
           ScaledDist = ifelse(ScaledDist > 1, 1, ScaledDist),
           Density = ifelse(NestSize == "Small", "High", "Low"), # Create a density colimn
           WorkerType = "RandSim") %>% # Create a worker type column
    ungroup() %>% # Ungroup the data set
    select(NestSize, Nest, RunNumber, ScaledX, ScaledY, Bin, ScaledDist, Density, WorkerType) %>% # Select the desired columns. 
    left_join(CornerFullSim) %>% # Join with data set that associates nest sections with the presence of a corner or not
    distinct()
}

# Run the shortest distance from the bin function for the Netlogo simulation data set
DistanceCoordsFunctionNetlogo(NetlogoBinnedFull)

####################################################################################################################
# SHORTEST DISTANCE TO NEST SECTIONS
# The following code is used to calculate the shortest distance each worker is from all nest sections but their own
# To do this, we use a reference data set that defines the shortest distance from the beginning of each bin to the nest entrance
####################################################################################################################

# NOTE: Most of below object assignments are redundant from the distance to the brood entrance calculations, but this allows you to run the script independently 

# TUBE NEST SHAPE 
# First we make data sets of distance segments either between a bin and the entrance or between two sections
# This is done outside of the actual function below to take less computing time
TubeDistances <- DistBinsFull %>% # Filtering out the Tube reference distances
  filter(Nest == "Tube")

# Filtering out only distances for bin 1
DistBin1 <- TubeDistances %>% 
  filter(Bin == 1) %>%
  mutate(BinX = 3.75, # Create x reference column
         BinY = 0, # Create y reference column
         Dist1 = Distance) %>% # Set reference distance
  select(Colony, Nest ,Dist1) # Select desired columns

# Filtering and creating reference distances for nest sections (bins) 2 - 8
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
    Distance7_0 = Dist7, # Distance from bin 7 to 0
    Distance3_0 = Dist3, # Distance from bin 4 to 0 
    Distance8_7 = abs(Dist8 - Dist7), # Distance from bin 8 to 6 or 4 to 2 or 4 to 5  
    Distance8_1 = abs(Dist8 - Dist1), # Distance from bin 8 to 1
    Distance7_1 = abs(Dist7 - Dist1), # Distance from bin 7 to 1
    Distance3_1 = abs(Dist3 - Dist1), # Distance from bin 4 to 1
    Distance2_1 = abs(Dist2 - Dist1), # Distance from bin 3 to 1
    Distance8_2 = abs(Dist8 - Dist2), # Distance from bin 8 to 2
    Distance7_2 = abs(Dist7 - Dist2), # Distance from bin 7 to 2
    Distance7_3 = abs(Dist7 - Dist3), # Distance from bin 7 to 3
    Distance8_4 = abs(Dist8 - Dist3), # Distance from bin 8 to 4
    Distance7_4 = abs(Dist7 - Dist3 - Distance8_7) # Distance from bin 7 to 4
  )

# Bin 4 and 7 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left,
# which we account for this with two sets of possible reference coordinates

# Creating the second bin 4 reference y coordinate
Bin4 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 3) %>%
  group_by(Colony) %>%
  mutate(BinY4 = BinY) %>%
  select(Colony, BinY4)

# Creating the second bin 7 reference y coordinate
Bin6 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 8) %>%
  group_by(Colony) %>%
  mutate(BinY7 = BinY) %>%
  select(Colony, BinY7)

# Creating the second bin 4 reference x coordinate
Bin7 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 4) %>%
  group_by(Colony) %>%
  mutate(BinX7 = BinX) %>%
  select(Colony, BinX7)

# Creating the second bin 6 reference x coordinate
Bin4.2 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 6) %>%
  group_by(Colony) %>%
  mutate(BinX6 = BinX,
         Distance4.2 = Distance) %>%
  select(Colony, BinX6)

# Creating the second bin 7 reference x coordinate
Bin7.2 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 7) %>%
  group_by(Colony) %>%
  mutate(BinX7.2 = BinX) %>%
  select(Colony, BinX7.2)

# Creating the second bin 3 reference x coordinate 
Bin3 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 3) %>%
  group_by(Colony) %>%
  mutate(BinX3 = BinX,
         Distance6 = Distance) %>%
  select(Colony, BinX3)

# Creating the alternative x and y coordinates for bin 3 reference towards the entrance
Bin2 <-DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 2) %>%
  group_by(Colony) %>%
  mutate(BinY3.2 = BinX,
         BinX3.2 = BinY,
         Distance3 = Distance) %>%
  select(Colony, BinY3.2, BinX3.2, Distance3)

# Creating columns for bins 2 and 3 to create reference distances in the function below
Bin2X <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 2) %>%
  group_by(Colony) %>%
  mutate(BinX.2 = BinX) %>%
  select(Colony, BinX.2)

# Combining the two different sets of reference coordinates
Bin4.6 <- full_join(Bin4, Bin6) %>%
  full_join(Bin7) %>%
  full_join(Bin4.2) %>%
  full_join(Bin7.2) %>%
  full_join(Bin3) %>%
  full_join(Bin2X) %>%
  full_join(Bin2)

# Function that uses the data sets of worker or simulation distances to the entrance
# The function uses these distances and either adds or subtracts reference distances to obtain the shortest distance to each bin
MinDistanceTube <- function(data_table) {
  # Filtering tube nest data and joining reference coordinate and distance columns
  DistanceToBins <<- data_table %>% 
    filter(Nest == "Tube") %>%
    left_join(DistBinsFull) %>%
    left_join(Bin4.6) %>%
    left_join(DistBin1_8Full) %>%
    # Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinX3 - BinX.2), # Shortest distance from Bin 3 to 4
           Distance3_6 = abs(Distance8_7 + Distance3_4), # Shortest distance from Bin 3 to 6
           Distance7_5 = abs(Distance7_0 - Distance3_0 - Distance3_6), # Shortest distance from Bin 7 to 5
           # Creating the second x and y coordinate used to find shortest distances to each Bin
           # These reference coordinates representing distances away from the entrance
           # Second x reference coordinates 
           BinX2 = ifelse(
             # If bin 1, use the same x reference
             Bin == 1, BinX, 
             # If bin 2, BinX2 = the Bin 2 x reference - the distance from bin 3 to 4
             ifelse(Bin == 2, BinX - Distance3_4, 
                    # If bin 3, use the same x reference 
                    ifelse(Bin == 3, BinX,
                           # If bin 4, BinX2 = the bin 4 x reference + the distance from the back of the nest to bin 7
                           ifelse(Bin == 4, BinX + Distance8_7, 
                                  # If bin 5, BinX2 = the Bin 5 x reference + the distance from bin 3 to 4
                                  ifelse(Bin == 5, BinX + Distance3_4,
                                         # If bin 6, BinX2 = the bin 6 x reference + the distance from the back of the nest to bin 7
                                         ifelse(Bin == 6, BinX + Distance8_7,
                                                # If bin 7, use the same x reference
                                                ifelse(Bin == 7, BinX,
                                                       # If bin 8, use the same x reference
                                                       ifelse(Bin == 8, BinX,
                                                              # All others would be 0, but they don't exist
                                                              0)))))))),
           # Second y reference coordinates
           BinY2 = ifelse(
             # If bin 1, BinY2 = the bin 1 y reference - the distance from bin 3 to 4 - the width of the bin
             # The width of tube nest bins are 0.2723 * the length of the bin
             # Since the distance from bin 3 to 4 is the length of the tube nest Bin, multiplying this by 0.2723 gives us the width. 
             Bin == 1, BinY + Distance3_4 - (0.2723 * Distance3_4),
             # If bin 2, BinY2 = the Bin 2 y reference + the width of tube nest Bins
             ifelse(Bin == 2, BinY + (0.2723 * Distance3_4),
                    # If bin 3, BinY2 = the bin 3 y reference + the distance from Bin 8 to Bin 6
                    ifelse(Bin == 3, BinY + Distance8_7,
                           # For bin 4, 5, and 6, use the same y reference
                           ifelse(Bin == 4, BinY,
                                  ifelse(Bin == 5, BinY,
                                         ifelse(Bin == 6, BinY,
                                                # For Bin 7, BinY2 = the bin 7 y reference - the distance from the back of the nest to bin 7
                                                ifelse(Bin == 7, BinY - Distance8_7,
                                                       # If bin 8, use the same y reference
                                                       ifelse(Bin == 8, BinY,
                                                              0)))))))),
           # Creating the x and y distances from each worker to all bin x and y references
           # Distances from the worker to the occupied bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, # Distance from each worker x coordinate to the bin x reference for bin closest to the entrance in a straight line
           DistanceY = ScaledY - BinY, # Distance from each worker y coordinate to the bin y reference for bin closest to the entrance in a straight line
           DistanceY4 = ScaledY - BinY4, # Distance from each worker y coordinate to the second bin 3 y coordinate creating a direct path to bin 2
           DistanceY6 = ScaledY - BinY7, # Distance from each worker y coordinate to the second bin 6 y coordinate creating a direct path to bin 8
           DistanceX7 = ScaledX - BinX7, # Distance from each worker x coordinate to the second bin 6 and 7 x coordinate creating a direct path to bin 4 
           DistanceX4.2 = ScaledX - BinX6, # Distance from each worker x coordinate to the second bin 4 x coordinate creating a direct path to bin 6
           DistanceX6 = ScaledX - BinX7.2, # Distance from each worker x coordinate to the bin 7 x coordinate creating a direct path to bin 7
           DistanceX3 = ScaledX - BinX3, # Distance from each worker x coordinate to the second bin 3 x coordinate, creating a direct path for bins 4-6 to bin 3
           DistanceX3.2 = ScaledX - BinX3.2, # Distance from each worker x coordinate to the second bin 3 x coordinate, creating a direct path to bin 4 for bin 3 coordinates
           DistanceY3.2 = ScaledY - BinY3.2, # Distance from each worker y coordinate to the second bin 3 y coordinate, creating a direct path to bin 4 for bin 3 coordinates
           # The below script is for calculating the shortest distance from each worker to each bin
           # This uses the Pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           # A set of shortest distances are calculated for bins towards and away from the nest entrance from the reference of each worker
           # Shortest distances from each worker to Bins towards the nest entrance
           # If the worker is in bin 4 but to the left of the x reference, use the first formula, else use the second
           PythagDist = ifelse(Bin == 4 & Nest == "Tube" & ScaledX < BinX,
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               # If the worker is in bin 7 but above the y reference, use the first formula, else use the second
                               ifelse(Bin == 7 & Nest == "Tube" & ScaledY > BinY,
                                      sqrt((DistanceX7^2) + ((DistanceY^2))),
                                      # If the worker is in bin 3 but below the y reference, use the first formula, else use the second
                                      ifelse(Bin == 3 & Nest == "Tube" & ScaledY < BinY4,
                                             sqrt((DistanceX3.2^2) + (DistanceY3.2^2)),
                                             sqrt((DistanceX^2) + (DistanceY^2))))),
           
           PythagDist4 = sqrt((DistanceX7^2) + (DistanceY^2)), # Creates the shortest distance to bin 4
           PythagDist3 = sqrt((DistanceX3^2) + (DistanceY^2)), # Creates the shortest distance to bin 3
           # Distances from the worker to the occupied Bins edge towards the back of the nest 
           DistanceX2 = ScaledX - BinX2, # Alternative x axis distance towards creating the shortest distance to the closest bin facing away from the nest
           DistanceY2 = ScaledY - BinY2, # Alternative y axis distance towards creating the shortest distance to the closest bin facing away from the nest
           PythagDist2 = sqrt((DistanceX2^2) + (DistanceY2^2)), # Shortest distance to the closest bin facing away from the nest 
           PythagDist8 = sqrt((DistanceX^2) + (DistanceY6^2)), # Shortest distance from bin 6 to bin 8 
           PythagDist4.2 = sqrt((DistanceX4.2^2) + (DistanceY2^2)), # Shortest distance from bin 4 to bin 6
           PythagDist7 = sqrt((DistanceX6^2) + (DistanceY2^2)), # Shortest distance to bin 7
           # Calculating the average distance of each worker to all bins but their own
           MeanDist = 
             # Distances from workers in bin 1 to all other bins divided by 7
             ifelse(Bin == 1, (((ifelse(ScaledY > BinY3.2,
                                        PythagDist2 + Distance8_2,
                                        Distance8_1 + PythagDist2)) + # To bin 8
                                  (ifelse(ScaledY > BinY3.2, 
                                          Distance7_2 + PythagDist2,
                                          Distance7_1 + PythagDist2)) + # To bin 7
                                  (ifelse(ScaledY > BinY3.2,
                                          Distance3_4 + PythagDist2,
                                          Distance3_1 + (Distance3_4 * 2) + PythagDist2)) + # To bin 6
                                  (ifelse(ScaledY > BinY3.2,
                                          Distance3_4 + Distance8_7 + PythagDist2,
                                          Distance3_1 + Distance8_7 + PythagDist2)) + # To bin 5
                                  (ifelse(ScaledY > BinY3.2,
                                          Distance3_4 + PythagDist2,
                                          Distance3_1 + PythagDist2)) + # To bin 4
                                  (ifelse(ScaledY > BinY3.2,
                                          PythagDist2,
                                          Distance2_1 + PythagDist2)) + # To bin 3
                                  PythagDist2) / 7), # To bin 2
                    # Distances from workers in bin 2 to all other bins divided by 7
                    ifelse(Bin == 2, (((Distance8_2 + PythagDist2) + # To bin 8
                                         ((Distance7_2 + PythagDist2)) + # To bin 7
                                         (Distance7_2 - Distance7_5 + PythagDist2) + # To bin 6
                                         (Distance8_7 * 2 + PythagDist2) + # To bin 5
                                         (Distance8_7 + PythagDist2) + # To bin 4
                                         (PythagDist2) + # To bin 3
                                         (PythagDist)) / 7), # To bin 1
                           # Distances from workers in Bin 3 to all other bins divided by 7
                           ifelse(Bin == 3, (((Distance8_4 + PythagDist2) + # To bin 8
                                                (Distance7_4 + PythagDist2) + # To bin 7
                                                (Distance8_7 + Distance3_4 + PythagDist2) + # To bin 6
                                                (Distance8_7 + PythagDist2) + # To bin 5
                                                (PythagDist2) + # To bin 4
                                                (PythagDist) + # To bin 2
                                                (ifelse(ScaledY < BinY4, PythagDist,
                                                        PythagDist + Distance2_1))) / 7), # To bin 1
                                  # Distances from workers in bin 4 to all other bins divided by 7
                                  ifelse(Bin == 4, (((Distance8_7 + PythagDist7) + # To bin 8
                                                       (PythagDist7) + # To bin 7
                                                       (PythagDist4.2) + # To bin 6
                                                       (PythagDist2) + # To bin 5
                                                       (PythagDist) + # To bin 3
                                                       (ifelse(ScaledX < BinX, 
                                                               PythagDist,
                                                               PythagDist + Distance8_7)) + # To bin 2
                                                       (ifelse(ScaledX < BinX, 
                                                               PythagDist + Distance3_1,
                                                               PythagDist + Distance3_4 + Distance3_1))) / 7), # To Bin 1
                                         # Distances from workers in bin 5 to all other bins divided by 7
                                         ifelse(Bin == 5, (((Distance8_7 +  PythagDist7) + # To bin 8
                                                              (PythagDist7) + # To bin 7
                                                              (PythagDist2) + # To bin 6
                                                              (PythagDist) + # To bin 4
                                                              (PythagDist3) + # To bin 3
                                                              (PythagDist3 + Distance8_7) + # To bin 2
                                                              (PythagDist3 + Distance3_1)) / 7), # To bin 1
                                                # Distances from workers in bin 6 to all other bins divided by 7
                                                ifelse(Bin == 6, (((ifelse(ScaledX > BinX,
                                                                           PythagDist8, 
                                                                           PythagDist2 + Distance8_7)) + # To bin 8
                                                                     (PythagDist2) + # To bin 7
                                                                     (PythagDist) + # To bin 5
                                                                     (PythagDist4) + # To bin 4
                                                                     (PythagDist3) + # To bin 3
                                                                     (PythagDist3 + Distance8_7) + # To bin 2
                                                                     (PythagDist3 + Distance3_4 + Distance8_7 + Distance3_1)) / 7), # To bin 1
                                                       # Distances from workers in bin 7 to all other bins divided by 7
                                                       ifelse(Bin == 7,(((PythagDist2) + # To bin 8
                                                                           (PythagDist) + # To bin 6
                                                                           (ifelse(ScaledY > BinY,
                                                                                   PythagDist,
                                                                                   PythagDist + Distance7_5)) + # To bin 5
                                                                           (ifelse(ScaledY > BinY,
                                                                                   PythagDist,
                                                                                   PythagDist + Distance7_4)) + # To bin 4
                                                                           (ifelse(ScaledY > BinY,
                                                                                   PythagDist,
                                                                                   PythagDist + Distance7_3)) + # To bin 3
                                                                           (PythagDist + Distance7_2) + # To bin 2
                                                                           (PythagDist + Distance7_1)) / 7), # To bin 1
                                                              # Distances from workers in bin 8 to all other bins divided by 7
                                                              ifelse(Bin == 8, (((PythagDist) + # To bin 7
                                                                                   (PythagDist + Distance8_7) + # To bin 6
                                                                                   (PythagDist + Distance7_5 + Distance8_7) + # To bin 5
                                                                                   (PythagDist + Distance7_4) + # To bin 4
                                                                                   (PythagDist + Distance7_4) + # To bin 3
                                                                                   (PythagDist + Distance7_4 + Distance8_7) + # To bin 2
                                                                                   (PythagDist + Distance7_4 + Distance3_1)) / 7), # To bin 1
                                                                     0 # Else a zero value, which shouldn't exist
                                                              )))))))),
           ScaledDistMean = MeanDist / MaxDist, # Scaling worker average distance to all other bins
           # In case ScaledDistMean is greater than 1, change to 1 because the values shouldn't exist, just a precautionary measure
           ScaledDistMean = ifelse(ScaledDistMean > 1, 1, ScaledDistMean)) %>% 
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Corner, ScaledDistMean, Density) # Select desired columns
}

# Run the tube nest shape mean shortest distance to all other nest sections for the WorkerDistScaledRD1_RD2 data set
MinDistanceTube(WorkerDistScaledRD1_RD2)

# CIRCLE NEST SHAPE
# Creating the reference x and y coordinates for the sequence and joining
# First we get the circle bin reference coordinates
# Then we mirror the reference x coordinate by the center of the nest, which is 7.5cm long, therefore the center is (7.5 / 2)
# Finally we round the x and y values to the 0.01 to make joining possible later
BinCoordCircRadRedTest <- NestBinCircle %>%
  mutate(XSub = BinX - (7.5 / 2),
         BinX2 = (7.5 / 2) - XSub) %>%
  mutate_at(vars(Colony, BinID, BinX, BinX2, BinY, BinY2), funs(round(., 2))) %>%
  select(-c(Bin))

# Create a sequence of numbers along each bin, for every 0.05 cm
# The function calls a row bind that is applied to the rows of circle nest worker coordinate data set created above 
SequencesDistances <- do.call(rbind, apply(BinCoordCircRadRedTest, 1, function(x) 
  # This creates a new data frame thatcreates a sequence of numbers (every 0.05cm) between Scaled X2 and Scaled X taken from the above data set,
  # Colony is also taken from above as a reference key for future joining
  data.frame(Colony = x[1], BinX2 = x[8], BinX = x[4], Sequence = seq(x[8], x[4], by = 0.05))))

# Common columns in both the SequencesDistances and BinCoordCircRadRedTest data sets need to be the same data type (numeric) for joining
# SequencesDistances
SequencesDistances$Colony <- as.numeric(SequencesDistances$Colony)
SequencesDistances$BinX2 <- as.numeric(SequencesDistances$BinX2)
SequencesDistances$BinX <- as.numeric(SequencesDistances$BinX)

# BinCoordCircRadRedTest
BinCoordCircRadRedTest$Colony <- as.numeric(BinCoordCircRadRedTest$Colony)
BinCoordCircRadRedTest$BinX2 <- as.numeric(BinCoordCircRadRedTest$BinX2)
BinCoordCircRadRedTest$BinX <- as.numeric(BinCoordCircRadRedTest$BinX)

# Joining the final reference coordinate data set for worker shortest distance to all nest sections but their own in circle nests
BinCoordCircRadRedTestGood <- full_join(BinCoordCircRadRedTest, SequencesDistances) %>% distinct()

# Calculating the average worker shortest distances to all other nest sections 
# This function uses the worker distance to the nest entrance data set, the sequences data set created above, and two reference data sets loaded above
# The first has the x and y reference coordinates from which we use to calculate shortest distances
# The second provides the max shortest distance possible from the back of the nest to the front of the nest, allowing us to scale the shortest distances
WorkerDistScaledRD1_RD2Seq <- WorkerDistScaledRD1_RD2 %>%
  # First we filter circle nest coordinates from the worker distance to the nest entrance data set
  filter(Nest == "Circle") %>%
  group_by(Colony) %>% # Group by the Colony column
  left_join(BinCoordCircRadRedTestGood) %>% # Join the sequences data set 
  group_by(Colony, Bin, BinID) %>% # Group by the Colony, Bin and BinID columns
  # Calculating the x and y distances from each worker to every sequence point
  # There are two y reference points, because workers in their Bins have both y values above and below their Bin
  # Using a conditional ifelse(), we can use either the first or second y reference to calculate these y axis distances 
  # This prevents duplicate distances, and overall the approach works best for how the reference datasheet is formatted
  mutate(DistanceXSeq = ScaledX - Sequence, # Distance along the x axis
         DistanceYSeq = ScaledY - BinY, # First distance along the y axis (towards the entrance)
         DistanceYSeq2 = ScaledY - BinY2, # Second distance along the y axis (towards the back of the nest)
         # Conditional statement that has Pythagorean distances calculated for Bins below the one the worker resides in first and then calculates those distances above 
         # Absolute values of both x and y distances are used in each formula to prevent negative distances
         PythagDistSeq = ifelse(BinID - Bin < 2, 
                                sqrt(abs(DistanceXSeq^2) + abs(DistanceYSeq^2)), 
                                sqrt(abs(DistanceXSeq^2) + abs(DistanceYSeq2^2))),
         # Finding the smallest distance along each sequence
         MinDist = min(PythagDistSeq)) %>%
  ungroup() %>% # Ungroup the data set
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, Corner, MinDist) %>% # Select desired columns
  distinct() %>% # Remove duplicates to avoid skewing the average calculations below
  left_join(DistBinsFull) %>% # Join the distance reference data set used above throughout
  group_by(Colony, Bin) %>% # Group by the Colony and Bin columns
  # Calculated the worker scaled average distances to all other nest sections but their own
  # The average is calculated first
  # Then these averages are scaled using the max possible distance in each nest
  # Finally, as a check, all values greater than 1 are made 1. This is impossible, but keeps consistency with the methods used for the tube nest shape
  mutate(MeanDist = mean(MinDist), 
         ScaledDistMean = MeanDist / MaxDist,
         ScaledDistMean = ifelse(ScaledDistMean > 1, 1, ScaledDistMean)) %>%
  ungroup() %>% # Ungroup the data set
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, Corner, ScaledDistMean) %>% # Select desired columns
  distinct() # Remove any duplicates

# Join the final worker mean scaled distance to all nest sections but their own data set
WorkerDistScaledMeanDist <- full_join(DistanceToBins, WorkerDistScaledRD1_RD2Seq) 

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
  mutate(BinX = 3.75, # Create x reference column
         BinY = 0, # Create y reference column
         Dist1 = Distance) %>% # Set reference distance
  select(Colony, Nest ,Dist1) # Select desired columns

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
    Distance7_0 = Dist7, # Distance from bin 7 to 0
    Distance3_0 = Dist3, # Distance from bin 4 to 0 
    Distance8_7 = abs(Dist8 - Dist7), # Distance from bin 8 to 6 or 4 to 2 or 4 to 5  
    Distance8_1 = abs(Dist8 - Dist1), # Distance from bin 8 to 1
    Distance7_1 = abs(Dist7 - Dist1), # Distance from bin 7 to 1
    Distance3_1 = abs(Dist3 - Dist1), # Distance from bin 4 to 1
    Distance2_1 = abs(Dist2 - Dist1), # Distance from bin 3 to 1
    Distance8_2 = abs(Dist8 - Dist2), # Distance from bin 8 to 2
    Distance7_2 = abs(Dist7 - Dist2), # Distance from bin 7 to 2
    Distance7_3 = abs(Dist7 - Dist3), # Distance from bin 7 to 3
    Distance8_4 = abs(Dist8 - Dist3), # Distance from bin 8 to 4
    Distance7_4 = abs(Dist7 - Dist3 - Distance8_7) # Distance from bin 7 to 4
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
# The code finds each mobile colony member's distance to the brood center, then scales it
# The scale is such that the longest shortest distance in the tube nest shape (entrance-to-back) is 1
# The code first finds these distances for the tube nest, using the references calculated above
# The code then finds these distances for the circle nest, which is just the Pythagorean distance between the points
####################################################################################################################

# WORKERS
WorkerToBroodFunction <- function(data.table) {
  # Create the alternative bin 4 reference y coordinate
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% # This first section of script makes "Bin4" a column that will make the shortest distance point the same as bin 3, which corrects the issue stated above
    group_by(Colony) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(Colony, BinY4 ,Distance4)
  
  # TUBE NEST
  DistanceToBroodWorkersTube <- data.table %>% # Full worker data set
    filter(Nest == "Tube") %>% # Filter out tube nest shape
    left_join(MeanBroodCoordFullTube) %>% # Join the tube nest brood centroid data set
    left_join(DistBin1_8Full) %>% # Join the tube nest reference distance data set
    left_join(Bin4) %>%
    drop_na() %>% # Drop any NAs
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    # Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), # Shortest distance from bins 3 to 4
           # Creating the x and y distances from each individual to all Bin x and y references
           # Distances from the individual to the occupied Bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference bin 4 coordinate
           # Calculating the shortest distance from each individual to each nest section
           # This uses the Pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           # A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each individual
           # Shortest distances from each individual to bins towards the nest entrance
           # If the individual is in bin 4 but to the left of the x reference, use the first formula, else use the second
           # This is only done for nest sections towards the entrance, the rest are done using reference coordinates
           PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               sqrt((DistanceX^2) + (DistanceY^2))),
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
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2)) + Distance2_1)),
                                         ifelse(Bin == 4,
                                                ifelse(ScaledX < BinXRef4_6,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1)),
                                                # Bin 5
                                                ifelse(Bin == 5,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       # Bin 6
                                                       ifelse(Bin == 6,
                                                              ifelse(BroodY > BinYRef2,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                              # Bin 7
                                                              ifelse(Bin == 7,
                                                                     ifelse(ScaledY > BinYRef7,
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
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
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
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
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
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
                                                                                          PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                   ifelse(Bin == 8,
                                                                                          PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          NA)))))))),
                                  # Distances from workers to brood centroids in bin 4
                                  ifelse(BroodBin == 4,
                                         # Bin 1
                                         ifelse(Bin == 1,
                                                ifelse(BroodX < BinXRef4_6,
                                                       ifelse(ScaledY > BinYRef2, 
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
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
                                                                                                 PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)), 
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
                                                                                                        PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
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
                                                                                                               PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
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
                                                                                          sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                          # Bin 5
                                                                                          ifelse(Bin == 5,
                                                                                                 sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                 # Bin 6
                                                                                                 ifelse(Bin == 6,
                                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                        # Bin 7
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
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, Corner, Bin, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # CIRCLE NEST
  DistanceToBroodWorkersCircle <- data.table %>% # Full worker data set 
    filter(Nest == "Circle") %>% # Filtering out circle nest coordinates
    left_join(MeanBroodCoordFullCircle) %>% # Joining the circle brood centroid data set
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)), # Shortest distance from each individual to the brood center
           ToBrood = BroodDist / MaxDist) %>% # Select the desired columns 
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Bin, Density, Corner, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # JOINING AND CREATING THE FINAL DATASETS
  BroodCentDistWorkersRD1_RD2 <<- full_join(DistanceToBroodWorkersTube, DistanceToBroodWorkersCircle)
  
  # Sub-setting the high density treatment worker data
  BroodCentDistWorkersRD1 <<- BroodCentDistWorkersRD1_RD2 %>%
    filter(Colony < 11)
  
  # Sub-setting the low density treatment worker data
  BroodCentDistWorkersRD2 <<- BroodCentDistWorkersRD1_RD2 %>%
    filter(Colony > 10)
}

# Run the worker distance to the brood center function for the FullDataCoordWorkersRD1_RD2 data set
WorkerToBroodFunction(FullDataCoordWorkersRD1_RD2)

# QUEENS
QueenToBroodFunction <- function(data.table){
  # Create the alternative bin 4 reference y coordinate
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% # This first section of script makes "Bin4" a column that will make the shortest distance point the same as bin 3, which corrects the issue stated above
    group_by(Colony) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(Colony, BinY4 ,Distance4)
  
  # TUBE NEST
  DistanceToBroodQueenTube <- data.table %>% # Full worker data set
    filter(Nest == "Tube") %>% # Filter out tube nest shape
    left_join(MeanBroodCoordFullTube) %>% # Join the tube nest brood centroid data set
    left_join(DistBin1_8Full) %>% # Join the tube nest reference distance data set
    left_join(Bin4) %>%
    drop_na() %>% # Drop any NAs
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    # Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), # Shortest distance from bins 3 to 4
           # Creating the x and y distances from each individual to all Bin x and y references
           # Distances from the individual to the occupied Bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference bin 4 coordinate
           # Calculating the shortest distance from each individual to each nest section
           # This uses the Pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           # A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each individual
           # Shortest distances from each individual to bins towards the nest entrance
           # If the individual is in bin 4 but to the left of the x reference, use the first formula, else use the second
           # This is only done for nest sections towards the entrance, the rest are done using reference coordinates
           PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               sqrt((DistanceX^2) + (DistanceY^2))),
           # Calculating the distance of each individual to the brood center
           # If the individual and brood centroid can be connected with a straight line, then this is done with the Pythagorean theorem 
           # If they aren't, the shortest distance to the bin closest to each other is calculated for each, than a reference distance is added where needed
           BroodDist = 
             # Distances from queens to brood centroids in bin 1
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
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2)) + Distance2_1)),
                                         ifelse(Bin == 4,
                                                ifelse(ScaledX < BinXRef4_6,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1)),
                                                # Bin 5
                                                ifelse(Bin == 5,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       # Bin 6
                                                       ifelse(Bin == 6,
                                                              ifelse(BroodY > BinYRef2,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                              # Bin 7
                                                              ifelse(Bin == 7,
                                                                     ifelse(ScaledY > BinYRef7,
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
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
                    # Distances from queens to brood centroids in bin 2
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
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                            # Bin 8
                                                                            ifelse(Bin == 8,
                                                                                   PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   NA)))))))),
                           # Distances from queens to brood centroids in bin 3
                           ifelse(BroodBin == 3,
                                  # Bin 1
                                  ifelse(Bin == 1,
                                         ifelse(BroodY < BinYRef3,
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
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
                                                                                          PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                   ifelse(Bin == 8,
                                                                                          PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          NA)))))))),
                                  # Distances from queens to brood centroids in bin 4
                                  ifelse(BroodBin == 4,
                                         # Bin 1
                                         ifelse(Bin == 1,
                                                ifelse(BroodX < BinXRef4_6,
                                                       ifelse(ScaledY > BinYRef2, 
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
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
                                                                                                 PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)), 
                                                                                                 NA)))))))),
                                         # Distances from queens to brood centroids in bin 5
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
                                                                                                        PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                        NA)))))))),
                                                # Distances from queens to brood centroids in bin 6
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
                                                                                                               PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                               NA)))))))),
                                                       # Distances from queens to brood centroids in bin 8, note that there weren't any brood centroids in bin 7, so it is skipped here
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
                                                                                          sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                          # Bin 5
                                                                                          ifelse(Bin == 5,
                                                                                                 sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                 # Bin 6
                                                                                                 ifelse(Bin == 6,
                                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                        # Bin 7
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
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, Corner, Bin, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # CIRCLE NEST
  DistanceToBroodQueenCircle <- data.table %>% # Full worker data set 
    filter(Nest == "Circle") %>% # Filtering out circle nest coordinates
    left_join(MeanBroodCoordFullCircle) %>% # Joining the circle brood centroid data set
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)), # Shortest distance from each individual to the brood center
           ToBrood = BroodDist / MaxDist) %>% # Select the desired columns 
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, Corner, Bin, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # JOINING AND CREATING THE FINAL DATASETS
  BroodCentDistQueensRD1_RD2 <<- full_join(DistanceToBroodQueenTube, DistanceToBroodQueenCircle)
  
  # Sub-setting the high density treatment worker data
  BroodCentDistQueensRD1 <<- BroodCentDistQueensRD1_RD2 %>%
    filter(Colony < 11)
  
  # Sub-setting the low density treatment worker data
  BroodCentDistQueensRD2 <<- BroodCentDistQueensRD1_RD2 %>%
    filter(Colony > 10)
}

# Run the queen distance to the brood center function for the FullDataCoordQueenRD1_RD2 data set
QueenToBroodFunction(FullDataCoordQueenRD1_RD2)

# ALATES
AlateToBroodFunction <- function(data.table){
  # Create the alternative bin 4 reference y coordinate
  Bin4 <- DistBinsFull %>% # Bin 4 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left
    filter(Nest == "Tube" & Bin == 3) %>% # This first section of script makes "Bin4" a column that will make the shortest distance point the same as bin 3, which corrects the issue stated above
    group_by(Colony) %>% 
    mutate(BinY4 = BinY,
           Distance4 = Distance) %>% # The distance that is added to the Pythagorean distance is also different
    select(Colony, BinY4 ,Distance4)
  
  # TUBE NEST
  DistanceToBroodAlateTube <- data.table %>% # Full worker data set
    filter(Nest == "Tube") %>% # Filter out tube nest shape
    left_join(MeanBroodCoordFullTube) %>% # Join the tube nest brood centroid data set
    left_join(DistBin1_8Full) %>% # Join the tube nest reference distance data set
    left_join(Bin4) %>%
    drop_na() %>% # Drop any NAs
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    # Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), # Shortest distance from bins 3 to 4
           # Creating the x and y distances from each individual to all Bin x and y references
           # Distances from the individual to the occupied Bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, # Distance from each individual's x coordinate to the x reference bin coordinate
           DistanceY = ScaledY - BinY, # Distance from each individual's y coordinate to the y reference bin coordinate
           DistanceY4 = ScaledY - BinY4, # Distance from each individual's y coordinate to the y reference bin 4 coordinate
           # Calculating the shortest distance from each individual to each nest section
           # This uses the Pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           # A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each individual
           # Shortest distances from each individual to bins towards the nest entrance
           # If the individual is in bin 4 but to the left of the x reference, use the first formula, else use the second
           # This is only done for nest sections towards the entrance, the rest are done using reference coordinates
           PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               sqrt((DistanceX^2) + (DistanceY^2))),
           # Calculating the distance of each individual to the brood center
           # If the individual and brood centroid can be connected with a straight line, then this is done with the Pythagorean theorem 
           # If they aren't, the shortest distance to the bin closest to each other is calculated for each, than a reference distance is added where needed
           BroodDist = 
             # Distances from alates to brood centroids in bin 1
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
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2)) + Distance2_1)),
                                         ifelse(Bin == 4,
                                                ifelse(ScaledX < BinXRef4_6,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1)),
                                                # Bin 5
                                                ifelse(Bin == 5,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       # Bin 6
                                                       ifelse(Bin == 6,
                                                              ifelse(BroodY > BinYRef2,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                              # Bin 7
                                                              ifelse(Bin == 7,
                                                                     ifelse(ScaledY > BinYRef7,
                                                                            ifelse(BroodY > BinYRef2,
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
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
                    # Distances from alates to brood centroids in bin 2
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
                                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                            # Bin 8
                                                                            ifelse(Bin == 8,
                                                                                   PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                   NA)))))))),
                           # Distances from alates to brood centroids in bin 3
                           ifelse(BroodBin == 3,
                                  # Bin 1
                                  ifelse(Bin == 1,
                                         ifelse(BroodY < BinYRef3,
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
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
                                                                                          PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                   ifelse(Bin == 8,
                                                                                          PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                          NA)))))))),
                                  # Distances from alates to brood centroids in bin 4
                                  ifelse(BroodBin == 4,
                                         # Bin 1
                                         ifelse(Bin == 1,
                                                ifelse(BroodX < BinXRef4_6,
                                                       ifelse(ScaledY > BinYRef2, 
                                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
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
                                                                                                 PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)), 
                                                                                                 NA)))))))),
                                         # Distances from alates to brood centroids in bin 5
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
                                                                                                        PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                        NA)))))))),
                                                # Distances from alates to brood centroids in bin 6
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
                                                                                                               PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                               NA)))))))),
                                                       # Distances from alates to brood centroids in bin 8, note that there weren't any brood centroids in bin 7, so it is skipped here
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
                                                                                          sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                          # Bin 5
                                                                                          ifelse(Bin == 5,
                                                                                                 sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                 # Bin 6
                                                                                                 ifelse(Bin == 6,
                                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                        # Bin 7
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
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, Corner, Bin, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # CIRCLE NEST
  DistanceToBroodAlateCircle <- data.table %>% # Full worker data set 
    filter(Nest == "Circle") %>% # Filtering out circle nest coordinates
    left_join(MeanBroodCoordFullCircle) %>% # Joining the circle brood centroid data set
    group_by(Colony, Day) %>% # Group by the Colony and Day columns
    mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)), # Shortest distance from each individual to the brood center
           ToBrood = BroodDist / MaxDist) %>% # Select the desired columns 
    left_join(CornerFull) %>% # Join with data set that associates nest sections with the presence of a corner or not.
    select(Colony, Nest, Day, ScaledX, ScaledY, Density, Corner, Bin, ToBrood) %>% # Select the desired columns
    distinct() %>% # Remove any duplicates (shouldn't exist)
    drop_na() # Remove NAs
  
  # JOINING AND CREATING THE FINAL DATASETS
  BroodCentDistAlatesRD2 <<- full_join(DistanceToBroodAlateTube, DistanceToBroodAlateCircle)
}

# Run the alate distance to the brood center function for the FullDataCoordAlate data set
AlateToBroodFunction(FullDataCoordAlate)
