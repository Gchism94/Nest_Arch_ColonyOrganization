#DISTANCE TO BROOD CENTER CALCULATION
#NOTE that observations for brood don't line up perfectly with the workers, queens, and alates
#This means that brood centers do not exist for those instances 
#When this was found, a different day was used to make the sample size 10, but another observation for workers was not created 
#The following calculates colony member distances from the brood center
#Worker to brood distance RD1
#High density treatment

MeanBroodCoordProps <- AntPropFullBroodRD1_RD2 %>%
  group_by(Colony, Nest, Day) %>%
  mutate(MaxBrood = max(PropBrood), MaxBin = ifelse(PropBrood == MaxBrood, Bin, NA)) %>%
  ungroup() %>%
  drop_na() %>%
  select(Colony, Nest, Density, Day, Bin)

MeanBroodCoordFull <- left_join(MeanBroodCoordProps, BroodDistScaledRD1_RD2) %>%
  drop_na() %>%
  group_by(Colony, Nest, Day) %>%
  mutate(BroodX = mean(ScaledX), BroodY = mean(ScaledY)) %>%
  rename(BroodBin = Bin) %>%
  select(Colony, Nest, Day, BroodBin, BroodX, BroodY, Density) %>%
  distinct()

TubeDistances <- DistBinsFull %>% #Filtering out the Tube reference distances
  filter(Nest == "Tube")

#Filtering out only distances for nest section 1
DistBin1 <- TubeDistances %>% 
  filter(Bin == 1) %>%
  mutate(Dist1 = Distance) %>% #Set reference distance
  select(Colony, Nest ,Dist1) #Select desired columns

#Filtering and creating reference distances for nest sections 2 - 8
DistBin2 <- TubeDistances %>%
  filter(Bin == 2) %>%
  mutate(Dist2 = Distance) %>%
  select(Colony, Nest, Dist2)

DistBin3 <- TubeDistances %>%
  filter(Bin == 3) %>%
  mutate(Dist3 = Distance) %>%
  select(Colony, Nest, Dist3)

DistBin4 <- TubeDistances %>%
  filter(Bin == 4) %>%
  mutate(Dist4 = Distance) %>%
  select(Colony, Nest, Dist4)

DistBin5 <- TubeDistances %>%
  filter(Bin == 5) %>%
  mutate(Dist5 = Distance) %>%
  select(Colony, Nest, Dist5)

DistBin6 <- TubeDistances %>%
  filter(Bin == 6) %>%
  mutate(Dist6 = Distance) %>%
  select(Colony, Nest, Dist6)

DistBin7 <- TubeDistances %>%
  filter(Bin == 7) %>%
  mutate(Dist7 = Distance) %>%
  select(Colony, Nest, Dist7)

DistBin8 <- TubeDistances %>%
  filter(Bin == 8) %>%
  mutate(Dist8 = Distance) %>%
  select(Colony, Nest, Dist8)

#Create a combined dataset of reference coordinates and distances
DistBin1_8Full <- full_join(DistBin1, DistBin2) %>%
  full_join(DistBin3) %>%
  full_join(DistBin4) %>%
  full_join(DistBin5) %>%
  full_join(DistBin6) %>%
  full_join(DistBin7) %>%
  full_join(DistBin8) %>%
  mutate(
    Distance7_0 = Dist7, #Distance from bin 7 to 0
    Distance3_0 = Dist3, #Distance from bin 4 to 0 
    Distance8_7 = abs(Dist8 - Dist7), #Distance from bin 8 to 6 or 4 to 2 or 4 to 5  
    Distance8_1 = abs(Dist8 - Dist1), #Distance from bin 8 to 1
    Distance7_1 = abs(Dist7 - Dist1), #Distance from bin 7 to 1
    Distance3_1 = abs(Dist3 - Dist1), #Distance from bin 4 to 1
    Distance2_1 = abs(Dist2 - Dist1), #Distance from bin 3 to 1
    Distance8_2 = abs(Dist8 - Dist2), #Distance from bin 8 to 2
    Distance7_2 = abs(Dist7 - Dist2), #Distance from bin 7 to 2
    Distance7_3 = abs(Dist7 - Dist3), #Distance from bin 7 to 3
    Distance8_4 = abs(Dist8 - Dist3), #Distance from bin 8 to 4
    Distance7_4 = abs(Dist7 - Dist3 - Distance8_7) #Distance from bin 7 to 4
  )
#Bin 4 and 7 has two possible shortest distances, one on the right past the farthest x coordinate of bin 3, and one to the left,
#which we account for this with two sets of possible reference coordinates

#Creating the second nest section 4 reference y coordinate
Bin4 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin==3) %>%
  group_by(Colony) %>%
  mutate(BinY4 = BinY) %>%
  select(Colony, BinY4)

#Creating the second nest section 7 reference y coordinate
Bin6 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 8) %>%
  group_by(Colony)%>%
  mutate(BinY7 = BinY) %>%
  select(Colony, BinY7)

#Creating the second nest section 4 reference x coordinate
Bin7 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 4) %>%
  group_by(Colony) %>%
  mutate(BinX7 = BinX) %>%
  select(Colony, BinX7)

#Creating the second nest section 6 reference x coordinate
Bin4.2 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 6) %>%
  group_by(Colony) %>%
  mutate(BinX6 = BinX,
         Distance6 = Distance) %>%
  select(Colony, BinX6)

#Creating the second nest section 7 reference x coordinate
Bin7.2 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 7) %>%
  group_by(Colony) %>%
  mutate(BinX7.2 = BinX) %>%
  select(Colony, BinX7.2)

#Creating the second nest section 6 reference x coordinate 
Bin3 <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 3) %>%
  group_by(Colony) %>%
  mutate(BinX3 = BinX,
         Distance6 = Distance) %>%
  select(Colony, BinX3)

#Creating columns for nest sections 2 and 3 to create reference distances in the function below
Bin2X <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 2) %>%
  group_by(Colony) %>%
  mutate(BinX.2 = BinX) %>%
  select(Colony, BinX.2)

Bin3X <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 3) %>%
  group_by(Colony) %>%
  mutate(BinX3 = BinX) %>%
  select(Colony, BinX3)

Bin1X <- DistBinsFull %>%
  filter(Nest == "Tube" & Bin == 1) %>%
  group_by(Colony) %>%
  mutate(BinX.1 = BinX) %>%
  select(Colony, BinX.1)

#Combining the two diffrent sets of reference coordinates
#they are for separate purposes so they will remain separate
Bin4.6 <- full_join(Bin4, Bin6) %>%
  full_join(Bin7) %>%
  full_join(Bin4.2) %>%
  full_join(Bin7.2) %>%
  full_join(Bin3)

Bin2.3 <- full_join(Bin2X, Bin3X) %>%
  full_join(Bin1X)

#Function that uses the datasets of worker or simulation distances to the entrance
#The function uses these distances and either adds or subtracts reference distances to obtain the shortest distance to each nest section

DistBinsTubeRefFull <- DistBinsFull %>%
    filter(Nest == "Tube") %>%
    select(-c(Distance, Xmax, Ymax, TubeRatio))

  DistBinsTubeRef1 <- DistBinsTubeRefFull %>%
    filter(Nest == "Tube" & Bin == 1) %>%
    rename(BinXRef1 = BinX, BinYRef1 = BinY) %>%
    select(-c(Bin))

  DistBinsTubeRef2 <- DistBinsTubeRefFull %>%
    filter(Nest == "Tube" & Bin == 2) %>%
    rename(BinXRef2 = BinX, BinYRef2 = BinY) %>%
    select(-c(Bin))

  DistBinsTubeRef3 <- DistBinsTubeRefFull %>%
    filter(Nest == "Tube" & Bin == 3) %>%
    rename(BinXRef3 = BinX, BinYRef3 = BinY)%>%
    select(-c(Bin))

  DistBinsTubeRef4_6 <- DistBinsTubeRefFull %>%
    filter(Nest == "Tube" & Bin == 4) %>%
    rename(BinXRef4_6 = BinX, BinYRef4_6 = BinY)%>%
    select(-c(Bin))

  DistBinsTubeRef7 <- DistBinsTubeRefFull %>%
    filter(Nest == "Tube" & Bin == 7) %>%
    rename(BinXRef7 = BinX, BinYRef7 = BinY)%>%
    select(-c(Bin))

  DistBinsTubeRef8 <- DistBinsTubeRefFull %>%
      filter(Nest == "Tube" & Bin == 8) %>%
      rename(BinXRef8 = BinX, BinYRef8 = BinY)%>%
    select(-c(Bin))

DistBinsTubeRef <- left_join(DistBinsTubeRefFull,DistBinsTubeRef1) %>%
    left_join(DistBinsTubeRef2) %>%
    left_join(DistBinsTubeRef3) %>%
    left_join(DistBinsTubeRef4_6) %>%
    left_join(DistBinsTubeRef7) %>%
    left_join(DistBinsTubeRef8)

  #Filtering tube nest data and joining reference coordinate and distance columns
MeanBroodCoordFullTube <- MeanBroodCoordFull %>% 
    filter(Nest == "Tube") %>%
    left_join(DistBinsTubeRef)


DistanceToBroodWorkersTube <- WorkerDistScaledRD1_RD2 %>% 
  filter(Nest == "Tube") %>%
  left_join(MeanBroodCoordFullTube) %>%
  left_join(DistBin1_8Full) %>%
  drop_na() %>%
  group_by(Colony, Day) %>%
    #Creating columns of reference distances from both the created distances and reference coordinates above 
    mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), #Shortest distance from Bin 3 to 4
           #Creating the x and y distances from each worker to all Bin x and y references
           #Distances from the worker to the occupied Bins edge towards the nest entrance
           DistanceX = ScaledX - BinX, #Distance from each worker x coordinate to the Bin x reference towards the entrance
           DistanceY = ScaledY - BinY, #Distance from each worker y coordinate to the Bin y reference towards the entrance
           DistanceY4 = ScaledY - BinYRef3, 
           #Calculating the shortest distance from each worker to each nest section
           #This uses the pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
           #A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each worker
           #Shortest distances from each worker to Bins towards the nest entrance
           #If the worker is in Bin 4 but to the left of the x reference, use the first formula, else use the second
           PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                               sqrt((DistanceX^2) + (DistanceY4^2)),
                               sqrt((DistanceX^2) + (DistanceY^2))),
           #Calculating the average distance of each worker to all Bins
           BroodDist = 
             #Distances from workers in Bin 1 to all other Bins divided by 7
             ifelse(BroodBin == 1,
                    ifelse(Bin == 1, sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
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
                                                ifelse(Bin == 5,
                                                       ifelse(BroodY > BinYRef2,
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                              PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       ifelse(Bin == 6,
                                                              ifelse(BroodY > BinYRef2,
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                     PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                       ifelse(Bin == 7,
                                                              ifelse(ScaledY > BinYRef7,
                                                                     ifelse(BroodY > BinYRef2,
                                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                            PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                                     ifelse(BroodY > BinYRef2,
                                                                            PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                            PythagDist + Distance7_1 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)))),
                                                              ifelse(Bin == 8,
                                                                     ifelse(BroodY > BinYRef2,
                                                                            PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                            PythagDist + Distance8_2 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2))),
                                                                     NA)
                                                              ))))))),
           ifelse(BroodBin == 2,
                  ifelse(Bin == 1, 
                         ifelse(ScaledY > BinYRef2,
                                sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                         ifelse(Bin == 2,
                                sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                ifelse(Bin == 3,
                                       ifelse(ScaledY < BinYRef3,
                                              sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                       ifelse(Bin == 4,
                                              ifelse(ScaledX < BinXRef4_6,
                                                     PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                     PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                     ifelse(Bin == 5,
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            ifelse(Bin == 6,
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            ifelse(Bin == 7,
                                                                   ifelse(ScaledY > BinYRef7,
                                                                          PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                          PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                   ifelse(Bin == 8,
                                                                          PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                          NA)))))))),
                  ifelse(BroodBin == 3,
                         ifelse(Bin == 1,
                                ifelse(BroodY < BinYRef3,
                                       ifelse(ScaledY > BinYRef2,
                                              sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                       ifelse(ScaledY > BinYRef2,
                                              sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                              sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
                                       ifelse(Bin == 2,
                                              ifelse(BroodY < BinYRef3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                              ifelse(Bin == 3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     ifelse(Bin == 4,
                                                            ifelse(ScaledX < BinXRef4_6,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                            ifelse(Bin == 5,
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 6,
                                                                          PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 7,
                                                                          ifelse(ScaledY > BinYRef7,
                                                                                 PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                 PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                          ifelse(Bin == 8,
                                                                                 PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                 NA)))))))),
                         ifelse(BroodBin == 4,
                                ifelse(Bin == 1,
                                       ifelse(
                                         BroodX < BinXRef4_6,
                                         ifelse(ScaledY > BinYRef2, 
                                                sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
                                         ifelse(ScaledY > BinYRef2,
                                                sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1)),
                                       ifelse(Bin == 2,
                                              ifelse(
                                                BroodX < BinXRef4_6,
                                                sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4),
                                              ifelse(Bin == 3,
                                                     ifelse(BroodX < BinXRef4_6,
                                                            sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                            sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                     ifelse(Bin == 4,
                                                            sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                            ifelse(Bin == 5,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   ifelse(Bin == 6,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                            ifelse(Bin == 7,
                                                                   ifelse(ScaledY > BinYRef7,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                   ifelse(Bin == 8,
                                                                          PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)), #OKAY
                                                                          NA)))))))),
                                ifelse(BroodBin == 5,
                                       ifelse(Bin == 1,
                                              ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                              ifelse(Bin == 2,
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                     ifelse(Bin == 3,
                                                            sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                            ifelse(Bin == 4,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   ifelse(Bin == 5,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 6,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   ifelse(Bin == 7,
                                                                          ifelse(ScaledY > BinYRef7,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                          ifelse(Bin == 8,
                                                                                 PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                 NA)))))))),
                                       ifelse(BroodBin == 6,
                                              ifelse(Bin == 1,
                                                     ifelse(ScaledY > BinYRef2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                     ifelse(Bin == 2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            ifelse(Bin == 3,
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 4,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 5,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 6,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 7,
                                                                                 ifelse(ScaledY > BinYRef7,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                 ifelse(Bin == 8,
                                                                                        PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                        NA)))))))),
                                       ifelse(BroodBin == 8,
                                              ifelse(Bin == 1,
                                                ifelse(ScaledY > BinYRef2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                            sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_1),
                                                     ifelse(Bin == 2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                            ifelse(Bin == 3,
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_4,
                                                                   ifelse(Bin == 4,
                                                                          sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                          ifelse(Bin == 5,
                                                                                 sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                 ifelse(Bin == 6,
                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                          ifelse(Bin == 7,
                                                                                 ifelse(ScaledY < BinYRef8,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        sqrt(((ScaledX - BinXRef8)^2) + ((ScaledY - BinYRef8)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2))),
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
  group_by(Colony) %>%
  mutate(ToBrood = BroodDist / MaxDist,
         ToBrood = ifelse(ToBrood > 1, 1, ToBrood)) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, ScaledDist, Density, Corner, Bin, ColonyMember, ToBrood) %>%
  distinct()
  drop_na()

DistanceToBroodQueenTube <- QueenDistScaledRD1_RD2 %>% 
  filter(Nest == "Tube") %>%
  left_join(MeanBroodCoordFullTube) %>%
  left_join(DistBin1_8Full) %>%
  drop_na() %>%
  group_by(Colony, Day) %>%
  #Creating columns of reference distances from both the created distances and reference coordinates above 
  mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), #Shortest distance from Bin 3 to 4
         #Creating the x and y distances from each worker to all Bin x and y references
         #Distances from the worker to the occupied Bins edge towards the nest entrance
         DistanceX = ScaledX - BinX, #Distance from each worker x coordinate to the Bin x reference towards the entrance
         DistanceY = ScaledY - BinY, #Distance from each worker y coordinate to the Bin y reference towards the entrance
         DistanceY4 = ScaledY - BinYRef3, 
         #Calculating the shortest distance from each worker to each nest section
         #This uses the pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
         #A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each worker
         #Shortest distances from each worker to Bins towards the nest entrance
         #If the worker is in Bin 4 but to the left of the x reference, use the first formula, else use the second
         PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                             sqrt((DistanceX^2) + (DistanceY4^2)),
                             sqrt((DistanceX^2) + (DistanceY^2))),
         #Calculating the average distance of each worker to all Bins
         BroodDist = 
           #Distances from workers in Bin 1 to all other Bins divided by 7
           ifelse(BroodBin == 1,
                  ifelse(Bin == 1, sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
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
                                              ifelse(Bin == 5,
                                                     ifelse(BroodY > BinYRef2,
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                     ifelse(Bin == 6,
                                                            ifelse(BroodY > BinYRef2,
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                            ifelse(Bin == 7,
                                                                   ifelse(ScaledY > BinYRef7,
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance7_1 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)))),
                                                                   ifelse(Bin == 8,
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2))),
                                                                          NA)
                                                            ))))))),
                  ifelse(BroodBin == 2,
                         ifelse(Bin == 1, 
                                ifelse(ScaledY > BinYRef2,
                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                ifelse(Bin == 2,
                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                       ifelse(Bin == 3,
                                              ifelse(ScaledY < BinYRef3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                              ifelse(Bin == 4,
                                                     ifelse(ScaledX < BinXRef4_6,
                                                            PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                     ifelse(Bin == 5,
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            ifelse(Bin == 6,
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                   ifelse(Bin == 7,
                                                                          ifelse(ScaledY > BinYRef7,
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                          ifelse(Bin == 8,
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 NA)))))))),
                         ifelse(BroodBin == 3,
                                ifelse(Bin == 1,
                                       ifelse(BroodY < BinYRef3,
                                              ifelse(ScaledY > BinYRef2,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                              ifelse(ScaledY > BinYRef2,
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                     sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
                                       ifelse(Bin == 2,
                                              ifelse(BroodY < BinYRef3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                              ifelse(Bin == 3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     ifelse(Bin == 4,
                                                            ifelse(ScaledX < BinXRef4_6,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                            ifelse(Bin == 5,
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 6,
                                                                          PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                          ifelse(Bin == 7,
                                                                                 ifelse(ScaledY > BinYRef7,
                                                                                        PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                        PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                 ifelse(Bin == 8,
                                                                                        PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                        NA)))))))),
                                ifelse(BroodBin == 4,
                                       ifelse(Bin == 1,
                                              ifelse(
                                                BroodX < BinXRef4_6,
                                                ifelse(ScaledY > BinYRef2, 
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1)),
                                              ifelse(Bin == 2,
                                                     ifelse(
                                                       BroodX < BinXRef4_6,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4),
                                                     ifelse(Bin == 3,
                                                            ifelse(BroodX < BinXRef4_6,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                            ifelse(Bin == 4,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   ifelse(Bin == 5,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 6,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 7,
                                                                                        ifelse(ScaledY > BinYRef7,
                                                                                               sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                               PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                        ifelse(Bin == 8,
                                                                                               PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                               NA)))))))),
                                       ifelse(BroodBin == 5,
                                              ifelse(Bin == 1,
                                                     ifelse(ScaledY > BinYRef2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                     ifelse(Bin == 2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            ifelse(Bin == 3,
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 4,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 5,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 6,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        ifelse(Bin == 7,
                                                                                               ifelse(ScaledY > BinYRef7,
                                                                                                      sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                      PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                               ifelse(Bin == 8,
                                                                                                      PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                      NA)))))))),
                                              ifelse(BroodBin == 6,
                                                     ifelse(Bin == 1,
                                                            ifelse(ScaledY > BinYRef2,
                                                                   sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                   sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                            ifelse(Bin == 2,
                                                                   sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                   ifelse(Bin == 3,
                                                                          sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                          ifelse(Bin == 4,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 5,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        ifelse(Bin == 6,
                                                                                               sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                               ifelse(Bin == 7,
                                                                                                      ifelse(ScaledY > BinYRef7,
                                                                                                             sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                             PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                                      ifelse(Bin == 8,
                                                                                                             PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                             NA)))))))),
                                                     ifelse(BroodBin == 8,
                                                            ifelse(Bin == 1,
                                                                   ifelse(ScaledY > BinYRef2,
                                                                          sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                          sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_1),
                                                                   ifelse(Bin == 2,
                                                                          sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                          ifelse(Bin == 3,
                                                                                 sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_4,
                                                                                 ifelse(Bin == 4,
                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                        ifelse(Bin == 5,
                                                                                               sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                               ifelse(Bin == 6,
                                                                                                      sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                      ifelse(Bin == 7,
                                                                                                             ifelse(ScaledY < BinYRef8,
                                                                                                                    sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                                    sqrt(((ScaledX - BinXRef8)^2) + ((ScaledY - BinYRef8)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2))),
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
  group_by(Colony) %>%
  mutate(ToBrood = BroodDist / MaxDist,
         ToBrood = ifelse(ToBrood > 1, 1, ToBrood)) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, ScaledDist, Density, Corner, Bin, ToBrood) %>%
  distinct() %>%
  drop_na()

DistanceToBroodAlateTube <- AlateDistScaledRD2 %>% 
  filter(Nest == "Tube") %>%
  left_join(MeanBroodCoordFullTube) %>%
  left_join(DistBin1_8Full) %>%
  drop_na() %>%
  group_by(Colony, Day) %>%
  #Creating columns of reference distances from both the created distances and reference coordinates above 
  mutate(Distance3_4 = abs(BinXRef3 - BinXRef4_6), #Shortest distance from Bin 3 to 4
         #Creating the x and y distances from each worker to all Bin x and y references
         #Distances from the worker to the occupied Bins edge towards the nest entrance
         DistanceX = ScaledX - BinX, #Distance from each worker x coordinate to the Bin x reference towards the entrance
         DistanceY = ScaledY - BinY, #Distance from each worker y coordinate to the Bin y reference towards the entrance
         DistanceY4 = ScaledY - BinYRef3, 
         #Calculating the shortest distance from each worker to each nest section
         #This uses the pythagorean theorem, which finds the hypotenuse through the formula sqrt((DistanceX^2) + (DistanceY^2))
         #A set of shortest distances are calculated for Bins towards and away from the nest entrance from the reference of each worker
         #Shortest distances from each worker to Bins towards the nest entrance
         #If the worker is in Bin 4 but to the left of the x reference, use the first formula, else use the second
         PythagDist = ifelse(Bin == 4 & ScaledX < BinX, 
                             sqrt((DistanceX^2) + (DistanceY4^2)),
                             sqrt((DistanceX^2) + (DistanceY^2))),
         #Calculating the average distance of each worker to all Bins
         BroodDist = 
           #Distances from workers in Bin 1 to all other Bins divided by 7
           ifelse(BroodBin == 1,
                  ifelse(Bin == 1, sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
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
                                              ifelse(Bin == 5,
                                                     ifelse(BroodY > BinYRef2,
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                     ifelse(Bin == 6,
                                                            ifelse(BroodY > BinYRef2,
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                            ifelse(Bin == 7,
                                                                   ifelse(ScaledY > BinYRef7,
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)) + Distance2_1),
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance7_1 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2)))),
                                                                   ifelse(Bin == 8,
                                                                          ifelse(BroodY > BinYRef2,
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef2 - BroodX)^2) + ((BinYRef2 - BroodY)^2))),
                                                                          NA)
                                                            ))))))),
                  ifelse(BroodBin == 2,
                         ifelse(Bin == 1, 
                                ifelse(ScaledY > BinYRef2,
                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                ifelse(Bin == 2,
                                       sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                       ifelse(Bin == 3,
                                              ifelse(ScaledY < BinYRef3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                              ifelse(Bin == 4,
                                                     ifelse(ScaledX < BinXRef4_6,
                                                            PythagDist + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                     ifelse(Bin == 5,
                                                            PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                            ifelse(Bin == 6,
                                                                   PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                   ifelse(Bin == 7,
                                                                          ifelse(ScaledY > BinYRef7,
                                                                                 PythagDist + Distance3_4 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 PythagDist + Distance7_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2))),
                                                                          ifelse(Bin == 8,
                                                                                 PythagDist + Distance8_2 + sqrt(((BinXRef3 - BroodX)^2) + ((BinYRef3 - BroodY)^2)),
                                                                                 NA)))))))),
                         ifelse(BroodBin == 3,
                                ifelse(Bin == 1,
                                       ifelse(BroodY < BinYRef3,
                                              ifelse(ScaledY > BinYRef2,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef2)^2) + ((BroodY - BinYRef2)^2))),
                                              ifelse(ScaledY > BinYRef2,
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                     sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1)),
                                       ifelse(Bin == 2,
                                              ifelse(BroodY < BinYRef3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2))),
                                              ifelse(Bin == 3,
                                                     sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                     ifelse(Bin == 4,
                                                            ifelse(ScaledX < BinXRef4_6,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                            ifelse(Bin == 5,
                                                                   PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 6,
                                                                          PythagDist + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                          ifelse(Bin == 7,
                                                                                 ifelse(ScaledY > BinYRef7,
                                                                                        PythagDist + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                        PythagDist + Distance7_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2))),
                                                                                 ifelse(Bin == 8,
                                                                                        PythagDist + Distance8_4 + sqrt(((BinXRef4_6 - BroodX)^2) + ((BinYRef4_6 - BroodY)^2)),
                                                                                        NA)))))))),
                                ifelse(BroodBin == 4,
                                       ifelse(Bin == 1,
                                              ifelse(
                                                BroodX < BinXRef4_6,
                                                ifelse(ScaledY > BinYRef2, 
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)) + Distance2_1),
                                                ifelse(ScaledY > BinYRef2,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                       sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1)),
                                              ifelse(Bin == 2,
                                                     ifelse(
                                                       BroodX < BinXRef4_6,
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef3)^2) + ((BroodY - BinYRef3)^2)),
                                                       sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4),
                                                     ifelse(Bin == 3,
                                                            ifelse(BroodX < BinXRef4_6,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2))),
                                                            ifelse(Bin == 4,
                                                                   sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                   ifelse(Bin == 5,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 6,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 7,
                                                                                        ifelse(ScaledY > BinYRef7,
                                                                                               sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                               PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                        ifelse(Bin == 8,
                                                                                               PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                               NA)))))))),
                                       ifelse(BroodBin == 5,
                                              ifelse(Bin == 1,
                                                     ifelse(ScaledY > BinYRef2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                     ifelse(Bin == 2,
                                                            sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                            ifelse(Bin == 3,
                                                                   sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                   ifelse(Bin == 4,
                                                                          sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                          ifelse(Bin == 5,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 6,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        ifelse(Bin == 7,
                                                                                               ifelse(ScaledY > BinYRef7,
                                                                                                      sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                      PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                               ifelse(Bin == 8,
                                                                                                      PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                      NA)))))))),
                                              ifelse(BroodBin == 6,
                                                     ifelse(Bin == 1,
                                                            ifelse(ScaledY > BinYRef2,
                                                                   sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                   sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_1),
                                                            ifelse(Bin == 2,
                                                                   sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)) + Distance3_4,
                                                                   ifelse(Bin == 3,
                                                                          sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef4_6)^2) + ((BroodY - BinYRef4_6)^2)),
                                                                          ifelse(Bin == 4,
                                                                                 sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                 ifelse(Bin == 5,
                                                                                        sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                        ifelse(Bin == 6,
                                                                                               sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                               ifelse(Bin == 7,
                                                                                                      ifelse(ScaledY > BinYRef7,
                                                                                                             sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                             PythagDist + sqrt(((BroodX - BinXRef7)^2) + ((BroodY - BinYRef7)^2))),
                                                                                                      ifelse(Bin == 8,
                                                                                                             PythagDist + Distance8_7 + sqrt(((BinXRef7 - BroodX)^2) + ((BinYRef7 - BroodY)^2)),
                                                                                                             NA)))))))),
                                                     ifelse(BroodBin == 8,
                                                            ifelse(Bin == 1,
                                                                   ifelse(ScaledY > BinYRef2,
                                                                          sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                          sqrt(((ScaledX - BinXRef2)^2) + ((ScaledY - BinYRef2)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_1),
                                                                   ifelse(Bin == 2,
                                                                          sqrt(((ScaledX - BinXRef3)^2) + ((ScaledY - BinYRef3)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_2,
                                                                          ifelse(Bin == 3,
                                                                                 sqrt(((ScaledX - BinXRef4_6)^2) + ((ScaledY - BinYRef4_6)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_4,
                                                                                 ifelse(Bin == 4,
                                                                                        sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                        ifelse(Bin == 5,
                                                                                               sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                               ifelse(Bin == 6,
                                                                                                      sqrt(((ScaledX - BinXRef7)^2) + ((ScaledY - BinYRef7)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2)) + Distance8_7,
                                                                                                      ifelse(Bin == 7,
                                                                                                             ifelse(ScaledY < BinYRef8,
                                                                                                                    sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
                                                                                                                    sqrt(((ScaledX - BinXRef8)^2) + ((ScaledY - BinYRef8)^2)) + sqrt(((BroodX - BinXRef8)^2) + ((BroodY - BinYRef8)^2))),
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
  group_by(Colony) %>%
  mutate(ToBrood = BroodDist / MaxDist,
         ToBrood = ifelse(ToBrood > 1, 1, ToBrood)) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, ScaledDist, Density, Corner, Bin, ToBrood, Sex, Ratio) %>%
  distinct() %>%
  drop_na()

MeanBroodCoordFullCircle <- MeanBroodCoordFull %>% 
  filter(Nest == "Circle") %>%
  left_join(DistBinsFull)

DistanceToBroodWorkersCircle <- WorkerDistScaledRD1_RD2 %>% 
  filter(Nest == "Circle") %>%
  left_join(MeanBroodCoordFullCircle) %>%
  group_by(Colony, Day) %>% 
  mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
         ToBrood = BroodDist / MaxDist) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, Corner, ColonyMember, ToBrood) %>%
  distinct() %>%
  drop_na()

DistanceToBroodQueenCircle <- QueenDistScaledRD1_RD2 %>% 
  filter(Nest == "Circle") %>%
  left_join(MeanBroodCoordFullCircle) %>%
  group_by(Colony, Day) %>% 
  mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
         ToBrood = BroodDist / MaxDist) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, Bin, ScaledDist, Density, Corner, ToBrood) %>%
  distinct() %>%
  drop_na()

DistanceToBroodAlateCircle <- AlateDistScaledRD2 %>% 
  filter(Nest == "Circle") %>%
  left_join(MeanBroodCoordFullCircle) %>%
  group_by(Colony, Day) %>% 
  mutate(BroodDist = sqrt(((ScaledX - BroodX)^2) + ((ScaledY - BroodY)^2)),
         ToBrood = BroodDist / MaxDist) %>%
  select(Colony, Nest, Day, ScaledX, ScaledY, ScaledDist, Density, Corner, Bin, ToBrood, Sex, Ratio) %>%
  distinct() %>%
  drop_na()

BroodCentDistWorkersRD1_RD2 <- full_join(DistanceToBroodWorkersTube, DistanceToBroodWorkersCircle)

BroodCentDistWorkersRD1 <- BroodCentDistWorkersRD1_RD2 %>%
  filter(Colony < 11)

BroodCentDistWorkersRD2 <- BroodCentDistWorkersRD1_RD2 %>%
  filter(Colony > 10)

BroodCentDistQueensRD1_RD2 <- full_join(DistanceToBroodQueenTube, DistanceToBroodQueenCircle)

BroodCentDistQueensRD1 <- BroodCentDistQueensRD1_RD2 %>%
  filter(Colony < 11)

BroodCentDistQueensRD2 <- BroodCentDistQueensRD1_RD2 %>%
  filter(Colony > 10)

BroodCentDistAlatesRD2 <- full_join(DistanceToBroodAlateTube, DistanceToBroodAlateCircle)
