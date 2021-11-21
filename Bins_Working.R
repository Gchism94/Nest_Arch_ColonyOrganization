#############################################
## Autor: Greg CHISM
## Date: Aug 2021
## email: gchism@email.arizona.edu
## Project: Nest shape influences colony organization in ants
## Title: Nest section bin functions 
#############################################

#BIN ASSIGNMENT FUNCTION
#The code below bins x and y coordinate colony data into eight even area nest sections
#Note that the code is specifically used below for the low density treatment colonies (11-20)
#This is because the high density treatment colonies (1-10) were done through excel but can be run through this script as well
#The code for Netlogo simulated x and y coordinate results is separate and below
#To do this, a reference dataset of bin coordinates is used and coordinates are run through a series of conditional statements
#Where each conditional statement checks whether the coordinate is in one of eight bins sequentially
#Note that you have to change the value of "Colony" before each iteration.
#A for loop could have accomplished to avoid changing the value, but computational power was a concern
CoordBinned<-function(data_table){
  #Selecting the colony bin references
  BinCoordAssign <- BinCoordFull %>%
    filter(Colony == "12")
  
  #Binning the tube colony member coordinates
  #Filtering out the tube data for the selected colony
  Colony12RD2BinnedTube <- data_table %>%
    filter(Colony == "12" & Nest == "Tube") %>%
    #Ifelse conditional statements for each bin
    mutate(Bin =
             #The bins determine whether the x and y coordinates are within a set of relevant bounds for each bin
             #The order is always x coordinates first, then y
             if_else(
               ScaledX >= BinCoordAssign$ScaledX[2] & 
                 ScaledX <= BinCoordAssign$ScaledX[1] &
                 ScaledY <= BinCoordAssign$ScaledY[2], 1,
               #Bin 2
               if_else(ScaledX >= BinCoordAssign$ScaledX[3] & 
                         ScaledX <= BinCoordAssign$ScaledX[2] &
                         ScaledY >= BinCoordAssign$ScaledY[3] &
                         ScaledY <= BinCoordAssign$ScaledY[2], 2,
                       #Bin 3
                       if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                 ScaledX <= BinCoordAssign$ScaledX[3] &
                                 ScaledY >= BinCoordAssign$ScaledY[3] &
                                 ScaledY <= BinCoordAssign$ScaledY[4], 3,
                               #Bin 4
                               if_else(ScaledX >= BinCoordAssign$ScaledX[4] & 
                                         ScaledX <= BinCoordAssign$ScaledX[5] &
                                         ScaledY >= BinCoordAssign$ScaledY[4] &
                                         ScaledY <= BinCoordAssign$ScaledY[5], 4,
                                       #Bin 5
                                       if_else(ScaledX >= BinCoordAssign$ScaledX[5] & 
                                                 ScaledX <= BinCoordAssign$ScaledX[6] &
                                                 ScaledY >= BinCoordAssign$ScaledY[4] &
                                                 ScaledY <= BinCoordAssign$ScaledY[5], 5,
                                               #Bin 6
                                               if_else(ScaledX >= BinCoordAssign$ScaledX[6] & 
                                                         ScaledX <= BinCoordAssign$ScaledX[7] &
                                                         ScaledY >= BinCoordAssign$ScaledY[7] &
                                                         ScaledY <= BinCoordAssign$ScaledY[6], 6,
                                                       #Bin 7
                                                       if_else(ScaledX >= BinCoordAssign$ScaledX[8] & 
                                                                 ScaledX <= BinCoordAssign$ScaledX[7] &
                                                                 ScaledY <= BinCoordAssign$ScaledY[7], 7,
                                                               #Bin 8
                                                               if_else(ScaledX >= BinCoordAssign$ScaledX[9] & 
                                                                         ScaledX <= BinCoordAssign$ScaledX[8] &
                                                                         ScaledY <= BinCoordAssign$ScaledY[9], 8, 0
                                                                       )))))))))
  #Binning circle nest coordinates 
  #Filtering out the circle data for the selected colony
  #The code is the same as above, except its only y coordinates since xlims don't matter for the binning to work 
  Colony12RD2BinnedCircle <- data_table %>%
    filter(Colony == "12" & Nest == "Circle") %>%
    mutate(Bin =
             #Bin 1
             if_else(ScaledY <= BinCoordAssign$ScaledY[10], 1,
                     #Bin 2
                     if_else(ScaledY >= BinCoordAssign$ScaledY[10] & 
                               ScaledY <= BinCoordAssign$ScaledY[11], 2,
                             #Bin 3
                             if_else(ScaledY >= BinCoordAssign$ScaledY[11] & 
                                       ScaledY <= BinCoordAssign$ScaledY[12], 3,
                                     #Bin 4
                                     if_else(ScaledY >= BinCoordAssign$ScaledY[12] & 
                                               ScaledY <= BinCoordAssign$ScaledY[13], 4,
                                             #Bin 5
                                             if_else(ScaledY >= BinCoordAssign$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssign$ScaledY[14], 5,
                                                     #Bin 6
                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssign$ScaledY[15], 6,
                                                             #Bin 7
                                                             if_else(ScaledY >= BinCoordAssign$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssign$ScaledY[16], 7,
                                                                     #Bin 8
                                                                     if_else(ScaledY >= BinCoordAssign$ScaledY[16], 8, 0
                                                                             )))))))))
  Colony12Binned <<- full_join(Colony12RD2BinnedTube, Colony12RD2BinnedCircle)
}

CoordBinned(WorkerTestFullWorking) #Colony member dataset

#The following can be used to check how many coordinates have zero bin values, meaning they were not placed in a bin
#These coordinates need to be removed as they are most likely an error (a coordinate where no ant exists)
#This also checks whether the correct colony assignment was used throughout
#Each issue was examined in Excel 


#e.g.Colony1
Colony1RD2Binned %>%
  filter(Bin == "0") %>%
  group_by(Bin) %>%
  summarise(n = n())

#Combine all datasets into the final one for workers, colonies 1-10
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
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
  #Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(Bin != 0) 

#Combine all datasets into the final one for workers, colonies 11-20
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordWorkersRD2 <- full_join(Colony11Binned) %>%
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
  filter(Bin != 0) 

#Combine all datasets into the final one for brood, colonies 1-10
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
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
  #Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(Bin != 0) 

#Combine all datasets into the final one for brood, colonies 11-20
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordBroodRD2 <- full_join(Colony11Binned) %>%
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
  filter(Bin != 0) 
 
#Combine all datasets into the final one for queens, colonies 1-10
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
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
  filter(Bin != 0) 

#Combine all datasets into the final one for queens, colonies 11-20
#This is because functions in other scripts keep colonies 1-10 (high nest density) and 11-20 (low nest density) separate
FullDataCoordQueenRD2 <- full_join(Colony11Binned) %>%
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
  filter(Bin != 0) 

#Combine all datasets into the final one for alates
#Note that alates are only found in the low density treatment so there is only one script to combine alates
FullDataCoordAlates <- full_join(Colony11Binned) %>%
  full_join(Colony12Binned) %>%
  full_join(Colony13Binned) %>%
  full_join(Colony14Binned) %>%
  full_join(Colony15Binned) %>%
  full_join(Colony18Binned) %>%
  full_join(Colony19Binned) %>%
  #Remove all of the zeros, since some coordinates exist that may be accidental or error. 
  filter(Bin != 0) 

#NETLOGO SIMULATION RESULTS BIN FUNCTION
#The below code is the same as for the experimental coordinates
#There is also a bin coordinate reference dataset used here 
#Except for a renaming line that changes "x" and "y" columns to "ScaledX" and "ScaledY" in order to work with the function
#This also makes later joining real data and simulated results easier 
#There is a separate set of code for each combination of nest shape (Tube, Circle) and size (Small, Large)

BinCoordNetlogo <- read_excel("~/Desktop/NetLogoGC/ArchitectureMoveModelFull_GC_11.Aug.2021_TestExp_Table.xlsx", 
                              sheet = "Bins") %>%
  mutate(ScaledX = ScaledX * 0.1,
         ScaledY = ScaledY * 0.1)
write.csv(BinCoordNetlogo, "BinCoordNetlogo.csv", row.names = F)
CoordBinnedNetlogo <- function(data_table){
  #Separating out the 
  BinCoordAssignSmall <- BinCoordNetlogo %>%
    filter(NestSize == "Small")
  
  BinCoordAssignLarge <- BinCoordNetlogo %>%
    filter(NestSize == "Large")
  
  #Binning small tube nest coordinates 
  NetlogoBinnedTubeSmall <- data_table %>%
    filter(NestSize == "Small" & Nest == "Tube") %>%
    mutate(Bin =
             #Bin 1
             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[2] & 
                       ScaledX <= BinCoordAssignSmall$ScaledX[1] &
                       ScaledY <= BinCoordAssignSmall$ScaledY[2], 1,
                     #Bin 2
                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[3] & 
                               ScaledX <= BinCoordAssignSmall$ScaledX[2] &
                               ScaledY >= BinCoordAssignSmall$ScaledY[3] &
                               ScaledY <= BinCoordAssignSmall$ScaledY[2], 2,
                             #Bin 3
                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[4] & 
                                       ScaledX <= BinCoordAssignSmall$ScaledX[3] &
                                       ScaledY >= BinCoordAssignSmall$ScaledY[3] &
                                       ScaledY <= BinCoordAssignSmall$ScaledY[4], 3,
                                     #Bin 4
                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[4] & 
                                               ScaledX <= BinCoordAssignSmall$ScaledX[5] &
                                               ScaledY >= BinCoordAssignSmall$ScaledY[4] &
                                               ScaledY <= BinCoordAssignSmall$ScaledY[5], 4,
                                             #Bin 5
                                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[5] & 
                                                       ScaledX <= BinCoordAssignSmall$ScaledX[6] &
                                                       ScaledY >= BinCoordAssignSmall$ScaledY[4] &
                                                       ScaledY <= BinCoordAssignSmall$ScaledY[5], 5,
                                                     #Bin 6
                                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[6] & 
                                                               ScaledX <= BinCoordAssignSmall$ScaledX[7] &
                                                               ScaledY >= BinCoordAssignSmall$ScaledY[7] &
                                                               ScaledY <= BinCoordAssignSmall$ScaledY[6], 6,
                                                             #Bin 7
                                                             if_else(ScaledX >= BinCoordAssignSmall$ScaledX[8] & 
                                                                       ScaledX <= BinCoordAssignSmall$ScaledX[7] &
                                                                       ScaledY <= BinCoordAssignSmall$ScaledY[7], 7,
                                                                     #Bin 8
                                                                     if_else(ScaledX >= BinCoordAssignSmall$ScaledX[9] & 
                                                                               ScaledX <= BinCoordAssignSmall$ScaledX[8] &
                                                                               ScaledY <= BinCoordAssignSmall$ScaledY[9], 8, 0
                                                                             )))))))))
  #Binning large tube nest coordinates 
  NetlogoBinnedTubeLarge <- data_table %>%
    filter(NestSize == "Large" & Nest == "Tube") %>%
    mutate(Bin =
             #Bin 1
             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[2] & 
                       ScaledX <= BinCoordAssignLarge$ScaledX[1] &
                       ScaledY <= BinCoordAssignLarge$ScaledY[2], 1,
                     #Bin 2
                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[3] & 
                               ScaledX <= BinCoordAssignLarge$ScaledX[2] &
                               ScaledY >= BinCoordAssignLarge$ScaledY[3] &
                               ScaledY <= BinCoordAssignLarge$ScaledY[2], 2,
                             #Bin 3
                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[4] & 
                                       ScaledX <= BinCoordAssignLarge$ScaledX[3] &
                                       ScaledY >= BinCoordAssignLarge$ScaledY[3] &
                                       ScaledY <= BinCoordAssignLarge$ScaledY[4], 3,
                                     #Bin 4
                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[4] & 
                                               ScaledX <= BinCoordAssignLarge$ScaledX[5] &
                                               ScaledY >= BinCoordAssignLarge$ScaledY[4] &
                                               ScaledY <= BinCoordAssignLarge$ScaledY[5], 4,
                                             #Bin 5
                                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[5] & 
                                                       ScaledX <= BinCoordAssignLarge$ScaledX[6] &
                                                       ScaledY >= BinCoordAssignLarge$ScaledY[4] &
                                                       ScaledY <= BinCoordAssignLarge$ScaledY[5], 5,
                                                     #Bin 6
                                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[6] & 
                                                               ScaledX <= BinCoordAssignLarge$ScaledX[7] &
                                                               ScaledY >= BinCoordAssignLarge$ScaledY[7] &
                                                               ScaledY <= BinCoordAssignLarge$ScaledY[6], 6,
                                                             #Bin 7
                                                             if_else(ScaledX >= BinCoordAssignLarge$ScaledX[8] & 
                                                                       ScaledX <= BinCoordAssignLarge$ScaledX[7] &
                                                                       ScaledY <= BinCoordAssignLarge$ScaledY[7], 7,
                                                                     #Bin 8
                                                                     if_else(ScaledX >= BinCoordAssignLarge$ScaledX[9] & 
                                                                               ScaledX <= BinCoordAssignLarge$ScaledX[8] &
                                                                               ScaledY <= BinCoordAssignLarge$ScaledY[9], 8, 0
                                                                     )))))))))
  #Binning small circle nest coordinates 
  NetlogoBinnedCircleSmall <- data_table%>%
    filter(NestSize == "Small" & Nest == "Circle") %>%
    mutate(Bin =
             #Bin 1
             if_else(ScaledY <= BinCoordAssignSmall$ScaledY[10], 1,
                     #Bin 2
                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[10] & 
                               ScaledY <= BinCoordAssignSmall$ScaledY[11], 2,
                             #Bin 3
                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[11] & 
                                       ScaledY <= BinCoordAssignSmall$ScaledY[12], 3,
                                     #Bin 4
                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[12] & 
                                               ScaledY <= BinCoordAssignSmall$ScaledY[13], 4,
                                             #Bin 5
                                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssignSmall$ScaledY[14], 5,
                                                     #Bin 6
                                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssignSmall$ScaledY[15], 6,
                                                             #Bin 7
                                                             if_else(ScaledY > BinCoordAssignSmall$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssignSmall$ScaledY[16], 7,
                                                                     #Bin 8
                                                                     if_else(ScaledY > BinCoordAssignSmall$ScaledY[16], 8, 0
                                                                             )))))))))
  #Binning large circle nest coordinates 
  NetlogoBinnedCircleLarge <- data_table %>%
    filter(NestSize == "Large" & Nest == "Circle") %>%
    mutate(Bin =
             #Bin 1
             if_else(ScaledY <= BinCoordAssignLarge$ScaledY[10], 1,
                     #Bin 2
                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[10] & 
                               ScaledY <= BinCoordAssignLarge$ScaledY[11], 2,
                             #Bin 3
                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[11] & 
                                       ScaledY <= BinCoordAssignLarge$ScaledY[12], 3,
                                     #Bin 4
                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[12] & 
                                               ScaledY <= BinCoordAssignLarge$ScaledY[13], 4,
                                             #Bin 5
                                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[13] & 
                                                       ScaledY <= BinCoordAssignLarge$ScaledY[14], 5,
                                                     #Bin 6
                                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[14] & 
                                                               ScaledY <= BinCoordAssignLarge$ScaledY[15], 6,
                                                             #Bin 7
                                                             if_else(ScaledY >= BinCoordAssignLarge$ScaledY[15] & 
                                                                       ScaledY <= BinCoordAssignLarge$ScaledY[16], 7,
                                                                     #Bin 8
                                                                     if_else(ScaledY >= BinCoordAssignLarge$ScaledY[16], 8, 0
                                                                     )))))))))
  
  NetlogoBinnedFull <<- full_join(NetlogoBinnedTubeSmall, NetlogoBinnedTubeLarge) %>%
    full_join(NetlogoBinnedCircleSmall) %>% 
    full_join(NetlogoBinnedCircleLarge)
}

CoordBinnedNetlogo(NetlogoTestFull) #NetlogoDataset

