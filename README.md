# README for Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space

***

## Overview
Data and R script used for the manuscript: Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space

***

## Purpose of the study
### Investigating how nest shape influences how _Temnothorax rugatulus_ colonies spatially organize in their nests. This includes physical location of colony members and their distances from the entrance, mobile colony member distance to the brood center, worker distance to the physical center of the nest, and comparing worker distributions with those predicted by a random walk model. 

## Dependencies 
##### Scripts for this manuscript should be executed in the following order: 
1. Stat_boxplot_custom.R: Custom boxplot function, replacing geom_boxplot() in ggplot2. The function extends the whisker range to the data range
2. Bins_Working.r: Functions to bin colony member and Netlogo simulation coordinates into nest sections
3. DistanceFunctions.R: Functions to find each colony member and Netlogo simulation coordinates shortest distances in the nest
4. SFZFunctions.R: Functions towards marker worker site fidelity (spatial fidelity & occurrence zones)
5. NestArchFuncts.Analyses.R: All analyses and plots
##### Several packages are required, however all are loaded through the package "pacman", so be certain to install this package before running any other code.
##### See the following documentation for further information on the "pacman" package: https://www.rdocumentation.org/packages/pacman/versions/0.5.1 

***

## Structure of the data
### EMPIRICAL DATA
#### WORKERS: FullDataCoordWorkers.csv, FullDataCoordWorkersRD2.csv
###### Raw experimental data with worker x and y position in nests
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Day: The experimental day that the observation was collected on
* ScaledX: X-axis coordinate, scaled from original (px) to (cm) in the software Fiji (Schindelin et al., 2012)
* ScaledY: Y-axis coordinate, scaled from original (px) to (cm) in the software Fiji 
* ColorID: The unique color marking assigned to an individual worker's head, thorax, abdomen1, abdomen2 (i.e., Yellow, White, Green, Green = Y,W,G,G) 
* Density: The density treatment (High / Low) 

***

#### BROOD / QUEENS: FullDataCoordBrood.csv, FullDataCoordBroodRD2.csv; FullDataCoordQueen.csv, FullDataCoordQueenRD2.csv
###### Raw experimental data with brood (OR) queen x and y position in nests
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Day: The experimental day that the observation was collected on
* ScaledX: X-axis coordinate, scaled from original (px) to (cm) in the software Fiji (Schindelin et al., 2012)
* ScaledY: Y-axis coordinate, scaled from original (px) to (cm) in the software Fiji 
* Density: The density treatment (High / Low) 

***

#### ALATES: FullDataCoordAlate.csv
###### Raw experimental data with alate (winged reproductive individuals) x and y position in nests
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Day: The experimental day that the observation was collected on
* ScaledX: X-axis coordinate, scaled from original (px) to (cm) in the software Fiji (Schindelin et al., 2012)
* ScaledY: Y-axis coordinate, scaled from original (px) to (cm) in the software Fiji 
* SexID: The unique sex assignment and number given to an individual alate: Sex, SexNumber, TotalNumber (i.e., the first male alate observation that came after three queen alates making it the fourth total observation = M,1,4) 

***

#### NETLOGO SIMULATIONS: ArchitectureMoveModelFull.csv
###### Raw netlogo simulation data with agent x and y positions in nests
* RunNumber: The simulation number - 1 to 4000 - there are 1000 simulations for each combination of nest shape and size
* NestSize: The size of the nest area that agents were allowed to move throughout (Small / Large)
* Nest: The nest shape treatment (Tube / Circle)
* TimeStep: The duration of each simulation (should be 50000)
* xcor: a list of every agent x coordinate position at the end of the simulation
* ycor: a list of every agent y coordinate position at the end of the simulation

*** 

### REFERENCE DATA (Bins_Working.R)
#### NEST BINS: Empirical
##### BinsNullFull.csv
###### Null data sheet with eight bins for tube and circle nests in every colony
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Bin: Nest section identifier (1-8)

***

##### BinCoordFull.csv
###### Reference binning coordinates to group empirical coordinates into nest sections
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* CoordID: The unique coordinate identifier within each colony and nest combination
* ScaledX: X-axis coordinate, scaled from original (px) to (cm) in the software Fiji (Schindelin et al., 2012)
* ScaledY: Y-axis coordinate, scaled from original (px) to (cm) in the software Fiji 

***

#### NEST BINS: Netlogo Simulations
##### BinsNullNetlogo.csv
###### Null data sheet with eight bins for tube and circle nests in each simulation treatment
* Nest: The nest shape treatment (Tube / Circle)
* NestSize: The size treatment for simulations (Small / Large)
* Bin: Nest section identifier (1-8)

***

##### BinCoordNetlogo.csv 
###### Reference binning coordinates to group Netlogo simulation coordinates into nest sections
* Nest: The nest shape treatment (Tube / Circle)
* NestSize: The size treatment for simulations (Small / Large)
* ScaledX: X-axis coordinate
* ScaledY: Y-axis coordinate
* CoordID: The unique coordinate identifier within each colony and nest combination

***

#### CORNERS: Empirical
##### CornerFull.csv
###### Whether a nest section has a corner or not
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Bin: Nest section identifier (1-8)
* Corner: Presence of a corner (Y / N)

***

#### CORNERS: Empirical
##### CornerFullSim.csv
* Nest: The nest shape treatment (Tube / Circle)
* Bin: Nest section identifier (1-8)
* Corner: Presence of a corner (Y / N)

***

### REFERENCE DATA (DistanceFunctions.R) 
#### DISTANCES IN THE NEST: Empirical
##### DistBinsFull.csv
###### Reference coordinates for the entrance of nest sections (Bin) front-to-back and shortest distance to the entrance from each nest section entrance
* Colony: Unique experimental colony identifiers
* Nest: The nest shape treatment (Tube / Circle)
* Distance: Reference shortest distance from a nest section to the entrance
* Bin: Nest section identifier (1-8)
* BinX: X-axis reference coodinate for a nest section entrance
* BinY: Y-axis reference coodinate for a nest section entrance
* Xmax: Max X-axis coordinate possible within the nest
* Ymax: Max Y-axis coordinate possible within the nest
* MaxDist: Max possible shortest distance from the nest entrance
* TubeRatio: Ratio of shortest distance to the nest entrance in circle nest / tube nest 

***

#### DISTANCES IN THE NEST: Netlogo Simulations
##### DistBinsFullNetlogo.csv
###### Reference coordinates for the entrance of nest sections (Bin) front-to-back and shortest distance to the entrance from each nest section entrance
* NestSize: The size treatment for simulations (Small / Large)
* Nest: The nest shape treatment (Tube / Circle)
* Distance: Reference shortest distance from a nest section to the entrance
* Bin: Nest section identifier (1-8)
* BinX: X-axis reference coodinate for a nest section entrance
* BinY: Y-axis reference coodinate for a nest section entrance
* Xmax: Max X-axis coordinate possible within the nest
* Ymax: Max Y-axis coordinate possible within the nest
* MaxDist: Max possible shortest distance from the nest entrance
* TubeRatio: Ratio of shortest distance to the nest entrance in circle nest / tube nest 

*** 

### REFERENCE DATA (SFZFunctions.R)
#### WORKER SITE FIDELITY (SPATIAL FIDELITY & OCCURRENCE ZONE SIZES), ALSO RELATING SIZES TO DISTANCES IN THE NEST
##### ColorRefFull.csv
###### Reference of all possible unique color identifiers paint marked workers
* Colony: Unique experimental colony identifiers
* Head: Head color mark
* Thorax: Thorax color mark
* Abd1: Left side abdomen mark
* Abd2: Right side abdomen mark

##### NestAreaFull.csv
###### Reference for colony size (number of workers in the colony) and nest area
* Colony: Unique experimental colony identifiers
* Number.ants: Number of workers in the colony after painting
* Diameter: The diameter of the circle nest 
* Area: The area of the nest 

***

##### ScalingCircleSFZ.csv
###### Reference to scale the radius of circle nests to make coordinates representing fidelity zone bins
* Colony: Unique experimental colony identifiers
* Scaling: The scaling factor that is applied to the radius of each circle nest

***
