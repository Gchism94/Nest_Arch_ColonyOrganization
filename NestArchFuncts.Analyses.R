####################################################################################################################
## Author: GREG CHISM
## Date: APR 2022
## email: gchism@arizona.edu
## Project: Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space
## Title: Calculating colony member proportions in nest sections, importing and preparing Netlogo simulation data, all paper analyses and figures 
####################################################################################################################

# This code is to replicate all analyses and figures for the "Nest shape influences colony organization in ants: spatial distribution and connectedness of colony members differs from that predicted by random movement and is affected by available space"

# READ ME: 
# YOU MUST RUN THE FOLLOWING R SCRIPTS IN ORDER
# (1) Bins_Working.R - Binning empirical and Netlogo simulation data and then calculating the proportions of each in the binned sections
# (2) DistanceFunctions.R - Distance to the nest entrance for all data, mobile colony member distance to the brood center, worker distance to all nest sections but their own
# (3) SFZFunctions.R - Paint-marked worker site fidelity zones

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
  labs(color = "KNN") +
  facet_wrap(~ Nest) +
  guides(color = guide_colorsteps(barheight = unit(0.5, "cm"),
                                  barwidth = unit(7.5, "cm"),
                                  even.steps = FALSE,
                                  frame.colour = "black")) +
  xlim(0, 7.5) +
  ylim(-0.25, 5) 

# BROOD 
BroodDensityColony <- ggplot(data = BroodDistScaledRD1_RD2 %>% filter(Colony == 11),
                                  aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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

AlateQueenDensity <- ggplot(data = AlateDistScaledRD2 %>% filter(Colony == 11),
                       aes(ScaledX, ScaledY)) +
  geom_pointdensity(alpha = 0.75) +
  geom_point(data = QueenDistScaledRD1_RD2 %>% filter(Colony == 11),
             aes(ScaledX, ScaledY), alpha = 0.75, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Alates & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig2.pdf", plot = DensityPlotColony, width = 9.375, height = 6.77, units = "in")

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
WorkerSim1TubePropAnalysis <- WorkerSim1PropAnalysisSmall %>% filter(Nest == "Tube")

# Circle nest shape
WorkerSim1CirclePropAnalysis <- WorkerSim1PropAnalysisSmall %>% filter(Nest == "Circle")

# Low nest density
WorkerSim2PropAnalysisLarge <- AntPropFullSim %>% filter(NestSize == "Large")

# Tube nest shape
WorkerSim2TubePropAnalysisRD2 <- WorkerSim2PropAnalysisLarge %>% filter(Nest == "Tube")

# Circle nest shape
WorkerSim2CirclePropAnalysisRD2 <- WorkerSim2PropAnalysisLarge %>% filter(Nest == "Circle")

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

# SIMULATION & EXPERIMENTAL DISTRIBUTION COMPARISIONS - CRAMER VON MISES TESTS WITH BENJAMINI-HOCHBERG FDR POST-HOC P-VALUE CORRECTIONS
# High density tube nest comparison
# Test statistic
cvm_EmpTube1.SimTube1.1 <- as.data.frame(cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)[1]) %>%
  mutate(Test = "EmpTube1.SimTube1") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvm_EmpTube1.SimTube1.2 <- as.data.frame(cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)[2]) %>%
  mutate(Test = "EmpTube1.SimTube1") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvm_EmpTube1.SimTube1 <- full_join(cvm_EmpTube1.SimTube1.1, cvm_EmpTube1.SimTube1.2)

# Low density tube nest comparison
# Test stat
cvm_EmpTube2.SimTube2.1 <- as.data.frame(cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "EmpTube2.SimTube2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvm_EmpTube2.SimTube2.2 <- as.data.frame(cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "EmpTube2.SimTube2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvm_EmpTube2.SimTube2 <- full_join(cvm_EmpTube2.SimTube2.1, cvm_EmpTube2.SimTube2.2)

# High density circle nest comparison
# Test statistic
cvm_EmpCircle1.SimCircle1.1 <- as.data.frame(cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)[1]) %>%
  mutate(Test = "EmpCircle1.SimCircle1") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvm_EmpCircle1.SimCircle1.2 <- as.data.frame(cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)[2]) %>%
  mutate(Test = "EmpCircle1.SimCircle1") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvm_EmpCircle1.SimCircle1 <- full_join(cvm_EmpCircle1.SimCircle1.1, cvm_EmpCircle1.SimCircle1.2)

# Low density circle nest comparison
# Test statistic
cvm_EmpCircle2.SimCircle2.1 <- as.data.frame(cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "EmpCircle2.SimCircle2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvm_EmpCircle2.SimCircle2.2 <- as.data.frame(cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "EmpCircle2.SimCircle2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvm_EmpCircle2.SimCircle2 <- full_join(cvm_EmpCircle2.SimCircle2.1, cvm_EmpCircle2.SimCircle2.2)

# High density tube nest empirical workers v. low density tube nest simulations
# Test statistic
cvmEmpTube1.SimTube2.1 <- as.data.frame(cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "EmpTube1.SimTube2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmEmpTube1.SimTube2.2 <- as.data.frame(cvm_test(WorkerTubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "EmpTube1.SimTube2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmEmpTube1.SimTube2 <- full_join(cvmEmpTube1.SimTube2.1, cvmEmpTube1.SimTube2.2)

# Low density tube nest empirical workers v. high density tube nest simulations
# Test statistic
cvmEmpTube2.SimTube1.1 <- as.data.frame(cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)[1]) %>%
  mutate(Test = "EmpTube2.SimTube1") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmEmpTube2.SimTube1.2 <- as.data.frame(cvm_test(WorkerTubePropAnalysisRD2$PropWorker, WorkerSim1TubePropAnalysis$PropWorker)[2]) %>%
  mutate(Test = "EmpTube2.SimTube1") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmEmpTube2.SimTube1 <- full_join(cvmEmpTube2.SimTube1.1, cvmEmpTube2.SimTube1.2)

# High density circle nest empirical workers v. low density circle nest simulations
# Test statistic
cvmEmpCircle1.SimCircle2.1 <- as.data.frame(cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "EmpCircle1.SimCircle2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmEmpCircle1.SimCircle2.2 <- as.data.frame(cvm_test(WorkerCirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "EmpCircle1.SimCircle2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmEmpCircle1.SimCircle2 <- full_join(cvmEmpCircle1.SimCircle2.1, cvmEmpCircle1.SimCircle2.2)

# Low density circle nest empirical workers v. high density circle nest simulations
# Test statistic
cvmEmpCircle2.SimCircle1.1 <- as.data.frame(cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)[1]) %>%
  mutate(Test = "EmpCircle2.SimCircle1") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmEmpCircle2.SimCircle1.2 <- as.data.frame(cvm_test(WorkerCirclePropAnalysisRD2$PropWorker, WorkerSim1CirclePropAnalysis$PropWorker)[2]) %>%
  mutate(Test = "EmpCircle2.SimCircle1") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmEmpCircle2.SimCircle1 <- full_join(cvmEmpCircle2.SimCircle1.1, cvmEmpCircle2.SimCircle1.2)

# COMPARING SIMULATIONS
# High density vs low density tube nest
# Note that because this statistic is calculated through bootstrap resamplings there is a slightly different p-value each time, but it is always > 0.95
# Test statistic
cvmSimTube1.SimTube2.1 <- as.data.frame(cvm_test(WorkerSim1TubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "SimTube1.SimTube2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmSimTube1.SimTube2.2 <- as.data.frame(cvm_test(WorkerSim1TubePropAnalysis$PropWorker, WorkerSim2TubePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "SimTube1.SimTube2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmSimTube1.SimTube2 <- full_join(cvmSimTube1.SimTube2.1, cvmSimTube1.SimTube2.2)

# High density vs low density circle nest
# Test Statistic
cvmSimCircle1.SimCircle2.1 <- as.data.frame(cvm_test(WorkerSim1CirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[1]) %>%
  mutate(Test = "SimCircle1.SimCircle2") %>%
  rename("Test_stat" = 1) %>%
  remove_rownames()

# P-value
cvmSimCircle1.SimCircle2.2 <- as.data.frame(cvm_test(WorkerSim1CirclePropAnalysis$PropWorker, WorkerSim2CirclePropAnalysisRD2$PropWorker)[2]) %>%
  mutate(Test = "SimCircle1.SimCircle2") %>%
  rename("P_value" = 1) %>%
  remove_rownames()

# Final data frame with test statistic and p-value
cvmSimCircle1.SimCircle2 <- full_join(cvmSimCircle1.SimCircle2.1, cvmSimCircle1.SimCircle2.2)

# Full set of CVM p-values
FullcvmDists <- full_join(cvm_EmpTube1.SimTube1, cvm_EmpTube2.SimTube2) %>%
  full_join(cvm_EmpCircle1.SimCircle1) %>% 
  full_join(cvm_EmpCircle2.SimCircle2) %>% 
  full_join(cvmEmpTube1.SimTube2) %>% 
  full_join(cvmEmpTube2.SimTube1) %>% 
  full_join(cvmEmpCircle1.SimCircle2) %>% 
  full_join(cvmEmpCircle2.SimCircle1) %>% 
  full_join(cvmSimTube1.SimTube2) %>% 
  full_join(cvmSimCircle1.SimCircle2)

# Benjamini-Hochberg method for correcting False Discovery Rates in multiple comparison testing p-values
BH.P.adjust <- as.data.frame(p.adjust(FullcvmDists$P_value, method = "BH")) %>%
  rename(BH.P.adjust = 1)

# Final adjusted p-values for multiple cvm comparisons 
FullcvmDists.adj <- cbind(FullcvmDists, BH.P.adjust)

# Creates a clean HTML file of the data frame
formattable(FullcvmDists.adj)

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
AntPropFullWorkersSimRD1_RD2 <- full_join(WorkerSimPropFullRD1, WorkerSimPropFullRD2) # Join all worker proportions data sets

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
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.87, vjust = 0.5),
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
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Worker type", color = "black")) +
  scale_fill_manual(values = c("red", "white"),
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.8)

# Compiling the empirical worker density in nest sections box plots 
WorkerPropPlot <- ggarrange(WorkerPlotRD1Tube, WorkerPlotRD2Tube,
                            labels = c("(a)", "(b)"),
                            font.label = list(size = 18, face = "plain"),
                            label.x = 0.9,
                            label.y = 1,
                            ncol = 2, nrow = 1,
                            common.legend = FALSE)

# Annotating the compiled tube nest plot to include a title
WorkerPropPlotAnnot <- annotate_figure(WorkerPropPlot,
                                      top = text_grob("Tube nest", color = "black",
                                                      size = 18, x = 0.08, y = -0.8),
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "white", hjust = 0.875, vjust = 0.5),
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
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, color = "white", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.75),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Worker type", color = "black")) +
  scale_fill_manual(values = c("blue", "white"),
                    labels = c("Obsv", "RandWalk")) +
  ylim(0, 0.5)

# Compiling the two above circle nest box plots 
WorkerPropPlot2 <- ggarrange(WorkerPlotRD1Circle, WorkerPlotRD2Circle,
                           labels = c("(c)", "(d)"),
                           font.label = list(size = 18, color = "black", face = "plain"),
                           label.x = 0.9,
                           label.y = 1,
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled circle nest plot to include a title
WorkerPropPlotAnnot2 <- annotate_figure(WorkerPropPlot2,
                                        top = text_grob("Circle nest", color = "black",
                                                        size = 18, x = 0.08, y = -0.8),
                                        bottom = NULL,
                                        left = NULL,
                                        right = NULL
)

# Compiling the empirical and Netlogo simulated worker densities in nest sections box plots 
WorkerPropPlotAnnotFull <- ggarrange(WorkerPropPlotAnnot, WorkerPropPlotAnnot2,
                             ncol = 1, nrow = 2,
                             common.legend = FALSE)

# Annotating the compiled plot to include shared axes lables
WorkerPropPlot <- annotate_figure(WorkerPropPlotAnnotFull,
                         top = NULL,
                         bottom = text_grob("Nest section", color = "black",
                                            size = 18, x = 0.5),
                         left = text_grob("Proportions of workers", color = "black",
                                          size = 18, rot = 90),
                         right = NULL)

# Save plot as a PDF
ggsave(file = "Fig3.pdf", plot = WorkerPropPlot, width = 10.4, height = 10.4, units = "in")

# LINEAR REGRESSION: observed & simulated workers
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

# SUPPLEMENTARY: NEST SECTION THAT HOLDS THE MAXIMUM EMPIRICAL AND NETLOGO SIMULATED WORKER PROPORTION
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = -0.1, vjust = 0.5),
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
        axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(size = 18, color = "black", hjust = -0.1, vjust = 0.5),
        legend.key = element_blank(),
        legend.position = c(0.975, 0.95),
        legend.direction = "horizontal",
        legend.justification = c(1, 1),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1, 'cm')) +
  xlab(NULL) +
  ylab(NULL) +
  guides(color = guide_legend(title = "Nest",
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
                        font.label = list(size = 18, face = "plain"),
                        ncol = 2, nrow = 1,
                        common.legend = FALSE)

# Annotating the compiled GLM line plots to include shared axes labels
MaxWorkerProp <- annotate_figure(MaxLogPlot,
                                  top = NULL,
                                  bottom = text_grob("Nest section", color = "black",
                                                     size = 18),
                                  left = text_grob("Max worker proportion", color = "black",
                                                   size = 18, rot = 90),
                                  right = NULL
)

# Save plot as a PDF
ggsave(file = "Fig_A11.pdf", plot = MaxWorkerProp, width = 10.4, height = 7.3, units = "in")

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
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.87, vjust = 0.5),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
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
  theme(axis.text.y = element_text(size = 18, color = "white"),
        axis.text.x = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)


# Compiling the two above brood densities in nest sections boxplots 
BroodPropPlot <- ggarrange(BroodProp1, BroodProp2,
                           labels = c("(a)", "(b)"),
                           font.label = list(size = 18, face = "plain"),
                           label.x = 0.9,
                           label.y = 1,
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled boxplots to include shared axes lables
BroodPropPlotFull <- annotate_figure(BroodPropPlot,
                                   top = text_grob("Brood", color = "black", 
                                                   size = 18, 
                                                   x = 0.055, y = -0.6),
                                   bottom = NULL,
                                   left = text_grob("Proportions of brood", color = "black",
                                                    size = 18, rot = 90),
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
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
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
        axis.text.x = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 18, color = "white"),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

# Compiling the two above queen densities in nest sections boxplots 
QueenPropPlot <- ggarrange(QueenProp1, QueenProp2,
                           labels = c("(c)", "(d)"),
                           font.label = list(size = 18, face = "plain"),
                           label.x = 0.9,
                           label.y = 1.045,
                           ncol = 2, nrow = 1)

# Annotating the compiled boxplots to include common axes labels
QueenPropPlotFull <- annotate_figure(QueenPropPlot,
                top = text_grob("Queens", color = "black", 
                                size = 18, 
                                x = 0.0625, y = 0),
                bottom = text_grob("Nest section", color = "black",
                                   size = 18, 
                                   x = 0.5),
                left = text_grob("Proportions of queens", color = "black",
                                 size = 18, rot = 90),
                right = NULL
)

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
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, color = "black"),
        plot.title = element_text(size = 18, face = "bold", color = "black"),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 1)

# Arranging the above alate densities in nest sections boxplot 
AlatePropPlot <- ggarrange(AlatePropFig,
                           labels = c("(e)"),
                           font.label = list(size = 18, face = "plain"),
                           label.x = 0.9,
                           label.y = 1.0325,
                           ncol = 2, nrow = 1)

# Annotating the alate proportions boxplot to include a title
AlateProp <- annotate_figure(AlatePropPlot,
                             top = text_grob("Alates", color = "black", 
                                             size = 18, 
                                             x = 0.0625, y = 0),
                                   bottom = NULL,
                                   left = NULL,
                                   right = NULL
)


# Compiling the brood and queen densities in nest sections plots
BroodQueenAlateProp <- ggarrange(BroodPropPlotFull, QueenPropPlotFull,
                                 AlateProp, 
                            ncol = 1, nrow = 3,
                            common.legend = TRUE)


# Save plot as a PDF
ggsave(file = "Fig4.pdf", plot = BroodQueenAlateProp, width = 10.4, height = 19.6, units = "in")

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
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
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.0, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 2500)

# Compiling the above empirical worker distances to the nest entrance histograms
WorkerDistPlot <- ggarrange(WorkerDist1, WorkerDist2,
                          labels = c("(a)", "(b)"),
                          label.x = 0.9,
                          font.label = list(size = 18, face = "plain"),
                          ncol = 2, nrow = 1,
                          common.legend = FALSE)

# Annotating the compiled plot to add a title
WorkerDistPlotFull <- annotate_figure(WorkerDistPlot,
                top = text_grob("Workers", color = "black",
                                size = 18, x = 0.055, y = -0.6),
                bottom = NULL,
                left =  text_grob("Obsv worker count", color = "black",
                                  size = 18, rot = 90),
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
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
                          font.label = list(size = 18, face = "plain"),
                          ncol = 2, nrow = 1)

# Annotating the compiled plot to add a title
WorkerSimPlotFull <- annotate_figure(WorkerSimPlot,
                top = text_grob("Random walk", color = "black",
                                size = 18, x = 0.089, y = -1),
                bottom = NULL,
                left =  text_grob("Sim worker count", color = "black",
                                  size = 18, rot = 90),
                right = NULL
)

# Compiling the worker and Netlogo simulation histograms
FullDistPlot <- ggarrange(WorkerDistPlotFull, WorkerSimPlotFull,
                        ncol = 1, nrow = 2,
                        common.legend = FALSE)

# Annotating the compiled histograms to include common axes labels 
WorkerSimDistPlot <- annotate_figure(FullDistPlot,
                top = NULL,
                bottom = text_grob("Scaled distance to nest entrance", color = "black",
                                   size = 18),
                left = NULL,
                right = NULL
)

# Save plot as a PDF
ggsave(file = "Fig_A8.pdf", plot = WorkerSimDistPlot, width = 10.4, height = 10.4, units = "in")

# Combining the empirical and Netlogo simulated worker distances from the nest entrance data sets
AllDistScaledRD1_RD2 <- SimDistScaled %>%
  mutate(Density = ifelse(NestSize == "Small", "High", "Low")) %>%
  full_join(WorkerDistScaledRD1_RD2)

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
summary(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1|Colony), data = WorkerDistScaledRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1|Colony), data = WorkerDistScaledRD1_RD2))

# LINEAR REGRESSION: Empirical and Netlogo simulated worker scaled distances from the nest entrance
# RESPONSE VARIABLE 
# ScaledDist - Worker scaled distances from the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# WorkerType - Empirical or Netlogo random walk simulated worker (Worker / RandSim)
# Density - Nest density (High / Low)
# Corner - Presence of a corner in the nest section (Y / N)
summary(lm(ScaledDist ~ Nest * WorkerType + Density + Corner, AllDistScaledRD1_RD2))

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
        axis.text = element_text(size = 18,color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
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
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.0, 0.685),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 8000) +
  xlim(0, 1)

# Compiling the brood distance to the nest entrance histograms
BroodDistPlot <- ggarrange(BroodDist1, BroodDist2,
                           labels = c("(a)", "(b)"),
                           label.x = 0.9,
                           font.label = list(size = 18, face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled histograms to include a common x axis labels and  title
BroodFullDist <- annotate_figure(BroodDistPlot,
                                 top = text_grob("Brood", color = "black",
                                                 size = 18, x = 0.055, y = -0.6),
                                 bottom = NULL,
                                 left =  text_grob("Brood count", color = "black",
                                                   size = 18, x = 0.525, rot = 90),
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
summary(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1|Colony), data = BroodDistScaledRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1|Colony), data = BroodDistScaledRD1_RD2))

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
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
                           font.label = list(size = 18, face = "plain"),
                           ncol = 2, nrow = 1)

# Annotating the compiled histograms to include a common x axis labels and  title
QueenFullDist <- annotate_figure(QueenDistPlot,
                               top = text_grob("Queens", color = "black",
                                               size = 18, x = 0.06, y = -1),
                               bottom = NULL,
                               left = text_grob("Queen count", color = "black",
                                                         size = 18, x = 0.525, rot = 90),
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
summary(lmer(ScaledDist ~ Nest * Density + Corner + Day + (1|Colony), data = QueenDistScaledRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest * Density + Day + Corner + (1|Colony), data = QueenDistScaledRD1_RD2))

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) 

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
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
                           font.label = list(size = 18, face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled plots to include a common y axis and title
AlateFullDist <- annotate_figure(AlateDistPlot,
                                      top = text_grob("Alates", color = "black",
                                                      size = 18, x = 0.06, y = -0.75),
                                      bottom = NULL,
                                      left = NULL,
                                      right = NULL
)

# Compiling the brood and queen densities in nest sections plots
BroodQueenAlateDist <- ggarrange(BroodFullDist, QueenFullDist,
                                 AlateFullDist, 
                                 ncol = 1, nrow = 3,
                                 common.legend = TRUE)

# Save plot as a PDF
ggsave(file = "Fig_A9.pdf", plot = BroodQueenAlateDist, width = 10.4, height = 19.6, units = "in")

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
summary(lmer(ScaledDist ~ Nest + Sex + Ratio + Day + Corner + (1|Colony), data = AlateDistScaledRD2Plot))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ScaledDist ~ Nest + Sex + Ratio + Day + Corner + (1|Colony), data = AlateDistScaledRD2Plot))

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) + 
  xlim(0, 0.85) + 
  ylim(0, 4000)

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
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1.0, 0.7),
        legend.direction = "horizontal",
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) + 
  xlim(0, 0.85) + 
  ylim(0, 4000)

#Compiling worker scaled distances from brood center plots
WorkerBroodDistPlot <- ggarrange(WorkerBroodDist1, WorkerBroodDist2,
                                 labels = c("(a)", "(b)"),
                                 label.x = 0.9,
                                 font.label = list(size = 18, face = "plain"),
                                 ncol = 2, nrow = 1,
                                 common.legend = FALSE)

#Annotating the compiled plot to include a common y axis and title
WorkerBroodFullDist <- annotate_figure(WorkerBroodDistPlot,
                                       top = text_grob("Workers", color = "black",
                                                       size = 18, x = 0.06, y = -0.675),
                                       bottom = NULL,
                                       left = text_grob("Worker count", color = "black",
                                                  size = 18, rot = 90),
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
summary(lmer(ToBrood ~ Nest * Density + Day + Corner + (1|Colony), data = BroodCentDistWorkersRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest * Density + Day + Corner + (1|Colony), data = BroodCentDistWorkersRD1_RD2))

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlim(0, 0.85) + 
  ylim(0, 200)

# Low density treatment
QueenBroodDist2 <- ggplot(BroodCentDistQueensRD2 %>% arrange(Nest),
                          aes(ToBrood, fill = Nest)) + 
  geom_histogram(position = "identity", 
                 alpha = 0.7, 
                 binwidth = 0.0416666) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlim(0, 0.85) + 
  ylim(0, 200)

# Compiling queens scaled distances to the brood center plots
QueenBroodDistPlot <- ggarrange(QueenBroodDist1, QueenBroodDist2,
                                labels = c("(c)", "(d)"),
                                label.x = 0.9,
                                font.label = list(size = 18, face = "plain"),
                                ncol = 2, nrow = 1,
                                common.legend = FALSE)

# Annotating the compiled plots to include a common y axis and title
QueenBroodFullDist <- annotate_figure(QueenBroodDistPlot,
                                      top = text_grob("Queens", color = "black",
                                                      size = 18, x = 0.06, y = -0.75),
                                      bottom = NULL,
                                      left = text_grob("Queen count", color = "black",
                                                  size = 18, rot = 90),
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
summary(lmer(ToBrood ~ Nest * Density + Day + Corner + (1|Colony), data = BroodCentDistQueensRD1_RD2))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest * Density + Day + Corner + (1|Colony), data = BroodCentDistQueensRD1_RD2))

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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  guides(fill = guide_legend(title = "Nest")) +
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
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.key = element_blank(),
        legend.justification = c(0.5, 1),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  xlab("Alate sex") +
  ylab("Scaled distance to brood center")

# Compiling the alate distance to the the brood center plots 
AlateBroodDistPlot <- ggarrange(AlateBroodPlot, AlateBroodPlot2,
                           labels = c("(e)", "(f)"),
                           label.x = 0.9,
                           label.y = 0.965, 
                           font.label = list(size = 18, face = "plain"),
                           ncol = 2, nrow = 1,
                           common.legend = FALSE)

# Annotating the compiled plots to include a common y axis and title
AlateBroodFullDist <- annotate_figure(AlateBroodDistPlot,
                                      top = text_grob("Alates", color = "black",
                                                      size = 18, x = 0.06, y = -0.75),
                                      bottom = NULL,
                                      left = NULL,
                                      right = NULL
)

# Compiling the brood and queen densities in nest sections plots
BroodQueenAlateBroodDist <- ggarrange(WorkerBroodFullDist, QueenBroodFullDist,
                                 AlateBroodFullDist, 
                                 ncol = 1, nrow = 3,
                                 common.legend = TRUE)

# Save plot as a PDF
ggsave(file = "Fig5.pdf", plot = BroodQueenAlateBroodDist, width = 10.4, height = 19.6, units = "in")

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
summary(lmer(ToBrood ~ Nest + Sex + Day + Corner + Ratio + (1|Colony), data = BroodCentDistAlatesRD2Plot))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(ToBrood ~ Nest + Sex + Day + Corner + Ratio + (1|Colony), data = BroodCentDistAlatesRD2Plot))

####################################################################################################################
# PLOTS AND ANALYSES: Scaled spatial fidelity and occurrence zones
# The scripts below are to analyze and visualize: 
# Scaled worker spatial fidelity and occurrence zone sizes, and consider how they relate to nest shape
####################################################################################################################

# SPATIAL FIDELITY ZONE SIZE AND NEST SHAPE
# BOXPLOTS
# High density treatment
FidZone.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                      filter(Colony < 11) %>% arrange(Nest), 
                    aes(x = Nest, y = SFZ, fill = Nest),
                    position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.85,
                      width = 0.5,
                      lwd = 1.25) +  
  theme_transparent() +  
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 4) 

# Low density treatment
FidZone.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                      filter(Colony > 10) %>% arrange(Nest), 
                    aes(x = Nest, y = SFZ, fill = Nest),
                    position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.85,
                      width = 0.5,
                      lwd = 1.25) +  
  theme_transparent() +  
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  ylim(0, 4) 

# OCCURRENCE ZONE SIZE AND NEST SHAPE
# BOXPLOTS
# High density treatment
Occur.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                    filter(Colony < 11) %>% arrange(Nest), 
                  aes(x = Nest, y = Occur, fill = Nest),
                  position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.85,
                      width = 0.5,
                      lwd = 1.25) +  
  theme_transparent() +  
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) 

# Low density treatment
Occur.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                    filter(Colony > 10) %>% arrange(Nest), 
                  aes(x = Nest, y = Occur, fill = Nest),
                  position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.85,
                      width = 0.5,
                      lwd = 1.25) +  
  theme_transparent() +  
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) 


# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(SFZ ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# Occur - Worker occurrence zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(Occur ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

####################################################################################################################
# PLOTS AND ANALYSES: True spatial fidelity and occurrence zones (cm^2)
# The scripts below are to analyze and visualize: 
# True worker spatial fidelity and occurrence zone sizes (cm^2), and consider how they relate to nest shape
####################################################################################################################

# SPATIAL FIDELITY ZONE SIZE AND NEST SHAPE
# High density treatment
scaleFUNArea <- function(x) sprintf("%.1f", x)
FidZoneArea.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                          filter(Colony < 11) %>% arrange(Nest), 
                        aes(x = Nest, y = SFZ_Area, fill = Nest),
                        position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65,
                      width = 0.15,
                      lwd = 1.25,
                      fatten = 2) +  
  stat_histinterval(slab_alpha = 0.65,
                    slab_color = "black",
                    slab_size = 1.25,
                    point_interval = "median_qi",
                    point_size = 0,
                    interval_alpha = 0,
                    scale = 0.5,
                    justification = -0.2, # Separate the interval from the bottom of the histogram
                    outline_bars = TRUE) + # Creating bars around each histogram bin
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = -7.75),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  coord_flip() +
  scale_y_continuous(labels = scaleFUNArea, limits = c(0, 1))  

# Low density treatment
FidZoneArea.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                          filter(Colony > 10) %>% arrange(Nest), 
                        aes(x = Nest, y = SFZ_Area, fill = Nest),
                        position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65,
                      width = 0.15,
                      lwd = 1.25,
                      fatten = 2) +  
  stat_histinterval(slab_alpha = 0.65,
                    slab_color = "black",
                    slab_size = 1.25,
                    point_interval = "median_qi",
                    point_size = 0,
                    interval_alpha = 0,
                    scale = 0.5,
                    justification = -0.2, # Separate the interval from the bottom of the histogram
                    outline_bars = TRUE) + # Creating bars around each histogram bin
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, color = "black", hjust = 0.875, vjust = -7.75),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = "none") +
  guides(fill = guide_legend(title = "Nest", color = "black")) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  coord_flip() +
  scale_y_continuous(labels = scaleFUNArea, limits = c(0, 1.5))  

#Compiling worker fidelity zone size plots
FidZoneAreaPlot <- ggarrange(FidZoneArea.1, FidZoneArea.2,
                             labels = c("(a)", "(b)"),
                             label.x = 0.9,
                             label.y = 0.9,
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 1,
                             common.legend = FALSE)

# Annotating the compiled plots to include a common y-axis
FidZoneAreaPlotFull <- annotate_figure(FidZoneAreaPlot,
                                       top = NULL,
                                       bottom = text_grob(expression(paste('Fidelity zone size ('*cm^2*')')), color = "black",
                                                          size = 18),
                                       left = NULL,
                                       right = NULL
)

# OCCURRENCE ZONE SIZE AND NEST SHAPE
# BOXPLOTS
# High density treatment
OccurArea.1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>%
                        filter(Colony < 11) %>% arrange(Nest), 
                      aes(x = Nest, y = Occur_Area, fill = Nest),
                      position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65,
                      width = 0.15,
                      lwd = 1.25,
                      fatten = 2) +  
  stat_histinterval(slab_alpha = 0.65,
                    slab_color = "black",
                    slab_size = 1.25,
                    point_interval = "median_qi",
                    point_size = 0,
                    interval_alpha = 0,
                    scale = 0.5,
                    justification = -0.2, # Separate the interval from the bottom of the histogram
                    outline_bars = TRUE) + # Creating bars around each histogram bin
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none",
        legend.key = element_blank()) +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  coord_flip() +
  scale_y_continuous(labels = scaleFUNArea, limits = c(0, 2.5))  

# Low density treatment
OccurArea.2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZFull %>%
                        filter(Colony > 10) %>% arrange(Nest), 
                      aes(x = Nest, y = Occur_Area, fill = Nest),
                      position = position_dodge(2)) + 
  stat_boxplot_custom(qs = c(0, 0.25, 0.5, 0.75, 1.00),
                      aes(fill = Nest), 
                      color = "grey25", 
                      alpha = 0.65,
                      width = 0.15,
                      lwd = 1.25,
                      fatten = 2) +  
  stat_histinterval(slab_alpha = 0.65,
                    slab_color = "black",
                    slab_size = 1.25,
                    point_interval = "median_qi",
                    point_size = 0,
                    interval_alpha = 0,
                    scale = 0.5,
                    justification = -0.2, # Separate the interval from the bottom of the histogram
                    outline_bars = TRUE) + # Creating bars around each histogram bin
  xlab(NULL) + 
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "white"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Tube", "Circle"), 
                    name = "Nest",
                    values = c("red", "blue")) +
  coord_flip() +
  scale_y_continuous(labels = scaleFUNArea, limits = c(0, 4))  

# Compiling worker occurrence zone size plots
OccurZoneAreaPlot <- ggarrange(OccurArea.1, OccurArea.2,
                               labels = c("(c)", "(d)"),
                               label.x = 0.9,
                               label.y = 0.8,
                               font.label = list(size = 18, face = "plain"),
                               ncol = 2, nrow = 1,
                               common.legend = FALSE)

# Annotating the compiled plots to include a common y-axis
OccurZoneAreaPlotFull <- annotate_figure(OccurZoneAreaPlot,
                                         top = NULL,
                                         bottom = text_grob(expression(paste('Occurrence zone size ('*cm^2*')')), color = "black",
                                                            size = 18),
                                         left = NULL,
                                         right = NULL
)

# Compile the spatial fidelity and occurrence zone plots
FidOccurZoneAreaPlot <- ggarrange(FidZoneAreaPlotFull, OccurZoneAreaPlotFull,
                                  ncol = 1, nrow = 2,
                                  common.legend = TRUE)

# Save plot as a PDF
ggsave(file = "Fig6.pdf", plot = FidOccurZoneAreaPlot, width = 10.4, height = 10.4, units = "in")

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ_Area - True worker fidelity zone size (cm^2), zones have at least 7 observations and be at least 15% of the total observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(SFZ_Area ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ_Area ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# Occur_Area - True worker occurrence zone size (cm^2), zones have at least 7 observations
# FIXED EFFECTS 
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(Occur_Area ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur_Area ~ Nest * Density + (1|Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes (cm^2) & colony size
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes(cm^2), and consider how they relate colony size
####################################################################################################################

# Spatial fidelity zone size and number of observations for an individual
SFZCol1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "High") %>% arrange(Nest), 
                  aes(x = Number.ants, y = SFZ_Area, 
                      color = Nest,
                      linetype = Nest,
                      shape = Nest)) +
  xlab(NULL) +
  ylab(expression(paste('Fidelity zone size ('*cm^2*')'))) +
  ggtitle("High density") +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 5, height = 0.05) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 1))  

# Low density treatment
SFZCol2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "Low") %>% arrange(Nest),
                  aes(x = Number.ants, y = SFZ_Area, 
                      color = Nest,
                      linetype = Nest,
                      shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 5, height = 0.05) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 1.75))  

# Occurrence zone size and number of observations for an individual
OccurCol1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "High") %>% arrange(Nest), 
                    aes(x = Number.ants, y = Occur_Area, 
                        color = Nest,
                        linetype = Nest,
                        shape = Nest)) +
  xlab(NULL) +
  ylab(expression(paste('Occurrence zone size ('*cm^2*')'))) +
  ggtitle("High density") +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 5, height = 0.05) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black") +
  theme_pubclean() + 
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, color = "white", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 3))  

# Low density treatment
OccurCol2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "Low") %>% arrange(Nest),
                    aes(x = Number.ants, y = Occur_Area, 
                        color = Nest,
                        linetype = Nest,
                        shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 5, height = 0.05) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "white", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 4)) 

# Compiling worker true occurrence zone size vs. colony size plots
SFZOccurZoneColPlot <- ggarrange(SFZCol1, SFZCol2, OccurCol1, OccurCol2,
                                 labels = c("(a)", "(b)", "(c)", "(d)"),
                                 label.x = 0.9,
                                 label.y = 0.995,
                                 font.label = list(size = 18, face = "plain"),
                                 ncol = 2, nrow = 2,
                                 common.legend = TRUE)

# Annotate the compiled plots to include a common x-axis
FidOccurColSizePlot <- annotate_figure(SFZOccurZoneColPlot,
                                       top = NULL,
                                       bottom = text_grob("Number of workers", color = "black",
                                                          size = 18, x = 0.525),
                                       left = NULL,
                                       right = NULL
)

# Save plot as a PDF
ggsave(file = "Fig7.pdf", plot = FidOccurColSizePlot, width = 10.4, height = 10.4, units = "in")

# LINEAR MIXED EFFECTS MODEL: Worker fidelity zone size (scaled) and the number of workers in a colony (colony size)
# RESPONSE VARIABLE
# SFZ - Worker fidelity zone size 
# FIXED EFFECTS 
# Number.ants - The number of workers in a colony
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(SFZ_Area ~ Number.ants * Nest * Density + (1|Colony), WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ_Area ~ Number.ants * Nest * Density + (1|Colony), WorkerDistScaledRD1_RD2SFZWorking))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and the number of workers in a colony (colony size)
# RESPONSE VARIABLE
# Occur - Worker occurrence zone size 
# FIXED EFFECTS 
# Number.ants - The number of workers in a colony
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(Occur_Area ~ Number.ants * Nest * Density + (1|Colony), WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur_Area ~ Number.ants * Nest * Density + (1|Colony), WorkerDistScaledRD1_RD2SFZWorking))


####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes & distances from the nest entrance
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes, and consider how they relate to both nest shape and mean distance from the nest entrance
####################################################################################################################

# FIDELITY ZONE SIZE AND DISTANCE TO THE NEST ENTRANCE
# LINE PLOTS
# High density treatment
scaleFUNDist <- function(x) sprintf("%.1f", x)
SFZDist1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "High") %>% arrange(Nest), 
                        aes(x = MeanScaledDist, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("High density") +
  geom_point(size = 6, alpha = 0.33) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black", alpha = 0) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1.1)) +  
  ylim(0, 4)

# Low density treatment
SFZDist2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "Low") %>% arrange(Nest),
                        aes(x = MeanScaledDist, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  geom_point(key_glyph = large_points, size = 6, alpha = 0.33) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black", alpha = 0) +
  theme_pubclean() +  
  theme(axis.text = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1.1)) +  
  ylim(0, 4)

SFZ_grob.1 = ggplotGrob(FidZone.1)

SFZ_grob.2 = ggplotGrob(FidZone.2)

fidymin <- min(WorkerDistScaledRD1_RD2SFZWorking$SFZ); fidymax <- max(WorkerDistScaledRD1_RD2SFZWorking$SFZ)

SFZDist.Box1 <- SFZDist1 + annotation_custom(grob = SFZ_grob.1, xmin = 0.85, xmax = 1.15,
                                                       ymin = fidymin, ymax = fidymax
)

SFZDist.Box2 <- SFZDist2 + annotation_custom(grob = SFZ_grob.2, xmin = 0.85, xmax = 1.15,
                                                       ymin = fidymin, ymax = fidymax
)

# Compiling worker fidelity zone size v. worker scaled distance to the brood center plots
SFZDistPlot <- ggarrange(SFZDist.Box1, SFZDist.Box2,
                              labels = c("(a)", "(b)"),
                              label.x = 0.9,
                              font.label = list(size = 20, face = "plain"),
                              ncol = 2, nrow = 1,
                              common.legend = TRUE)

# Annotating the compiled plots to include a common y-axis
SFZFullDistPlot <- annotate_figure(SFZDistPlot,
                                        top = NULL,
                                        bottom = NULL,
                                        left = text_grob("Scaled fidelity zone size", color = "black",
                                                         size = 20, rot = 90),
                                        right = NULL
)

# OCCURRENCE ZONE SIZE AND DISTANCE TO THE BROOD CENTER
# LINE PLOTS
# High density treatment
OccurDist1 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "High") %>% arrange(Nest),
                          aes(x = MeanScaledDist, y = Occur,
                              color = Nest,
                              linetype = Nest,
                              shape = Nest)) +
  geom_point(size = 6, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 2, color = "black") +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1.1)) +  
  scale_y_continuous(labels = scaleFUN, limits = c(1, 10))

# Low density treatment
OccurDist2 <- ggplot(data = WorkerDistScaledRD1_RD2SFZWorking %>% filter(Density == "Low") %>% arrange(Nest),
                          aes(x = MeanScaledDist, y = Occur, 
                              color = Nest, 
                              linetype = Nest,
                              shape = Nest)) +
  geom_point(size = 6, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 2, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1.1)) +  
  scale_y_continuous(labels = scaleFUN, limits = c(1, 10))

Occur_grob.1 = ggplotGrob(Occur.1)

Occur_grob.2 = ggplotGrob(Occur.1)

occurymin <- min(WorkerDistScaledRD1_RD2SFZWorking$Occur); occurymax <- max(WorkerDistScaledRD1_RD2SFZWorking$Occur)

OccurDist.Box1 <- OccurDist1 + annotation_custom(grob = Occur_grob.1, xmin = 0.9, xmax = 1.15,
                                                           ymin = occurymin, ymax = occurymax
)

OccurDist.Box2 <- OccurDist2 + annotation_custom(grob = Occur_grob.1, xmin = 0.9, xmax = 1.15,
                                                           ymin = occurymin, ymax = occurymax
)

# Compile the spatial fidelity and occurrence zone v. worker scaled distance to the brood center plots and include a common legend
OccurDistPlot <- ggarrange(OccurDist.Box1, OccurDist.Box2,
                                labels = c("(c)", "(d)"),
                                label.x = 0.9,
                                font.label = list(size = 20, face = "plain"),
                                ncol = 2, nrow = 1,
                                common.legend = FALSE)

# Annotate the compiled plots to include a common x-axis
OccurFullDistPlot <- annotate_figure(OccurDistPlot,
                                          top = NULL,
                                          bottom = text_grob("Average scaled distance to entrance", color = "black",
                                                             size = 20, x = 0.525),
                                          left = text_grob("Scaled occurrence zone size", color = "black",
                                                           size = 20, rot = 90),
                                          right = NULL
)


FidOccurDistPlot <- ggarrange(SFZFullDistPlot, OccurFullDistPlot,
                                   ncol = 1, nrow = 2,
                                   common.legend = FALSE)

# Save plot as a PDF
ggsave(file = "Fig_A10.pdf", plot = FidOccurDistPlot, width = 10.4, height = 10.4, units = "in")

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# MeanScaledDist - each worker's average distance to the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(SFZ ~ MeanScaledDist * Nest * Density + (1 | Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ MeanScaledDist * Nest * Density + (1 | Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations 
# FIXED EFFECTS 
# MeanScaledDist - each worker's average distance to the nest entrance (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(Occur ~ MeanScaledDist * Nest * Density + (1 | Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ MeanScaledDist * Nest * Density + (1 | Colony), data = WorkerDistScaledRD1_RD2SFZWorking))

####################################################################################################################
# PLOTS AND ANALYSES: Spatial fidelity and occurrence zone sizes & distances from the nest entrance and brood center
# The scripts below are to analyze and visualize: 
# Worker spatial fidelity and occurrence zone sizes, and consider how they relate to both nest shape and mean distance from the brood center
####################################################################################################################

# FIDELITY ZONE SIZE AND DISTANCE TO THE BROOD CENTER
# LINE PLOTS
# High density treatment
SFZBroodDist1 <- ggplot(data = BroodCentDistWorkersSFZ %>% filter(Density == "High") %>% arrange(Nest), 
                        aes(x = MeanToBrood, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("High density") +
  geom_point(size = 6, alpha = 0.33) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black", alpha = 0) +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 0.9)) +  
  ylim(0, 4)

# Low density treatment
SFZBroodDist2 <- ggplot(data = BroodCentDistWorkersSFZ %>% filter(Density == "Low") %>% arrange(Nest),
                        aes(x = MeanToBrood, y = SFZ, 
                            color = Nest,
                            linetype = Nest,
                            shape = Nest)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  geom_point(key_glyph = large_points, size = 6, alpha = 0.33) +
  geom_line(stat = "smooth", method = lm, se = FALSE, size = 2, color = "black", alpha = 0) +
  theme_pubclean() +  
  theme(axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, color = "black", hjust = 0.875, vjust = 0.5),
        legend.key = element_blank(),
        legend.justification = c(1, -0.7),
        legend.position = c(1, 0.71),
        legend.direction = "horizontal",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key.size = unit(1, 'cm')) +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  guides(shape = guide_legend(override.aes = list(alpha = 0.75))) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1)) +  
  ylim(0, 4)

SFZ_grob.1 = ggplotGrob(FidZone.1)

SFZ_grob.2 = ggplotGrob(FidZone.2)

fidxmin <- min(BroodCentDistWorkersSFZ$MeanToBrood); fidxmax <- max(BroodCentDistWorkersSFZ$MeanToBrood)

fidymin <- min(BroodCentDistWorkersSFZ$SFZ); fidymax <- max(BroodCentDistWorkersSFZ$SFZ)

SFZBroodDist.Box1 <- SFZBroodDist1 + annotation_custom(grob = SFZ_grob.1, xmin = 0.6, xmax = 0.95,
                                                       ymin = fidymin, ymax = fidymax
)

SFZBroodDist.Box2 <- SFZBroodDist2 + annotation_custom(grob = SFZ_grob.2, xmin = 0.65, xmax = 1,
                                                       ymin = fidymin, ymax = fidymax
)

# Compiling worker fidelity zone size v. worker scaled distance to the brood center plots
SFZBroodDistPlot <- ggarrange(SFZBroodDist.Box1, SFZBroodDist.Box2,
                              labels = c("(a)", "(b)"),
                              label.x = 0.9,
                              font.label = list(size = 20, face = "plain"),
                              ncol = 2, nrow = 1,
                              common.legend = TRUE)

# Annotating the compiled plots to include a common y-axis
SFZFullBroodDistPlot <- annotate_figure(SFZBroodDistPlot,
                                        top = NULL,
                                        bottom = NULL,
                                        left = text_grob("Scaled fidelity zone size", color = "black",
                                                         size = 20, rot = 90),
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
  geom_point(size = 6, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 2, color = "black") +
  ggtitle("High density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 0.9)) +  
  scale_y_continuous(labels = scaleFUN, limits = c(1, 10))

# Low density treatment
OccurBroodDist2 <- ggplot(data = BroodCentDistWorkersSFZ %>% filter(Density == "Low") %>% arrange(Nest),
                          aes(x = MeanToBrood, y = Occur, 
                              color = Nest, 
                              linetype = Nest,
                              shape = Nest)) +
  geom_point(size = 6, alpha = 0.33) +
  geom_smooth(method = 'lm', se = FALSE, size = 2, color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Low density") +
  theme_pubclean() +  
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 20, face = "bold", color = "white", hjust = 0.75, vjust = 0.5),
        legend.position = "none") +
  labs(color = "Nest", linetype = "Nest", shape = "Nest") +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red"),
                     labels = c("Circle", "Tube")) +
  scale_x_continuous(labels = scaleFUNDist, limits = c(0, 1)) +  
  scale_y_continuous(labels = scaleFUN, limits = c(1, 10))

Occur_grob.1 = ggplotGrob(Occur.1)

Occur_grob.2 = ggplotGrob(Occur.1)

occurxmin <- min(BroodCentDistWorkersSFZ$MeanToBrood); occurxmax <- max(BroodCentDistWorkersSFZ$MeanToBrood)

occurymin <- min(BroodCentDistWorkersSFZ$Occur); occurymax <- max(BroodCentDistWorkersSFZ$Occur)

OccurBroodDist.Box1 <- OccurBroodDist1 + annotation_custom(grob = Occur_grob.1, xmin = 0.6, xmax = 0.95,
                                                           ymin = occurymin, ymax = occurymax
)

OccurBroodDist.Box2 <- OccurBroodDist2 + annotation_custom(grob = Occur_grob.1, xmin = 0.65, xmax = 1,
                                                           ymin = occurymin, ymax = occurymax
)

# Compile the spatial fidelity and occurrence zone v. worker scaled distance to the brood center plots and include a common legend
OccurBroodDistPlot <- ggarrange(OccurBroodDist.Box1, OccurBroodDist.Box2,
                                labels = c("(c)", "(d)"),
                                label.x = 0.9,
                                font.label = list(size = 20, face = "plain"),
                                ncol = 2, nrow = 1,
                                common.legend = FALSE)

# Annotate the compiled plots to include a common x-axis
OccurFullBroodDistPlot <- annotate_figure(OccurBroodDistPlot,
                                          top = NULL,
                                          bottom = text_grob("Average scaled distance to brood center", color = "black",
                                                             size = 20, x = 0.525),
                                          left = text_grob("Scaled occurrence zone size", color = "black",
                                                           size = 20, rot = 90),
                                          right = NULL
)


FidOccurBroodDistPlot <- ggarrange(SFZFullBroodDistPlot, OccurFullBroodDistPlot,
                              ncol = 1, nrow = 2,
                              common.legend = FALSE)


# Save plot as a PDF
ggsave(file = "Fig8.pdf", plot = FidOccurBroodDistPlot, width = 10.4, height = 10.4, units = "in")

# LINEAR MIXED EFFECTS MODEL: Worker spatial fidelity zone size (scaled) and nest shape
# RESPONSE VARIABLE
# SFZ - Worker spatial fidelity zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations and at least 15% of total observations
# FIXED EFFECTS 
# BroodDist - each worker's average distance to the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(SFZ ~ MeanToBrood * Nest * Density + (1 | Colony), data = BroodCentDistWorkersSFZ))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(SFZ ~ MeanToBrood * Nest * Density + (1 | Colony), data = BroodCentDistWorkersSFZ))

# LINEAR MIXED EFFECTS MODEL: Worker occurrence zone size (scaled) and nest shape
# RESPONSE VARIABLE
# Occur - Worker occurrence zone size (0 - 1, where 1 is the entire area of the nest), zones have at least 3 observations 
# FIXED EFFECTS 
# BroodDist - each worker's average distance to the brood center (0 - 1, where 1 is the longest, shortest distance from the nest entrance)
# Nest - Nest shape (Tube / Circle)
# Density - Nest density (High / Low)
# RANDOM EFFECTS
# (1|Colony) - Colony identification 
summary(lmer(Occur ~ MeanToBrood * Nest * Density + (1 | Colony), data = BroodCentDistWorkersSFZ))

# Marginal and conditional R-squared, showing the influence of the random effect on the model
r.squaredGLMM(lmer(Occur ~ MeanToBrood * Nest * Density + (1 | Colony), data = BroodCentDistWorkersSFZ))

# ANOVA to test whether individual ants hold spatial fidelity in relation to distance from the brood center over all observations
summary(aov(ToBrood ~ Nest * Density + AntIDColNest, data = BroodCentDistWorkersSFZFull))

####################################################################################################################
# SUPPLEMENTARY PLOTS
# SUPPLEMENTAL SFZ & OCCUR VS FREQUENCY PLOTS
####################################################################################################################
# Spatial fidelity zone size and number of observations for an individual
# NOTE The figures are produced by changing the number of observations required to > 2 in the FidelityZones function in the SFZFunctions.R script
SFZFreq <- ggplot(data = FidelityZonesDataRD1_RD2Supp, aes(y = SFZ, x = Freq)) +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 0.15, height = 0.005) +
  xlab(NULL) +
  ylab("Scaled fidelity zone size") +
  theme_pubclean() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20 ,colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black")) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red")) 

# Occurrence zone size and number of observations for an individual
OccurFreq <- ggplot(data = FidelityZonesDataRD1_RD2Supp, aes(y = Occur, x = Freq)) +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 0.15, height = 0.005) +
  xlab(NULL) +
  ylab("Scaled occurrence zone size") +
  theme_pubclean() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20 ,colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black")) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red")) +
  scale_y_continuous(labels = scaleFUN, limits = c(0, 10))   

# True site fidelity (cm^2)
SFZFreq_Area <- ggplot(data = FidelityZonesDataRD1_RD2Supp, aes(y = SFZ_Area, x = Freq)) +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 0.15, height = 0.005) +
  xlab(NULL) +
  ylab(expression(paste('Fidelity zone size ('*cm^2*')'))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20 ,colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black")) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red")) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 3)) 

# Occurrence zone size and number of observations for an individual
OccurFreq_Area <- ggplot(data = FidelityZonesDataRD1_RD2Supp, aes(y = Occur_Area, x = Freq)) +
  geom_jitter(size = 6, alpha = 0.5, aes(color = Nest, shape = Nest), width = 0.15, height = 0.005) +
  xlab(NULL) +
  ylab(expression(paste('Occurrence zone size ('*cm^2*')'))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20 ,colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 20, colour = "black")) +
  scale_color_manual(breaks = c("Circle", "Tube"), 
                     name = "Nest",
                     values = c("blue", "red")) +
  scale_y_continuous(labels = scaleFUNDist, limits = c(0, 4)) 

# Compiling worker site fidelity and observations plots
Fid.OccurZoneFreqPlot <- ggarrange(SFZFreq, OccurFreq, SFZFreq_Area, OccurFreq_Area,
                                   labels = c("(a)", "(b)", "(c)", "(d)"),
                                   label.x = 0.9,
                                   label.y = 1.025,
                                   font.label = list(size = 20, face = "plain"),
                                   ncol = 2, nrow = 2,
                                   common.legend = TRUE)



# Annotate the compiled plots to include a common x-axis
FidOccurFreqFull <- annotate_figure(Fid.OccurZoneFreqPlot,
                                    top = NULL,
                                    bottom = text_grob("Number of observations", color = "black",
                                                       size = 20, x = 0.525),
                                    left = NULL,
                                    right = NULL
)

# Save plot as a PDF
ggsave(file = "Fig_A7.pdf", plot = FidOccurFreqFull, width = 10.4, height = 10.4, units = "in")

####################################################################################################################
# SUPPLEMENTARY PLOTS
# SUPPLEMENTAL COLONY DENSITY PLOTS
####################################################################################################################

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
             aes(ScaledX, ScaledY), alpha = 0.33, color = "red2", size = 2, shape = 17) +
  scale_color_viridis() +
  coord_fixed() +
  theme_pubclean() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Brood & Queens") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony1 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A12.pdf", plot = DensityPlotColony1, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony2 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A13.pdf", plot = DensityPlotColony2, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which can be used below or later
DensityPlotColony3 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A14.pdf", plot = DensityPlotColony3, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony4 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A15.pdf", plot = DensityPlotColony4, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony5 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A16.pdf", plot = DensityPlotColony5, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony6 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A17.pdf", plot = DensityPlotColony6, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony7 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A18.pdf", plot = DensityPlotColony7, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony8 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A19.pdf", plot = DensityPlotColony8, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony9 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A20.pdf", plot = DensityPlotColony9, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony10 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A21.pdf", plot = DensityPlotColony10, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony11 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A22.pdf", plot = DensityPlotColony11, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony12 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A23.pdf", plot = DensityPlotColony12, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony13 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A24.pdf", plot = DensityPlotColony13, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony14 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A25.pdf", plot = DensityPlotColony14, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony15 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A26.pdf", plot = DensityPlotColony15, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony16 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A27.pdf", plot = DensityPlotColony16, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony17 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A28.pdf", plot = DensityPlotColony17, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony18 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A29.pdf", plot = DensityPlotColony18, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony19 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A30.pdf", plot = DensityPlotColony19, width = 9.375, height = 6.77, units = "in")

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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
        legend.text = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(size = 18, color = "black", vjust = -6.75, hjust = 0.01)) +
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
                             font.label = list(size = 18, face = "plain"),
                             ncol = 2, nrow = 2,
                             widths = 0.05,
                             vjust = 1.85,
                             label.x = 0.85)

# Assigning the compiled plot to an object, which will be used below 
DensityPlotColony20 <- FullDensityPlot

# Save plot as a PDF
ggsave(file = "Fig_A31.pdf", plot = DensityPlotColony20, width = 9.375, height = 6.77, units = "in")
