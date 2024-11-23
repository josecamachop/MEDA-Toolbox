ï»¿MEDA Toolbox for its use in MATLAB.

Codified by: JosÃ© Camacho PÃ¡ez.
GUI by: Elena JimÃ©nez MaÃ±as
Version: 2.1
Last modification of this document: 14/Jul/13

References:

 Data set: 

      R.S. Gray, D. Kotz, C. Newport, N. Dubrovsky, A. Fiske, J. Liu, C. Masone, S. McGrath, Y.
	Yuan, Outdoor experimental comparison of four ad hoc routing algorithms, in: Proceedings of
	the 7th ACM international symposium on Modeling, analysis and simulation of wireless and
	mobile systems, MSWiM â04, 2004, pp. 220-229.


This folder contains the data for an EDA example based on PLS. The data set 
was downloaded from http://crawdad.cs.dartmouth.edu/keyword-MANET.html. It consists of an outdoor 
experiment for the comparison of four different routing algorithms in a mobile ad.hoc network (MANET) 
formed by 33 laptops in movement. The evaluated algorithms are: Any Path Routing without Loops (APRL), 
Ad hoc On-demand Distance Vector (AODV), On-Demand Multicast Routing Protocol (ODMRP) and System and
Traffic-dependent Adaptive Routing Algorithm (STARA). A set of statistics listed below are computed 
from the original data at regular intervals of time, yielding a total of 100 intervals. 

Indices Descriptors
1 PD Average distance between laptops
2 mM Minimum value for max. distances
3 Mm Minimum value for min. distances
4 cX X centroid X
5 cY Y centroid Y
6 cZ Z centroid Z
7 n1 Amount of laptops with a distance to the centroid lower than 1/32 of the max. distance
8 n2 Amount of laptops with a distance to the centroid between 1/32 and 2/32 of the max. distance
9 n3 Amount of laptops with a distance to the centroid between 2/32 and 3/32 of the max. distance
10 n4 Amount of laptops with a distance to the centroid higher than 3/32 of the max. distance
11 nTI Number of TIN
12 nTO Number ofTOUT
13 nSI Number ofSIN
14 nSO Number of SOUT
15 vTI Volumeof TIN
16 vTO VolumeofTOUT
17 vSI VolumeofSIN
18 vSO VolumeofSOUT	

These descriptors are used for predicting the differences between the four routing algorithms.

The objective of this experiment is to illustrate that the MEDA framework can be used to
understand the data.


Items in the folder:

- manet.mat: Data prepared for analysis. The file includes three MATLAB variables:

	- x: matrix 70x18 with the data. 

	- y: matrix 70x4, with the PLS clasification.

	- clases: classes of the observations.

	- laby: labels of the observations.

	- lavel_v: labels of the variables.
