MEDA Toolbox for its use in MATLAB.

Coded by: José Camacho Páez, Alejandro Pérez-Villegas
GUI by: Elena Jiménez Mañas, Rafael Rodríguez Gómez
Last modification of this document: 05/Apr/16

References:

 Data set: 

    D.L. Selwood, D.J. Livingstone, J.C.W. Comley, A.B. ODowd, A.T. Hudson, 
       P. Jackson, K.S. Jandu, V.S. Rose, J.N. Stables Structure-Activity 
       Relationships of Antifiral Antimycin Analogues: A Multivariate 
       Pattern Recognition Study, Journal of Medical Chemistry 33 (1990) 
       136142.  

 Analysis:

   J. Camacho. Exploratory Data Analysis using latent subspace models.
       INTECH. ISBN 978-953-51-0438-4. Pages 63 - 90. 2012

This folder contains the data and MATLAB scripts for an EDA example based on PLS. The data set 
was downloaded from http://michem.disat.unimib.it/chm/download/datasets.htm. It consists of 31 
antifilarial antimycin A1 analogues for which 53 physicochemical descriptors were calculated 
for Quantitative Structure-Activity Relationship (QSAR) modelling. The set of descriptors is 
listed below:

Indices Descriptors
1:10 	ATCH1 ATCH2 ATCH3 ATCH4 ATCH5 ATCH6 ATCH7 ATCH8 ATCH9 ATCH10
11:20 	DIPV_X DIPV_Y DIPV_Z DIPMOM ESDL1 ESDL2 ESDL3 ESDL4 ESDL5 ESDL6
21:30 	ESDL7 ESDL8 ESDL9 ESDL10 NSDL1 NSDL2 NSDL3 NSDL4 NSDL5 NSDL6
31:40 	NSDL7 NSDL8 NSDL9 NSDL10 VDWVOL SURF_A MOFI_X MOFI_Y MOFI_Z PEAX_X
41:50 	PEAX_Y PEAX_Z MOL_WT S8_1DX S8_1DY S8_1DZ S8_1CX S8_1CY S8_1CZ LOGP
51:53 	M_PNT SUM_F SUM_R	

These descriptors are used for predicting in vitro antifilarial activity 
(-LOGEC50). This data set has been employed for testing variable selection methods.

The objective of this experiment is to illustrate that the MEDA framework can be used to
understand the data and as a valid means to perform non-automated variable selection.


Items in the folder:

- prepare_data.m: Script to prepare the data for analysis.

- selwood.mat: Data prepared for analysis. The file includes three MATLAB variables:

	- x: matrix 31x54 with the data. 

	- Obs: labels of the observations.

	- Vars: labels of the variables.

- run.m: Script to perform the EDA using the MEDA Toolbox commands. To run the script, the 
	current directory should be the one where the script is stored.

- html: web version of the example
 
