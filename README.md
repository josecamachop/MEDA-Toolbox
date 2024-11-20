Multivariate Exploratory Data Analysis (MEDA) Toolbox for its use in MATLAB 

Contact person: Jose Camacho (josecamacho@ugr.es)
Last modification of this document: 20/Nov/24

Installation

1 - Download and extract (or git clone) in a directory of your choice <directory_path>

2 - Add to the MATLAB path the following directory with its subdirectories:

	- <directory_path>\toolbox

For step 2, in the MATLAB command line, you may use commands 'addpath' and 'genpath', e.g.: 

	>> addpath(genpath('<directory_path>\toolbox')) 

and then 'savepath'. You may alternatively use the 'pathtool' from the command line or the menu.

Please, acknowledge the use of this software by referencing: "Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E. Multivariate Exploratory Data Analysis (MEDA) Toolbox. Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49-57, available at https://github.com/CoDaSLab/MEDA-Toolbox". Please check the documentation of the routines for more related references. 

Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have suggestions, please contact josecamacho@ugr.es

Items in the folder:

- GUIDELINES.md: Guidelines for the use of the MEDA Toolbox (Please, read first)

- LICENSE.txt: License for reuse

- toolbox: latent variable models and auxility methods, including pcaEig (PCA), kernelpls & simpls (PLS), asca & apca & vasca (ANOVA+PCA) 

   	- anova: routines for factorization of data coming from an experimental design
   	- auxiliary: general routines for hanlding input/output data
   	- bigData: extended routines for massive observations
   	- dataSim: routines for data simulation
   	- dataViz: routines for data visualization, including line plots, scatter plots, and map plots
   	- GUI: Graphical User Interface routines
   	- internal: routines not meant to be directly used, but rather called from other routines
   	- modelSel: routines for model selection, including scree plots and crossvalidation, among others
   	- modelVal: double cross-validation routines
   	- modelViz: routines for the visualization of latent variable models, including score and loading plots, biplots and 		other tools
   	- sparse: sparse versions of latent variable models, inclding group-wise models

- examples: Examples of Exploratory Data Analysis, including data sets and MATLAB scripts based on the toolbox
  
- techReports: Technical Reports that make use of the tolbox

We would like to thank the direct or indirect contribution of several colleagues:

- E. Szymanska, G.H. Tinnevelt and T.P.J. Offermans for the Sparse Partial Least Squares (SPLS) routine.

- G. Zwanenburg, H.C.J. Hoefsloot, J.A. Westerhuis, J.J. Jansen and A.K. Smilde for the original ANOVA Simultaneous Component Analysis (ASCA) routine.

- E. Saccenti for the Horn's Parallel Analysis to determine the number of Principal Components.

- R. Vitale for the Dray's method and permutation testing method to determine the number of Principal Components.

Copyright (C) 2024  Universidad de Granada
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
