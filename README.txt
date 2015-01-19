MEDA Toolbox for its use in MATLAB.

Coded by: José Camacho Páez
GUI by: Elena Jiménez Mañas
Last modification of this document: 05/Aug/14

Installation

	- Extract the rar file in a directory of your choice <directory_path>

	- Add to the MATLAB path the following directories (use command addpath, e.g. addpath '<path>'):
		- <directory_path>
		- <directory_path>/BigData
		- <directory_path>/GUI


Please, acknowledge the use of this software by referecing it: "MEDA Toolbox, available at https://github.com/josecamachop/MEDA-Toolbox/archive/master.zip. José Camacho and Elena Jiménez, 
Network Engineering and Security Group, University of Granada, Spain." 


Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have suggestions, please contact josecamacho@ugr.es


Items in the folder:

- GUIDELINES.txt: Guidelines for the use of the EDA Toolbox (Please, read first)

- toolbox routines:

	- projection models: pca_pp.m (PCA), kernel_pls.m (PLS)

	- exploratory tools: SVIplot.m (SVI plots), scores_pca.m & scores_pls.m (Socre plots), loadings_pca.m & loadings_pls.m (Loading plots), 
		meda.m, meda_pca.m & meda_pls.m (MEDA), omeda.m, omeda_pca.m & omeda_pls.m (oMEDA) 

	- statistical process control tools: sqresiduals_pca.m & sqresiduals_pls.m (Squared Residuals plots), leverage_pca.m & leverage_pls.m (Leverage or Hotelling T2 plots)

	- tools to select the number of LVs: var_pca.m & var_pls.m (Variance plots) crossval2D_pca & crossval2D_pls (Cross-validation rutines)

	- graphical tools: plot_scatter.m, plot_vec.m, plot_map.m 

	- auxiliary routines: preprocess2D.m, seriation.m, spe_lim.m, hot_lim.m   

	- data simulation tools: ADICOV.m  

- BigData: Big Data routines.

- GUI: Graphical User Interface routines.

- Examples: Examples of Exploratory Data Analysis, including data sets and MATLAB scripts based on the toolbox.


Copyright (C) 2014  José Camacho Páez
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.