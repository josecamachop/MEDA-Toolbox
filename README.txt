MEDA Toolbox for its use in MATLAB. 

Contact person: José Camacho Páez (josecamacho@ugr.es)
Last modification of this document: 14/Oct/19

Installation

	- Extract the rar file in a directory of your choice <directory_path>

	- Add to the MATLAB path the following directories (use command addpath, e.g. addpath '<path>'):
		- <directory_path>
		- <directory_path>/BigData
		- <directory_path>/GUI


Please, acknowledge the use of this software by refercing it: "Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E. Multivariate 
Exploratory Data Analysis (MEDA) Toolbox. Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49-57, available at 
https://github.com/josecamachop/MEDA-Toolbox" Also, please check the documentation of the routines used for more related references. 

Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have 
suggestions, please contact josecamacho@ugr.es

We would like to thanks the direct or indirect contribution of several colleagues:

- E. Szymanska, G.H. Tinnevelt and T.P.J. Offermans for the Sparse Partial Least Squares (SPLS) routine.
- G. Zwanenburg, H.C.J. Hoefsloot, J.A. Westerhuis, J.J. Jansen and A.K. Smilde for the ANOVA Simultaneous Component Analysis (ASCA) routine.
- E. Saccenti for the Horn's Parallel Analysis to determine the number of Principal Components.
- R. Vitale for the Dray's method and permutation testing method to determine the number of Principal Components.

Items in the folder:

- GUIDELINES.txt: Guidelines for the use of the MEDA Toolbox (Please, read first)

- toolbox routines:

	- projection models: pca_pp.m (PCA), kernel_pls.m & simpls.m (PLS), parglm & asca & apca (ANOVA+PCA) gia & gpca (GPCA), sparsepls2 (SPLS), gpls, gasca, SDI, parglmVS (VASCA) 

	- exploratory & visualization tools: SVIplot.m (SVI plots), scores_pca.m & scores_pls.m (Socre plots), loadings_pca.m & loadings_pls.m (Loading plots), 
		meda_pca.m & meda_pls.m (MEDA), omeda_pca.m & omeda_pls.m (oMEDA), mspc_pca.m & mspc_pls.m (MSPC), leverages_pca.m & leverages_pls.m (leverages of variables)
		
	- exploratory tools without visualization: meda.m, omeda.m, mspc.m

	- tools to select the number of LVs: var_pca.m & var_pls.m (Variance plots) crossval_pca, ckf, crossval_pls, dcrossval_pls (Cross-validation routines), PAtest, permsvd, dray

	- tools to select the number of LVs & sparsity in SPLS: crossval_spls, dcrossval_spls, crossval_spls_da, dcrossval_spls_da

	- tools to select the number of LVs & sparsity in GPLS: crossval_gpls, dcrossval_gpls

	- graphical tools: plot_scatter.m, plot_vec.m, plot_map.m 

	- auxiliary routines: preprocess2D.m, preprocess2Dapp.m, seriation.m, spe_lim.m, hot_lim.m   

	- data simulation tools: ADICOV.m, simuleMV.m  

- BigData: Big Data routines.

- GUI: Graphical User Interface routines.

- Examples: Examples of Exploratory Data Analysis, including data sets and MATLAB scripts based on the toolbox.

Copyright (C) 2022  Universidad de Granada
Copyright (C) 2022  José Camacho Páez
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.