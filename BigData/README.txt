Large Data Extension for the MEDA Toolbox for its use in MATLAB.

Coded by: José Camacho Páez (josecamacho@ugr.es)
Last modification of this document: 05/Aug/14

Please, acknowledge the use of this software by referecing it: "MEDA Toolbox, available at https://github.com/josecamachop/MEDA-Toolbox/archive/master.zip. José Camacho and Elena Jiménez, 
Network Engineering and Security Group, University of Granada, Spain."

Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have suggestions, please contact josecamacho@ugr.es


Items in the folder:

- Big Data toolbox routines:

	- projection models: Lpca.m (PCA), Lpls.m (PLS)

	- exploratory tools: scores_Lpca.m & scores_Lpls.m (Socre plots), loadings_Lpca.m & loadings_Lpls.m (Loading plots), 
		meda_Lpca.m & meda_Lpls.m (MEDA), omeda_Lpca.m & omeda_Lpls.m (oMEDA) 

	- statistical process control tools: sqresiduals_Lpca.m & sqresiduals_Lpls.m (Squared Residuals plots)

	- tools to select the number of LVs: var_Lpca.m & var_Lpls.m (Variance plots) 

	- graphical tools: plot_Lscatter.m, plot_Lvec.m

	- auxiliary routines: preprocess2Di.m, centN.m, histN.m, psc.m 

	- Lmodel management routines: Lmodel_ini.m, update_ewma.m, update_iterative.m

	- file management routines: add_data.m, add_indices.m, cfilesys.m, read_data.m, read_indices.m, VCfile.m 
  

You may find an example of use in the Examples\Networkmetric\KDD directory.

Copyright (C) 2014  José Camacho Páez
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.