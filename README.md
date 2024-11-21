<h1>Multivariate Exploratory Data Analysis (MEDA) Toolbox for its use in MATLAB</h1> 

Contact person: Jose Camacho (josecamacho@ugr.es) <br>
Last modification of this document: 20/Nov/24

The Multivariate Exploratory Data Analysis (MEDA) Toolbox in Matlab is a set of multivariate analysis tools for the exploration of data sets. 
There are several alternative tools in the market for that purpose, both commercial and free. The PLS_Toolbox from Eigenvector Inc. 
is a very nice example. The MEDA Toolbox is not intended to replace or being a competitor of any of these toolkits. Rather, the MEDA
Toolbox is a complementary tool which includes several of our recent contributions to the field. Thus, traditional exploratory plots 
based on Principal Component Analysis (PCA) or Partial Least Squares (PLS), such as score, loading and residual plots, are combined 
with new methods: MEDA, oMEDA, SVI plots, ADICOV, EKF & CKF crossvalidation, CSP, GPCA, SPLS crossvalidation, ASCA and VASCA .... 

The MEDA Toolbox can be used to analyze normal size data sets (several hundreds of observations times several hundreds of variables)
There is also an extension of the toolbox for large data sets, with millions of items, under folder BigData. 

The MEDA Toolbox is expected to work on Octave.

Please, acknowledge the use of this software by referencing: "Camacho, J., Pérez, A., Rodríguez, R., Jiménez-Mañas, E. Multivariate Exploratory Data Analysis (MEDA) Toolbox. Chemometrics and Intelligent Laboratory Systems, 2015, 143: 49-57, available at https://github.com/CoDaSLab/MEDA-Toolbox". Please check the documentation of the routines for more related references. 

Please, note that the software is provided "as is" and we do not accept any responsibility or liability. Should you find any bug or have suggestions, please contact josecamacho@ugr.es

<h2>Installation</h2>

1 - Download and extract (or git clone) in a directory of your choice <directory_path>

2 - Add to the MATLAB path the following directory with its subdirectories:

	- <directory_path>\toolbox

For step 2, in the MATLAB command line, you may use commands 'addpath' and 'genpath', e.g.: 

	>> addpath(genpath('<directory_path>\toolbox')) 

and then 'savepath'. You may alternatively use the 'pathtool' from the command line or the menu.

<h2>Working with the MEDA Toolbox</h2>
 
There are two ways to work with the MEDA toolbox: using the GUI (starting users) and using the commands (expert user). The GUI is self 
explanatory. It incorporates help messages to guide its use. To launch the GUI, type "MEDA" in the command line of Matlab after the 
installation. 

The commands provide of much more functionality. Each command includes helping information, that can be seen by typing "help <command>" 
in the command line of Matlab. This helping information includes examples. Also, in the Examples directory, several real data examples 
are included. It is suggested to have a look at these examples.

There are two types of commands or routines. Those that do not provide visualization (e.g. meda, mspc, omeda, ...) are optimized for data 
computing. Those that issue a plot (e.g medaPca, medaPls, scoresPca, ...) are optimized for interactive scripting


<h2>Items in the folder</h2>

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
   	- modelViz: routines for the visualization of latent variable models, including score and loading plots, biplots and other tools
   	- sparse: sparse versions of latent variable models, inclding group-wise models

- examples: Examples of Exploratory Data Analysis, including data sets and MATLAB scripts based on the toolbox
  
- techReports: Technical Reports that make use of the tolbox


<h2>The Tools</h2>

An introduction to the exploratory tools in the MEDA toolbox can be found in:

- Loading plots, Score plots, Residual plots, Leverage plots ==> Kim H. Esbensen. Multivariate Data Analysis: in practice. Camo. ISBN: 82-993330-3-2.

- MEDA ==> Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp. 8-18.  

- oMEDA ==> Journal of Chemometrics, 2011, 25 (11), pp. 592-600.

- SVI plots ==> Chemometrics and Intelligent Laboratory Systems 100, 2010, pp. 48-56.

- ADICOV ==> Chemometrics and Intelligent Laboratory Systems 105(2), 2011, pp. 171-180.

- EKF ==> Chemometrics and Intelligent Laboratory Systems 131, 2014, pp. 37-50

- CKF ==> Journal of Chemometrics, 29(8), 2015, pp. 467-478.

- CSP ==> Chemometrics and Intelligent Laboratory Systems, 135, 2014, pp. 110-125.

- GPCA ==> Journal of Computational and Graphical Statistics , 2017, 26 (3): 501-512.

- simuleMV ==> Chemometrics and Intelligent Laboratory Systems, 2017, 160: 40-51.

- GPLS ==> Journal of Chemometrics, 2018, 32: 1-11.

- GASCA ==> Metabolomics. 2018; 14(6): 73.

- XCAN ==> Chemometrics and Intelligent Laboratory Systems, 2020, 203: 104038.

- SDI ==> Chemometrics and Intelligent Laboratory Systems, 2020, 206: 104160.

- VASCA ==> Bioinformatics, 2023, 39 (1): btac795.

Loading and Score plots in the toolbox are 2D scatter plots in the sub-space of the PCA or PLS model. Residual plots are line/bar plots
of the squared residuals of the variables/observations. Leverage plots are line/bar plots of the squared scores of the variables/observations.

MEDA plots are map plots that contain the relationship between each pair of variables. They look similar to correlation plots, but 
they present better properties for exploration.

oMEDA plots are bar plots that compute the differences between groups of observations in the latent subspace. They are similar to 
contribution plots but again with better properties.

SVI plots are useful to understand how a variable is related to the rest of the data set. The reading of the reference above is recommended.

ADICOV is a method to approximate a specific data distribution imposing a certain correlation structure. It is very useful for realistic data
simulation.

EKF and CKF are different versions of cross-validation in PCA, useful in some applications to infer an adequate number of PCs.

CSP are compressed versions of score plots, where millions of observations can be visualized thanks to clustering techniques.

GPCA is a PCA modification where loadings are constrain to one group of variables.

simuleMV is a simulation engine for multivariate data.

GPLS and GASCA are the GPCA extension to PLS and ASCA, respectively.

SDI is a procedure to find interesting subspaces in a PLS-DA model.

VASCA is a variable-selection version of ASCA.

<h2>ACKs</h2>

We would like to thank the direct or indirect contribution of several colleagues:

- E. Szymanska, G.H. Tinnevelt and T.P.J. Offermans for the Sparse Partial Least Squares (SPLS) routine.

- G. Zwanenburg, H.C.J. Hoefsloot, J.A. Westerhuis, J.J. Jansen and A.K. Smilde for the original ANOVA Simultaneous Component Analysis (ASCA) routine.

- E. Saccenti for the Horn's Parallel Analysis to determine the number of Principal Components.

- R. Vitale for the Dray's method and permutation testing method to determine the number of Principal Components.


<h2>Copyright</h2>

Copyright (C) 2024  Universidad de Granada
 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
