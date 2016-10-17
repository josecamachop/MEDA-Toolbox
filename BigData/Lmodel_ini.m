function Lmodel = Lmodel_ini

% Large model inicialization
%
% Lmodel = Lmodel_ini % complete call
%
%
% OUTPUTS:
%
% Lmodel.type: (1x1) PCA (1, default) o PLS (2)
%
% Lmodel.update: (1x1) EWMA (1) or ITERATIVE (2, default)
%
% Lmodel.lvs: (1x1) number of latent variables (e.g. lvs = 1:2 selects the
%   first two LVs). By default, Lmodel.lvs = 1:rank(xcs)
%
% Lmodel.N: (1x1) number of effective observations in the model.
%
% Lmodel.prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering 
%       2: auto-scaling (default)
%
% Lmodel.av: (1xM) sample average according to the preprocessing method.
%
% Lmodel.sc: (1xM) sample scale according to the preprocessing method.
%
% Lmodel.weight: (1xM) weight applied after the preprocessing method.
%
% Lmodel.XX: (MxM) sample cross-product matrix of X.
%
% Lmodel.prepy: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (default) 
%
% Lmodel.avy: (1xM) sample average according to the preprocessing method.
%
% Lmodel.scy: (1xM) sample scale according to the preprocessing method.
%
% Lmodel.weighty: (1xM) weight applied after the preprocessing method.
%
% Lmodel.XY: (MxO) sample cross-product matrix of X and Y.
%
% Lmodel.YY: (OxO) sample cross-product matrix of Y.
%
% Lmodel.nc: (1x1) number of clusters in the model.
%
% Lmodel.centr: (NxM) centroids of the clusters of observations
%
% Lmodel.multr: (Nx1) multiplicity of each cluster.
%
% Lmodel.class: (Nx1) class associated to each cluster.
%
% Lmodel.updated: (Nx1) specifies whether a data point is new.
%
% Lmodel.obs_l: {Nx1} label of each cluster.
%
% Lmodel.var_l: {Nx1} label of each variable.
%
% Lmodel.mat: (MxA) projection matrix for distance computation.
%
% Lmodel.index_fich: {Nx1} file system with the original observations in
%   each cluster for ITERATIVE models.
%
% Lmodel.path: (str) path to the file system for ITERATIVE models.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Sep/15.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

Lmodel.type = 1;
Lmodel.update = 2;
Lmodel.lvs = 0;
Lmodel.N = 0;
Lmodel.prep = 2;
Lmodel.av = 0;
Lmodel.sc = 0;
Lmodel.weight = 0;
Lmodel.XX = 0;
Lmodel.prepy =2;
Lmodel.avy = 0;
Lmodel.scy = 0;
Lmodel.weighty = 0;
Lmodel.XY = 0;
Lmodel.YY = 0;
Lmodel.nc = 0;
Lmodel.centr = [];
Lmodel.multr = [];
Lmodel.class = [];
Lmodel.updated = [];
Lmodel.obs_l = {};
Lmodel.vars_l = {};
Lmodel.mat = [];
Lmodel.index_fich = {};
Lmodel.path='';