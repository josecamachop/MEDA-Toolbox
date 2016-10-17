function ok = check_Lmodel(Lmodel)

% Check Lmodel integrity
%
% ok = check_Lmodel(Lmodel) % complete call
%
%
% INPUT:
%
% Lmodel.type: (1x1) PCA (1) o PLS (2)
%
% Lmodel.update: (1x1) EWMA (1) or ITERATIVE (2)
%
% Lmodel.lvs: (1x1) number of latent variables (e.g. lvs = 1:2 selects the
%   first two LVs). By default, Lmodel.lvs = 1:rank(xcs)
%
% Lmodel.N: (1x1) number of effective observations in the model.
%
% Lmodel.prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
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
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
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
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 17/Oct/16.
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
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
M = size(Lmodel.XX, 2);

% Convert column arrays to row arrays
if size(Lmodel.lvs,2) == 1, Lmodel.lvs = Lmodel.lvs'; end;
if size(Lmodel.av,2) == 1, Lmodel.av = Lmodel.av'; end;
if size(Lmodel.sc,2) == 1, Lmodel.sc = Lmodel.sc'; end;
if size(Lmodel.weight,2) == 1, Lmodel.weight = Lmodel.weight'; end;
if size(Lmodel.avy,2) == 1, Lmodel.avy = Lmodel.avy'; end;
if size(Lmodel.scy,2) == 1, Lmodel.scy = Lmodel.scy'; end;
if size(Lmodel.weighty,2) == 1, Lmodel.weighty = Lmodel.weighty'; end;

% Preprocessing
Lmodel.lvs = unique(Lmodel.lvs);
A = length(Lmodel.lvs);

% Validate dimensions of input data
assert (M>0, 'Dimension Error: Lmodel.XX with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.XX), [M M]), 'Dimension Error: Lmodel.XX must be M-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (A>0, 'Dimension Error: Lmodel.lvs with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.lvs), [1 A]), 'Dimension Error: Lmodel.lvs must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.prep), [1 1]), 'Dimension Error: Lmodel.prep must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.prepy), [1 1]), 'Dimension Error: Lmodel.prep must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(Lmodel.type~=1 & Lmodel.type~=2)), 'Value Error: Lmodel.type must contain 1 or 2. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(Lmodel.update~=1 & Lmodel.update~=2)), 'Value Error: Lmodel.update must contain 1 or 2. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(Lmodel.lvs<0)), 'Value Error: Lmodel.lvs must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(Lmodel.lvs), Lmodel.lvs), 'Value Error: Lmodel.lvs must contain integers. Type ''help %s'' for more info.', routine(1).name);

ok = true;