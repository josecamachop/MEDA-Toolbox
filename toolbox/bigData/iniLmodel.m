function Lmodel = iniLmodel(X,Y,varargin)

% Large model inicialization
%
% iniLmodel % minimum call
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting
%
% Y: [NxO] billinear data set of predicted variables
%
% Optional INPUTS (parameter):
%
% 'ObsLabel': {Nx1} label of each observation.
%
% 'VarLabel': {Mx1} label of each variable.
%
% OUTPUTS:
%
% Lmodel.centr: [NxM] centroids of the clusters of observations.
%
% Lmodel.centrY: [NxL] responses of centroids of the clusters of
% observations.
%
% Lmodel.nc: [1x1] number of clusters in the model.
%
% Lmodel.multr: [ncx1] multiplicity of each cluster.
%
% Lmodel.class: [ncx1] class associated to each cluster.
%
% Lmodel.vclass: [Mx1] class associated to each variable.
%
% Lmodel.N: [1x1] number of effective observations in the model.
%
% Lmodel.type: [1x1] 'PCA', 'PLS' or 'ASCA'
%
% Lmodel.update: [1x1] EWMA (1) or ITERATIVE (2)
%
% Lmodel.XX: [MxM] sample cross-product matrix of X.
%
% Lmodel.lvs: [1x1] number of latent variables (e.g. lvs = 1:2 selects the
%   first two LVs). By default, Lmodel.lvs = 1:rank(xcs)
%
% Lmodel.prep: [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
%
% Lmodel.av: [1xM] sample average according to the preprocessing method.
%
% Lmodel.sc: [1xM] sample scale according to the preprocessing method.
%
% Lmodel.weight: [1xM] weight applied after the preprocessing method.
%
% Lmodel.updated: [ncx1] specifies whether a data point is new.
%
% Lmodel.obsl: {ncx1} label of each cluster.
%
% Lmodel.varl: {ncx1} label of each variable.
%
% Lmodel.mat: [MxA] projection matrix for distance computation.
%
% Lmodel.prepy: [1x1] preprocesing of the data
%       0: no preprocessing.
%       1: mean centering (default) 
%       2: auto-scaling (centers and scales data so that each variable 
%           has variance 1)
%
% Lmodel.avy: [1xM] sample average according to the preprocessing method.
%
% Lmodel.scy: [1xM] sample scale according to the preprocessing method.
%
% Lmodel.weighty: [1xM] weight applied after the preprocessing method.
%
% Lmodel.XY: [MxO] sample cross-product matrix of X and Y.
%
% Lmodel.YY: [OxO] sample cross-product matrix of Y.
%
% Lmodel.indexfich: {ncx1} file system with the original observations in
%   each cluster for ITERATIVE models.
%
% Lmodel.path: (str) path to the file system for ITERATIVE models.
%
%
% EXAMPLE OF USE:
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 21/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
% 
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


if nargin < 1, X = []; end;
N = size(X, 1);
M = size(X, 2);
if nargin < 2, Y = []; end;

p = inputParser;
    if N>0 
        obsl = cellstr(num2str((1:N)')); 
    else
        obsl = {}; 
    end;
addParameter(p,'ObsLabel',obsl);   
    if M>0
        varl = cellstr(num2str((1:M)')); 
    else
        varl = {};
    end
addParameter(p,'VarLabel',varl);    
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
obsl = p.Results.ObsLabel;
varl = p.Results.VarLabel;


Lmodel.centr = X;
Lmodel.centrY = Y;
Lmodel.obsl = obsl;
Lmodel.varl =  varl;

Lmodel.multr = []; 
Lmodel.class = []; 
Lmodel.vclass = []; 
Lmodel.updated = [];
Lmodel.lvs = [];
Lmodel.prep = [];
Lmodel.sc = [];
Lmodel.weight = [];
Lmodel.prepy = [];
Lmodel.avy = [];
Lmodel.scy = [];
Lmodel.weighty = [];
Lmodel.XX = [];
Lmodel.XY = [];
Lmodel.YY = [];
Lmodel.mat = [];
Lmodel.indexfich = {};
Lmodel.path = '';


[kk,Lmodel] = checkLmodel(Lmodel);