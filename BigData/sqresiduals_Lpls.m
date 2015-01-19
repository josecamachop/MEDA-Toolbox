function E = sqresiduals_Lpca(Lmodel,pcs,Ltest,opt,label)

% Compute and plot squared residuals in PCA for large data.
%
% sqresiduals_pca(cal,pcs) % minimum call
% sqresiduals_pca(cal,pcs,test,prep,opt,label) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.centr: (LxM) centroids of the clusters of observations
%       Lmodel.multr: (Lx1) multiplicity of each cluster.
%       Lmodel.class: (Lx1) class associated to each cluster.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs)
%
% Ltest: (struct Lmodel) model with test data:
%       Ltest.XX: (MxM) X-block cross-product matrix.
%       Ltest.centr: (NxM) centroids of the clusters of observations
%       Ltest.multr: (Nx1) multiplicity of each cluster.
%       Ltest.class: (Nx1) class associated to each cluster.
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: Squared residuals in the observations (default)
%       2: Squared residuals in the variables 
%
% label: name of the observations (opt 0 o 1, dimension ((L+N)x1) or 
%   variables (opt 2, dimension (Mx1)) (numbers are used by default), eg.
%   num2str((1:L+N))')'
%
%
% OUTPUTS:
%
% E: squared residuals (opt 0 o 1, dimension {1x(L+N)} or opt 2, dimension 
%   {1x(M)})
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 07/May/13.
%
% Copyright (C) 2014  José Camacho Páez
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

%% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3 || isempty(Ltest), x = Lmodel.centr; test = []; else x = [Lmodel.centr;Ltest.centr]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(pcs);
if sp < 2, error('Error in the dimension of the arguments.'); end;
if nargin < 4, opt = 1; end;
if nargin < 5, 
    switch opt,
        case 2
            label=num2str((1:s(2))');
            a = 'v'*ones(size(cal,2),1);
            label = [a label];
        otherwise
            label=num2str([1:size(cal,1) 1:size(test,1)]');
            a = ['c'*ones(size(cal,1),1);'t'*ones(size(test,1),1)];
            label = [a label];
    end;
end

%% Main code

Lmodel.lv = max(pcs);
P = Lpca(Lmodel);
T = Lmodel.centr*P;

if exist('test')&~isempty(test),
    TT = Ltest.centr*P;
    classes = [Lmodel.class;Ltest.class];
    mult = [Lmodel.multr;Ltest.multr];
else
    TT = [];
    classes = Lmodel.class;
    mult = [Lmodel.multr];
end


switch opt,
    case 2
        E = sum((Lmodel.centr - T*P').^2,1);
        if ~isempty(test)
            E = [E;sum((Ltest.centr - TT*P').^2,1)]';
        end
    otherwise
        E = sum(([Lmodel.centr;Ltest.centr] - [T;TT]*P').^2,2);
end;

plot_vec(E,label,'Squared Residuals');
        