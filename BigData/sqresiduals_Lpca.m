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
%   variables (opt 2, dimension (Mx1)) (numbers are used by default for
%   variables and multiplicity for observations), eg.   num2str((1:L+N))')'
%
%
% OUTPUTS:
%
% E: squared residuals (opt 0 o 1, dimension {1x(L+N)} or opt 2, dimension 
%   {1x(M)})
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Sep/15.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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
    label = [];
end

%% Main code

Lmodel.lv = max(pcs);
P = Lpca(Lmodel);
P = P(:,pcs);
T = Lmodel.centr*P;

if exist('Ltest')&~isempty(Ltest),
    TT = Ltest.centr*P;
    mult = [Lmodel.multr;Ltest.multr];
else
    TT = [];
    classes = Lmodel.class;
    mult = [Lmodel.multr];
end

maxv = [0 1 10 100 1000 Inf];
switch opt,
    case 2
        E = diag(Lmodel.XX - P*P'*Lmodel.XX*P*P');
        if ~isempty(Ltest)
            E = [E diag(Ltest.XX - P*P'*Ltest.XX*P*P')];
        end
        plot_vec(E,label,'Squared Residuals');
    case 1
        for n_clus = 1:s(1), % Cálculo de la covarianza
            if Lmodel.multr(n_clus)>1,
                if isempty(Lmodel.index_fich)
                    x = Lmodel.centr(n_clus,:)*Lmodel.multr(n_clus);
                    XX = x'*x; % approximate
                else
                    XX = VCfile(Lmodel.index_fich{n_clus},s(2),0,Lmodel.path);
                end
            else
                XX = Lmodel.centr(n_clus,:)'*Lmodel.centr(n_clus,:);
            end
            E(n_clus) = sum(diag(XX - P*P'*XX*P*P'));
        end
        mult=Lmodel.multr;
        if ~isempty(Ltest)
            mult = [mult;Ltest.multr];
        for n_clus2 = 1:size(Ltest.centr,1), % Cálculo de la covarianza
            if Ltest.multn(n_clus2)>1,
                if isempty(Ltest.index_fich)
                    x = Ltest.centr(n_clus,:)*Ltest.multr(n_clus);
                    XX = x'*x; % approximate
                else
                    XX = VCfile(Ltest.index_fich{n_clus},s(2),0,Ltest.path);
                end
            else
                XX = Ltest.centr(n_clus,:)'*Ltest.centr(n_clus,:);
            end
            E(n_clus+n_clus2) = sum(diag(XX - P*P'*XX*P*P'));
        end
        end
        plot_Lvec(E,mult,label,'Squared Residuals',[],[10 100 1000]);
    otherwise 
        E = sum((Lmodel.centr - T*P').^2,2);
        mult=Lmodel.multr;
        if ~isempty(Ltest)
            mult = [mult;Ltest.multr];
            E = [E ; sum((Ltest.centr - TT*P').^2,2)];
        end
        plot_Lvec(E,mult,label,'Squared Residuals',[],[10 100 1000]);
end;


        