function T2 = leverage_Lpca(Lmodel,pcs,Ltest,opt,label)

% Compute and plot leverages in PCA for large data.
%
% leverage_Lpca(Lmodel,pcs) % minimum call
% leverage_Lpca(Lmodel,pcs,Ltest,opt,label) % complete call
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
% T2: leverage (opt 0,1 or 3, dimension {1x(L+N)} or opt 2 or 4, dimension 
%   {1x(M)})
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 24/Sep/15.
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

T2 = [];

if ~isempty(pcs)
    Lmodel.lv = max(pcs);
    P = Lpca(Lmodel);
    P = P(:,pcs);
    T = Lmodel.centr*P;


    if exist('Ltest')&~isempty(Ltest),
        TT = Ltest.centr*P;
        mult = [Lmodel.multr;Ltest.multr];
    else
        TT = [];
        mult = [Lmodel.multr];
    end
    
    [Ts,notused,dtT2] = preprocess2D(T,2);

    switch opt,
        case 2
            T2 = diag(P*P');
        otherwise
            T2c = diag(Ts*Ts');
            if ~isempty(test)
                T2 = diag(TT*diag(1./(dtT2.^2))*TT');
            else
                T2 = T2c;
            end
    end;
    
    switch opt,
        case 1
            plot_Lvec(T2,mult,label,'D-statistic',[],[10 100 1000]);
        case 2
            plot_vec(T2,label,'D-statistic');
        case 3
            plot_Lvec(T2,mult,label,'D-statistic',(ones(size(T2,1),1)*[hot_lim(length(pcs),length(T2c),0.05) hot_lim(length(pcs),length(T2c),0.01)])',[10 100 1000]);
        case 4
             plot_Lvec(T2,mult,label,'D-statistic',(ones(size(T2,1),1)*[hot_lim(length(pcs),length(T2c),0.05,1) hot_lim(length(pcs),length(T2c),0.01,1)])',[10 100 1000]);
    end
end


        