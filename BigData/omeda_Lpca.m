function omeda_vec = omeda_Lpca(Lmodel,pcs,Ltest,dummy,opt,label)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA) for PCA. The original paper is Journal of Chemometrics, 
% DOI: 10.1002/cem.1405. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda_Lpca(Lmodel,pcs,Ltest,dummy) % minimum call
% omeda_vec = omeda_Lpca(Lmodel,pcs,Ltest,dummy,opt,label) %complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs)
%
% Ltest: (struct Lmodel) model with test data:
%       Ltest.XX: (MxM) X-block cross-product matrix.
%       Ltest.centr: (NxM) centroids of the clusters of observations
%       Ltest.multr: (Nx1) multiplicity of each cluster.
%
% dummy: (Nx1) dummy variable containing 1 for the observations in the
%   first group (in test) for the comparison performed in oMEDA, -1 for the 
%   observations in the second group, and 0 for the rest of observations.
%   Also, weights can be introduced.
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: plot oMEDA vector (default)
%       2: plot oMEDA vector and significance limits
%       3: plot oMEDA vector normalized by significance limits
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
%
% OUTPUTS:
%
% omeda_vec: (Mx1) oMEDA vector.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 04/Jul/13.
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
    
%

%% Parameters checking

if nargin < 4, error('Error in the number of arguments.'); end;
s = size(Lmodel.XX);
if ndims(dummy)==2 & find(size(dummy)==max(size(dummy)))==2, dummy = dummy'; end
if s(1) ~= s(2) || ndims(Lmodel.XX)~=2, error('Error in the dimension of the arguments.'); end;
st = size(Ltest.centr);
if st(2)~=s(2) || size(dummy,1)~=st(1), error('Error in the dimension of the arguments.'); end;
if nargin < 5, opt = 1; end;
if nargin < 6 || isempty(label)
    label=[]; 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

Lmodel.lv = max(pcs);
P = Lpca(Lmodel);

omeda_vec = omeda(Ltest.centr,Ltest.multr.*dummy,P);

%% Show results

if opt == 1,
    plot_vec(omeda_vec,label,'d^2_A');
elseif opt == 2 | opt == 3,
    
    calr = Lmodel.centr*P*P';
    l = length(find(dummy));
    vcal = l*sum(calr.^2,1)/size(Lmodel.centr,1);
    lim = 0.9*vcal + 0.1*sum(vcal)/length(vcal); % heuristic
    
    if opt==2
        plot_vec(omeda_vec,label,'d^2_A',[lim;lim]);
    else       
        omeda_vec = omeda_vec./(lim');
        plot_vec(omeda_vec,label,'d^2_A',[ones(1,s(2));-ones(1,s(2))]);
    end
        
end
    


        