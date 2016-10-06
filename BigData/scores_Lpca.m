function [T,TT] = scores_Lpca(Lmodel,pcs,Ltest,opt,label)

% Compute and plot scores in PCA for large data.
%
% scores_Lpca(Lmodel,pcs) % minimum call
% scores_Lpca(Lmodel,pcs,Ltest,opt,label,classes) % complete call
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
%       0: no plots.
%       1: score plot 
%       2: score plot with empty marks
%       3: 2D compressed score plot, with the multiplicity info in the
%           markers
%       4: 2D compressed score plot (by default), with the multiplicity
%           info in the size of the markers
%       5: 3D compressed score plot, with the multiplicity information in
%           the Z axis
%       6: 2D hexagonal binning plot, with the multiplicity info in the 
%           size of the markers and class in the Z axis.
%
% label: ((L+N)x1) labels of the observations. Some possible inputs are:
%   - num2str((1:L+N)') or [] for observation numbers (by default)
%   - num2str(Lmodel_it.class) for observation classes
%   - ' ' to avoid labels.
%
%
% OUTPUTS:
%
% T: (LxA) calibration scores.
%
% TT: (NxA) test scores.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/Apr/14.
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
    

%% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if nargin < 3 || isempty(Ltest), x = Lmodel.centr; test = []; else x = [Lmodel.centr;Ltest.centr]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(pcs);
if sp < 2, error('Error in the dimension of the arguments.'); end;
if nargin < 4, opt = 4; end;
if nargin < 5 || isempty(label)
    label=num2str((1:s(1))'); 
% elseif ~isequal(label,' '),
%     if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
%     if size(label,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end

%% Main code

Lmodel.lv = max(pcs);
P = Lpca(Lmodel);
T = Lmodel.centr*P;

if exist('Ltest')&~isempty(Ltest)&~isempty(Ltest.centr),
    TT = Ltest.centr*P;
    classes = [Lmodel.class;Ltest.class];
    mult = [Lmodel.multr;Ltest.multr];
else
    TT = [];
    classes = Lmodel.class;
    mult = [Lmodel.multr];
end

if opt,
    T = [T;TT];
    for i=1:length(pcs)-1,
        for j=i+1:length(pcs),
            plot_Lscatter([T(:,pcs(i)),T(:,pcs(j))],label,classes,{sprintf('PC %d',pcs(i)),sprintf('PC %d',pcs(j))},opt-1,mult,10.^(floor(log10(max(mult)))/3:floor(log10(max(mult)))/3:floor(log10(max(mult)))));
        end      
    end
end



