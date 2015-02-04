
function [T,TT] = scores_pca(cal,pcs,test,prep,opt,label,classes)

% Compute and plot scores in PCA.
%
% scores_pca(cal,pcs) % minimum call
% scores_pca(cal,pcs,test,prep,opt,label,classes) % complete call
%
% INPUTS:
%
% cal: (LxM) billinear data set for model fitting.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs)
%
% test: (NxM) billinear data set for test. These data are preprocessed in
%   the same way than calibration data.
%
% prep: (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: score plot (default)
%       2: score plot with empty marks
%
% label: ((L+N)x1) name of the variables (numbers are used by default), eg.
%   num2str((1:L+N))')', use ' ' to avoid labels.
%
% classes: ((L+N)x1) vector with the assignment of the variables to classes, 
%   numbered from 1 onwards (1 class by default), eg. ones(L+N),1)
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
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 02/Feb/15.
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
if nargin < 3, x = cal; test = []; else x = [cal;test]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;

[kk, index] = unique(pcs, 'first');
pcs = pcs(sort(index));

if nargin < 4, prep = 2; end;
if nargin < 5, opt = 1; end;
if nargin < 6 || isempty(label)
    label = [];
    %label=num2str((1:s(1))'); 
elseif ~isequal(label,' '),
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end
if nargin < 7 || isempty(classes)
    classes = [];
else
    if ndims(classes)==2 & find(size(classes)==max(size(classes)))==2, classes = classes'; end
    if size(classes,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end

%% Main code

[calp,m,dt] = preprocess2D(cal,prep);
[P,T] = pca_pp(calp,max(pcs));

if exist('test')&~isempty(test),
    testp = (test - ones(size(test,1),1)*m)./(ones(size(test,1),1)*dt);
    TT = testp*P;
else
    TT = [];
end

if opt,
    T = [T;TT];
    if length(pcs) == 1,
        plot_vec(T(:,pcs), label, sprintf('PC %d',pcs), [], 0, 'r');
    end
    for i=1:length(pcs)-1,
        for j=i+1:length(pcs),
            plot_scatter([T(:,pcs(i)),T(:,pcs(j))],label,classes,{sprintf('PC %d',pcs(i)),sprintf('PC %d',pcs(j))}',opt-1);
        end      
    end
end
        