
function P = loadings_Lpca(Lmodel,pcs,opt,label,classes)

% Compute and plot loadings in PCA for large data.
%
% loadings_Lpca(Lmodel,pcs) % minimum call
% loadings_Lpca(Lmodel,pcs,opt,label,classes) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%
% pcs: (1xA) Principal Components considered (e.g. pcs = 1:2 selects the
%   first two pcs)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: loading plot of P (default)
%       2: loading plot of P with empty marks
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
% classes: (Mx1) vector with the assignment of the variables to classes, 
%   numbered from 1 onwards (1 class by default), eg. ones(M,1)
%
%
% OUTPUTS:
%
% P: (MxA) scores.
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
s = size(Lmodel.XX);
if nargin < 3, opt = 1; end;
if nargin < 4 || isempty(label)
    label=num2str((1:s(2))'); 
else
    if ndims(label)==2 & find(size(label)==max(size(label)))==2, label = label'; end
    if size(label,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end
if nargin < 5 || isempty(classes)
    classes = ones(s(2),1); 
else
    if ndims(classes)==2 & find(size(classes)==max(size(classes)))==2, classes = classes'; end
    if size(classes,1)~=s(2), error('Error in the dimension of the arguments.'); end;
end

%% Main code

Lmodel.lv = max(pcs);
P = Lpca(Lmodel);

if opt,
    for i=1:length(pcs)-1,
        for j=i+1:length(pcs),
            plot_scatter([P(:,pcs(i)),P(:,pcs(j))],label,classes,{sprintf('PC %d',pcs(i)),sprintf('PC %d',pcs(j))},opt-1);     
        end
    end
end
        