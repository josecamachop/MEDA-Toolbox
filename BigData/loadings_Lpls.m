
function [P,Q,R] = loadings_Lpls(Lmodel,lvs,opt,label,classes)

% Compute and plot loadings in PLS for large data.
%
% loadings_Lpls(Lmodel,lvs) % minimum call
% loadings_Lpls(Lmodel,lvs,opt,label,classes) % complete call
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% opt: (1x1) options for data plotting.
%       0: no plots.
%       1: loading plot of P (default)
%       2: loading plot of P with empty marks
%       3: loading plot of W*inv(P'*W)
%       4: loading plot of W*inv(P'*W) with empty marks
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
% P: (MxA) x-block scores.
% Q: (OxA) y-block scores.
% R: (MxA) W*inv(P'*W).
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 06/May/13.
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

Lmodel.lv = max(lvs);
[beta,W,P,Q,R] = Lpls(Lmodel);

if opt,
    for i=1:length(lvs)-1,
        for j=i+1:length(lvs),
            if opt ==1 || opt ==2,
                plot_scatter([P(:,lvs(i)),P(:,lvs(j))],label,classes,{sprintf('LV %d',lvs(i)),sprintf('LV %d',lvs(j))},opt-1);   
            elseif opt == 3 || opt ==4,
                plot_scatter([R(:,lvs(i)),R(:,lvs(j))],label,classes,{sprintf('LV %d',lvs(i)),sprintf('LV %d',lvs(j))},opt-3);                
            end    
        end
    end
end
        