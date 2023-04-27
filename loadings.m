
function fig_h =loadings(model,opt,tit,label,classes,blur)

% Compute and plot loadings.
%
% fig_h =loadings_pca(model) % minimum call
% fig_h =loadings_pca(model,opt,tit,label,classes,blur) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores. 
%
% opt: [1X1]
%       0: scatter plot of pairs of PCs (by default)
%       1: bar plot of each single PC
%
% tit: (str) title for the plots. Empty by default;
%
% label: [Mx1] name of the variables (numbers are used by default)
%
% classes: [Mx1] groups for different visualization (a single group 
%   by default)
%
% blur: [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
%
% OUTPUTS:
%
% fig_h: set of figure handles
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,8);
% [~,~,model] = pca_pp(X,1:2);
%
% P = loadings(model);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/Apr/2023
%
% Copyright (C) 2023  University of Granada, Granada
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(model.scores, 1);
M = size(model.loads, 1);

if nargin < 2 || isempty(opt), opt = 0; end; 
if nargin < 3, tit = ''; end 
if nargin < 4 || isempty(label), label = 1:M; end
if nargin < 5 || isempty(classes), classes = ones(M,1); end
if nargin < 6 || isempty(blur),    blur    = 1;       end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
assert (length(opt)==1, 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 4th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: 5th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~=0 & opt~=1)), 'Value Error: 2nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

P = model.loads;


%% Show results

fig_h = [];
if length(model.lvs) == 1 || opt
    for i=1:length(model.lvs)
        fig_h = [fig_h plot_vec(P(:,i), label, classes, {'',sprintf('Loadings PC %d',model.lvs(i))})];
        title(tit);
    end
else
    for i=1:length(model.lvs)-1
        for j=i+1:length(model.lvs)
            fig_h = [fig_h plot_scatter([P(:,i),P(:,j)], label, classes, {sprintf('Loadings PC %d',model.lvs(i)),sprintf('Loadings PC %d',model.lvs(j))}',[],[],[],[],blur)];
            title(tit);
        end
    end
end
        