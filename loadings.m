
function P = loadings(model,opt,tit,label,classes)

% Compute and plot loadings.
%
% P = loadings_pca(model) % minimum call
% P = loadings_pca(model,opt,tit,label,classes) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   av: [1xM] centering parameters. 
%   sc: [1xM] scaling parameters. 
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores. 
%
% opt: (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: no plots
%           1: plot loadings
%       b:
%           0: scatter plot of pairs of PCs
%           1: bar plot of each single PC
%   By deafult, opt = '10'. If less than 2 digits are specified, least 
%   significant digit is set to 0, i.e. opt = 1 means a=1 and b=0. If a=0, 
%   then b is ignored.
%
% tit: (str) title for the plots. Empty by default;
%
% label: [Mx1] name of the variables (numbers are used by default)
%
% classes: [Mx1] groups for different visualization (a single group 
%   by default)
%
%
% OUTPUTS:
%
% P: [MxA] loadings
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,8);
%
% model.lvs = 1:3;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% [model.loads,model.scores] = pca_pp(Xcs,model.lvs);
%
% P = loadings(model);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 25/Apr/2018
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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
[M,A] = size(model.loads);
if nargin < 2 || isempty(opt), opt = '10'; end; 
if nargin < 3, tit = ''; end 
if nargin < 4 || isempty(label), label = [1:M]; end
if nargin < 5 || isempty(classes), classes = ones(M,1); end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;


% Validate dimensions of input data
assert (ischar(opt) && length(opt)==2, 'Dimension Error: 2nd argument must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 4th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: 5th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 2nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

P = model.loads;


%% Show results

if opt(1) == '1',
    
    if length(model.lvs) == 1 || opt(2) == '1',
        for i=1:length(model.lvs),
        	plot_vec(P(:,i), label, classes, {'',sprintf('Loadings PC %d',model.lvs(i))});
            title(tit);
        end
    else
        for i=1:length(model.lvs)-1,
            for j=i+1:length(model.lvs),
                plot_scatter([P(:,i),P(:,j)], label, classes, {sprintf('Loadings PC %d',model.lvs(i)),sprintf('Loadings PC %d',model.lvs(j))}');
                title(tit);
            end      
        end
    end
end
        