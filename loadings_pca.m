
function P = loadings_pca(x,pcs,prep,opt,label,classes)

% Compute and plot loadings in PCA
%
% P = loadings_pca(x) % minimum call
% P = loadings_pca(x,pcs,prep,opt,label,classes) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default) 
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
% EXAMPLE OF USE: Scatter plot of random scores
%
% X = simuleMV(20,10,8);
% P = loadings_pca(X,1:3);
%
%
% EXAMPLE OF USE: Line plot of random scores
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% P = loadings_pca(X,1:3,[],11);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 19/Apr/2016
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 1:rank(x); end;
if nargin < 3 || isempty(prep), prep = 2; end;
if nargin < 4 || isempty(opt), opt = '10'; end; 
if nargin < 5 || isempty(label), label = [1:M]; end
if nargin < 6 || isempty(classes), classes = ones(M,1); end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 2nd argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==2, 'Dimension Error: 4th argument must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 5th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: 6th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 4th argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,prep);
P = pca_pp(xcs,pcs);


%% Show results

if opt(1) == '1',
    
    if length(pcs) == 1 || opt(2) == '1',
        for i=1:length(pcs),
                plot_vec(P(:,i), label, classes, {'',sprintf('Loadings PC %d',pcs(i))});
        end
    else
        for i=1:length(pcs)-1,
            for j=i+1:length(pcs),
                plot_scatter([P(:,i),P(:,j)], label, classes, {sprintf('Loadings PC %d',pcs(i)),sprintf('Loadings PC %d',pcs(j))}');
            end      
        end
    end
end
        