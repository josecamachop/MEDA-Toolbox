
function [P,W,Q] = loadings_pls(x,y,lvs,prepx,prepy,opt,label,classes)

% Compute and plot loadings in PLS
%
% P = loadings_pls(x,y) % minimum call
% [P,W,Q] = loadings_pls(x,y,lvs,prepx,prepy,opt,label,classes) % complete call
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 1:rank(x)
%
% prepx: [1x1] preprocesing of the x-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% prepy: [1x1] preprocesing of the y-block
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)   
%
% opt: (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot loadings
%       b:
%           0: scatter plot of pairs of LVs
%           1: bar plot of each single LV
%       c:
%           0: plot weights
%           1: plot X-block loadings
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
%
% label: [Mx1] name of the variables (numbers are used by default)
%
% classes: [Mx1] groups for different visualization (a single group 
%   by default)
%
%
% OUTPUTS:
%
% P: [MxA] X-block loadings
%
% W: [MxA] X-block weights
%
% Q: [OxA] Y-block loadings
%
%
% EXAMPLE OF USE: Random loadings: bar and scatter plot of loadings
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% loadings_pls(X,Y,1);
% [P,W,Q] = loadings_pls(X,Y,1:3);
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);
if nargin < 3 || isempty(lvs), lvs = 1:rank(x); end;
if nargin < 4 || isempty(prepx), prepx = 2; end;
if nargin < 5 || isempty(prepy), prepy = 2; end;
if nargin < 6 || isempty(opt), opt = '100'; end; 
if nargin < 7 || isempty(label), label = [1:M]; end
if nargin < 8 || isempty(classes), classes = ones(M,1); end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert int arrays to str
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 3rd argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 6th argument must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 7th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [M 1]), 'Dimension Error: 8th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 6th argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,prepx);
ycs = preprocess2D(y,prepy);

[beta,W,P,Q] = kernel_pls(xcs'*xcs,xcs'*ycs,lvs);


%% Show results

if opt(1) == '1',
    
    if opt(3) == '0',
        Pt = W;
        text = 'Weights';
    else
        Pt = P;
        text = 'X-block loadings';
    end
    
    if length(lvs) == 1 || opt(2) == '1',
        for i=1:length(lvs),
            plot_vec(Pt(:,i), label, classes, {'',sprintf('%s LV %d',text,lvs(i))});
        end
    else
        for i=1:length(lvs)-1,
            for j=i+1:length(lvs),
                plot_scatter([Pt(:,i),Pt(:,j)], label, classes, {sprintf('%s LV %d',text,lvs(i)),sprintf('%s LV %d',text,lvs(j))}');
            end      
        end
    end
end
        