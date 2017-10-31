
function [meda_map,ind,ord] = meda_pls(x,y,lvs,prepx,prepy,thres,opt,label,vars)

% Missing data methods for exploratory data analysis in PCA. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. 
%
% meda_map = meda_pls(x,y) % minimum call
% [meda_map,ind,ord] = meda_pls(x,y,lvs,prepx,prepy,thres,opt,label,vars) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
% y: [NxO] billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. pcs = 1:2 selects the
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
% thres: [1x1] threshold (0,1] for discretization and discarding (0.1 by default)
%
% opt: (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot MEDA matrix
%       b:
%           0: no seriated
%           1: seriated
%       c:
%           0: no discard
%           1: discard variables with diagonal below threshold 
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
%
% label: [Mx1] name of the variables (numbers are used by default)
%
% vars: [Sx1] Subset of variables to plot (1:M by default)
%
%
% OUTPUTS:
%
% meda_map: [MxM] non-seriated MEDA matrix.
%
% ind: [M2x1] indices of seriated variables over the input threshold (M2 < M).
%
% ord: [M2x1] order of seriated variables.
%
%
%
% EXAMPLE OF USE: Seriation and discarding uninformative variables
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% lvs = 1:3;
% map = meda_pls(X,Y,lvs,2,2,0.3,'111');
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/May/17
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(x, 1);
M = size(x, 2);
O = size(y, 2);
if nargin < 3 || isempty(lvs), lvs = 1:rank(x); end;
if nargin < 4 || isempty(prepx), prepx = 2; end;
if nargin < 5 || isempty(prepy), prepy = 2; end;
if nargin < 6 || isempty(thres), thres = 0.1; end; 
if nargin < 7 || isempty(opt), opt = '100'; end; 
if nargin < 8 || isempty(label), label = 1:M; end
if nargin < 9 || isempty(vars), vars = 1:M; end;

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;
if size(vars,1)  == 1, vars = vars'; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
A = length(lvs);
if isstruct(opt) % opt backward compatibility
    opt = opt.plot + 10*opt.seriated + 100*opt.discard;
end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

% Validate dimensions of input data
assert (A>0, 'Dimension Error: 3rd argument with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepx), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prepy), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(thres), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 7th argument must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 8th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(size(vars) > [M 1])), 'Dimension Error: 9th argument must be at most M-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert (thres>0 && thres<=1, 'Value Error: 6th argument must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 7th argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(vars<=0)) && isequal(fix(vars), vars) && isempty(find(vars>M)), 'Value Error: 9th argument must contain positive integers below or equal to M. Type ''help %s'' for more info.', routine(1).name);


%% Main code

x = x(:,vars);
label = label(vars);

x2 = preprocess2D(x,prepx);
y2 = preprocess2D(y,prepy);

[beta,W,P,Q,R] = kernel_pls(x2'*x2,x2'*y2,lvs);

meda_map = meda(x2'*x2,R,P);

if opt(3) == '1'
    Dmap = diag(meda_map);
    ind = find(Dmap > thres);
else
    ind = 1:length(vars);
end

if opt(2) == '1',
    [map, ord] = seriation(meda_map(ind,ind));
else
    ord = 1:length(vars);
end
    

%% Show results

if opt(1) == '1',
    
    map = meda_map;

    if opt(3) == '1',
        ind2 = ind;
    else
        ind2 = 1:length(vars);
    end
    
    map = map(ind2,ind2);
    label = label(ind2);
    
    if opt(2) == '1',
        ord2 = ord;
    else
        ord2 = 1:length(vars);
    end
    
    map = map(ord2,ord2);
    label = label(ord2);
    
    plot_map(map,label);
    
end  