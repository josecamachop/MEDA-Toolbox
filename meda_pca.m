
function [meda_map,meda_dis,ord] = meda_pca(x,pcs,prep,thres,opt,label,vars)

% Missing data methods for exploratory data analysis in PCA. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. 
%
% meda_map = meda_pca(x) % minimum call
% [meda_map,meda_dis,ord] = meda_pca(x,pcs,prep,thres,opt,label,vars) % complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set 
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 1:rank(xcs)
%
% prep: [1x1] preprocesing of the data
%       0: no preprocessing
%       1: mean centering
%       2: autoscaling (default)  
%
% thres: [1x1] threshold (0,1] for discretization and discarding (0.1 by default)
%
% opt: (str) binary code of the form 'cba' for:
%       a:
%           0: no plots
%           1: plot MEDA matrix (default)
%       b:
%           0: no seriated
%           1: seriated (default)
%       c:
%           0: no discard
%           1: discard 0 variance variables (default)
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
% meda_dis: [MxM] discretized MEDA matrix.
%
% ord: [Mx1] order of seriated variables.
%
%
% EXAMPLE OF USE: Seriation and discarding uninformative variables
%
% X = real(ADICOV(randn(10,10).^19,randn(100,10),10));
% pcs = 1:3;
% map = meda_pca(X,pcs,2,0.3,'111');
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 29/Mar/16
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(x, 1);
M = size(x, 2);
if nargin < 2 || isempty(pcs), pcs = 1:rank(x); end;
if nargin < 3 || isempty(prep), prep = 2; end;
if nargin < 4 || isempty(thres), thres = 0.1; end; 
if nargin < 5 || isempty(opt), opt = '111'; end; 
if nargin < 6 || isempty(label), label = 1:M; end
if nargin < 7 || isempty(vars), vars = 1:M; end;

% Convert row arrays to column arrays
if size(label,1) == 1, label = label'; end;
if size(vars,1)  == 1, vars = vars'; end;

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
A = length(pcs);
if isstruct(opt) % opt backward compatibility
    opt = opt.plot + 10*opt.seriated + 100*opt.discard;
end

% Convert int arrays to str
if isnumeric(opt), opt=fliplr(num2str(opt,'%.3d')); end

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: 2nd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: 3rd argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(thres), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
assert (ischar(opt), 'Dimension Error: 5th argument must be a string. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(label), [M 1]), 'Dimension Error: 6th argument must be M-by-1. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(size(vars) > [M 1])), 'Dimension Error: 7th argument must be at most M-by-1. Type ''help %s'' for more info.', routine.name);

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: 2nd argument must contain positive integers. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(pcs>rank(x))), 'Value Error: 2nd argument must contain values below the rank of the data. Type ''help %s'' for more info.', routine.name);
assert (thres>0 && thres<=1, 'Value Error: 4th argument must be in (0,1]. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(vars<=0)) && isequal(fix(vars), vars) && isempty(find(vars>M)), 'Value Error: 7th argument must contain positive integers below or equal to M. Type ''help %s'' for more info.', routine.name);


%% Main code

x = x(:,vars);
label = label(vars);

x2 = preprocess2D(x,prep);

P = pca_pp(x2,pcs);
        
meda_map = meda(x2'*x2,P,P);

if nargout > 1
    meda_dis = meda_map; % discretize
    ind = find(meda_dis(:)<thres);
    meda_dis(ind) = 0;
end

if nargout > 2 || opt(2),
    [map1, ord] = seriation(meda_map);
end

    
%% Show results

if opt(1) == '1',
    
    map = meda_map;
 
    if opt(2) == '1',
        ord2 = ord;
    else
        ord2 = 1:M;
    end
    
    map = map(ord2,ord2);
    label = label(ord2);
    
    if opt(3) == '1',
        Dmap = diag(map);
        ind = find(Dmap > thres);
    else
        ind = 1:M;
    end
    
    map = map(ind,ind);
    label = label(ind);
    
    plot_map(map,label);
    
end

        