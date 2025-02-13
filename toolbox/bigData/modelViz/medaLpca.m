function [medamap,ind,ord,Lmodel] = medaLpca(Lmodel,varargin)

% Missing data methods for exploratory data analysis in PCA. The original
% paper is Chemometrics and Intelligent Laboratory Systems 103(1), 2010, pp.
% 8-18. This algorithm follows the suggested computation by Arteaga, which
% makes use of the covariance matrices.
%
% [medamap,ind,ord,Lmodel] = medaLpca(Lmodel) % minimum call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.lvs: [1x1] number of PCs.
%
% Optional INPUTS (parameters):
%
% 'Threshold': [1x1] threshold (0,1] for discretization and discarding (0.1 by default) 
%
% 'Option': (str or num) options for data plotting: binary code of the form 'abc' for:
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
% 'VarIndex': [Sx1] Subset of variables to plot (1:M by default)
%
%
% OUTPUTS:
%
% medamap: [MxM] non-seriated MEDA matrix.
%
% ind: [M2x1] indices of seriated variables over the input threshold (M2 < M).
%
% ord: [M2x1] order of seriated variables.
%
% Lmodel: (struct Lmodel) model after integrity checking.
%
%
% EXAMPLE OF USE:
%
% X = simuleMV(20,10,'LevelCorr',8);
% Lmodel = iniLmodel(X);
% Lmodel.lvs = 1:3;
% map = medaLpca(Lmodel,'Threshold', 0.3,'Option','111');
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 22/Nov/2024
%
% Copyright (C) 2024  University of Granada, Granada
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
[ok, Lmodel] = checkLmodel(Lmodel);
M = size(Lmodel.centr,2);
% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Threshold',0.1);
addParameter(p,'Option','100');
addParameter(p,'VarIndex',1:M);  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
thres = p.Results.Threshold;
opt = p.Results.Option;
vars = p.Results.VarIndex;

% Convert row arrays to column arrays
if size(vars,1)  == 1, vars = vars'; end;

% opt backward compatibility
if isstruct(opt) 
    opt = opt.plot + 10*opt.seriated + 100*opt.discard;
end

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end

% Validate dimensions of input data
assert (isequal(size(thres), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 3rd argument must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(size(vars) > [M 1])), 'Dimension Error: 4th argument must be at most M-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (thres>0 && thres<=1, 'Value Error: 2nd argument must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(vars<=0)) && isequal(fix(vars), vars) && isempty(find(vars>M)), 'Value Error: 4th argument must contain positive integers below or equal to M. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel.XX = Lmodel.XX(vars,vars);
s = size(Lmodel.XX);

Lmodel = Lpca(Lmodel);
P = Lmodel.loads;

medamap = meda(Lmodel.XX,P);

if nargout > 1 || opt(3) == '1'
    Dmap = diag(medamap);
    ind = find(Dmap > thres);
end

if nargout > 2 || opt(2) == '1'
    [map, ord] = seriation(medamap(ind,ind));
end
   

%% Show results

if opt(1) == '1'
    
    map = medamap;

    if opt(3) == '1'
        ind2 = ind;
    else
        ind2 = 1:length(vars);
    end
    
    map = map(ind2,ind2);
    label = Lmodel.varl(ind2);
    
    if opt(2) == '1'
        ord2 = ord;
    else
        ord2 = 1:length(vars);
    end
    
    map = map(ord2,ord2);
    label = label(ord2);
    
    plotMap(map,'VarsLabel',label);
    
end  

        