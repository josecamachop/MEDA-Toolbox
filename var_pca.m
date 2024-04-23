
function [x_var,cumpress] = var_pca(x,varargin)

% Variability captured in terms of the number of PCs. It includes the ckf
% algorithm.
%
% var_pca(x,pcs) % minimum call
% var_pca(x,'Pcs',pcs,'Preprocessing',prep,'Option',opt) %complete call
%
%
% INPUTS:
%
% x: [NxM] billinear data set for model fitting
%
%
% Optional INPUTS (parameters):
%
% 'Pcs': [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs). By default, pcs = 0:rank(x). The value for 0 PCs is
%   added at the begining if not specified.
%
% 'Preprocessing': [1x1] preprocesing
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)  
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: no plots
%           1: plot residual variance
%       b:
%           0: Residual Variance in X and ckf 
%           1: Residual Variance in X 
%   By deafult, opt = '10'. If less than 2 digits are specified, least 
%   significant digit is set to 0, i.e. opt = 1 means a=1 and b=0. If a=0, 
%   then b is ignored.
%
%
% OUTPUTS:
%
% x_var: [Ax1] Percentage of captured variance of X.
%
% cumpress: [Ax1] ckf curve.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% pcs = 0:10;
% x_var = var_pca(X,'Pcs',pcs);
%
%
% codified by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/2024
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
N = size(x, 1);
M = size(x, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Pcs',0:rank(x));   
addParameter(p,'Preprocessing',2);   
addParameter(p,'Option',10);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
pcs = p.Results.Pcs;
opt = p.Results.Option;
prep = p.Results.Preprocessing;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'0'); end

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end;

% Preprocessing
pcs = unique(pcs);
A = length(pcs);

% Validate dimensions of input data
assert (A>0, 'Dimension Error: parameter ''Pcs'' with non valid content. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (ischar(opt) && length(opt)==2, 'Dimension Error: parameter ''Option'' must be a string or num of 2 bits. Type ''help %s'' for more info.', routine(1).name);

% Preprocessing
pcs = unique([0 pcs]);

% Validate values of input data
assert (isempty(find(pcs<0)), 'Value Error: parameter ''Pcs'' must not contain negative values. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' contain integers. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

xcs = preprocess2D(x,'Preprocessing',prep); 

[P,T] = pca_pp(xcs,'Pcs',1:max(pcs));
pcs(find(pcs>size(P,2))) = [];

totalVx = sum(sum(xcs.^2));
x_var = ones(length(pcs),1);

for i = 1:length(pcs)
    x_var(i) = x_var(i) - sum(eig(T(:,1:pcs(i))'*T(:,1:pcs(i))))/totalVx;
end

cumpress = zeros(length(pcs),1);
if nargout>1 || opt(2) == '0'
    for i = 1:length(pcs)
         c = ckf(xcs,T(:,1:pcs(i)),P(:,1:pcs(i)),'Option',0);
         cumpress(i) = c(end);
    end
end

    
%% Show results

if opt(1) == '1'
    if opt(2) == '1'
        plot_vec(x_var,'EleLabel',pcs,'XYLabel',{'#PCs','% Residual Variance'},'Option','01');
    else
        plot_vec([x_var cumpress/cumpress(1)],'EleLabel',pcs,'XYLabel',{'#PCs','% Residual Variance'},'Option','01','VecLabel',{'X','ckf'});
        legend('show');
    end
end

        