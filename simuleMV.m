function X = simuleMV(obs,vars,corr,lcorr)

% Simulation of MV data with ADICOV (Approximation of a DIstribution for a 
% given COVariance, Chemometrics and Intelligent Laboratory Systems 105(2), 
% 2011, pp. 171-180.
%
% X = simuleMV(obs,vars) % minimum call
% X = simuleMV(obs,vars,corr,lcorr)% complete call
%
%
% INPUTS:
%
% obs: [1x1] number of observations (rows) in the output.
%
% vars: [1x1] number of variables (columns) in the output.
%
% corr: [vars x vars] pre-specified correlation. NANs in this input are 
%   replaced by random values. All random ( nan(vars) ) used by default.
%
% lcorr: [1x1] level of correlation among variables, in [0,100] (10 by 
%   default). The correlation used is equal to: corr.^(2*lcorr+1).
%
%
% OUTPUTS:
%
% X: [obsxvars] data matrix generated.
%
%
% EXAMPLE OF USE: To obtain a matrix 100x10 with random covariance matrix, 
% use the following call:
%
% X = simuleMV(100,10);
% meda_pca(X); % visualization (auto-scaled data)
%
%
% EXAMPLE OF USE: Matrix 100x10 with the first five variables with 
% correlation 1:
%
% corr = nan(10);
% corr(1:5,1:5) = 1;
% X = simuleMV(100,10,corr);
% meda_pca(X); % visualization (auto-scaled data)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 25/May/16.
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 3 || isempty(corr), corr = nan(vars); end;
ind = find(~isnan(corr));
if nargin < 4 || isempty(lcorr), lcorr = 20; end;

% Validate dimensions of input data
assert (isequal(size(obs), [1 1]), 'Dimension Error: 1st argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(vars), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(corr), [vars vars]), 'Dimension Error: 3rd argument must be vars-by-vars. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lcorr), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (obs>0, 'Value Error: 1st argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(obs), obs), 'Value Error: 1st argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
assert (vars>0, 'Value Error: 2nd argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(vars), vars), 'Value Error: 2nd argument must be an integer. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(ind), 
    assert (isempty(find(corr(ind)<-1)), 'Value Error: 3rd argument elements must be equal or above -1. Type ''help %s'' for more info.', routine(1).name);
    assert (isempty(find(corr(ind)>1)), 'Value Error: 3rd argument elements must be equal or below 1. Type ''help %s'' for more info.', routine(1).name);
end
assert (lcorr>=0, 'Value Error: 4th argument must be equal or above 0. Type ''help %s'' for more info.', routine(1).name);
assert (lcorr<=100, 'Value Error: 4th argument must be equal or below 100. Type ''help %s'' for more info.', routine(1).name);


%% Main code

COV = randn(vars).^(2*lcorr+1);
COV = COV./max(max(abs(COV)));
COV(ind) = corr(ind);
COV = COV - diag(diag(COV)) + eye(vars);
X = real(ADICOV(COV,randn(obs,vars),vars));
