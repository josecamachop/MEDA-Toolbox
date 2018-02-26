function testcs = preprocess2Dapp(test,average,scale,weight)

% Preprocess 2-way data with previously computed average and scale: 
%   testcs = weight.*(test - average)./scale.
%
% testcs = preprocess2Dapp(test,average)         % minimum call
% testcs = preprocess2Dapp(test,average,scale,weight)     % complete call
%
% INPUTS:
%
% test: [NxM] billinear data set
%
% average: [1xM] average to subtract from test.
%
% scale: [1xM] scale to divide test. A vector or ones is used by default.
%
% weight: [1xM] weight applied after preprocessing scaling. Set to 1 by
%   default.
%
% OUTPUTS:
%
% test: [NxM] preprocessed data.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 19/Feb/2018
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez, 
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
N = size(test, 1);
M = size(test, 2);
if nargin < 3 || isempty(scale), scale = ones(1,M); end;
if nargin < 4 || isempty(weight), weight = ones(1,M); end;

% Convert column arrays to row arrays
if size(average,2) == 1, average = average'; end;
if size(scale,2) == 1, scale = scale'; end;
if size(weight,2) == 1, weight = weight'; end;

% Validate dimensions of input data
assert (isequal(size(average), [1 M]), 'Dimension Error: 2nd argument must be 1-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(scale), [1 M]), 'Dimension Error: 3rd argument must be 1-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(weight), [1 M]), 'Dimension Error: 4th argument must be 1-by-M. Type ''help %s'' for more info.', routine(1).name); 


%% Main code

testcs = (ones(N,1)*weight).*(test - ones(N,1)*average)./(ones(N,1)*scale);
