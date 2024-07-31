function lim = hot_lim(npc,nob,p_value,varargin)

% Control limit for D statistic. This routine needs the Statistics and
% Machine Learning Toolbox in Matlab.
%
% lim = hot_lim(npc,nob,p_value)       % minimum call
% lim = hot_lim(npc,nob,p_value,'Phase',phase)       % complete call
%
%
% INPUTS:
%
% npc: [1x1] Number of PCs
%
% nob: [1x1] Number of observations
%
% p_value: [1x1] p-value of the test, in (0,1]
% 
% Optional INPUTS (parameter):
%
% 'Phase': [1x1] SPC phase:
%   - 1: Phase I
%   - 2: Phase II (by default)
%
%
% OUTPUTS:
%
% lim: [1x1] control limit at a 1-p_value confidence level.
%
%
% EXAMPLE OF USE: For 2 PCs, 100 observations, the 99% confidence limit in
% phase II follows:
%
% lim = hot_lim(2,100,0.01,'Phase',2)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/2024
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Phase',2);  
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
phase = p.Results.Phase;


% Convert char values in numerical 
if ischar(phase), 
    if isequal(phase,'I'), 
        phase = 1; 
    elseif isequal(phase,'II'), 
        phase = 2;
    end;
end

% Validate dimensions of input data
assert (isequal(size(npc), [1 1]), 'Dimension Error: parameter ''npc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nob), [1 1]), 'Dimension Error: parameter ''nob'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(p_value), [1 1]), 'Dimension Error: parameter ''p_vlaue'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(phase), [1 1]), 'Dimension Error: parameter ''Phase'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (npc>0 && isequal(fix(npc), npc), 'Value Error: parameter ''npc'' must be a positive integer. Type ''help %s'' for more info.', routine(1).name);
assert (nob>0 && isequal(fix(nob), nob), 'Value Error: parameter ''nob'' must be a positive integer. Type ''help %s'' for more info.', routine(1).name);
assert (p_value>=0 && p_value<1, 'Value Error: parameter ''p_vlaue'' must be in (0,1]. Type ''help %s'' for more info.', routine(1).name);
assert (phase==1 || phase==2, 'Value Error: parameter ''Phase'' must be 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

if phase ==2
    b=nob-npc;
    if b < 0.1, %  b must be positive, and does not make a difference below 0.1
        b = 0.1;
    end
    lim = (npc*(nob*nob-1)/(nob*(b)))*finv(1-p_value,npc,b);
else
    b=(nob-npc-1)/2;
    if b < 0.1 %  b must be positive, and does not make a difference below 0.1
        b = 0.1;
    end
    lim = (nob-1)^2/nob*betainv(1-p_value,npc/2,b);
end


