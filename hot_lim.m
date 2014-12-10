function lim = hot_lim(pc,nob,p_value)

% Control limit for D statistic.
%
% lim = hot_lim(pc,nob,p_value)       % complete call
%
% INPUTS:
%
% pc: (1x1) Number of PCs.
%
% nob: (1x1) Number of observations.
%
% p_value: (1x1) p-value of the test.
%
%
% OUTPUTS:
%
% lim: (1x1) control limit at a 1-p_value confidence level.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 04/Jan/13.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  José Camacho Páez
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

% Parameters checking

if nargin < 3, error('Error in the number of arguments.'); end;
if pc<1, error('Incorrect content of pc.'); end;
if nob<1, error('Incorrect content of nob.'); end;
if (p_value<0||p_value>1), error('Incorrect value of p_value.'); end;

% Computation

lim = (pc*(nob*nob-1)/(nob*(nob-pc)))*finv(1-p_value,pc,nob-pc);


