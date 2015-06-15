function lim = hot_lim(pc,nob,p_value,phase)

% Control limit for D statistic.
%
% lim = hot_lim(pc,nob,p_value)       % minimum call
% lim = hot_lim(pc,nob,p_value,phase)       % complete call
%
% INPUTS:
%
% pc: (1x1) Number of PCs.
%
% nob: (1x1) Number of observations.
%
% p_value: (1x1) p-value of the test.
%
% phase: (1x1) SPC phase:
%   - 1: Phase I
%   - 2: Phase II (by default)
%
%
% OUTPUTS:
%
% lim: (1x1) control limit at a 1-p_value confidence level.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 09/Jul/15.
%
% Copyright (C) 2015  University of Granada, Granada
% Copyright (C) 2015  Jose Camacho Paez
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

if nargin < 4, phase = 2; end;
if (phase<1||phase>2), error('Incorrect value of phase.'); end;

% Computation

if phase ==2,
    lim = (pc*(nob*nob-1)/(nob*(nob-pc)))*finv(1-p_value,pc,nob-pc);
else
    lim = (nob-1)^2/nob*icdf('Beta',1-p_value,pc/2,(nob-pc-1)/2);
end


