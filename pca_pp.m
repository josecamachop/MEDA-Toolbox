function [p,t] = pca_pp(x,pc)

% Principal Component Analysis.
%
% [p,t] = pca_pp(x,pc)     % complete call
%
% INPUTS:
%
% x: (NxM) Two-way batch data matrix, N(observations) x M(variables)
%
% pc: number of principal components.
%
%
% OUTPUTS:
%
% p: (M x pc) matrix of loadings.
%
% t: (N x pc) matrix of scores.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es).
% version: 2.1
% last modification: 03/Jul/14.
%
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

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of prep.'); end;
dmin = min(s);
if pc>dmin, pc=dmin; end;

% Computation

[u,d,p]=svd(x,0);
t = u*d;
p = p(:,1:pc);
t = t(:,1:pc);

        



