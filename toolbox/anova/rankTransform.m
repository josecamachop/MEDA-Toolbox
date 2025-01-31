function Xtrans = rank_transform(X)

% Transformation on rankings in the columns.  
%
% Xtrans = rank_transform(X)     % minimum call
%
%
% INPUTS:
%
% X: [NxM] data set 
%
%
% OUTPUTS:
%
% Xtrans: [NxM] output data set 
%
%
% EXAMPLE OF USE: Random
%
% x = randn(5)
% rank_transform(x)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 31/Jan/2025
%
% Copyright (C) 2025  University of Granada, Granada
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

Xtrans = X;
for i=1:size(X,2)
    uX = unique(X(:,i));
    uX = sort(uX);
    jj=0;
    for j=1:length(uX)
        ind2 = find(X(:,i)==uX(j));
        Xtrans(ind2,i) = jj+mean(1:length(ind2));
        jj = jj+length(ind2);
    end
end