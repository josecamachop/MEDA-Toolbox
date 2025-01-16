function map = okabeIto(nElems)

% Colormap for color blindness
%
% okabeIto(nElems) % minimum call
%
%
% See also: plotVec, plotScatter
%
%
% INPUTS:
%
% nElems: [1x1] number of color items. 
%
%
% OUTPUTS:
%
% map: (nElems x 3) colormap.
%
%
% EXAMPLE OF USE: Compare okabeIto with hsv
%
% okabeIto(5)
% hsv(5)
%
%
% coded by: Michael Sorochan Armstrong ()
%           Jose Camacho (josecamacho@ugr.es)
% last modification: 13/Jan/2025
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

% Choosing the color
okabeIto = [0.1,0.1,0.1;
    0.902,0.624,0;
    0.337,0.706,0.914;
    0,0.620,0.451;
    0.941,0.894,0.259;
    0,0.447,0.698;
    0.835,0.369,0;
    0.8,0.475,0.655];

if nElems > size(okabeIto,1) % if there are a lot of input colors, repmat to extend the palette.
    map = repmat(okabeIto, ceil(nElems/size(okabeIto, 1)), 1);
else
    %map = okabeIto(1:round(size(okabeIto,1)/nElems):end,:);
    map = okabeIto(1:nElems,:);
end

end