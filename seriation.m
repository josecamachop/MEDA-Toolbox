
function [mapo,ord] = seriation(mapi)

% Seriation (ordination) of a covariance-like matrix.
%
% [mapo, ord] = seriation(mapi) % complete call
%
%
% INPUTS:
%
% mapi: [MxM] symmetric input matrix. 
%
%
% OUTPUTS:
%
% mapo: [MxM] symmetric output matrix.
%
% ord: [Mx1] seriated indices.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Xcs = preprocess2D(X,'Preprocessing',2);
% mapo = seriation(cov(Xcs));
% plot_map(mapo);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 23/Apr/2024.
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
M = size(mapi, 1);

% Validate dimensions of input data
assert (isequal(size(mapi), [M M]), 'Dimension Error: parameter ''mapi'' must be M-by-M. Type ''help %s'' for more info.', routine(1).name);


%% Main code

mapo = mapi;
for i=1:M
    fragment{i} = i;
end
mapoa = abs(mapo);
mapoa(1:(M+1):end) = -Inf;

finish = false;
iter = 0;
while ~finish
   i = find(mapoa(:)==max(mapoa(:)),1);
   ci = 1+fix((i-1)/M);
   fi = i-M*fix((i-1)/M);
   
   fci = 0;
   ffi = 0;
   for i=1:length(fragment)
       if fragment{i}(1)==ci
           fragc = fragment{i};
           fci = i;
       elseif fragment{i}(end)==ci
           fragc = fliplr(fragment{i});
           fci = i;
       elseif fragment{i}(1)==fi
           fragf = fliplr(fragment{i});
           ffi = i;
       elseif fragment{i}(end)==fi
           fragf = fragment{i};
           ffi = i;
       end
   end
   
   if fci && ffi
       if fci<ffi
            fragment = fragment([1:fci-1 fci+1:ffi-1 ffi+1:end]);
       else
            fragment = fragment([1:ffi-1 ffi+1:fci-1 fci+1:end]);
       end
       fragment{end+1} = [fragf fragc];       
   end
   
   mapoa(ci,fi) = -Inf;
   mapoa(fi,ci) = -Inf;
   
   if length(fragment)==1
       finish = true;
   end
   
   iter = iter+1;
   
   if iter > 1e4
       finish = true;
       frag2 = [];
       for i=1:length(fragment)
           frag2 = [frag2 fragment{i}];
       end
       fragment{1} = frag2;
   end
   
end

mapo = mapo(fragment{1},fragment{1});
ord = fragment{1};