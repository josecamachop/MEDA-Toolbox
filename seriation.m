
function [mapo,ord] = seriation(mapi)

% Seriation (ordination) of a covariance-like matrix.
%
% mapo = seriation(mapi) % complete call
%
%
% INPUTS:
%
% mapi: (MxM) symmetric input matrix. 
%
% OUTPUTS:
%
% mapo: (MxM) symmetric output matrix.
%
% ord: (1xM) seriated indices.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es).
% version: 0.0
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

%% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
s = size(mapi);
if s(1)~=s(2), error('Error in the dimension of the arguments.'); end;


%% Main code

mapo = mapi;
for i=1:s(1),
    fragment{i} = i;
end
mapoa = abs(mapo);
mapoa(1:(s(1)+1):end) = -Inf;

finish = false;
while ~finish,
   i = find(mapoa(:)==max(mapoa(:)),1);
   ci = 1+fix((i-1)/s(1));
   fi = i-s(1)*fix((i-1)/s(1));
   
   fci = 0;
   ffi = 0;
   for i=1:length(fragment),
       if fragment{i}(1)==ci,
           fragc = fragment{i};
           fci = i;
       elseif fragment{i}(end)==ci,
           fragc = fliplr(fragment{i});
           fci = i;
       elseif fragment{i}(1)==fi,
           fragf = fliplr(fragment{i});
           ffi = i;
       elseif fragment{i}(end)==fi,
           fragf = fragment{i};
           ffi = i;
       end
   end
   
   if fci & ffi,
       if fci<ffi,
            fragment = fragment([1:fci-1 fci+1:ffi-1 ffi+1:end]);
       else
            fragment = fragment([1:ffi-1 ffi+1:fci-1 fci+1:end]);
       end
       fragment{end+1} = [fragf fragc];       
   end
   
   mapoa(ci,fi) = -Inf;
   mapoa(fi,ci) = -Inf;
   
   if length(fragment)==1,
       finish = true;
   end
   
end

mapo = mapo(fragment{1},fragment{1});
ord = fragment{1};