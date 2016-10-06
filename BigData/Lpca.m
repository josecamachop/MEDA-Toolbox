function [P,sdT] = Lpca(Lmodel)
  
% PCA for large data. 
%
% [P,S] = Lpca(Lmodel) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.lv: (1x1) number of PCs A.
%
%
% OUTPUTS:
%
% P: (MxA) matrix of loadings in the PCA model.
%
% sdT: (1xA) standard deviations of the scores.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Apr/16.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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
    
%

[P,d2] = eig(Lmodel.XX);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
P = P(:,inddd(1:Lmodel.lv));

sdT = real(sqrt(dd(1:Lmodel.lv)/(Lmodel.N-1)));