function [P,sdT] = Lpca(Lmodel)
  
% PCA for large data. 
%
% [P,S] = Lpca(Lmodel) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PCA..
%
%
% OUTPUTS:
%
% P: (MxA) matrix of loadings in the PCA model.
%
% sdT: (1xA) standard deviations of the scores.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,8);
% Lmodel = Lmodel_ini;
% Lmodel.XX = X'*X;
% Lmodel.lvs = 0:10;
% [P,sdT] = Lpca(Lmodel);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 17/Oct/2016
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
    
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

check_Lmodel(Lmodel);


%% Main code

[P,d2] = eig(Lmodel.XX);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
P = P(:,inddd(1:max(Lmodel.lvs)));

sdT = real(sqrt(dd(1:max(Lmodel.lvs))/(Lmodel.N-1)));