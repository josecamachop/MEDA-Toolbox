function Lmodel = Lpls(Lmodel)
    
% PLS for large data. 
%
% Lmodel = Lpls(Lmodel) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PLS
%   model:
%       Lmodel.XX: [MxM] X-block cross-product matrix.
%       Lmodel.XY: [MxO] cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lv: [1x1] number of LVs A.
%
%
% OUTPUTS:
%
% Lmodel: (struct Lmodel) output model with coefficients
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Lmodel = iniLmodel(X,Y);
% Lmodel.lvs = 0:10;
% Lmodel = Lpls(Lmodel)
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Nov/2024
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

% Validate values of input data
[ok,Lmodel] = checkLmodel(Lmodel);


%% Main code

model = kernelpls(Lmodel.XX,Lmodel.XY,'LVs',1:max(Lmodel.lvs));
R = model.altweights;
beta = model.beta;
        
[V2,d2] = eig(R'*Lmodel.XX*R);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
        
model.sdT = real(sqrt(dd/(Lmodel.N-1)));

fnames = fieldnames(model);
for n = 1:length(fnames)
    Lmodel = setfield(Lmodel,fnames{n},getfield(model,fnames{n}));
end
