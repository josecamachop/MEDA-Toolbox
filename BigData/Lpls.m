function [beta,W,P,Q,R,sdT,Lmodel] = Lpls(Lmodel)
    
% PLS for large data. 
%
% [beta,W,P,Q,R,sdT,Lmodel] = Lpls(Lmodel) % complete call
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
% beta: [Mx1] matrix with regression coefficients.
%
% W: [MxA] matrix of weights in the PLS model.
%
% P: [MxA] matrix of loadings of the x-block in the PLS model.
%
% Q: [OxA] matrix of loadings of the y-block in the PLS model.
%
% R: [OxA] equals to W·inv(P'·W).
%
% sdT: [1xA] standard deviations of the scores.
%
% Lmodel: (struct Lmodel) model after integrity checking.
%
%
% EXAMPLE OF USE: Random data
%
% X = simuleMV(20,10,8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Lmodel = Lmodel_ini(X,Y);
% Lmodel.lvs = 0:10;
% [beta,W,P,Q,R,sdT] = Lpls(Lmodel)
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 21/May/17.
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
[ok,Lmodel] = check_Lmodel(Lmodel);


%% Main code

[beta,W,P,Q,R] = kernel_pls(Lmodel.XX,Lmodel.XY,1:max(Lmodel.lvs));
W = W(:,Lmodel.lvs);
P = P(:,Lmodel.lvs);
Q = Q(:,Lmodel.lvs);
R = R(:,Lmodel.lvs);
beta=R*Q';
        
[V2,d2] = eig(R'*Lmodel.XX*R);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
        
sdT = real(sqrt(dd/(Lmodel.N-1)));
