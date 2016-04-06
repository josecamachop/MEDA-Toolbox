function [beta,W,P,Q,R,sdT] = Lpls(Lmodel)
    
% PLS for large data. 
%
% [beta,W,P,Q,R,S] = Lpls(Lmodel) % complete call
%
%
% INPUTS:
%
% Lmodel: (struct Lmodel) model with the information to compute the PLS
%   model:
%       Lmodel.XX: (MxM) X-block cross-product matrix.
%       Lmodel.XY: (MxL) cross-product matrix between the x-block and the
%           y-block.
%       Lmodel.lv: (1x1) number of LVs A.
%
%
% OUTPUTS:
%
% beta: (Mx1) matrix with regression coefficients.
%
% W: (MxA) matrix of weights in the PLS model.
%
% P: (MxA) matrix of loadings of the x-block in the PLS model.
%
% Q: (LxA) matrix of loadings of the y-block in the PLS model.
%
% R: (LxA) equals to W·inv(P'·W).
%
% sdT: (1xA) standard deviations of the scores.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 05/Apr/16.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

[beta,W,P,Q,R] = kernel_pls(Lmodel.XX,Lmodel.XY,1:Lmodel.lv);
        
[V2,d2] = eig(R'*Lmodel.XX*R);
dd = diag(d2);
[dd,inddd]=sort(dd,'descend');
        
sdT = real(sqrt(dd/(Lmodel.N-1)));
