
function omeda_vec = omeda(test,dummy,R,Q)

% Observation-based Missing data methods for Exploratory Data Analysis 
% (oMEDA). The original paper is Journal of Chemometrics, 2011, 25 
% (11): 592-600. This algorithm follows the direct computation for
% Known Data Regression (KDR) missing data imputation.
%
% omeda_vec = omeda(test,dummy,R) % minimum call
% omeda_vec = omeda(test,dummy,R,Q) % complete call
%
%
% INPUTS:
%
% test: (NxM) data set with the observations to be compared. These data 
%   should be preprocessed in the same way than calibration data for R and 
%   Q.
%
% dummy: (Nx1) dummy variable containing 1 for the observations in the
%   first group for the comparison performed in oMEDA, -1 for the 
%   observations in the second group, and 0 for the rest of observations.
%
% R: (MxLVs) Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (test = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = test*W*inv(P'*W)*Q), this matrix is W*inv(P'*W).  
%
% Q: (MxLVs) Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (test = T*P'), this is the matrix of 
%   loadings P. For PLS (Y = test*W*inv(P'*W)*Q), this matrix is also P. 
%   (Q=R by default)
%
%
% OUTPUTS:
%
% omeda_vec: (Mx1) oMEDA vector.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
%
% Copyright (C) 2014  University of Granada, Granada
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

if nargin < 3, error('Error in the number of arguments.'); end;
s = size(test);
if ndims(dummy)==2 & find(size(dummy)==max(size(dummy)))==2, dummy = dummy'; end
if size(R,2)<1 || size(R,1)~=s(2) || size(dummy,1)~=s(1), error('Error in the dimension of the arguments.'); end;
if nargin < 4, Q = R; end;
if size(R) ~= size(Q), error('Error in the dimension of the arguments.'); end;

%% Main code

ind=find(dummy>0);
dummy(ind) = dummy(ind)/max((dummy(ind)));
ind=find(dummy<0);
dummy(ind) = -dummy(ind)/min((dummy(ind)));

xA = test*R*Q';
sumA = xA'*dummy;

omeda_vec = (2*(test'*dummy).*abs(sumA) - sumA.*abs(sumA))./sqrt(dummy'*dummy);

    


        