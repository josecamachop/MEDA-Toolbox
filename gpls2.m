function [beta,W,P,Q,R,bel] = gpls2(xcs,ycs,lvs,mtype)

% Group-wise Partial Least Squares for prediction, algorithm 2. This routine
% incorporates the map estimation, groups identification with GIA and model
% calibration with GPLS. The original paper is Camacho, J., 
% Saccenti, E. Group-wise Partial Least Squares Regression. Submitted to
% Chemometrics and Intelligent Laboratory Systems, 2016.
%
% beta = gpls2(xcs,ycs)     % minimum call
% [beta,W,P,Q,R,bel] = gpls2(xcs,ycs,lvs,mtype)    % complete call
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% ycs: [NxO] preprocessed billinear data set of predicted variables
%
% lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:rank(xcs)
%
% mtype: [1x1] type of correlation map used (3 by default)
%   1: Common correlation matrix.
%   2: XYYX normalized.
%   3: MEDA map.
%   4: oMEDA map.
%
%
% OUTPUTS:
%
% beta: [MxO] matrix of regression coefficients: W*inv(P'*W)*Q'
%
% W: [MxA] matrix of weights
%
% P: [MxA] matrix of x-loadings
%
% Q: [OxA] matrix of y-loadings
%
% R: [MxA] matrix of modified weights: W*inv(P'*W)
%
% bel: [Ax1] correspondence between LVs and States.
%
%
% EXAMPLE OF USE: Random data:
%
% Y = randn(100,2);
% X(:,1:2) = Y + 0.1*randn(100,2);
% X(:,3:10) = simuleMV(100,8,6);
% lvs = 1:2;
% Xcs = preprocess2D(X,2);
% Ycs = preprocess2D(Y,2);
% [beta,W,P,Q,R,bel] = gpls2(Xcs,Ycs,lvs);
% 
% for i=lvs,
%   plot_vec(R(:,i),[],[],{'',sprintf('Regression coefficients LV %d',i)});
% end
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 16/Nov/16.
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
N = size(xcs, 1);
M = size(xcs, 2);
O = size(ycs, 2);
if nargin < 3 || isempty(lvs), lvs = 0:rank(xcs); end;
if nargin < 4 || isempty(mtype), mtype = 3; end;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>M)) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(ycs), [N O]), 'Dimension Error: 2nd argument must be N-by-O. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: 3rd argument must be 1-by-A. Type ''help %s'' for more info.', routine.name);
assert (isequal(size(mtype), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: 3rd argument must contain positive integers. Type ''help %s'' for more info.', routine.name);
assert (isempty(find(mtype~=1 & mtype~=2 & mtype~=3 & mtype~=4)), 'Value Error: 4th argument must contain an integer from 1 to 4. Type ''help %s'' for more info.', routine(1).name);


%% Main code

switch mtype,
    case 1,
        dtx = sqrt(sum(xcs.^2));
        map = (xcs'*xcs)./(dtx'*dtx);
    case 2,
        map = xcs'*ycs*ycs'*xcs;
		map = map/max(max(abs(map)));
    case 3,
        map = meda_pls(xcs,ycs,lvs,0,0,0.1,'000');
    case 4,
        map = ones(M);
		o = zeros(M,1);
		for j=1:O,
			o2 = omeda_pls(xcs,ycs,lvs,xcs,ycs(:,j),0,0,000);
			o = o + o2./max(abs(o2));
        end

        o2 = (o*o');
		map = map.*o2/max(max(abs(o2)));
        
end


lambda=1;
states={};
while length(states)< max(lvs),
    lambda = lambda*0.99;
    [bel2,states] = gia(map,lambda,1);
end
 
[beta,W,P,Q,R,bel] = gpls(xcs,ycs,states,lvs);

