
function [Dst,Qst] = mspc(testcs,varargin)

% Multivariate Statistical Process Control statistics
%
% Dst = mspc(testcs) % minimum call (only Q-st)
%
%
% INPUTS:
%
% testcs: [NxM] preprocessed billinear data set with the observations to be 
%   monitored.
%
%
% Optional INPUTS (parameters):
%
% 'InvCovarT': [AxA] inverse of covariance matrix of T, where T are the 
%   calibration scores.
%
% 'InSubspace': [MxA] Matrix to perform the projection from the original to the  
%   latent subspace. For PCA (testcs = T*P'), this is the matrix of loadings 
%   P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is W*inv(P'*W). For the 
%   original space (default) the identity matrix is used.  
%
% 'OutSubspace': [MxA] Matrix to perform the projection from the latent subspace to 
%   the original space. For PCA (testcs = T*P'), this is the matrix of 
%   loadings P. For PLS (Y = testcs*W*inv(P'*W)*Q), this matrix is also P. 
%   For the original space the identity matrix is used. Q=R is used by 
%   default. 
%
%
% OUTPUTS:
%
% Dst: [Nx1] D-statistic or Hotelling T2
%
% Qst: [Nx1] Q-statistic
%
%
%
% EXAMPLE OF USE: PCA-based MSPC on NOC test data and anomalies.
%
% nobs = 100;
% nvars = 10;
% nPCs = 1;
% X = simuleMV(nobs,nvars,'LevelCorr',6);
% [Xcs, m, sc] = preprocess2D(X,'Preprocessing',2);
% 
% pcs = 1:nPCs;
% model = pcaEig(Xcs,'PCs',pcs);
% p = model.loads;
% t = model.scores;
% e = Xcs - t*p';
% UCLd = hotLim(nPCs,nobs,0.05,'Phase',2);
% UCLq = speLim(e,0.01);
% 
% nobst = 10;
% test = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',cov(X)*(nobst-1));
% test(6:10,:) = 3*test(6:10,:);
% testcs = preprocess2Dapp(test,m,'Scale',sc);
% 
% [Dst,Qst] = mspc(testcs,'InvCovarT',inv(cov(t)),'InSubspace',p);
% 
% plotScatter([Dst,Qst],'ObsClass',[ones(5,1);2*ones(5,1)],'XYLabel',{'D-st','Q-st'},'LimCont',{UCLd,UCLq}); 

% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 20/Nov/2024
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
M = size(testcs, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'InvCovarT',[]); 
A = size('InvCovarT', 1);
addParameter(p,'InSubspace',zeros(M,A));
addParameter(p,'OutSubspace',[]);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
invCT = p.Results.InvCovarT;
A = size(invCT, 1);
R = p.Results.InSubspace;
Q = p.Results.OutSubspace;
if isempty(R), zeros(M,A); end;
if isempty(Q), Q = R; end;


% Validate dimensions of input data
assert (isequal(size(invCT), [A A]), 'Dimension Error: parameter ''InCovarT'' must be LVs-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(R), [M A]), 'Dimension Error: parameter ''InSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Q), [M A]), 'Dimension Error: parameter ''OutSubspace'' must be M-by-LVs. Type ''help %s'' for more info.', routine(1).name);


%% Main code

t = testcs * R;
e = testcs - t * Q';

Dst = sum((t * invCT) .* t,2);
Qst = sum(e.^2,2);


    


        