
function T2 = leverage_pls(cal,y,lvs,test,prepx,prepy,opt,label)

% Compute and plot leverages in PCA.
%
% leverage_pls(cal,y,lvs) % minimum call
% leverage_pls(cal,y,lvs,test,prepx,prepy,opt,label) % complete call
%
% INPUTS:
%
% cal: (LxM) billinear data set for model fitting.
%
% y: (LxO) billinear data set of predicted variables
%
% lvs: (1xA) Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two lvs)
%
% test: (NxM) data set with the observations to be compared. These data 
%   are preprocessed in the same way than calibration data.
%
% prepx: (1x1) preprocesing of the x-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% prepy: (1x1) preprocesing of the y-block
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)
%
% opt: (1x1) options for data plotting.
%       0: no plots
%       1: Squared residuals in the observations (default)
%       2: Squared residuals in the variables 
%       3: Squared residuals in the observations with control limits
%
% label: name of the observations (opt 1, dimension ((L+N)x1) or 
%   variables (opt 2, dimension (Mx1)) (numbers are used by default), eg.
%   num2str((1:L+N))')'
%
%
% OUTPUTS:
%
% T2: leverage (opt 0,1 or 3, dimension {1x(L+N)} or opt 2 or 4, dimension 
%   {1x(M)})
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 03/Jul/14.
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

%% Parameters checking

if nargin < 3, error('Error in the number of arguments.'); end;
if nargin < 4, x = cal; else x = [cal;test]; end;
s = size(x);
if s(1) < 1 || s(2) < 1 || ndims(x)~=2, error('Error in the dimension of the arguments.'); end;
sp = length(lvs);
if nargin < 5, prepx = 2; end;
if nargin < 6, prepy = 2; end; 
if nargin < 7, opt = 1; end;
if nargin < 8, label = []; end

%% Main code

[calp,m,dt] = preprocess2D(cal,prepx);
yp = preprocess2D(y,prepy);

[beta,W,P] = kernel_pls(calp'*calp,calp'*yp,max(lvs));
W2 = W*inv(P'*W);
W2 = W2(:,lvs);
P = P(:,lvs);
T = calp*W2;

if exist('test')&~isempty(test),
    testp = (test - ones(size(test,1),1)*m)./(ones(size(test,1),1)*dt);
    TT = testp*P;
else
    testp = [];
    TT = [];
end

[Ts,notused,dtT2] = preprocess2D(T,2);

switch opt,
    case 2
        T2 = diag(P*P');
    otherwise
        T2 = diag(Ts*Ts');
        if ~isempty(test)
            T2 = [T2;diag(TT*diag(1./(dtT2.^2))*TT')];
        end
end;

if opt, 
    if opt<3,
        plot_vec(T2,label,'Leverage');
    else
        plot_vec(T2,label,'Leverage',(ones(size(T2,1),1)*[hot_lim(length(lvs),length(T2),0.05) hot_lim(length(lvs),length(T2),0.01)])');
    end
end