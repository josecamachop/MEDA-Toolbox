function PEVpq = razorPlot(X,Gram,K,varargin)

% Razor plot to select the number of components and sparsity in sparse 
% Principal Component Analysis (sPCA) following Camacho et al. "All sparse 
% PCA models are wrong, but some are useful. Part III: model interpretation", 
% Chemometrics and Intelligent Laboratory Systems.
%
% PEVpq = razorPlot(X,Gram,K)     % minimum call
%
%
% INPUTS:
%
% X: [NxM] preprocessed billinear data set 
%
% Gram: [NxN] Gram matrix from X (either X or Gram have to be inputted) 
%
% K: [1x1] Number of sparse components.
%
%
% Optional INPUTS (parameters):
%
% 'Tolerance': [1x1] tolerance value. By default, 1e-15.
%
% 'MaxIters': [1x1] maximum iterations. By default, 1e3.
%
% 'Threshold': [1x1] threshold to use a truncated plot. By default, 0.05.
%
% 'Reference': [1x1] reference variance. By default, 1.
%
%
% OUTPUTS:
%
% PEVpq: [MxM...xMxK] percentage of explained variance of the model
% variants.
%
%
% EXAMPLE OF USE: Pitprops
%
% var_l = {'topdiam' 'length'  'moist' 'testsg' 'ovensg' 'ringtop' 'ringbut' 'bowmax' 'bowdist' 'whorls'  'clear'  'knots' 'diaknot'};
%
% XX=[ 1.000  0.954  0.364  0.342 -0.129   0.313   0.496  0.424   0.592  0.545  0.084 -0.019   0.134
%      0.954  1.000  0.297  0.284 -0.118   0.291   0.503  0.419   0.648  0.569  0.076 -0.036   0.144
%      0.364  0.297  1.000  0.882 -0.148   0.153  -0.029 -0.054   0.125 -0.081  0.162  0.220   0.126
%      0.342  0.284  0.882  1.000  0.220   0.381   0.174 -0.059   0.137 -0.014  0.097  0.169   0.015
%     -0.129 -0.118 -0.148  0.220  1.000   0.364   0.296  0.004  -0.039  0.037 -0.091 -0.145  -0.208
%      0.313  0.291  0.153  0.381  0.364   1.000   0.813  0.090   0.211  0.274 -0.036  0.024  -0.329
%      0.496  0.503 -0.029  0.174  0.296   0.813   1.000  0.372   0.465  0.679 -0.113 -0.232  -0.424
%      0.424  0.419 -0.054 -0.059  0.004   0.090   0.372  1.000   0.482  0.557  0.061 -0.357  -0.202
%      0.592  0.648  0.125  0.137 -0.039   0.211   0.465  0.482   1.000  0.526  0.085 -0.127  -0.076
%      0.545  0.569 -0.081 -0.014  0.037   0.274   0.679  0.557   0.526  1.000 -0.319 -0.368  -0.291
%      0.084  0.076  0.162  0.097 -0.091  -0.036  -0.113  0.061   0.085 -0.319  1.000  0.029   0.007
%     -0.019 -0.036  0.220  0.169 -0.145   0.024  -0.232 -0.357  -0.127 -0.368  0.029  1.000   0.184
%      0.134  0.144  0.126  0.015 -0.208  -0.329  -0.424 -0.202  -0.076 -0.291  0.007  0.184   1.000];
%
% ev = eig(XX);
% reference = sum(ev(end:-1:end-5)); % Use 6PCs PCA as a reference 
% razorPlot([], XX, 6, 'Reference', reference);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 1/Feb/2025
%
% Copyright (C) 2025  University of Granada, Granada
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
assert (nargin >= 3, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
    
if isempty(X)
  % Infer X from X'*X
  [Vg Dg] = eig(Gram);
  X = Vg*sqrt(abs(Dg))*Vg';
end

[N M] = size(X);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Tolerance',1e-15);
addParameter(p,'MaxIters',1e3);
addParameter(p,'Threshold',0.05);
addParameter(p,'Reference',1);          
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
tol = p.Results.Tolerance;
max_iter = p.Results.MaxIters;
thres = p.Results.Threshold;
reference = p.Results.Reference/sum(sum(X.^2)); 

%% PEV vs sparsity: SPCA-Z multi-component, truncated search

clc
PEVpq = [];
fp = [];
pcs = 1:K;
tic
[fp, PEVpq] =  computeModels(X,pcs,tol,max_iter,thres,reference);
fp = shiftdim(fp,length(pcs));
PEVpq = shiftdim(PEVpq,length(pcs));

total_time=toc;
disp(sprintf('Finished, time elapsed %d',total_time))


%% Plot the truncated razor plot

ufp = unique(fp);
PEVfp = [];
for i=2:length(ufp)
    ind = find(fp == ufp(i));
    mind = find(PEVpq(ind)==max(PEVpq(ind)),1);
    PEVfp(i-1) = PEVpq(ind(mind));
end

val = num2cell(ufp(2:end));
val{end+1} = 'Ref';
f = plotVec([PEVfp reference],'ObsClass',[2*ones(1,length(PEVfp)) 1]);
legend('off')
ylabel('PEV')
xlabel('f')
a=get(f,'Children');
set(a,'XTickLabel',val);
set(a,'XTick',1:length(val));
set(a,'XTickLabelRotation',45);

    
%% Plot the truncated razor plot per component

ufp = unique(fp);
PEVfp2D = [];
for i=2:length(ufp)
    for j = pcs
        ind = find(fp(j,:) == ufp(i));
        mind = find(PEVpq(j,ind)==max(PEVpq(j,ind)),1);
        if ~isempty(mind)
            PEVfp2D(j,i-1) = PEVpq(j,ind(mind));
        end
    end
end

figure
surf((((ones(length(pcs),1)*ufp(2:end)')))',(pcs'*ones(1,length(ufp)-1))',(PEVfp2D)')
hold on
pcolor((((ones(length(pcs),1)*ufp(2:end)')))',(pcs'*ones(1,length(ufp)-1))',(PEVfp2D)')
axis([ufp(2) ufp(end) pcs(1) pcs(end)])
colorbar
ylabel('# Components')
xlabel('f')
zlabel('PEV')


%% Recursive function 

function [fp, PEVpq, flag] =  computeModels(X,pcs,tol,max_iter,thres,reference,K,vec)

% Set default values
routine=dbstack;
assert (nargin >= 6, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
 
M = size(X,2);
flag = false;

if nargin < 7
    K = 1;
    prev = M;
else
    prev = vec(K-1);
end

if K <= length(pcs)
    fp = zeros([M*ones(1,length(pcs)-K+1),length(pcs)]);
    PEVpq = zeros([M*ones(1,length(pcs)-K+1),length(pcs)]);
end

for j=1:prev
    vec(K) = j;
    if K <= length(pcs)
        [f, p, flag] = computeModels(X,pcs,tol,max_iter,thres,reference,K+1,vec);
        fp(j,:) = f(:);
        PEVpq(j,:) = p(:);
        if flag
            return
        end
    else
        for i=1:length(pcs)
            [p,q] = spcaZou(X,[],pcs(i),-vec(1:pcs(i)),'Tolerance',tol,'MaxIters',max_iter);
            
            fp(i) = length(find(p.^2)) - length(find(sum(p.^2,1)));
            PEVpq(i) = 1 - sum(sum((X - X*p*inv(q'*p)*q').^2))/sum(sum(X.^2));
            
            if (PEVpq(i)/reference) > (1 - thres)
                flag = true;
                return
            end
        end
        
        % disp(sprintf('Model with %s non-zero elements, time elapsed %d',num2str(vec),toc))
    end
end
        