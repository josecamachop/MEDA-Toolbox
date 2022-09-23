function [bpvals, pboot] = pbootasca(X,F,ascao,nfact,nboot)

% Bootstraping in ASCA models.
%
% bpvals = pbootasca(X,F,ascao,nfact) % minimum call
% [bpvals,pboot] = pbootasca(X,F,ascao,nfact,nboot) complete
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, where columns correspond to factors and rows to
% levels.
%
% ascao: (structure) structure that contains scores, loadings, singular
% values and projections of the factors and interactions in ASCA.
%
% nfact: [1x1] factor where bootstrapping is performed
%
% nboot: [1x1] number of runs (1000 by default)
%
%
% OUTPUTS:
%
% bpvals: p-value per variables
%
% pboot: bootstrapped loadings
%
%
% EXAMPLE OF USE: Random data, three variables with information on the factor:
% This example takes long to compute, you may reduce the number of
% variables or permutations.
%
% n_obs = 40;
% n_vars = 400;
%
% class = (randn(n_obs,1)>0)+1;
% X = simuleMV(n_obs,n_vars,8);
% X(class==2,1:3) = X(class==2,1:3) + 10;
%
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
%
% ascao = asca(parglmo); % With variable selection through bootstrapping
% bpvals = pbootasca(X, class, ascao, 1);
% 
%
% h = figure; hold on
% plot([1 n_vars],[parglmo.p parglmo.p],'b-.')
% plot(parglmoVS.p(parglmoVS.ord_factors),'g-o')
% plot(bpvals(parglmoVS.ord_factors),'k-')
% plot([0,size(X,2)],[0.05 0.05],'r:')
% plot([0,size(X,2)],[0.01 0.01],'r--')
% legend('ASCA','VASCA','Bootstrap','alpha=0.05','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('p-values','FontSize',18)
% xlabel('Variables in selected order','FontSize',18)
%
%
% coded by: Rafa Vitale ()
%           José Camacho (josecamacho@ugr.es)
% last modification: 23/Sep/22
%
% Copyright (C) 2022  Raffael Vitale, Universidad de ....
% Copyright (C) 2022  José Camacho, Universidad de Granada
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
assert (nargin >= 4, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(X, 1);
M = size(X, 2);
if nargin < 5 || isempty(nboot), nboot = 1000; end;

% Validate dimensions of input data
assert (isequal(size(nfact), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nboot), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

lev = F(:,nfact);
p = ascao.factors{1}.loads;

for boot=1:nboot
    
    for nl=1:max(lev)
        
        ind=find(lev==nl);
        bootind=ceil(length(ind)*rand(length(ind),1));
        yboot(ind,:)=X(ind(bootind),:);
        
    end
    
    [~, parglmo] = parglm(yboot,F); 
    ascao = asca(parglmo); 
    pb = ascao.factors{nfact}.loads; 
    [~,pboot(boot,:,:)]=orth_proc(p,pb);
    clear yboot
    
end

for j=1:size(pboot,2)
    bpvals(j) = (min(length(find(pboot(:,j)<0)),length(find(pboot(:,j)>=0)))+1)/1001;
end


function [r,yrot]=orth_proc(x,y)

%Computes orthogonal procrustes rotation projecting matrix y onto the
%subspace spanned by matrix x.
% syntax: [r,yrot]=orth_proc(x,y)
% where r is the rotation matrix and yrot is the procrustes rotated y

[u,~,v]=svd(y'*x,0);
r=u*v';
yrot=y*r;