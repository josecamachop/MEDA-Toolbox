function parglmo = parglmVS(X, F, interactions, prep, n_perm, ts)

% Parallel GLM for ANOVA with Variable Selection to obtain multivariate factor and 
% interaction matrices in a crossed experimental design and permutation 
% test for significance. This is the basis of VASCA.
%
% parglmoVS = parglmVS(X, F)   % minimum call
% parglmoVS = parglmVS(X, F, interactions, prep, n_perm, ts)   % complete call
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
% interactions: [Ix2] matrix where rows contain the factors for which
% interactions are to be calculated.
%
% prep: [1x1] preprocesing:
%       1: mean preping
%       2: autoscaling (default)
%
% n_perm: [1x1] number of permutations (1000 by default).
%
% ts: [1x1] Use SSQ as test statistic (0, by default) or SSQ/TotalSSQ (1) 
%
%
% OUTPUTS:
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values and explained variance. 
%
%
% EXAMPLE OF USE: Random data, three variables with information on the factor:
%
% n_obs = 40;
% n_vars = 400;
%
% class = (randn(n_obs,1)>0)+1;
% X = simuleMV(n_obs,n_vars,8);
% X(class==2,1:3) = X(class==2,1:3) + 10;
%
% parglmo = parglm(X,class); % No variables selection 
% parglmoVS = parglmVS(X, class); % With variable selection
%
% h = figure; hold on
% plot([1 n_vars],[parglmo.p parglmo.p],'b-.')
% plot(parglmoVS.p(parglmoVS.ord_factors),'g-o')
% plot([0,size(X,2)],[0.05 0.05],'r:')
% plot([0,size(X,2)],[0.01 0.01],'r--')
% legend('ASCA','VASCA','alpha=0.05','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('p-values','FontSize',18)
% xlabel('Variables in selected order','FontSize',18)
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 10/Dec/21
%
% Copyright (C) 2021  José Camacho, University of Granada
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(X, 1);
M = size(X, 2);
if nargin < 3 || isempty(interactions), interactions = []; end;
if nargin < 4 || isempty(prep), prep = 2; end;
if nargin < 5 || isempty(n_perm), n_perm = 1000; end;
if nargin < 6 || isempty(ts), ts = 0; end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

size_data           = size(X);                   % size of the data matrix
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
SSQ_factors         = zeros(n_factors,size_data(2));      % sum of squares for factors
SSQ_interactions    = zeros(n_interactions,size_data(2));   % sum of squares for interactions
F_factors = zeros(n_perm+1,n_factors,size_data(2));         % F-value to make significance tests
F_interactions = zeros(n_perm+1,n_interactions,size_data(2)); % F-value to make significance tests

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,prep);
X = X./(ones(size(X,1),1)*dt);

SSQ_X = sum(X.^2);

% Make structure with unchanging 'variables'
parglmo.data           = X;
parglmo.design         = F;
parglmo.n_factors      = n_factors;
parglmo.n_interactions = n_interactions;

% Create Design Matrix
n = 1;
D = ones(size(X,1),1);

for f = 1 : n_factors
    uF = unique(F(:,f));
    parglmo.n_levels(f) = length(uF); 
    for i = 1:length(uF)-1
        D(find(F(:,f)==uF(i)),n+i) = 1;
    end
    parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
    D(find(F(:,f)==uF(end)),parglmo.factors{f}.Dvars) = -1;
    n = n + length(uF) - 1;
end

for i = 1 : n_interactions
    for j = parglmo.factors{interactions(i,1)}.Dvars
        for k = parglmo.factors{interactions(i,2)}.Dvars
            D(:,end+1) = D(:,j).* D(:,k);
        end
    end
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    n = size(D,2);
end
    
% GLM model calibration with LS, only fixed factors

B = pinv(D'*D)*D'*X;
X_residuals = X - D*B;
parglmo.D = D;
parglmo.B = B;

% Create Effect Matrices

parglmo.inter = D(:,1)*B(1,:);
SSQ_inter = sum(parglmo.inter.^2);

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(f,:) = sum(parglmo.factors{f}.matrix.^2);
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(i,:) = sum(parglmo.interactions{i}.matrix.^2);
end

SSQ_residuals = sum(X_residuals.^2);

if ts
    for f = 1 : n_factors
        F_factors(1,f,:) = SSQ_factors(f,:)./SSQ_X;
    end
    for i = 1 : n_interactions
        F_interactions(1,i,:) = SSQ_interactions(i,:)./SSQ_X;
    end
else
    for f = 1 : n_factors
        F_factors(1,f,:) = SSQ_factors(f,:);
    end
    for i = 1 : n_interactions
        F_interactions(1,i,:) = SSQ_interactions(i,:);
    end
end
    
parglmo.effects = 100*([SSQ_inter' SSQ_factors' SSQ_interactions' SSQ_residuals']./(SSQ_X'*ones(1,2+n_factors+n_interactions)));
parglmo.residuals = X_residuals;

disp('Percentage each effect contributes to the total sum of squares')
disp('Overall means: variables average')
disp(mean(parglmo.effects(:,1)))
disp('Factors: variables average')
disp(mean(parglmo.effects(:,1 + (1 : n_factors))))
if n_interactions>0
    disp('Interactions: variables average')
    disp(mean(parglmo.effects(:,1 + n_factors + (1 : n_interactions))))
end
disp('Residuals: variables average')
disp(mean(parglmo.effects(:,end)))

for j = 1 : n_perm
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
    
    B = pinv(D'*D)*D'*X(perms, :);
    
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors(f,:) = sum(factors{f}.matrix.^2);
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQ_interactions(i,:) = sum(interacts{i}.matrix.^2);
    end
    
    if ts
        for f = 1 : n_factors
            F_factors(1 + j,f,:) = SSQ_factors(f,:)./SSQ_X;
        end
        for i = 1 : n_interactions
            F_interactions(1 + j,i,:) = SSQ_interactions(i,:)./SSQ_X;
        end
    else
        for f = 1 : n_factors
            F_factors(1 + j,f,:) = SSQ_factors(f,:);
        end
        for i = 1 : n_interactions
            F_interactions(1 + j,i,:) = SSQ_interactions(i,:);
        end
    end
    
end        % permutations
    
% Order variables by relevance
ord_factors =  zeros(n_factors,size_data(2));
ord_interactions =  zeros(n_interactions,size_data(2));
for factor = 1 : n_factors
    [~,ord_factors(factor,:)] = sort(F_factors(1,factor,:),'descend');
    for j = 1 : n_perm
        F_factors(1 + j,factor,:) = sort(F_factors(1 + j,factor,:),'descend');
    end
end
for interaction = 1 : n_interactions
    [~,ord_interactions(interaction,:)] = sort(F_interactions(1,interaction,:),'descend');
    for j = 1 : n_perm
        F_interactions(1 + j,interaction,:) = sort(F_interactions(1 + j,interaction,:),'descend');
    end
end

parglmo.ord_factors = ord_factors;
if n_interactions>0
    parglmo.ord_interactions = ord_interactions;
end

% Calculate multivariate p-values
p_factor = zeros(n_factors,size_data(2));
p_interaction = zeros(n_interactions,size_data(2));
for factor = 1 : n_factors
    for var = 1 : size_data(2)
        p_factor(factor,ord_factors(factor,var)) = (size(find( sum(F_factors(2:end,factor,1:var),3) ...
            >= sum(F_factors(1, factor,ord_factors(factor,1:var)))) ,1) + 1)/(n_perm+1);
    end
end
for interaction = 1 : n_interactions
    for var = 1 : size_data(2)
        p_interaction(interaction,ord_interactions(interaction,var)) = (size(find( sum(F_interactions(2:end,interaction,1:var),3) ...
            >= sum(F_interactions(1, interaction, ord_interactions(interaction,1:var)))), 1) + 1)/(n_perm+1);
    end
end
        
disp('p-values factors: variables average')
disp (mean(p_factor'))
if n_interactions>0
    disp('p-values intractions:  variables average')
    disp(mean(p_interaction'))
end
parglmo.p = [p_factor' p_interaction'];



 

    