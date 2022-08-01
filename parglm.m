function [parglmo,T] = parglm(X, F, interactions, prep, n_perm, ordinal)

% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for significance. This
% approach permutes the raw values, which is sub-optimal in terms of power 
% according to Andreson and Ter Braak.
%
% parglmo = parglm(X, F)   % minimum call
% parglmo = parglm(X, F, interactions, prep, n_perm, ordinal)   % complete call
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
% ordinal: [1xF] whether factors are nominal or ordinal
%       0: nominal
%       1: ordinal
%
%
% OUTPUTS:
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values and explained variance. 
%
% T (table): ANOVA-like output.
%
%
% EXAMPLE OF USE: Random data, two significative factors, with 4 and 3
%   levels, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(randn(1,vars),reps,1);
%     end
% end
%
% parglmo = parglm(X, F);
%
%
% EXAMPLE OF USE: Random data, two factors, with 4 and 3 levels, but only
%   the first one is significative, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
%
% parglmo = parglm(X, F);
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 01/Aug/22
%
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(X, 1);
M = size(X, 2);
if nargin < 3 || isempty(interactions), interactions = []; end;
if nargin < 4 || isempty(prep), prep = 2; end;
if nargin < 5 || isempty(n_perm), n_perm = 1000; end;
if nargin < 6 || isempty(ordinal), ordinal = zeros(1,size(F,2)); end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
SSQ_factors         = zeros(1,n_factors,1);      % sum of squares for factors
SSQ_interactions    = zeros(1,n_interactions);   % sum of squares for interactions
F_factors         = zeros(n_perm + 1,n_factors,1);      % F for factors
F_interactions    = zeros(n_perm + 1,n_interactions);   % F for interactions
p_factor            = zeros(1, n_factors);       % p-values factors
p_interaction       = zeros(1, n_interactions);  % p-values interactions

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,prep);
X = X./(ones(size(X,1),1)*dt);

SSQ_X = sum(sum(X.^2));

% Make structure with general 'variables'
parglmo.data           = X;
parglmo.design         = F;
parglmo.n_factors      = n_factors;
parglmo.n_interactions = n_interactions;

% Create Design Matrix
n = 1;
D = ones(size(X,1),1);

for f = 1 : n_factors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),1);
        parglmo.factors{f}.Dvars = n+1;
        n = n + 1;
    else
        uF = unique(F(:,f));
        paranovao.n_levels(f) = length(uF);
        for i = 1:length(uF)-1
            D(find(F(:,f)==uF(i)),n+i) = 1;
        end
        parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
        D(find(F(:,f)==uF(end)),parglmo.factors{f}.Dvars) = -1;
        n = n + length(uF) - 1;
    end
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

% Degrees of freedom

Tdf = size(X,1);      
Rdf = Tdf-1;
for f = 1 : n_factors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = size(unique(D(:,parglmo.factors{f}.Dvars),'Rows'),1) -1;
    end
    Rdf = Rdf-df(f);
end
df_int = [];
for i = 1 : n_interactions
    df_int(i)=(df(interactions(i,1))+1)*(df(interactions(i,2))+1) - 1;
    Rdf = Rdf-df_int(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end;
    
% GLM model calibration with LS, only fixed factors

B = pinv(D'*D)*D'*X;
X_residuals = X - D*B;
parglmo.D = D;
parglmo.B = B;

% Create Effect Matrices

parglmo.inter = D(:,1)*B(1,:);
SSQ_inter = sum(sum(parglmo.inter.^2));

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(1,f) = sum(sum(parglmo.factors{f}.matrix.^2));
    F_factors(1,f) = (sum(sum(parglmo.factors{f}.matrix.^2))/df(f))/(sum(sum((X_residuals).^2))/Rdf);
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(1,i) = sum(sum(parglmo.interactions{i}.matrix.^2));
    F_interactions(1,i) = (sum(sum(parglmo.interactions{i}.matrix.^2))/df_int(i))/(sum(sum((X_residuals).^2))/Rdf);
end

SSQ_residuals = sum(sum(X_residuals.^2));
parglmo.effects = 100*([SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals]./SSQ_X);
parglmo.residuals = X_residuals;

% Interactions p-values are calculated through permutation of X - Xa - Xa
% Do the permutations (do this, for now, for two factors)
for j = 1 : n_perm
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
      
    B = pinv(D'*D)*D'*X(perms, :);
    X_residuals = X(perms, :) - D*B;
    
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        F_factors(1 + j,f) = (sum(sum(factors{f}.matrix.^2))/df(f))/(sum(sum((X_residuals).^2))/Rdf);
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        F_interactions(1 + j,i) = (sum(sum(interacts{i}.matrix.^2))/df_int(i))/(sum(sum((X_residuals).^2))/Rdf);
    end

end        % permutations

% Calculate p-values
% how many ssq's are larger than measurement ssq?
for factor = 1 : n_factors
    p_factor(factor) = (size(find(F_factors(2:n_perm + 1, factor) >= F_factors(1, factor)),1) + 1)/(n_perm);
end
for interaction = 1 : n_interactions
    p_interaction(interaction) = (size(find(F_interactions(2:n_perm + 1, interaction) ...
        >= F_interactions(1, interaction)),1) + 1)/(n_perm);
end
parglmo.p = [p_factor p_interaction];


%% ANOVA-like output table

name={'Mean'};
for factor = 1 : n_factors
    name{end+1} = sprintf('Factor %d',factor);
end
for interaction = 1 : n_interactions
    name{end+1} = sprintf('Interaction %d',interaction);
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
      
SSQ = [SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals SSQ_X];
par = [parglmo.effects 100];
DoF = [1 df df_int Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan F_factors(1,:) F_interactions(1,:) nan nan];
p_value = [nan p_factor p_interaction nan nan];

T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F','Pvalue'});


