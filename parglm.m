function [T, parglmo] = parglm(X, F, interactions, prep, n_perm, ts, ordinal, fmtc)

% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for multivariate 
% statistical significance. 
%
% Related routines: asca, apca, parglmVS, parglmMC, create_design
%
% T = parglm(X, F)   % minimum call
% [T, parglmo] = parglm(X, F, interactions, prep, n_perm, ts, ordinal, fmtc)   % complete call
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
% n_perm: [1x1] number of permutations (1000 by default)
%
% ts: [1x1] Use SSQ (0) or the F-value (otherwise, by default) as test statistic  
%
% ordinal: [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
% 
% fmtc: [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis
%       0: do not correct (default)
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR)
%       4: Q-value from Benjamini-Hochberg step-down
%
%
% OUTPUTS:
%
% T (table): ANOVA-like output table
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on fmtc) and explained variance 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line) 
%   Random data, two factors, with 4 and 3 levels, but only the first one 
%   is significative, and 4 replicates
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
% table = parglm(X, F)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two significative factors, with 4 and 3 levels, and 4 
%   replicates
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
%
% F = create_design(levels,reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     fi{i} = randn(1,vars);
% end
% for j = 1:length(levels{2}),
%     fj{j} = randn(1,vars);
% end
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(fi{i} + fj{j},reps,1);
%     end
% end
%
% table = parglm(X, F, [1 2])
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Random data, two factors with 4 and 3 levels, and 4 replicates, with 
%   significant interaction
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
% table = parglm(X, F, [1 2])
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 23/Nov/22
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
if nargin < 6 || isempty(ts), ts = 1; end;
if nargin < 7 || isempty(ordinal), ordinal = zeros(1,size(F,2)); end;
if nargin < 8 || isempty(fmtc), fmtc = 0; end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
if fmtc
    mtcc                = n_factors + n_interactions;        % correction for the number of tests
else
    mtcc = 1;
end
SSQ_factors         = zeros(n_perm*mtcc,n_factors,1);      % sum of squares for factors
SSQ_interactions    = zeros(n_perm*mtcc,n_interactions);   % sum of squares for interactions
F_factors         = zeros(n_perm*mtcc + 1,n_factors,1);      % F for factors
F_interactions    = zeros(n_perm*mtcc + 1,n_interactions);   % F for interactions
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
parglmo.prep           = prep;
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
    df_int(i) = (df(interactions(i,1)))*(df(interactions(i,2)));
    Rdf = Rdf-df_int(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end
    
% GLM model calibration with LS, only fixed factors
pD =  pinv(D'*D)*D';
B = pD*X;
X_residuals = X - D*B;
parglmo.D = D;
parglmo.B = B;

% Create Effect Matrices
parglmo.inter = D(:,1)*B(1,:);
SSQ_inter = sum(sum(parglmo.inter.^2));
SSQ_residuals = sum(sum(X_residuals.^2));

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(1,f) = sum(sum(parglmo.factors{f}.matrix.^2));
    F_factors(1,f) = (SSQ_factors(1,f)/df(f))/(SSQ_residuals/Rdf);
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(1,i) = sum(sum(parglmo.interactions{i}.matrix.^2));
    F_interactions(1,i) = (SSQ_interactions(1,i)/df_int(i))/(SSQ_residuals/Rdf);
end

parglmo.effects = 100*([SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals]./SSQ_X);
parglmo.residuals = X_residuals;

% Permutations
for j = 1 : n_perm*mtcc
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
      
    B = pD*X(perms, :);
    X_residuals = X(perms, :) - D*B;
    SSQ_residualsp = sum(sum(X_residuals.^2));
    
    % Factors
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors(1 + j,f) = sum(sum(factors{f}.matrix.^2));
        F_factors(1 + j,f) = (SSQ_factors(1 + j,f)/df(f))/(SSQ_residualsp/Rdf);
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);    
        SSQ_interactions(1 + j,i) = sum(sum(interacts{i}.matrix.^2));
        F_interactions(1 + j,i) = (SSQ_interactions(1 + j,i)/df_int(i))/(SSQ_residualsp/Rdf);
    end

end        

% Select test statistic
if ts
    ts_factors = F_factors;
    ts_interactions = F_interactions;
else
    ts_factors = SSQ_factors;
    ts_interactions = SSQ_interactions;
end
    
% Calculate p-values
for f = 1 : n_factors
    p_factor(f) = (size(find(ts_factors(2:(n_perm*mtcc + 1), f) >= ts_factors(1, f)),1) + 1)/(n_perm*mtcc+1);
end
for i = 1 : n_interactions
    p_interaction(i) = (size(find(ts_interactions(2:(n_perm*mtcc + 1), i) ...
        >= ts_interactions(1, i)),1) + 1)/(n_perm*mtcc+1);
end

% Multiple test correction for several factors/interactions
parglmo.p = [p_factor p_interaction]; 
if mtcc > 1
    switch fmtc
        case 1 % Bonferroni 
            parglmo.p = min(1,parglmo.p * mtcc); 

        case 2 % Holm/Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            for ind = 1 : mtcc 
                parglmo.p(indx(ind)) = min(1,parglmo.p(indx(ind)) * (mtcc-ind+1));
            end

        case 3 % Benjamini & Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            parglmo.p(indx(mtcc)) = parglmo.p(indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * mtcc/ind]);
            end

        case 4 % Q-value from Benjamini & Hochberg
            [~,indx] = sort(parglmo.p,'ascend');
            parglmo.p(indx(mtcc)) = parglmo.p(indx(mtcc));
            for ind = mtcc-1 : -1 : 1 
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * mtcc/ind,parglmo.p(indx(ind+1))]);
            end
    end
end


%% ANOVA-like output table

name={'Mean'};
for f = 1 : n_factors
    name{end+1} = sprintf('Factor %d',f);
end
for i = 1 : n_interactions
    name{end+1} = sprintf('Interaction %d',i);
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
      
SSQ = [SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals SSQ_X];
par = [parglmo.effects 100];
DoF = [1 df df_int Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan F_factors(1,:) F_interactions(1,:) nan nan];
p_value = [nan parglmo.p nan nan];

T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F','Pvalue'});


