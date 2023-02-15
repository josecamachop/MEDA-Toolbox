function [T, parglmo] = parglmMC(X, F, interactions, prep, n_perm, ts, ordinal, mtc, fmtc)

% Parallel General Linear Model to obtain factor and interaction matrices 
% in a crossed experimental design and permutation test for univariate 
% statistical significance through multiple-test correction methods that 
% allows variable selection
%
% Related routines: parglm, parglmVS, asca, apca, create_design
%
% T = parglmMC(X, F)   % minimum call
% [T, parglmoMC] = parglmMC(X, F, interactions, prep, n_perm, ts, ordinal, mtc, fmtc)   % complete call
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
%       0: nominal
%       1: ordinal
%
% mtc: [1x1] Multiple test correction
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR, by default)
%       -pi0: Q-value assuming 0<pi0<=1 the ratio of null variables   
% 
% fmtc: [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis.
%       0: do not correct (by default)
%       1: same as mtc 
%
% OUTPUTS:
%
% T (table): ANOVA-like output table
%
% parglmoMC (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on mtc and fmtc) and explained 
% variance  
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor. This example takes long to compute, you  may reduce the 
% number of variables or permutations.
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
% rng(S);
% [TMC, parglmoMC] = parglmMC(X, class); % With variable selection through multiple testing correction
%
% h = figure; hold on
% plot([1 n_vars],[parglmo.p parglmo.p],'b-.')
% plot(parglmoVS.p(parglmoVS.ord_factors),'g-o')
% plot(parglmoMC.p(parglmoMC.ord_factors),'k-')
% plot([0,size(X,2)],[0.05 0.05],'r:')
% plot([0,size(X,2)],[0.01 0.01],'r--')
% legend('ASCA','VASCA','BH-FDR','alpha=0.05','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('p-values','FontSize',18)
% xlabel('Variables in selected order','FontSize',18)
%
%
% coded by: Jos� Camacho (josecamacho@ugr.es)
% last modification: 14/Dec/22
%
% Copyright (C) 2022  Jos� Camacho, Universidad de Granada
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
if nargin < 8 || isempty(mtc), mtc = 3; end;
if nargin < 9 || isempty(fmtc), fmtc = 0; end;

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(mtc), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: 9th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

n_interactions      = size(interactions,1);              % number of interactions
n_factors           = size(F,2);                         % number of factors
mtcc                = n_factors + n_interactions;        % correction for the number of tests
ts_factors         = zeros(n_perm*M+1,n_factors,M);       % sum of squares for factors
ts_interactions    = zeros(n_perm*M+1,n_interactions,M);  % sum of squares for interactions
p_factor            = zeros(n_factors,M);                % p-values factors
p_interaction       = zeros(n_interactions,M);           % p-values interactions

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,prep);
X = X./(ones(size(X,1),1)*dt);

SSQ_X = sum(X.^2);

% Make structure with unchanging 'variables'
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
        parglmo.n_levels(f) = length(uF);
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
        df(f) = size(unique(D(:,parglmo.factors{f}.Dvars),'rows'),1) -1;
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
SSQ_inter = sum(parglmo.inter.^2);
SSQ_residuals = sum(X_residuals.^2);

SSQ_factors = [];
F_factors = [];
for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(f,:) = sum(parglmo.factors{f}.matrix.^2);
    F_factors(f,:) = (SSQ_factors(f,:)/df(f))./(SSQ_residuals/Rdf);
end

% Interactions
SSQ_interactions = [];
F_interactions = [];
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(i,:) = sum(parglmo.interactions{i}.matrix.^2);
    F_interactions(i,:) = (SSQ_interactions(i,:)/df_int(i))./(SSQ_residuals/Rdf);
end
    
if n_interactions
    parglmo.effects = 100*([SSQ_inter' SSQ_factors' SSQ_interactions' SSQ_residuals']./(SSQ_X'*ones(1,2+n_factors+n_interactions)));
else
    parglmo.effects = 100*([SSQ_inter' SSQ_factors' SSQ_residuals']./(SSQ_X'*ones(1,2+n_factors+n_interactions)));
end
parglmo.residuals = X_residuals;

% Select test statistic to order variables
if ts
    ts_factors(1,:,:) = F_factors;
    ts_interactions(1,:,:) = F_interactions;
else
    ts_factors(1,:,:) = SSQ_factors;
    ts_interactions(1,:,:) = SSQ_interactions;
end
    
% Permutations
for j = 1 : (n_perm * M) % Increase the number of permutations to perform MTC
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
    
    B = pD*X(perms, :);
    X_residuals = X(perms, :) - D*B;
    SSQ_residuals_p = sum(X_residuals.^2);
    
    % Factors
    SSQ_factors_p = [];
    F_factors_p = [];
    for f = 1 : n_factors
        factor_matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors_p(f,:) = sum(factor_matrix.^2);
        F_factors_p(f,:) = (SSQ_factors_p(f,:)/df(f))./(SSQ_residuals_p/Rdf);
    end
    
    % Interactions
    SSQ_interactions_p = [];
    F_interactions_p = [];
    for i = 1 : n_interactions
        interacts_matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQ_interactions_p(i,:) = sum(interacts_matrix.^2);
        F_interactions_p(i,:) = (SSQ_interactions_p(i,:)/df_int(i))./(SSQ_residuals_p/Rdf);
    end
   
    % Select test statistic to order variables
    if ts
        ts_factors(1 + j,:,:) = F_factors_p;
        ts_interactions(1 + j,:,:) = F_interactions_p;
    else
        ts_factors(1 + j,:,:) = SSQ_factors_p;
        ts_interactions(1 + j,:,:) = SSQ_interactions_p;
    end
end        
   

% Calculate univariate p-values and order variables by relevance
for f = 1 : n_factors
    for var = 1 : M
        p_factor(f,var) = (size(find(ts_factors(2:end,f,var) ...
            >= ts_factors(1,f,var)) ,1) + 1)/(n_perm*M + 1);
    end
    [~,ord_factors(f,:)] = sort(p_factor(f,:),'ascend');
end
for i = 1 : n_interactions
    for var = 1 : M
        p_interaction(i,var) = (size(find(ts_interactions(2:end,i,var) ...
            >= ts_interactions(1,i,var)) ,1) + 1)/(n_perm*M + 1);
    end
    [~,ord_interactions(i,:)] = sort(p_interaction(i,:),'ascend');
end

parglmo.ord_factors = ord_factors;
if n_interactions>0,
    parglmo.ord_interactions = ord_interactions;
end

parglmo.p = [p_factor' p_interaction'];

% Multiple test correction
switch mtc
    
    case 1 % Bonferroni
        if fmtc
            parglmo.p = min(1,parglmo.p * M * mtcc);
        else
            parglmo.p = min(1,parglmo.p * M);
        end
        
    case 2 % Holm/Hochberg
        
        if fmtc
            parglmo.p = parglmo.p;
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = 1 : length(indx)
                parglmo.p(indx(ind)) = min(1,parglmo.p(indx(ind)) * ((M*mtcc)-ind+1));
            end
        else
            for f = 1 : n_factors
                for var = 1 : M
                    p_factor(f,ord_factors(f,var)) = min(1,p_factor(f,ord_factors(f,var)) * (M-var+1));
                end
            end
            for i = 1 : n_interactions
                for var = 1 : M
                    p_interaction(i,ord_interactions(i,var)) = min(1,p_interaction(i,ord_interactions(i,var)) * (M-var+1));
                end
            end
            parglmo.p = [p_factor' p_interaction'];
        end
        
    case 3 % Benjamini & Hochberg
        
        if fmtc
            parglmo.p = parglmo.p;
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = length(indx)-1 : -1 : 1
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * (M*mtcc)/ind]);
            end
        else
            for f = 1 : n_factors
                for var = M-1 : -1 : 1
                    p_factor(f,ord_factors(f,var)) = min([1,p_factor(f,ord_factors(f,var)) * M/var]);
                end
            end
            for i = 1 : n_interactions
                for var = M-1 : -1 : 1
                    p_interaction(i,ord_interactions(i,var)) = min([1,p_interaction(i,ord_interactions(i,var)) * M/var]);
                end
            end
            parglmo.p = [p_factor' p_interaction'];
        end
        
end

if mtc <0 % Q-value     
    if fmtc
        parglmo.p = parglmo.p;
        [~,indx] = sort(parglmo.p(:),'ascend');
        parglmo.p(indx(end)) = parglmo.p(indx(end))*(-mtc);
        for ind = length(indx)-1 : -1 : 1
            parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * (-mtc*M*mtcc)/ind,parglmo.p(indx(ind+1))]);
        end
    else
        for f = 1 : n_factors
            p_factor(f,ord_factors(f,M)) = p_factor(f,ord_factors(f,M))*(-mtc);
            for var = M-1 : -1 : 1
                p_factor(f,ord_factors(f,var)) = min([1,p_factor(f,ord_factors(f,var)) * (-mtc) * M/var,p_factor(f,ord_factors(f,var+1))]);
            end
        end
        for i = 1 : n_interactions
            p_interaction(i,ord_interactions(i,M)) = p_interaction(i,ord_interactions(i,M))*(-mtc);
            for var = M-1 : -1 : 1
                p_interaction(i,ord_interactions(i,var)) = min([1,p_interaction(i,ord_interaction(i,var)) * (-mtc) * M/var,p_interaction(i,ord_interactions(i,var+1))]);
            end
        end
        parglmo.p = [p_factor' p_interaction'];
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
      
if n_interactions
    SSQ = sum([SSQ_inter' SSQ_factors' SSQ_interactions' SSQ_residuals' SSQ_X'],1);
else
    SSQ = sum([SSQ_inter' SSQ_factors' SSQ_residuals' SSQ_X'],1);
end
par = [mean(parglmo.effects) 100];
DoF = [1 df df_int Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan mean(F_factors,2)' mean(F_interactions,2)' nan nan];
p_value = [nan mean(parglmo.p) nan nan];

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    T.data = [name'; SSQ'; par'; DoF'; MSQ'; F'; p_value'];
    T.labels = {'Source','SumSq','AvPercSumSq','df','MeanSq','AvF','AvPvalue'};
else
    T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','AvPercSumSq','df','MeanSq','AvF','AvPvalue'});
end


