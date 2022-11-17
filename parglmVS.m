function [T, parglmo] = parglmVS(X, F, interactions, prep, n_perm, ts, ordinal, fmtc)

% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for incremental
% multivariate statistical significance that allows variable selection. This 
% is the basis of VASCA (Variable-selection ASCA).
%
% Related routines: parglm, parglmMC, asca, apca, create_design
%
% T = parglmVS(X, F)   % minimum call
% [T, parglmoVS] = parglmVS(X, F, interactions, prep, n_perm, ts, ordinal, fmtc)   % complete call
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
% fmtc: [1x1] whether to correct for multiple-tesis when multifactorial
% analysis or not.
%       0: do not correct (default)
%       1: correct
%
%
% OUTPUTS:
%
% T (table): ANOVA-like output table
%
% parglmoMV (structure): structure with the factor and interaction
% matrices, p-values and explained variance 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor.
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
% [TVS, parglmoVS] = parglmVS(X, class); % With variable selection
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
% last modification: 16/Nov/22
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

n_interactions      = size(interactions,1);              % number of interactions
n_factors           = size(F,2);                         % number of factors
if fmtc,
    mtcc                = n_factors + n_interactions;        % correction for the number of tests
else,
    mtcc = 1;
end
SSQ_factors         = zeros(n_perm*mtcc+1,n_factors,M);       % sum of squares for factors
SSQ_interactions    = zeros(n_perm*mtcc+1,n_interactions,M);  % sum of squares for interactions
SSQ_residuals       = zeros(n_perm*mtcc+1,M);                 % sum of squares for residuals
F_factors           = zeros(n_perm*mtcc+1,n_factors,M);       % F-value 
F_interactions      = zeros(n_perm*mtcc+1,n_interactions,M);  % F-value 
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
SSQ_inter = sum(parglmo.inter.^2);
SSQ_residuals(1,:) = sum(X_residuals.^2);

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(1,f,:) = sum(parglmo.factors{f}.matrix.^2);
    F_factors(1,f,:) = squeeze(SSQ_factors(1,f,:)/df(f))./(SSQ_residuals(1,:)/Rdf)';
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(1,i,:) = sum(parglmo.interactions{i}.matrix.^2);
    F_interactions(1,i,:) = squeeze(SSQ_interactions(1,i,:)/df_int(i))./(SSQ_residuals(1,:)/Rdf)';
end

if n_interactions
    parglmo.effects = 100*([SSQ_inter' permute(SSQ_factors(1,:,:),[3 2 1]) permute(SSQ_interactions(1,:,:),[3 2 1]) SSQ_residuals(1,:)']./(SSQ_X'*ones(1,2+n_factors+n_interactions)));
else
    parglmo.effects = 100*([SSQ_inter' permute(SSQ_factors(1,:,:),[3 2 1]) SSQ_residuals(1,:)']./(SSQ_X'*ones(1,2+n_factors+n_interactions)));
end
parglmo.residuals = X_residuals;

% Permutations
for j = 1 : n_perm*mtcc
    
    perms = randperm(size(X,1)); % permuted data (permute whole rows)
    
    B = pD*X(perms, :);
    X_residuals = X(perms, :) - D*B;
    SSQ_residuals(1 + j,:) = sum(X_residuals.^2);
    
    % Factors
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors(1 + j,f,:) = sum(factors{f}.matrix.^2);
        F_factors(1 + j,f,:) = squeeze(SSQ_factors(1 + j,f,:)/df(f))./(SSQ_residuals(1 + j,:)/Rdf)';
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQ_interactions(1 + j,i,:) = sum(interacts{i}.matrix.^2);
        F_interactions(1 + j,i,:) = squeeze(SSQ_interactions(1 + j,i,:)/df_int(i))./(SSQ_residuals(1 + j,:)/Rdf)';
    end

end       

% Select test statistic to order variables
if ts
    ts_factors = F_factors;
    ts_interactions = F_interactions;
else
    ts_factors = SSQ_factors;
    ts_interactions = SSQ_interactions;
end
    
% Order variables by relevance
ord_factors =  zeros(n_factors,M);
ord_interactions =  zeros(n_interactions,M);
for f = 1 : n_factors
    [~,ord_factors(f,:)] = sort(ts_factors(1,f,:),'descend');
    for var = 1 : M
        F_factors(1,f,ord_factors(f,var)) = (sum(SSQ_factors(1,f,ord_factors(f,1:var)),3)/df(f))./(sum(SSQ_residuals(1,ord_factors(f,1:var)),2)/Rdf);
    end    
        
    for j = 1 : n_perm*mtcc
        [~,ord] = sort(ts_factors(1 + j,f,:),'descend');
        SSQ_factors(1 + j,f,:) = SSQ_factors(1 + j,f,ord);
        for var = 1 : M
            F_factors(1 + j,f,var) = (sum(SSQ_factors(1 + j,f,1:var),3)/df(f))./(sum(SSQ_residuals(1 + j,ord(1:var)),2)/Rdf);
        end    
    end
end
for i = 1 : n_interactions
    [~,ord_interactions(i,:)] = sort(ts_interactions(1,i,:),'descend');
    for var = 1 : M
        F_interactions(1,i,ord_interactions(i,var)) = (sum(SSQ_interactions(1,i,ord_interactions(i,1:var)),3)/df_int(i))./(sum(SSQ_residuals(1,ord_interactions(i,1:var)),2)/Rdf);
    end  
    
    for j = 1 : n_perm*mtcc
        [~,ord] = sort(ts_interactions(1 + j,i,:),'descend');
        SSQ_interactions(1 + j,i,:) = SSQ_interactions(1 + j,i,ord);
        for var = 1 : M
            F_interactions(1 + j,i,var) = (sum(SSQ_interactions(1 + j,i,1:var),3)/df_int(i))./(sum(SSQ_residuals(1 + j,ord(1:var)),2)/Rdf);
        end    
    end
end

parglmo.ord_factors = ord_factors;
if n_interactions>0
    parglmo.ord_interactions = ord_interactions;
end

% Calculate multivariate p-values
for f = 1 : n_factors
    for var = 1 : M
        if ts
            p_factor(f,ord_factors(f,var)) = (size(find( F_factors(2:(n_perm*mtcc + 1),f,var) ...
                >= F_factors(1, f,ord_factors(f,var))), 1) + 1)/(n_perm*mtcc+1);
        else
            p_factor(f,ord_factors(f,var)) = (size(find( sum(SSQ_factors(2:(n_perm*mtcc + 1),f,1:var),3) ...
                >= sum(SSQ_factors(1, f,ord_factors(f,1:var)))), 1) + 1)/(n_perm*mtcc+1);
        end
    end
end
for i = 1 : n_interactions
    for var = 1 : M
        if ts
            p_interaction(i,ord_interactions(i,var)) = (size(find( F_interactions(2:(n_perm*mtcc + 1),i,var) ...
                >= F_interactions(1, i, ord_interactions(i,var))), 1) + 1)/(n_perm*mtcc+1);
        else
            p_interaction(i,ord_interactions(i,var)) = (size(find( sum(SSQ_interactions(2:(n_perm*mtcc + 1),i,1:var),3) ...
                >= sum(SSQ_interactions(1, i, ord_interactions(i,1:var)))), 1) + 1)/(n_perm*mtcc+1);
        end
    end
end
        
% Multiple test correction for several factors/interactions
p_factor = min(1,p_factor * mtcc); 
p_interaction = min(1,p_interaction * mtcc);

parglmo.p = [p_factor' p_interaction'];


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
    SSQ = sum([SSQ_inter' permute(SSQ_factors(1,:,:),[3 2 1]) permute(SSQ_interactions(1,:,:),[3 2 1]) SSQ_residuals(1,:)' SSQ_X'],1);
else
    SSQ = sum([SSQ_inter' permute(SSQ_factors(1,:,:),[3 2 1]) SSQ_residuals(1,:)' SSQ_X'],1);
end
par = [mean(parglmo.effects) 100];
DoF = [1 df df_int Rdf Tdf];
MSQ = SSQ./DoF;
F = [nan mean(F_factors(1,:,:),3) mean(F_interactions(1,:,:),3) nan nan];
p_value = [nan mean(p_factor,2)' mean(p_interaction,2)' nan nan];

T = table(name', SSQ', par', DoF', MSQ', F', p_value','VariableNames', {'Source','SumSq','AvPercSumSq','df','MeanSq','AvF','AvPvalue'});


 

    