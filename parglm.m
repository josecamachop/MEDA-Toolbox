function parglmo = parglm(X, F, interactions, center, n_perm)

% Parallel General Linear Model to obtain multivariate factor and interaction 
% matrices in a crossed experimental design and permutation test for significance. This
% approach permutes the raw values, which is sub-optimal in terms of power 
% according to Andreson and Ter Braak.
%
% parglmo = paranova(X, F)   % minimum call
% parglmo = paranova(X, F, interactions, center, n_perm)   % complete call
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
% center: [1x1] preprocesing:
%       1: mean centering
%       2: autoscaling (default)
%
% n_perm: [1x1] number of permutations (1000 by default).
%
%
% OUTPUTS:
%
% parglmo (structure): structure with the factor and interaction
% matrices, p-values and explained variance. 
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
%           Gooitzen Zwanenburg (G.Zwanenburg@uva.nl)
% last modification: 26/Apr/21
%
% Copyright (C) 2021  José Camacho, Universidad de Granada
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
if nargin < 4 || isempty(center), center = 2; end;
if nargin < 5 || isempty(n_perm), n_perm = 1000; end;

% Validate dimensions of input data
assert (isequal(size(center), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
SSQ_factors         = zeros(n_perm + 1,n_factors,1);      % sum of squares for factors
SSQ_interactions    = zeros(n_perm + 1,n_interactions);   % sum of squares for interactions
p_factor            = zeros(1, n_factors);       % p-values factors
p_interaction       = zeros(1, n_interactions);  % p-values interactions

% In column space
parglmo.factors                = cell(n_factors,1);
parglmo.interactions           = cell(n_interactions,1);

% preprocess the data
[Xs,m,dt] = preprocess2D(X,center);
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
    uF = unique(F(:,f));
    paranovao.n_levels(f) = length(uF); 
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
SSQ_inter = sum(sum(parglmo.inter.^2));

for f = 1 : n_factors
    parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
    SSQ_factors(1,f) = sum(sum(parglmo.factors{f}.matrix.^2));
end

% Interactions
for i = 1 : n_interactions
    parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
    SSQ_interactions(1,i) = sum(sum(parglmo.interactions{i}.matrix.^2));
end

SSQ_residuals = sum(sum(X_residuals.^2));
parglmo.effects = 100*([SSQ_inter SSQ_factors(1,:) SSQ_interactions(1,:) SSQ_residuals]./SSQ_X);
parglmo.residuals = X_residuals;

disp('Percentage each effect contributes to the total sum of squares')
disp('Overall means')
disp(parglmo.effects (1))
disp('Factors')
disp(parglmo.effects (1 + (1 : n_factors)))
if n_interactions>0
    disp('Interactions')
    disp(parglmo.effects (1 + n_factors + (1 : n_interactions)))
end
disp('Residuals')
disp(parglmo.effects (end))

% Interactions p-values are calculated through permutation of X - Xa - Xa
% Do the permutations (do this, for now, for two factors)
for j = 1 : n_perm
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
      
    B = pinv(D'*D)*D'*X(perms, :);
    
    for f = 1 : n_factors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQ_factors(1 + j,f) = sum(sum(factors{f}.matrix.^2));
    end
    
    % Interactions
    for i = 1 : n_interactions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQ_interactions(1 + j,i) = sum(sum(interacts{i}.matrix.^2));
    end

end        % permutations

% Calculate p-values
% how many ssq's are larger than measurement ssq?
for factor = 1 : n_factors
    p_factor(factor) = (size(find(SSQ_factors(2:n_perm + 1, factor) >= SSQ_factors(1, factor)),1) + 1)/(n_perm);
end
for interaction = 1 : n_interactions
    p_interaction(interaction) = (size(find(SSQ_interactions(2:n_perm + 1, interaction) ...
        >= SSQ_interactions(1, interaction)),1) + 1)/(n_perm);
end
disp('p-values factors:')
disp (p_factor)
if n_interactions>0
    disp('p-values intractions')
    disp(p_interaction)
end
parglmo.p = [p_factor p_interaction];
