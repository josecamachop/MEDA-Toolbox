function paranovao = paranova(X, F,varargin)

% Parallel ANOVA to obtain multivariate factor and interaction matrices in
% a crossed experimental design and permutation test for significance. This
% approach permutes the raw values, which is sub-optimal in terms of power 
% according to Andreson and Ter Braak. This routine is superseded by parglm.m
% (please, use the latter)
%
% paranovao = paranova(X, F)   % minimum call
% paranovao = paranova(X, F, 'Interactions',interactions, 'Preprocessing',center, 'Permutations',n_perm)   % complete call
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, where columns correspond to factors and rows to
% levels. Levels start at 1 and should be correlative.
%
% Optional INPUTS (parameters):
%
% 'Interactions': [Ix2] matrix where rows contain the factors for which
% interactions are to be calculated.
%
% 'Preprocessing': [1x1] preprocesing:
%       1: mean centering
%       2: autoscaling (default)
%
% 'Permutations': [1x1] number of permutations (1000 by default).
%
%
% OUTPUTS:
%
% paranovao (structure): structure with the factor and interaction
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
% F = create_design(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,8) + repmat(randn(1,vars),reps,1);
%     end
% end
% 
% paranovao = paranova(X, F);

%
% EXAMPLE OF USE: Random data, two factors, with 4 and 3 levels, but only
%   the first one is significative, and 4 replicates:
%
% reps = 4;
% vars = 400;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = create_design(levels,'replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,'LevelCorr',8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
% 
% paranovao = paranova(X, F);
%
%
% coded by: Gooitzen Zwanenburg (G.Zwanenburg@uva.nl)
%           José Camacho (josecamacho@ugr.es)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024  Gooitzen Zwanenburg, University of Amsterdam
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

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Interactions',[]); 
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);     
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
interactions = p.Results.Interactions;
center = p.Results.Preprocessing;
n_perm = p.Results.Permutations;

% Validate dimensions of input data
assert (isequal(size(center), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

size_data           = size(X);                   % size of the data matrix
n_levels            = max(F);                    % number of levels / factor
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
factors             = 1 : n_factors;             % factors to be evaluated
X_raw               = X;                         % Input data matrix
X_level_means       = cell(n_factors,1);         % level_means per factor
SSQ_factors         = zeros(n_perm + 1,n_factors,1);      % sum of squares for factors
X_interaction_means = cell(n_interactions);      % cell_means for interactions
SSQ_interactions    = zeros(n_perm + 1,n_interactions);   % sum of squares for interactions
X_permuted          = cell(n_factors,1);         % permuted data per factor
p_factor            = zeros(1, n_factors);       % p-values factors
p_interaction       = zeros(1, n_interactions);  % p-values interactions

% In column space
paranovao.factors                = cell(n_factors,1);
paranovao.interactions           = cell(n_interactions,1);

% center/standardize the data
if center == 1
    Mean = ones(size_data(1),1)*mean(X_raw);        % Overall mean
    X = (X_raw - Mean);                             % Center
    SSQ_mean = sum(sum(Mean.^2));                   % SSQ overall means
    SSQ_X = sum(sum(X_raw.^2));                     % Sum of squares data matrix
elseif center == 2
    stdm = std(X_raw);
    stdm(find(stdm==0)) = 1;
    Mean_std = ones(size_data(1),1)*mean(X_raw)./...
        (ones(size_data(1),1)*stdm);
    X_std = X_raw./(ones(size_data(1),1)*stdm);
    X = (X_std - Mean_std);                         % Standardize
    SSQ_mean = sum(sum(Mean_std.^2));               % SSQ overall means
    SSQ_X = sum(sum(X_std.^2));
end
X_residuals         = X;                            % initial matrix with residuals

% Make structure with unchanging 'variables'
paranovao.data           = X;
paranovao.design         = F;
paranovao.n_factors      = n_factors;
paranovao.n_levels       = n_levels;
paranovao.n_interactions = n_interactions;

% Collect level means for the factors indicated in the model
for factor = factors
    X_level_means{factor} = level_means(X, paranovao, factor);
    SSQ_factors(1,factor) = sum(sum(X_level_means{factor}.^2));
    X_residuals = X_residuals - X_level_means{factor};
    paranovao.factors{factor}.matrix = X_level_means{factor};
end

X_residuals_afterF = X_residuals;

% Interactions
for i = 1 : n_interactions
    factor_1 = interactions(i,1);
    factor_2 = interactions(i,2);
    % Calculate cell means
    X_interaction_means{i} = zeros(size(X));
    for level_factor_1 = 1 : n_levels(factor_1)          % levels for first factor
        for level_factor_2 = 1 : n_levels(factor_2)      % levels for second factor
            tmp = zeros(size_data(1),1);
            found = find((F(:,factor_2) == level_factor_2) & ...  % find rows
                (F(:,factor_1) == level_factor_1));
            m = mean(X_residuals(found,:),1);                         % average over cell
            tmp(found) = 1;
            X_interaction_means{i} = X_interaction_means{i} + tmp*m;
        end
    end
    paranovao.interactions{i}.matrix = X_interaction_means{i};
    SSQ_interactions(1,i) = sum(sum(X_interaction_means{i}.^2));
    X_residuals = X_residuals - X_interaction_means{i};
end

SSQ_residuals = sum(sum(X_residuals.^2));
perc_effects = effect_explains(SSQ_X, SSQ_mean, SSQ_factors(1,:), ...
    SSQ_interactions(1,:), SSQ_residuals);
paranovao.effects = perc_effects;
paranovao.residuals = X_residuals;

disp('Percentage each effect contributes to the total sum of squares')
disp('Overall means')
disp(perc_effects(1))
disp('Factors')
disp(perc_effects(1 + (1 : n_factors)))
if n_interactions>0
    disp('Interactions')
    disp(perc_effects(1 + n_factors + (1 : n_interactions)))
end
disp('Residuals')
disp(perc_effects(end))

% Interactions p-values are calculated through permutation of X - Xa - Xa
% Do the permutations (do this, for now, for two factors)
for j = 1 : n_perm
    
    perms = randperm(size(X,1)); % permuted data (permute whole data matrix)
        
    for factor = 1 : n_factors
        % Level means for each factor
        X_level_means{factor} = level_means(X(perms, :), paranovao, factor);
        SSQ_factors(1 + j, factor) = sum(sum( (X_level_means{factor}).^2));
    end
    
    % Permutations for interactions.
    for i = 1 : n_interactions
        factor_1 = interactions(i,1);
        factor_2 = interactions(i,2);
        
        X_interaction_means{i} = zeros(size(X));
        for level_factor_1 = 1 : n_levels(factor_1)                   % levels for first factor
            for level_factor_2 = 1 : n_levels(factor_2)               % levels for second factor
                tmp = zeros(size_data(1),1);
                found = find((F(:,factor_2) == level_factor_2) & ...  % find rows
                    (F(:,factor_1) == level_factor_1));
                m = mean(X_residuals_afterF(perms(found),:));                       % average over cell
                tmp(found) = 1;
                X_interaction_means{i} = X_interaction_means{i} + tmp*m;
            end
        end
        SSQ_interactions(1 + j, i) = sum(sum( (X_interaction_means{i}).^2));
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
paranovao.p = [p_factor p_interaction];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function perc_explained = pc_explains(sv)
        
        % Function to calculate the percentage of variance explained by each PC
        
        % Input
        % vector with singular values
        
        % Output
        % vector with percentages explained variance
        
        sv_squared     = sv.^2;
        total_variance = sum(sv_squared);
        perc_explained = (sv_squared/total_variance)*100;
    end

    function perc_explained_effect = effect_explains(ssq_X, ssq_mean, ssq_factors, ...
            ssq_interactions, ssq_residuals)
        
        % Function to calculate the percentage of variance explained by
        % each effect.
        
        % Input
        % sums of squares of the data matrix, the factors and the
        % interactions
        
        % Output
        % vector with percentages explained variance for each effect.
        
        ssq = [ssq_mean ssq_factors ssq_interactions ssq_residuals];
        perc_explained_effect = 100*(ssq./ssq_X);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function M = level_means(Y, Input, factor)
        % Calculate level means for factor in matrix Y
        size_Y = size(Y);
        M = zeros(size_Y);
        for level = 1 : Input.n_levels(factor)
            tmp = zeros(size_Y(1),1);
            found = find(Input.design(:,factor) == level);  % find rows that belong to level
            m = mean(Y(found,:));                      % calculate level mean
            tmp(found) = 1;                            % flag the rows found
            M = M + tmp*m;
        end
    end

end


 

    