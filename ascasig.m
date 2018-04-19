function p = ascasig(X, F, interactions, center, n_perm)

% Does permutation test on the factors of a crossed experimental design
%
% p = ascasig(X, F)   % minimum call
% p = ascasig(X, F, interactions, center, n_perm)   % complete call
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
% p: [1x1] p_value.
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
% p = ascasig(X, F);
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
% p = ascasig(X, F);
%
%
% coded by: Gooitzen Zwanenburg (G.Zwanenburg@uva.nl)
% last modification: 19/Apr/18.
%
% Copyright (C) 2018  Gooitzen Zwanenburg, University of Amsterdam
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

size_data           = size(X);                              % size of the data matrix
n_levels            = max(F);                               % number of levels / factor
n_interactions      = size(interactions,1);                 % number of interactions
n_factors           = size(F,2);                            % number of factors
factors             = 1 : n_factors;                        % factors to be evaluated
X_level_means       = cell(n_factors,1);                    % level_means per factor
X_permuted          = cell(n_factors,1);                    % permuted data per factor
X_interaction_means = cell(n_interactions);                 % cell_means X_int - X_a - X_b
ssq                 = zeros(n_perm + 1,n_factors);          % within sum of squares data matrix
ssq_interaction     = zeros(n_perm + 1,n_interactions);     % within sum of squares data matrix
p_factor            = zeros(1, n_factors);                  % p-values factors
p_interaction       = zeros(1, n_interactions);             % p-values interactions

% center/standardize the data
if center == 1
    X = (X - ones(size_data(1),1)*mean(X));      % Centre
elseif center == 2
    X = (X - ones(size_data(1),1)*mean(X))./...  % Standardize
        (ones(size_data(1),1)*std(X));
end

% Make structure with unchanging 'variables'
Fixed.F              = F;
Fixed.interactions   = interactions;
Fixed.n_factors      = n_factors;
Fixed.n_levels       = n_levels;
Fixed.n_interactions = n_interactions;

% Collect level means for the factors indicated in the model
for factor = factors
    X_level_means{factor} = level_means(X, Fixed, factor);
    
    % sum of squares of level averages of measured data
    ssq(1, factor) = sum(sum( (X_level_means{factor}).^2));
end

for i = 1 : n_interactions
    factor_1 = interactions(i,1);
    factor_2 = interactions(i,2);
    % X_residuals = X - Xa - Xb
    X_residuals = X - X_level_means{factor_1} - X_level_means{factor_2};
    % Calculate cell means
    
    X_interaction_means{i} = zeros(size(X));
    for level_factor_1 = 1 : n_levels(factor_1)          % levels for first factor
        for level_factor_2 = 1 : n_levels(factor_2)      % levels for second factor
            tmp = zeros(size_data(1),1);
            found = find((F(:,factor_2) == level_factor_2) & ...  % find rows
                (F(:,factor_1) == level_factor_1));
            m = mean(X_residuals(found,:));                         % average over cell
            tmp(found) = 1;
            X_interaction_means{i} = X_interaction_means{i} + tmp*m;
        end
    end
    ssq_interaction(1, i) = sum(sum( (X_interaction_means{i}).^2));
end


% Interactions p-values are calculated through permutation of X - Xa - Xa
% Do the permutations (do this, for now, for two factors)
for j = 1 : n_perm
    
    % Permutations of experimental factors
    % Restrict permutations to the factor considered
    reverse_factors = n_factors : -1 : 1;
    for factor = 1 : n_factors
        other_factor       = reverse_factors(factor);
        X_permuted{factor} = zeros(size(X));
        for level = 1 : n_levels(other_factor)
            found = find(F(:,other_factor) == level);
            n_rows = length(found);                               % number of rows for level
            perms  = randperm(n_rows);                            % random permutations of n_rows
            permuted_level = found(perms);
            X_permuted{factor}(found,:) = X(permuted_level,:);
        end
        
        % Level means for each factor
        X_level_means{factor} = level_means(X_permuted{factor}, Fixed, factor);
        
        ssq(1 + j, factor) = sum(sum( (X_level_means{factor}).^2));
    end
    
    % Permutations for interactions.
    for i = 1 : n_interactions
        factor_1 = interactions(i,1);
        factor_2 = interactions(i,2);
        perms    = randperm(size(X,1));
        X_r      = X(perms, :);          % permuted data (permute whole data matrix)
        
        % Calculate level means experimental factors for permuted data
        for factor = [factor_1 factor_2]
            X_level_means{factor} = level_means(X_r, Fixed, factor);
        end
        
        % X_residuals = Xr - Xa - Xb
        X_residuals = X_r - X_level_means{factor_1} - X_level_means{factor_2};
        X_interaction_means{i} = zeros(size(X));
        for level_factor_1 = 1 : n_levels(factor_1)                   % levels for first factor
            for level_factor_2 = 1 : n_levels(factor_2)               % levels for second factor
                tmp = zeros(size_data(1),1);
                found = find((F(:,factor_2) == level_factor_2) & ...  % find rows
                    (F(:,factor_1) == level_factor_1));
                m = mean(X_residuals(found,:));                       % average over cell
                tmp(found) = 1;
                X_interaction_means{i} = X_interaction_means{i} + tmp*m;
            end
        end
        ssq_interaction(1 + j, i) = sum(sum( (X_interaction_means{i}).^2));
    end
    
end        % permutations

% Calculate p-values
% how many ssq's are larger than measurement ssq?
for factor = 1 : n_factors
    p_factor(factor) = (size(find(ssq(2:n_perm + 1, factor) >= ssq(1, factor)),1) + 1)/(n_perm);
end
for interaction = 1 : n_interactions
    p_interaction(interaction) = (size(find(ssq_interaction(2:n_perm + 1, interaction) ...
        >= ssq_interaction(1, interaction)),1) + 1)/(n_perm);
end
disp('p-values factors:')
disp (p_factor)
disp('p-values intractions')
disp(p_interaction)
p = [p_factor p_interaction];


end            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = level_means(Y, Input, factor)
% Calculate level means for factor in matrix Y
size_Y = size(Y);
M = zeros(size_Y);
for level = 1 : Input.n_levels(factor)
    tmp = zeros(size_Y(1),1);
    found = find(Input.F(:,factor) == level);  % find rows that belong to level
    m = mean(Y(found,:));                      % calculate level mean
    tmp(found) = 1;                            % flag the rows found
    M = M + tmp*m;
end
end
 

    