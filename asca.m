function ascao = asca(paranovao, alpha)

% ASCA is a data analysis algorithm for designed experiments. It does a 
% principal component analysis on the level averages of each experimental 
% factor in a designed experiment with balanced data. Interactions between 
% two factors can also be calculated. The original paper for this software 
% is Zwanenburg, G, Hoefsloot, HCJ, Westerhuis, JA, Jansen, JJ, Smilde, AK.
% ANOVA–principal component analysis and ANOVA–simultaneous component 
% analysis: a comparison. Journal of Chemometrics, 2018, 25:561-567.
%
% ascao = asca(paranovao)   % minimum call
% ascao = asca(paranovao, alpha)   % complete call
%
%
% INPUTS:
%
% paranovao (structure): structure with the factor and interaction
% matrices, p-values and explained variance. Obtained with parallel anova
% (paranova)
%
% alpha: [1x1] significance level (0.05 by default)
%
%
% OUTPUTS:
%
% ascao (structure): structure that contains scores, loadings, singular
% values and projections of the factors and interactions.
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
% paranovao = paranova(X, F);
%
% ascao = asca(paranovao);
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
% paranovao = paranova(X, F);
%
% ascao = asca(paranovao);
%
%
% coded by: Gooitzen Zwanenburg (G.Zwanenburg@uva.nl)
%           José Camacho (josecamacho@ugr.es)
% last modification: 21/Mar/18
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 2 || isempty(alpha), alpha = 0.05; end;

% Validate dimensions of input data
assert (isequal(size(alpha), [1 1]), 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

ascao = paranovao;

% Structure with results
ascao.factors.scores               = cell(ascao.n_factors,1);
ascao.factors.loadings             = cell(ascao.n_factors,1);
ascao.factors.projected            = cell(ascao.n_factors,1);
ascao.factors.singular             = cell(ascao.n_factors,1);
ascao.interactions.scores          = cell(ascao.n_interactions);
ascao.interactions.loadings        = cell(ascao.n_interactions);
ascao.interactions.singular        = cell(ascao.n_interactions);


%Do PCA on level averages for each factor
for factor = 1 : ascao.n_factors
    [t,s,p] = do_pca(ascao.factors.means{factor});
    projected_data = (ascao.residuals*p);       % project residuals on loadings
    ascao.factors.scores{factor}    = t;
    ascao.factors.loadings{factor}  = p;
    ascao.factors.singular{factor}  = s;
    ascao.factors.projected{factor} = projected_data;
    ascao.factors.explained{factor} = pc_explains(s);
end

%Do PCA on interactions
for interaction = 1 : ascao.n_interactions
    [t,s,p] = do_pca(ascao.interactions.means{interaction});
    projected_data = (ascao.residuals*p);       % project residuals on loadings
    ascao.interactions.scores{interaction}    = t;
    ascao.interactions.loadings{interaction}  = p;
    ascao.interactions.singular{interaction}  = s;
    ascao.interactions.projected{interaction} = projected_data;
    ascao.interactions.explained{interaction} = pc_explains(s);
end

plot_asca(ascao);
plot_interactions(ascao, ascao.n_interactions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [scores, singular_values, loadings] = do_pca(D)
        % This function does the work: do PCA through singular value
        % decomposition on input matrix. Returns scores, loadings and
        % singular values.
        [u, sv, v] = svd(D);
        scores = u*sv;
        singular_values = diag(sv);
        loadings = v;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_asca(ascao)
        
        % Function to plot the scores, loadings and projections for each
        % factor.
        
        % Input
        % strucuture ascao
        
        % Output
        % none
        
        f_F              = ascao.design;
        f_n_levels       = max(f_F);                    % number of levels / factor
        f_n_factors      = size(f_F,2);                 % number of factors
        legend_labels    = cell(max(f_n_levels),1);
        
        % Plot results for experimental factors
        for f_factor = 1 : f_n_factors
            loadings       = ascao.factors.loadings{f_factor};
            scores         = ascao.factors.scores{f_factor};
            projections    = ascao.factors.projected{f_factor};
            plot_residuals = scores + projections;
            figure;
            bar(loadings(:,1:2));
            legend('First loading', 'Second loading');
            sload = ['Loadings factor: ' int2str(f_factor)];
            title(sload);
            figure;
            hold on;
            j = 1;                         % Count legends
            clear legend_labels;
            for ii = 1 : f_n_levels(f_factor)
                f_found = find(f_F(:,f_factor) == ii);
                if ii == 1
                    sym = '*';
                    co = 'b';
                    clr = [sym co];
                elseif ii == 2
                    sym = '*';
                    co = 'r';
                    clr = [sym co];
                elseif ii == 3
                    sym = '*';
                    co = 'g';
                    clr = [sym co];
                elseif ii == 4
                    sym = '*';
                    co = 'k';
                    clr = [sym co];
                elseif ii == 5
                    sym = '*';
                    co = 'y';
                    clr = [sym co];
                end
                legend_labels{j} = ['Score level ' int2str(ii)];
                legend_labels{j + 1} = ['Level ' int2str(ii)];
                j = j + 2;
                plot(scores(f_found,1), scores(f_found,2), clr, 'LineWidth', 3 );
                plot(plot_residuals(f_found,1), plot_residuals(f_found,2), clr);
                stitle = [ 'Score plot factor: ' int2str(f_factor)];
                xlabel(sprintf('PC 1 (%%%.1f)',ascao.factors.explained{f_factor}(1)))
                ylabel(sprintf('PC 2 (%%%.1f)',ascao.factors.explained{f_factor}(2)))
                title(stitle);
                legend(legend_labels);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_interactions(ascao, i_interactions)
        
        % Function to plot the scores, loadings and projections for each
        % interaction.
        
        % Input
        % strucuture ascao
        % number of interactions
        
        % Output
        % none
        
        % Plot results for experimental factors
        for i_inter = 1 : i_interactions
            loadings       = ascao.interactions.loadings{i_inter};
            scores         = ascao.interactions.scores{i_inter};
            projections    = ascao.interactions.projected{i_inter};
            figure;
            bar(loadings(:,1:2));
            legend('First loading', 'Second loading');
            sload = ['Loadings interaction: ' int2str(i_inter)];
            title(sload);
            figure;
            hold on
            plot(scores(:,1), scores(:,2), '*r', 'LineWidth', 3 );
            plot(projections(:,1), projections(:,2), 'ob');
            % We want the text shifted with respect to the points
            x_axis_size = xlim;
            y_axis_size = ylim;
            x_dist = abs(x_axis_size(1) - x_axis_size(2));
            y_dist = abs(y_axis_size(1) - y_axis_size(2));
            % Shift text 2% up and to the right
            text(scores(:,1) + 0.02*x_dist, scores(:,2) + 0.02*y_dist, ...
                int2str(ascao.design), 'FontSize', 8 );
            
            
            stitle = [ 'Score plot interaction: ' int2str(i_inter)];
            title(stitle);
        end
    end

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

end