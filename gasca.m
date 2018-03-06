function gascao = gasca(X, F, interactions, options)

% GASCA is a data analysis algorithm for designed experiments. It does a
% group-wise principal component analysis on the level averages of each
% experimental factor in a designed experiment with balanced data.
% Interactions between two factors can also be calculated. The original
% paper is Saccenti, E., Smilde, A.K. and Camacho, J. Group-wise ANOVA
% simultaneous component analysis for designed omics experiments.
% Submitted to Metabolomics, 2018.
%
% gascao = gasca(X, F)   % minimum call
% gascao = gasca(X, F, interactions, center)   % complete call
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
% options (structure):
%   options.center:         center == 1: center data;
%                           center == 2 standardize data (default)
%   options.corrtype:       to be specified when maptype corr is used.
%                           Availabale methods: those implmemented in function
%                           corr of Matlab (default:'Pearson')
%   options.significance:   if > 0 select correlation on the base of
%                           p-value. Correlation for which p-value >
%                           signficance are set to 0 (Default:0.01)
%   options.groupdis:       if true, a plot of the distribution of the
%                           groups is shown (default false)
%   options.ordered:        if true, loadings are reordered according to 
%                           groups (default true)
%
%
% OUTPUTS:
%
% gascao (structure): structure that contains scores, loadings, singular
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
% gascao = gasca(X, F);
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
% gascao = gasca(X, F);
%
%
% coded by: Edoardo Saccenti (edoardo.saccenti@wur.nl)
%           Jose Camacho Paez (josecamacho@ugr.es)
%           Gooitzen Zwanenburg (G.Zwanenburg@uva.nl)
% last modification: 6/Mar/18.
%
% Copyright (C) 2018  Edoardo Saccenti, Wageningen University
% Copyright (C) 2018  Jose Camacho Paez, University of Granada
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
if nargin < 4 || isempty(options) || ~isfield(options,'center'), options.center = 2; end;
if ~isfield(options,'corrtype'), options.corrtype = 'Pearson'; end;
if ~isfield(options,'significance'), options.significance = 0.01; end;
if ~isfield(options,'groupdis'), options.groupdis = false; end;
if ~isfield(options,'ordered'), options.ordered = true; end;

% Validate dimensions of input data
assert (isequal(size(options.center), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code

size_data           = size(X);                   % size of the data matrix
n_levels            = max(F);                    % number of levels / factor
n_interactions      = size(interactions,1);      % number of interactions
n_factors           = size(F,2);                 % number of factors
factors             = 1 : n_factors;             % factors to be evaluated
X_raw               = X;                         % Input data matrix
X_level_means       = cell(n_factors,1);         % level_means per factor
SSQ_factors         = zeros(1,n_factors,1);      % sum of squares for factors
X_interaction_means = cell(n_interactions);      % cell_means for interactions
SSQ_interactions    = zeros(1,n_interactions);   % sum of squares for interactions

% Structure with results
gascao.factors.scores               = cell(n_factors,1);
gascao.factors.loadings             = cell(n_factors,1);
gascao.factors.projected            = cell(n_factors,1);
gascao.factors.singular             = cell(n_factors,1);
gascao.factors.explained            = cell(n_factors,1);
gascao.interactions.scores          = cell(n_interactions);
gascao.interactions.loadings        = cell(n_interactions);
gascao.interactions.singular        = cell(n_interactions);
gascao.interactions.explained       = cell(n_interactions);

% In column space
gascao.factors.means                = cell(n_factors,1);
gascao.interactions.means           = cell(n_interactions);

% center/standardize the data
center = options.center;

if center == 1
    Mean = ones(size_data(1),1)*mean(X_raw);        % Overall mean
    X = (X_raw - Mean);                             % Center
    SSQ_mean = sum(sum(Mean.^2));                   % SSQ overall means
    SSQ_X = sum(sum(X_raw.^2));                     % Sum of squares data matrix
elseif center == 2
    Mean_std = ones(size_data(1),1)*mean(X_raw)./...
        (ones(size_data(1),1)*std(X_raw));
    X_std = X_raw./(ones(size_data(1),1)*std(X_raw));
    X = (X_std - Mean_std);                         % Standardize
    SSQ_mean = sum(sum(Mean_std.^2));               % SSQ overall means
    SSQ_X = sum(sum(X_std.^2));
end
X_residuals         = X;                            % initial matrix with residuals
gascao.data           = X;
gascao.design         = F;

% Collect level means for the factors indicated in the model
for factor = factors
    X_level_means{factor} = zeros(size_data);
    for level = 1 : n_levels(factor)
        tmp = zeros(size_data(1),1);
        found = find(F(:,factor) == level);      % find rows that belong to level
        m = mean(X(found,:));                    % calculate level mean
        tmp(found) = 1;                          % flag the rows found
        X_level_means{factor} = X_level_means{factor} + tmp*m;
    end
    SSQ_factors(factor) = sum(sum(X_level_means{factor}.^2));
    X_residuals = X_residuals - X_level_means{factor};
    gascao.factors.means{factor} = X_level_means{factor};
end

% Interactions
for i = 1 : n_interactions
    factor_1 = interactions(i,1);
    factor_2 = interactions(i,2);
    X_interaction_means{i} = zeros(size_data);
    for level_factor_1 = 1 : n_levels(factor_1)          % levels for first factor
        for level_factor_2 = 1 : n_levels(factor_2)      % levels for second factor
            tmp = zeros(size_data(1),1);
            found = find((F(:,factor_2) == level_factor_2) & ...  % find rows
                (F(:,factor_1) == level_factor_1));
            if size(found,1) == 1                        % only one subject/cell
                m = X(found,:);                          % average is row
            else
                m = mean(X(found,:));                    % average over cell
            end
            tmp(found) = 1;
            X_interaction_means{i} = X_interaction_means{i} + tmp*m;
        end
    end
    X_interaction_means{i} = X_interaction_means{i} - ...
        X_level_means{factor_1} - X_level_means{factor_2};
    SSQ_interactions(i) = sum(sum(X_interaction_means{i}.^2));
    gascao.interactions.means{i} = X_interaction_means{i};
    X_residuals = X_residuals - X_interaction_means{i};
end

SSQ_residuals = sum(sum(X_residuals.^2));
perc_effects = effect_explains(SSQ_X, SSQ_mean, SSQ_factors, ...
    SSQ_interactions, SSQ_residuals);
gascao.effects = perc_effects;
gascao.residuals = X_residuals;

%Do GPCA on level averages for each factor
for factor = 1 : n_factors
    figureindex = factor;
    [t,s,p,o] = do_gpca(X_level_means{factor},sprintf('factor %d',factor),options,figureindex);
    projected_data = (X_residuals*p);       % project residuals on loadings
    gascao.factors.scores{factor}    = t;
    gascao.factors.loadings{factor}  = p;
    gascao.factors.singular{factor}  = s;
    if options.ordered,
        gascao.factors.order{factor}  = o;
    else
        gascao.factors.order{factor}  = 1:size_data(2);
    end
    gascao.factors.projected{factor} = projected_data;
    gascao.factors.explained{factor} = pc_explains(s);
end

%Do GPCA on interactions
for interaction = 1 : n_interactions
    figureindex = 2*n_factors + 1;
    [t,s,p,o] = do_gpca(X_interaction_means{interaction},sprintf('interaction %d',interaction),options, figureindex);
    projected_data = (X_residuals*p);       % project residuals on loadings
    gascao.interactions.scores{interaction}    = t;
    gascao.interactions.loadings{interaction}  = p;
    gascao.interactions.singular{interaction}  = s;
    if options.ordered,
        gascao.interactions.order{factor}  = o;
    else
        gascao.interactions.order{factor}  = 1:size_data(2);
    end
    gascao.interactions.projected{interaction} = projected_data;
    gascao.interactions.explained{interaction} = pc_explains(s);
end

plot_gasca(gascao);
plot_interactions(gascao, n_interactions);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [scores, singular_values, loadings, ord] = do_gpca(D,label,options,figureindex)
        % This function does the work: do GPCA. Returns scores, loadings and
        % singular values.
        
        if options.groupdis,
            ExploreMap(D,label,options.corrtype,options.significance,figureindex);
        end
        
        [meda_map, pval] = corr(D,'type',options.corrtype);
        
        if  options.significance > 0
            meda_map(find(pval>options.significance)) = 0;
        end

        disp(['Visualizing relationship map for ', label,' ...']);

        [map, ord] = seriation(meda_map);
        
        plot_map(map);
        title(label)
        
        axis square;
        
        c = input(sprintf('Please, select c for %s: ',label));
        
        min_group_size = floor(sqrt(size(meda_map,2)));
        
        disp('Computing group selection...');
        
        [~,states,~] =  gia(meda_map,c,min_group_size);
        
        
        [loadings,scores] =  gpca(D,states,1:rank(D));
        singular_values(1) = trace(loadings(:,1)*scores(:,1)'*scores(:,1)*loadings(:,1)');
        for i=2:size(loadings, 2),
            singular_values(i) = trace(loadings(:,i)*scores(:,i)'*scores(:,i)*loadings(:,i)')-singular_values(i-1);
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_gasca(GASCA)
        
        % Function to plot the scores, loadings and projections for each
        % factor.
        
        % Input
        % strucuture GASCA
        
        % Output
        % none
        
        f_F              = GASCA.design;
        f_n_levels       = max(f_F);                    % number of levels / factor
        f_n_factors      = size(f_F,2);                 % number of factors
        legend_labels    = cell(max(f_n_levels),1);
        
        % Plot results for experimental factors
        for f_factor = 1 : f_n_factors
            loadings       = GASCA.factors.loadings{f_factor};
            ordering       = GASCA.factors.order{f_factor};
            scores         = GASCA.factors.scores{f_factor};
            projections    = GASCA.factors.projected{f_factor};
            plot_residuals = scores + projections;
            figure;
            if size(loadings,2)>1,
                bar(loadings(ordering,1:2));
                ortho = loadings(ordering,1)'*loadings(ordering,2) < 1e-5;
                legend('First loading', 'Second loading');
            else
                bar(loadings(ordering,1));
                legend('First loading');
            end
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
                if size(loadings,2)>1,
                    if ortho,
                        plot(scores(f_found,1), scores(f_found,2), clr, 'LineWidth', 3 );
                        plot(plot_residuals(f_found,1), plot_residuals(f_found,2), clr);
                        xlabel(sprintf('GPC 1 (%%%.1f)',GASCA.factors.explained{f_factor}(1)))
                        ylabel(sprintf('GPC 2 (%%%.1f)',GASCA.factors.explained{f_factor}(2)))
                    else
                        subplot(2,1,1)
                        plot(zeros(size(f_found)),scores(f_found,1), clr, 'LineWidth', 3 );
                        plot(f_found,plot_residuals(f_found,1), clr);
                        xlabel('Sample')
                        ylabel(sprintf('GPC 1 (%%%.1f)',GASCA.factors.explained{f_factor}(1)))
                        subplot(2,1,2)
                        plot(zeros(size(f_found)),scores(f_found,2), clr, 'LineWidth', 3 );
                        plot(f_found,plot_residuals(f_found,2), clr);
                        xlabel('Sample')
                        ylabel(sprintf('GPC 2 (%%%.1f)',GASCA.factors.explained{f_factor}(2)))
                    end
                else
                    plot(zeros(size(f_found)),scores(f_found,1), clr, 'LineWidth', 3 );
                    plot(f_found,plot_residuals(f_found,1), clr);
                    xlabel('Sample')
                    ylabel(sprintf('GPC 1 (%%%.1f)',GASCA.factors.explained{f_factor}(1)))
                end
                stitle = [ 'Score plot factor: ' int2str(f_factor)];
                title(stitle);
                legend(legend_labels,'Location','NorthEastOutside');
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_interactions(GASCA, i_interactions)
        
        % Function to plot the scores, loadings and projections for each
        % interaction.
        
        % Input
        % strucuture GASCA
        % number of interactions
        
        % Output
        % none
        
        % Plot results for experimental factors
        for i_inter = 1 : i_interactions
            loadings       = GASCA.interactions.loadings{i_inter};
            ordering       = GASCA.interactions.order{i_inter};
            scores         = GASCA.interactions.scores{i_inter};
            projections    = GASCA.interactions.projected{i_inter};
            figure;
            bar(loadings(ordering,1:2));
            legend({'First loading', 'Second loading'},'Location','NorthEastOutside');
            sload = ['Loadings interaction: ' int2str(i_inter)];
            title(sload);
            figure;
            hold on;
            plot(scores(:,1), scores(:,2), '*r', 'LineWidth', 3 );
            plot(projections(:,1), projections(:,2), 'ob');
            % We want the text shifted with respect to the points
            x_axis_size = xlim;
            y_axis_size = ylim;
            x_dist = abs(x_axis_size(1) - x_axis_size(2));
            y_dist = abs(y_axis_size(1) - y_axis_size(2));
            % Shift text 2% up and to the right
            text(scores(:,1) + 0.02*x_dist, scores(:,2) + 0.02*y_dist, ...
                int2str(GASCA.design), 'FontSize', 8 );
            
            
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function [h, states] = ExploreMap(S,label,corrtype,alpha,figureindex)
        
        disp(['Exploring the CORRELATION map with various values of gamma for ', label])
        
        min_group_size = floor(sqrt(size(S,2)));
        
        lenStates = [];
        
        [map4gia, pval] = corr(S,'type',corrtype);
        
        if ~isempty(alpha) < 1
            map4gia(pval>0.05) = 0;
        end
        
        C = [0.5:0.15:0.95];
        [~,~,stree] = gia(map4gia,C(1),min_group_size);
        
        for i = 1 : length(C);
            
            c = C(i);
            disp(sprintf('Running GIA with c = %2.2f',c));
            
            
            [~,states,~] = gia(map4gia,c,min_group_size,stree);
            
            lenStates(i) = length(states);
            
            NumStat = [];
            for jj = 1 : length(states);
                NumStat(jj) = length(states{1,jj});
            end
            
            NumStat1(i) = length(find(NumStat == 1));
            NumStat2(i) = length(find(NumStat == 2));
            NumStat3(i) = length(find(NumStat > 2));
            
            Averagenum(i) = median(NumStat(find(NumStat >= min_group_size)));
            
            if isnan(Averagenum(i)) == 1
                MINnum(i) = NaN;
            else
                MINnum(i) = min(NumStat(find(NumStat >= min_group_size)));
            end
            
            if isnan(Averagenum(i)) == 1
                MAXnum(i) = NaN;
            else
                MAXnum(i) = max(NumStat(find(NumStat >= min_group_size)));
            end
            
        end
        
        h = figure(figureindex);
        set(gca,'FontSize',15);
        hold on;
        
        plot(C,lenStates(1,:),'.-','MarkerSize',25);
        plot(C,Averagenum,'.-','MarkerSize',25);
        plot(C,MINnum,'.-','MarkerSize',25);
        plot(C,MAXnum,'.-','MarkerSize',25);
        
        
        xlabel('Threshold \gamma','FontSize',14);
        legend({'Total # Groups', sprintf('Median group size (>%i)',min_group_size-1),'Min group size','Max group size'},'Location','NorthEastOutside')
        grid on;
        box on;
        axis square;
        xlim([0.5,1]);
    end

end