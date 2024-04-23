function [PCmean, PCrep, powercurveo] = powercurve(X, F, varargin)

% ASCA Power Curves (PCs). We derive four different types of PCs organized 
% in two taxonomies. On the one hand, we distinguish Population PCs from 
% Sample PCs, where the former are derived from population parameters and 
% the latter from descriptive statistics of a specific sample of data. On 
% the other hand, we distinguish Relative PCs from Absolute PCs, where the 
% former represent statistical power in terms of the relative effect size 
% between structure and noise, and the latter in terms of the sample size. 
% All possible four combinations can be derived. 
%
% Related routines: parglm, asca, apca, parglmVS, parglmMC, create_design
%
% PCmean = powercurve(X, F)   % minimum call
% [PCmean, PCrep] = powercurve(X, F, 'Model',model, 'Type',type, 'Repetitions',n_rep, 'RandomGen',randg, 'RandomGenC',randgC, 'Theta',theta, 'ALpha',alpha,'Preprocessing',prep, 'Permutations',n_perm, 'Ts',ts, 'Ordinal',ordinal, 'Fmtc',fmtc, 'Coding',coding, nested, 'Replicates',replicates)   % complete call
%
%
% INPUTS:
%
% X: can be one of the following
%
%   - [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable, to generate Sample PCs
% 
%   - (struct) to generate Population PCs
%       X.N: number of rows
%       X.M: number of columns
%       X.k: standard deviation coefficients
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels
%
% Option INPUTS (parameters):
%
% 'Model': This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
%
% 'Type': [1x1] type of power curve
%   - 1: Relative power curves
%   - 2: Absoute power curves
%
% 'Repetitions': [1x1] number of repetitions to compute the power curves (1000 by default)
%
% 'RandomGen': (func) random generator (@randn by default) suggested alternatives
%   - @(N,M)simuleMV(N,M,8): multivariate correlated with level 8 (other values may be used)
%   - @rand: uniform (moderately non-normal)
%   - @(N,M)exprnd(1,N,M).^3: very non-normal 
%
% 'RandomGenC': (func) random generator in effect size coefficients (@()0.1*randn+1 by default)
%
% 'Theta': [1xT] For type equal to 1, theta controls the compromise of 
%   true significance vs random (0:0.1:1 by default). For type equal to 2, 
%   theta controls the number of replicates (1:10 by default)
%
% 'Alpha': [1x1] significance level (0.01 by defult)
%
% 'Preprocessing': [1x1] preprocesing:
%       0: no preprocessing 
%       1: mean-centering 
%       2: auto-scaling (default)
%
% 'Permutations': [1x1] number of permutations (1000 by default)
%
% 'Ts': [1x1] Use SSQ (0) or the F-value (otherwise, by default) as test statistic  
%       0: Sum-of-squares of the factor/interaction
%       1: F-ratio of the SS of the factor/interaction divided by the SS of 
%       the residuals (by default)
%       2: F-ratio following the factors/interactions hierarchy
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
% 
% 'Fmtc': [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis
%       0: do not correct (default)
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR)
%       4: Q-value from Benjamini-Hochberg step-down
%
% 'Coding': [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
% 'Nested': [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
% 'Replicates': [1x1] index of the factor with replicates (only used for type
% 3), 0 by default, meaning no factor with replicates
%
%
% OUTPUTS:
%
% PCmean: [1xT] Average PC for the n_rep repetitions
%
% PCrep: [n_repxT] PC for eacg repetition
%
% powercurveo (structure): structure with general information 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Simulated Power Curve (type 1) from effect coefficients, two factors 
%   with 4 and 3 levels, and 4 replicates, with interaction. The relative 
%   effect (standard deviation) is .1 and .2 for the first and second 
%   factor and .3 for the interaction, respectively. 
%
% reps = 4;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = create_design(levels,'Replicates',reps);
% 
% X.N = size(F,1);
% X.M = 400;
% X.k = [.1,.2,.3];
% 
% PCmean = powercurve(X, F, 'Model',{[1 2]},'Type',1,'Repetitions',200)
% legend('Factor A','Factor B','Interaction')
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 23/Apr/24
%
% Copyright (C) 2024  Universidad de Granada
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

n_factors = size(F,2);                 % number of factors

if isstruct(X)
    N = X.N;
    M = X.M;
    coeff = X.k;
else
    N = size(X, 1);
    M = size(X, 2);
end

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
    if isstruct(X), tip = 1; 
    else, tip = 2;
    end
addParameter(p,'Type',tip);   
addParameter(p,'Model','linear');
addParameter(p,'RandomGen',@randn);
addParameter(p,'Repetitions',1000);
addParameter(p,'RamdonGenC',@()0.1*randn+1);
    if tip == 2
        THeta = 1:10; 
    else
        THeta = 0:0.1:1; 
    end
addParameter(p,'Theta',THeta);
addParameter(p,'Alpha',0.01);
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);
addParameter(p,'Ts',1);
addParameter(p,'Ordinal',zeros(1,size(F,2)));
addParameter(p,'Fmtc',0);
addParameter(p,'Coding',zeros(1,size(F,2)));
addParameter(p,'Nested',[]);
addParameter(p,'Replicates',0);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
type = p.Results.Type;
n_rep = p.Results.Repetitions;
model = p.Results.Model;
randg = p.Results.RandomGen;
randgC = p.Results.RamdonGenC;
theta = p.Results.Theta;
alpha = p.Results.Alpha;
prep = p.Results.Preprocessing;
n_perm = p.Results.Permutations;
ts = p.Results.Ts;
ordinal = p.Results.Ordinal;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;
replicates = p.Results.Replicates;
theta = sort(theta,'ascend');


if isequal(model,'linear')
    interactions = [];
end  
    
f = 1:n_factors;
if ~isempty(nested), f(nested(:,2)) = []; end
if isequal(model,'interaction')
    interactions = allinter(f,2);
end    

if isequal(model,'full')
    interactions = allinter(f,length(f));
end    

if isnumeric(model) && isscalar(model) && model >= 2 && model <= n_factors
        interactions = allinter(f,model);
end    

if isnumeric(model) && ~isscalar(model)
    for i = 1:size(model,1)
        interactions{i} = model(i,:);
    end
end    

if iscell(model), interactions = model; end

% Validate dimensions of input data
assert (isequal(size(type), [1 1]), 'Dimension Error: parameter ''Type'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_rep), [1 1]), 'Dimension Error: parameter ''Repetitions'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(alpha), [1 1]), 'Dimension Error: parameter ''Alpha'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(n_perm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(replicates), [1 1]), 'Dimension Error: parameter ''Replicates'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
n_interactions      = length(interactions);      % number of interactions
n_factors           = size(F,2);                 % number of factors
if fmtc
    mtcc                = n_factors + n_interactions;        % correction for the number of tests
else
    mtcc = 1;
end

% Make structure with general 'variables'
powercurveo.data           = X;
powercurveo.prep           = prep;
powercurveo.design         = F;
powercurveo.n_factors      = n_factors;
powercurveo.n_interactions = n_interactions;
powercurveo.n_perm         = n_perm;
powercurveo.ts             = ts;
powercurveo.ordinal        = ordinal;
powercurveo.fmtc           = fmtc;
powercurveo.coding         = coding;
powercurveo.nested         = nested;
powercurveo.n_rep          = n_rep;
powercurveo.theta          = theta;
powercurveo.alpha          = alpha;
powercurveo.randg          = randg;
powercurveo.randgC         = randgC;
powercurveo.type           = type;
powercurveo.replicates     = replicates;

% Create Design Matrix

n = 1; % This is necessary for indirect orthogonalization, to avoid leakage of variance to the residuals
D = ones(size(X,1),1);

powercurveo.Dvars = [0];

for f = 1 : n_factors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),'Preprocessing',1);
        powercurveo.factors{f}.Dvars = n+1;
        n = n + 1;
        powercurveo.factors{f}.order = 1;
        powercurveo.Dvars(n+1) = 1;
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            powercurveo.factors{f}.factors = [];
            uF = unique(F(:,f));
            powercurveo.n_levels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            powercurveo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            powercurveo.factors{f}.order = 1;
            powercurveo.Dvars(powercurveo.factors{f}.Dvars) = powercurveo.factors{f}.order;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            powercurveo.factors{f}.factors = [ref powercurveo.factors{ref}.factors];
            urF = unique(F(:,ref));
            powercurveo.n_levels(f) = 0;
            powercurveo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                powercurveo.n_levels(f) = powercurveo.n_levels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                powercurveo.factors{f}.Dvars = [powercurveo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            powercurveo.factors{f}.order = powercurveo.factors{ref}.order + 1;
            powercurveo.Dvars(powercurveo.factors{f}.Dvars) = powercurveo.factors{f}.order;
        end
    end
end

for i = 1 : n_interactions
    Dout = computaDint(interactions{i},powercurveo.factors,D);
    D = [D Dout];
    powercurveo.interactions{i}.Dvars = n+1:size(D,2);
    powercurveo.interactions{i}.factors = interactions{i};
    n = size(D,2);
    powercurveo.interactions{i}.order = max(powercurveo.factors{interactions{i}(1)}.order,powercurveo.factors{interactions{i}(2)}.order) + 1;
    powercurveo.Dvars(powercurveo.interactions{i}.Dvars) = powercurveo.interactions{i}.order;
end

pD =  pinv(D'*D)*D';
    
% Degrees of freedom
Tdf = N;      
if prep
    mdf = 1;
else
    mdf = 0;
end
Rdf = Tdf-mdf;
for f = 1 : n_factors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(powercurveo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
df_int = [];
for i = 1 : n_interactions
    df_int(i) = prod(df(powercurveo.interactions{i}.factors));
    Rdf = Rdf-df_int(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end

if ~isstruct(X) % Sample PCs
    
    % preprocess the data
    [Xs,m,dt] = preprocess2D(X,'Preprocessing',prep);
    X = X./(ones(size(X,1),1)*dt);

    % Handle missing data 
    [r,c]=find(isnan(X));
    Xnan = X;
    ru = unique(r);
    for i=1:length(ru)
        ind = find(r==ru(i));
        ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
        for j=1:length(c(ind))
            ind3 = find(isnan(X(ind2,c(ind(j)))));
            if length(ind2)>length(ind3)
                X(r(ind(j)),c(ind(j))) = nanmean(X(ind2,c(ind(j)))); % use conditional mean replacement
            else
                X(r(ind(j)),c(ind(j))) = nanmean(X(:,c(ind(j)))); % use unconditional mean replacement if CMR not possible
            end
        end
    end
    SSQ_X = sum(sum(X.^2));

    % GLM model calibration with LS, only fixed factors
    B = pD*X;
    X_residuals = X - D*B;
    powercurveo.D = D;
    powercurveo.B = B;

    % Create Effect Matrices
    if prep
        parglmo.inter = D(:,1)*B(1,:);
    else
        parglmo.inter = 0;
    end   

    for f = 1 : n_factors
        powercurveo.factors{f}.matrix = D(:,powercurveo.factors{f}.Dvars)*B(powercurveo.factors{f}.Dvars,:);
    end
    
    for i = 1 : n_interactions
        powercurveo.interactions{i}.matrix = D(:,powercurveo.interactions{i}.Dvars)*B(powercurveo.interactions{i}.Dvars,:);
    end
    
    powercurveo.coeffs = ones(1,n_factors+n_interactions);
    powercurveo.rescoef = 1;
      
else % Population PCs
    powercurveo.coeffs = X.k;
    powercurveo.rescoef = 1;
end


%% Compute the Power Curve
 
eD = zeros(length(theta),length(powercurveo.coeffs),n_rep);

F2 = F;        
for i2=1:n_rep
    
    %disp(i2)
    
    rng(i2);
    
    if type == 1 % Relative PCs
        
        if isstruct(X) 
            for f = 1 : n_factors
                if ordinal(f)
                    powercurveo.factors{f}.matrix = randg(N,M); 
                    powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');   
                else 
                    powercurveo.factors{f}.matrix = randg(powercurveo.n_levels(f),M);
                    powercurveo.factors{f}.matrix = sqrt(powercurveo.n_levels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro'); 

                    if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
                        powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
                    else % nested
                        mati = powercurveo.factors{f}.matrix;
                        Fi = F(:,[powercurveo.factors{f}.factors f]);
                        Fi(:,2:end) = Fi(:,2:end)-1;
                        Li = powercurveo.n_levels([powercurveo.factors{f}.factors f]);
                        powercurveo.factors{f}.matrix = [];
                        for n = 1: N
                            powercurveo.factors{f}.matrix(n,:) = mati(Fi(n,:)*[1 Li(1:end-1)]',:);
                        end
                    end
                end            
            end

            for i = 1 : n_interactions
                mati = randg(prod(powercurveo.n_levels(powercurveo.interactions{i}.factors)),M);
                mati = sqrt(size(mati,1))*mati/norm(mati,'fro'); 
                Fi = F(:,powercurveo.interactions{i}.factors);
                Fi(:,2:end) = Fi(:,2:end)-1;
                Li = powercurveo.n_levels(powercurveo.interactions{i}.factors);
                powercurveo.interactions{i}.matrix = [];
                for n = 1: N
                    powercurveo.interactions{i}.matrix(n,:) = mati(Fi(n,:)*[1 Li(1:end-1)]',:);
                end
            end
        end

        Xstruct = zeros(N,M);
        for f = 1 : n_factors
            if ~ordinal(f)
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
            end
        end

        for i = 1 : n_interactions
            Xstruct = Xstruct + randgC() * powercurveo.coeffs(n_factors+i) * powercurveo.interactions{i}.matrix;
        end

        Xnoise = randg(N,M); 
        Xnoise = randgC() * powercurveo.rescoef * sqrt(N)*Xnoise/norm(Xnoise,'fro');
        
        for a = 1:length(theta)

            Xm = (1-theta(a))*Xnoise;
            Xm = Xm + theta(a)*Xstruct;

            for f = 1 : n_factors
                if ordinal(f)
                    Xm = Xm + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
                end
            end
            
            % Parallel GLM
            [T, parglmo] = parglm(Xm, F, 'Model',model, 'Preprocessing',prep, 'Permutations',n_perm, 'Ts',ts, 'Ordinal',ordinal, 'Fmtc',fmtc, 'Coding',coding, 'Nested',nested);
            
            powercurveo.T{i2,a} = T;
            
            for o = 1:length(powercurveo.coeffs)
                if parglmo.p(o) <= alpha
                    eD(a,o,i2) = 1;
                end
            end
        end
        
     else % Absolute PCs
         
        for a = 1:length(theta) 
         
            if replicates>0
                f = replicates;
                uF = unique(F2(:,f));

                F = F2;
                if theta(a) < length(uF)
                    F(find(F2(:,f)>theta(a)),:) = [];
                else
                    for t = 1:theta(a)-length(uF)
                        Fi = F2(find(F2(:,f)==uF(1)),:);
                        Fi(:,f) = max(uF) + t;
                        F = [F;Fi];
                    end
                end
                
                for f = 1 : n_factors
                    if ~ordinal(f)
                        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
                            uF = unique(F(:,f));
                            powercurveo.n_levels(f) = length(uF);
                        else % if nested
                            ind = find(nested(:,2)==f);
                            ref = nested(ind,1);
                            powercurveo.factors{f}.factors = [ref powercurveo.factors{ref}.factors];
                            urF = unique(F(:,ref));
                            powercurveo.n_levels(f) = 0;
                            for j = 1:length(urF)
                                rind = find(ismember(F(:,ref),urF(j)));
                                uF = unique(F(rind,f));
                                powercurveo.n_levels(f) = powercurveo.n_levels(f) + length(uF);
                            end
                        end
                    end
                end
                
            else
                F = repmat(F2,theta(a),1);
            end
            
            N = size(F,1);
        
            if isstruct(X)
                for f = 1 : n_factors
                    if ordinal(f)
                        powercurveo.factors{f}.matrix = randg(N,M);
                        powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
                    else
                        powercurveo.factors{f}.matrix = randg(powercurveo.n_levels(f),M);
                        powercurveo.factors{f}.matrix = sqrt(powercurveo.n_levels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
                        
                        if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
                            powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
                        else % nested
                            mati = powercurveo.factors{f}.matrix;
                            Fi = F(:,[powercurveo.factors{f}.factors f]);
                            Fi(:,2:end) = Fi(:,2:end)-1;
                            Li = powercurveo.n_levels([powercurveo.factors{f}.factors f]);
                            powercurveo.factors{f}.matrix = [];
                            for n = 1: N
                                powercurveo.factors{f}.matrix(n,:) = mati(Fi(n,:)*[1 Li(1:end-1)]',:);
                            end
                        end
                    end
                end
                
                for i = 1 : n_interactions
                    mati = randg(prod(powercurveo.n_levels(powercurveo.interactions{i}.factors)),M);
                    mati = sqrt(size(mati,1))*mati/norm(mati,'fro');
                    Fi = F(:,powercurveo.interactions{i}.factors);
                    Fi(:,2:end) = Fi(:,2:end)-1;
                    Li = powercurveo.n_levels(powercurveo.interactions{i}.factors);
                    powercurveo.interactions{i}.matrix = [];
                    for n = 1: N
                        powercurveo.interactions{i}.matrix(n,:) = mati(Fi(n,:)*[1 Li(1:end-1)]',:);
                    end
                end
            end

            Xstruct = zeros(N,M);
            for f = 1 : n_factors
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
            end
            
            for i = 1 : n_interactions
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(n_factors+i) * powercurveo.interactions{i}.matrix;
            end
            
            Xnoise = randg(N,M);
            Xnoise = randgC() * powercurveo.rescoef * sqrt(N)*Xnoise/norm(Xnoise,'fro');

            Xm = Xnoise + Xstruct;
            
            % Parallel GLM
            [T, parglmo] = parglm(Xm, F, 'Model',model, 'preprocessing',prep, 'Permutations',n_perm, 'Ts',ts, 'Ordinal',ordinal, 'Fmtc',fmtc, 'Coding',coding, 'Nested',nested);
            
            powercurveo.T{i2,a} = T;
            
            for o = 1:length(powercurveo.coeffs)
                if parglmo.p(o) <= alpha
                    eD(a,o,i2) = 1;
                end
            end
        end    
             
    end
    
end

PCmean = mean(eD,3);
PCrep = eD; 

%% Show results

figure;
SPC_bootstrap(theta,PCrep,n_rep,true,0.05,false);    
%for o = 1:length(powercurveo.coeffs), plot(theta,PCmean(:,o)); end;
if type == 1 % Relative PCs
    
    xlabel('Effect size (\theta)','FontSize', 16);
    
    if isstruct(X)
        title('Relative Population Curve','FontSize', 16);
    else
        title('Relative Sample Curve','FontSize', 16);
    end
else
    xlabel('Number of replicates (\eta)','FontSize', 16);
    
    if isstruct(X)
        title('Absolute Population Curve','FontSize', 16);
    else
        title('Absolute Sample Curve','FontSize', 16);
    end
end
ylabel('Power','FontSize', 16);

end

%% Auxiliary function for interactions

function interactions = allinter(factors,order)
    
    if order > 2
        interactions = allinter(factors,order-1);
        for i = 1:length(interactions)
            for j = factors(find(factors > max(interactions{i})))
                interactions{end+1} = [interactions{i} j];
            end
        end
    else
        interactions = {};
        for i = factors
            for j = factors(find(factors >i))
                interactions{end+1} = [i j];
            end
        end
    end
    
end
    
        
function Dout = computaDint(interactions,factors,D) % Compute coding matrix

    if length(interactions)>1
        deepD = computaDint(interactions(2:end),factors,D);
        Dout = [];
        for k = factors{interactions(1)}.Dvars
            for l = 1:size(deepD,2)
                Dout(:,end+1) = D(:,k).* deepD(:,l);
            end
        end
    else
        Dout = D(:,factors{interactions}.Dvars);
    end

end









