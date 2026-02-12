function [T, parglmo] = parglm(X, F, varargin)

% Parallel General Linear Model to factorize in multivariate factor and 
% interaction matrices in an experimental design and permutation test for 
% multivariate statistical significance. These represent the two first 
% steps in an ANOVA Simultaneous Component Analsysis (ASCA) 
%
% T = parglm(X, F)   % minimum call
%
% See also: asca, apca, parglmVS, parglmMC, createDesign
%
%
% INPUTS:
%
% X: [NxM] billinear data set for model fitting, where each row is a
% measurement, each column a variable
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels.
%
%
% Optional INPUTS (parameters):
%
% 'Model': This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
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
%       the residuals
%       2: F-ratio following the factors/interactions hierarchy (see
%       Montogomery, by default)
%
% 'Ordinal': [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
%
% 'Random': [1xF] whether factors are fixed or random
%       0: fixed (default)
%       1: random
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
% 'Type': Type of ANOVA factorization
%   'Simultaneous': All factors at once (by default, check %SS)
%   'Sequential': Sequential, marginalizing in order of variance
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
% F = createDesign(levels,'Replicates',reps);
% 
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     X(find(F(:,1) == levels{1}(i)),:) = simuleMV(length(find(F(:,1) == levels{1}(i))),vars,'LevelCorr',8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(i))),1);
% end
% X = X + 100*ones(size(F,1),1)*rand(1,vars);
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
% F = createDesign(levels,'Replicates',reps);
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
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,'LevelCorr',8) + repmat(fi{i} + fj{j},reps,1);
%     end
% end
% X = X + 100*ones(size(F,1),1)*rand(1,vars); 
% 
% table = parglm(X, F, 'Model',{[1 2]})
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
% F = createDesign(levels,'Replicates',reps);
%
% X = zeros(size(F,1),vars);
% for i = 1:length(levels{1}),
%     for j = 1:length(levels{2}),
%         X(find(F(:,1) == levels{1}(i) & F(:,2) == levels{2}(j)),:) = simuleMV(reps,vars,'LevelCorr',8) + repmat(randn(1,vars),reps,1);
%     end
% end
% X = X + 100*ones(size(F,1),1)*rand(1,vars); 
%
% table = parglm(X, F, 'Model',{[1 2]})
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 29/Jan/2026
% Dependencies: Matlab R2017b, MEDA v1.10
%
% Copyright (C) 2026  University of Granada, Granada
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

nFactors = size(F,2);                 % number of factors

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;

addParameter(p,'Model','linear'); 
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);  
addParameter(p,'Ts',2); 
addParameter(p,'Ordinal',zeros(1,size(F,2))); 
addParameter(p,'Random',zeros(1,size(F,2))); 
addParameter(p,'Fmtc',0); 
addParameter(p,'Coding',zeros(1,size(F,2))); 
addParameter(p,'Nested',[]); 
addParameter(p,'Type','Simultaneous'); 
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
model = p.Results.Model;
prep = p.Results.Preprocessing;
nPerm = p.Results.Permutations;
ts = p.Results.Ts;
ordinal = p.Results.Ordinal;
random = p.Results.Random;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;
type = p.Results.Type;

if isempty(model), model = 'linear'; end

if isequal(model,'linear')
    interactions = [];
end  
    
f = 1:nFactors;
if ~isempty(nested), f(nested(:,2)) = []; end
if isequal(model,'interaction')
    interactions = allinter(f,2);
end    

if isequal(model,'full')
    interactions = allinter(f,length(f));
end    

if isnumeric(model) && isscalar(model) && model >= 2 && model <= nFactors
        interactions = allinter(f,model);
end    

if isnumeric(model) && ~isscalar(model)
    for i = 1:size(model,1)
        interactions{i} = model(i,:);
    end
end    

if iscell(model), interactions = model; end

% Validate dimensions of input data
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nPerm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(random), [1 size(F,2)]), 'Dimension Error: parameter ''Random'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
nInteractions      = prod(size(interactions));              % number of interactions
nFactors           = size(F,2);                 % number of factors
if fmtc
    mtcc                = nFactors + nInteractions;        % correction for the number of tests
else
    mtcc = 1;
end
SSQfactors         = zeros(nPerm*mtcc+1,nFactors,1);      % sum of squares for factors
SSQinteractions    = zeros(nPerm*mtcc+1,nInteractions);   % sum of squares for interactions
Ffactors           = zeros(nPerm*mtcc+1,nFactors,1);      % F for factors
Finteractions      = zeros(nPerm*mtcc+1,nInteractions);   % F for interactions
pFactor            = zeros(1,nFactors);       % p-values factors
pInteraction       = zeros(1,nInteractions);  % p-values interactions

% In column space
parglmo.factors                = cell(nFactors,1);
parglmo.interactions           = cell(nInteractions,1);

% preprocess the data: data is only scaled at this point
[Xs,m,dt] = preprocess2D(X,'Preprocessing',prep);
X = X./(ones(size(X,1),1)*dt);

% Make structure with parameters
parglmo.data           = X;
parglmo.prep           = prep;
parglmo.scale           = dt;
parglmo.design         = F;
parglmo.nFactors      = nFactors;
parglmo.nInteractions = nInteractions;
parglmo.nPerm          = nPerm;
parglmo.ts              = ts;
parglmo.ordinal         = ordinal;
parglmo.random         = random;
parglmo.fmtc            = fmtc;
parglmo.coding          = coding;
parglmo.nested          = nested; 
parglmo.type          = type; 

% Create Design Matrix: mean centering added in the desing or not
if prep
    n = 1;
    D = ones(size(X,1),1);
else
    n = 0;
    D = [];
end

for f = 1 : nFactors
    parglmo.factors{f}.factors = [];
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),'Preprocessing',1); % Should we normalize a covariate? Not relevant while using GLM
        parglmo.factors{f}.Dvars = n+1;
        n = n + 1;
        parglmo.factors{f}.order = 1;
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            uF = unique(F(:,f));
            parglmo.nLevels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            parglmo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),parglmo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            parglmo.factors{f}.order = 1;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            parglmo.factors{f}.factors = [ref parglmo.factors{ref}.factors]; % Careful! nested should be in order in index f
            urF = unique(F(:,ref));
            parglmo.nLevels(f) = 0;
            parglmo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                parglmo.nLevels(f) = parglmo.nLevels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                parglmo.factors{f}.Dvars = [parglmo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            parglmo.factors{f}.order = parglmo.factors{ref}.order + 1;
        end
    end
end

for i = 1 : nInteractions
    Dout = computaDint(interactions{i},parglmo.factors,D);
    D = [D Dout];
    parglmo.interactions{i}.Dvars = n+1:size(D,2);
    parglmo.interactions{i}.factors = interactions{i};
    parglmo.interactions{i}.type = sum(random(interactions{i}) == 1) > 0;
    n = size(D,2);
    
    Mm=0;
    for j=1:length(interactions{i})
        Mm = max(Mm,parglmo.factors{j}.order);
    end

    for i2 = 1 : i-1
        [~,ia,ib] = intersect(parglmo.interactions{i2}.factors,parglmo.interactions{i}.factors);
        if (length(ia) == length(parglmo.interactions{i2}.factors) & length(parglmo.interactions{i}.factors) > length(parglmo.interactions{i2}.factors))
             Mm = max(Mm,parglmo.interactions{i2}.order);
        end
    end

    parglmo.interactions{i}.order = Mm + 1;
    parglmo.interactions{i}.refI = [];
end

% Degrees of freedom
Tdf = size(X,1); 
if prep
    mdf = 1;
else
    mdf = 0;
end
Rdf = Tdf-mdf;
for f = 1 : nFactors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(parglmo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
dfint = [];
for i = 1 : nInteractions
    dfint(i) = prod(df(parglmo.interactions{i}.factors));
    Rdf = Rdf-dfint(i);
end
if Rdf < 0
    disp('Warning: degrees of freedom exhausted');
    return
end

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
            X(r(ind(j)),c(ind(j))) = mean(X(ind2,c(ind(j))), 'omitnan'); % use conditional mean replacement
        else
            X(r(ind(j)),c(ind(j))) = mean(X(:,c(ind(j))), 'omitnan'); % use unconditional mean replacement if CMR not possible
        end
    end
end
parglmo.data = X;
parglmo.Xnan = Xnan;
parglmo.df = df;
parglmo.dfint = dfint;
parglmo.Tdf = Tdf;
parglmo.Rdf = Rdf;

SSQX = sum(sum(X.^2));

if type == "Sequential", rep = 2; else, rep = 1; end

% GLM model calibration with LS
while rep > 0
    pD =  pinv(D'*D)*D';
    B = pD*X;
    Xresiduals  = X - D*B;
    parglmo.D = D;
    parglmo.B = B;
    
    % Create Effect Matrices
    if prep
        parglmo.inter = D(:,1)*B(1,:);
        SSQinter = sum(sum(parglmo.inter.^2));
        if rep == 1 && SSQinter < 0.5*SSQX 
            disp('Warning: average with less than 50% of SS, consider not to mean center.');
        end
        SSQXc = sum(sum((X-parglmo.inter).^2));
    else
        parglmo.inter = 0;
        SSQinter = 0;
    end    
    SSQresiduals = sum(sum(Xresiduals .^2));
    
    for f = 1 : nFactors
        parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQfactors(1,f) = sum(sum(parglmo.factors{f}.matrix.^2)); % Note: we are not using Type III sum of squares, and probably we should, although we did not find any difference in our experiments
    end
    
    for i = 1 : nInteractions
        parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQinteractions(1,i) = sum(sum(parglmo.interactions{i}.matrix.^2));
    end
    
    if rep > 1 % Sequential
        if prep
            offset = 1;
            effects = [SSQinter SSQfactors(1,:) SSQinteractions(1,:)];
        else
            offset = 0; 
            effects = [SSQfactors(1,:) SSQinteractions(1,:)];
        end
        [~,ord] = sort(effects,'descend');
        DN = [];
        for t = 1:length(ord)
            if offset && ord(t) == 1
                if isempty(DN)
                    D2(:,1) = D(:,1);
                else
                    D2(:,1) = (eye(size(D,1)) - (DN*pinv(DN'*DN)*DN')) * D(:,1);
                end
                DN = [DN D2(:,1)];
            elseif ord(t) < nFactors+1+offset
                if isempty(DN)
                    D2(:,parglmo.factors{ord(t)-offset}.Dvars) = D(:,parglmo.factors{ord(t)-offset}.Dvars);
                else
                    D2(:,parglmo.factors{ord(t)-offset}.Dvars) = (eye(size(D,1)) - (DN*pinv(DN'*DN)*DN')) * D(:,parglmo.factors{ord(t)-offset}.Dvars);
                end
                DN = [DN D2(:,parglmo.factors{ord(t)-offset}.Dvars)];
            else
                if isempty(DN)
                    D2(:,parglmo.interactions{ord(t)-nFactors-offset}.Dvars) = D(:,parglmo.interactions{ord(t)-nFactors-offset}.Dvars);
                else
                    D2(:,parglmo.interactions{ord(t)-nFactors-offset}.Dvars) = (eye(size(D,1)) - (DN*pinv(DN'*DN)*DN')) * D(:,parglmo.interactions{ord(t)-nFactors-offset}.Dvars); 
                end
                DN = [DN D2(:,parglmo.interactions{ord(t)-nFactors-offset}.Dvars)];
            end
        end
        D = D2;
    end
    rep = rep-1;
end

MSq = SSQresiduals/Rdf;
for f = 1 : nFactors
    parglmo.factors{f}.refF = [];
    parglmo.factors{f}.refI = [];
    if ts==2 % reference MSE
        SSref = 0;
        Dfref = 0;
        for f2 = 1 : nFactors % Nested
            if ~isempty(find(f==parglmo.factors{f2}.factors & random(f2) == 1)) && (MSq < SSQfactors(1,f2)/df(f2)) % when a nested factor is random and significantly larger than the background noise
                SSref = SSref + SSQfactors(1,f2);
                Dfref = Dfref + df(f2);
                parglmo.factors{f}.refF = [parglmo.factors{f}.refF f2]; 
            end
        end
        for i = 1 : nInteractions
            if ~isempty(find(f==parglmo.interactions{i}.factors))
                rest = setdiff(parglmo.interactions{i}.factors,f);
                if (sum(random(rest) == 1) > 0) && (MSq < SSQinteractions(1,i)/dfint(i)) % when an interaction is random and larger than the background noise
                    SSref = SSref + SSQinteractions(1,i);
                    Dfref = Dfref + dfint(i);
                    parglmo.factors{f}.refI = [parglmo.factors{f}.refI i];
                end
            end
        end
        if SSref == 0
            Ffactors(1,f) = (SSQfactors(1,f)/df(f))/MSq;
        else
            Ffactors(1,f) = (SSQfactors(1,f)/df(f))/(SSref/Dfref);
        end
    else
        Ffactors(1,f) = (SSQfactors(1,f)/df(f))/MSq;
    end
end

for i = 1 : nInteractions
    if ts==2 % reference MSE
        SSref = 0;
        Dfref = 0;
        for i2 = 1 : nInteractions
            [~,ia,ib] = intersect(parglmo.interactions{i}.factors,parglmo.interactions{i2}.factors);
            if (length(ia) == length(parglmo.interactions{i}.factors) & length(parglmo.interactions{i2}.factors) > length(parglmo.interactions{i}.factors))
                rest = setdiff(parglmo.interactions{i2}.factors,parglmo.interactions{i}.factors);
                if (sum(random(parglmo.interactions{i2}.factors(rest)) == 1) > 0) && (MSq < SSQinteractions(1,i2)/dfint(i2)) % when the higher order interaction is random and larger than the background noise
                    SSref = SSref + SSQinteractions(1,i2);
                    Dfref = Dfref + dfint(i2);
                    parglmo.interactions{i}.refI = [parglmo.interactions{i}.refI i2];
                end
            end
        end
        if SSref == 0
            Finteractions(1,i) = (SSQinteractions(1,i)/dfint(i))/MSq;
        else
            Finteractions(1,i) = (SSQinteractions(1,i)/dfint(i))/(SSref/Dfref);
        end
    else
        Finteractions(1,i) = (SSQinteractions(1,i)/dfint(i))/MSq;
    end
end

parglmo.effects = 100*([SSQfactors(1,:) SSQinteractions(1,:) SSQresiduals]./(SSQXc));
parglmo.residuals = Xresiduals;

% Permutations
parfor j = 1 : nPerm*mtcc
   
    perms = randperm(size(Xnan,1)); % permuted data (complete rows)
    
    X = Xnan(perms, :);
    [r,c]=find(isnan(X));
    ru = unique(r);
    for i=1:length(ru)
        ind = find(r==ru(i));
        ind2 = find(sum((D-ones(size(D,1),1)*D(r(ind(1)),:)).^2,2)==0);
        for f=1:length(c(ind))
            ind3 = find(isnan(X(ind2,c(ind(f)))));
            if length(ind2)>length(ind3)
                X(r(ind(f)),c(ind(f))) = mean(X(ind2,c(ind(f))), 'omitnan'); % use conditional mean replacement
            else
                X(r(ind(f)),c(ind(f))) = mean(X(:,c(ind(f))), 'omitnan'); % use unconditional mean replacement if CMR not possible
            end
        end
    end
      
    B = pD*X;
    Xresiduals  = X - D*B;
    SSQresidualsp = sum(sum(Xresiduals .^2));
    
    % Factors
    factors = cell(1,nFactors);
    SSQf = zeros(1,nFactors);
    for f = 1 : nFactors
        factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQf(f) = sum(sum(factors{f}.matrix.^2));
    end
    SSQfactors(1 + j,:) = SSQf;
    
    % Interactions
    interacts = {};
    SSQi = zeros(1,nInteractions);
    for i = 1 : nInteractions
        interacts{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);    
        SSQi(i) = sum(sum(interacts{i}.matrix.^2));
    end
    SSQinteractions(1 + j,:) = SSQi;
    
    % F Factors
    Ff = zeros(1,nFactors);
    for f = 1 : nFactors
        if ts==2% pooled reference MSE
            SSref = 0;
            Dfref = 0;
            for f2 = 1 : nFactors
                if ~isempty(find(f2==parglmo.factors{f}.refF))
                    SSref = SSref + SSQf(f2);
                    Dfref = Dfref + df(f2);
                end
            end
            
            for i = 1 : nInteractions
                if ~isempty(find(i==parglmo.factors{f}.refI))
                    SSref = SSref + SSQi(i);
                    Dfref = Dfref + dfint(i);
                end
            end
            if SSref == 0
                Ff(f) = (SSQf(f)/df(f))/(SSQresidualsp/Rdf);
            else
                Ff(f) = (SSQf(f)/df(f))/(SSref/Dfref);
            end
        else
            Ff(f) = (SSQf(f)/df(f))/(SSQresidualsp/Rdf);
        end
    end
    Ffactors(1 + j,:) = Ff; 

    % F Interactions
    Fi = zeros(1,nInteractions);
    for i = 1 : nInteractions
        if ts==2 % reference MSE
            SSref = 0;
            Dfref = 0;
            for i2 = 1 : nInteractions
                if ~isempty(find(i2==parglmo.interactions{i}.refI))
                    SSref = SSref + SSQi(1,i2);
                    Dfref = Dfref + dfint(i2);
                end
            end
            if SSref == 0
                Fi(i) = (SSQi(i)/dfint(i))/(SSQresidualsp/Rdf);
            else
                Fi(i)= (SSQi(i)/dfint(i))/(SSref/Dfref);
            end
        else
            Fi(i) = (SSQi(i)/dfint(i))/(SSQresidualsp/Rdf);
        end
    end
    Finteractions(1 + j,:) = Fi; 

end        

% Select test statistic
if ts
    tsFactors = Ffactors;
    tsInteractions = Finteractions;
else
    tsFactors = SSQfactors;
    tsInteractions = SSQinteractions;
end
    
% Calculate p-values
for f = 1 : nFactors
    pFactor(f) = (size(find(tsFactors(2:(nPerm*mtcc + 1), f) >= tsFactors(1, f)),1) + 1)/(nPerm*mtcc+1);
end
for i = 1 : nInteractions
    pInteraction(i) = (size(find(tsInteractions(2:(nPerm*mtcc + 1), i) ...
        >= tsInteractions(1, i)),1) + 1)/(nPerm*mtcc+1);
end

% Multiple test correction for several factors/interactions
parglmo.p = [pFactor pInteraction]; 
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

name={};
for f = 1 : nFactors
    name{end+1} = sprintf('Factor %d',f);
end
for i = 1 : nInteractions
    name{end+1} = sprintf('Interaction %s',strrep(num2str(parglmo.interactions{i}.factors),'  ','-'));
end
name{end+1} = 'Residuals';
name{end+1} = 'Total';
   

SSQ = [SSQfactors(1,:) SSQinteractions(1,:) SSQresiduals SSQXc];
par = [parglmo.effects sum(parglmo.effects)];
DoF = [df dfint Rdf Tdf-mdf];
MSQ = SSQ./DoF;
F = [Ffactors(1,:) Finteractions(1,:) nan nan];
pValue = [parglmo.p nan nan];

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
  T.mat = [SSQ', par', DoF', MSQ', F', pValue'];
  T.var = {'SumSq', 'PercSumSq', 'df', 'MeanSq', 'F', 'Pvalue'};
  T.source = name';
else
  T = table(name', SSQ', par', DoF', MSQ', F', pValue','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F','Pvalue'});
end

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

