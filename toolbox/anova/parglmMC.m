function [T, parglmo, tsFactors, tsInteractions, SSQX, SSQinter, SSQFactorsT, SSQInteractionsT, SSQresidualsT] = parglmMC(X, F, varargin)

% Parallel General Linear Model to factorize in multivariate factor and 
% interaction matrices in an experimental design and permutation test for 
% univariate statistical significance with mutiple-test correction. This 
% routine provides the computation of FDR and q-values. 
%
% T = parglmMC(X, F)   % minimum call
%
% See also: parglm, parglmVS, asca, apca, createDesign
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
% 'Mtc': [1x1] Multiple test correction
%       1: Bonferroni 
%       2: Holm step-up or Hochberg step-down
%       3: Benjamini-Hochberg step-down (FDR)
%       -pi0: Q-value assuming 0<pi0<=1 the ratio of null variables (by default -1, assuming all variables are null, and so equivalen to the FDR more conservative)  
% 
% 'Fmtc': [1x1] correct for multiple-tesis when multifactorial (multi-way)
% analysis.
%       0: do not correct (by default)
%       1: same as mtc 
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
% parglmoMC (structure): structure with the factor and interaction
% matrices, p-values (corrected, depending on mtc and fmtc) and explained 
% variance  
%
% tsFactors: [nPerm*M+1 x nFactors x M] test statistic for factors
%
% tsInteractions: [nPerm*M+1 x nInteractions x M] test statistic for interactions
%
% SSQX, SSQInter, SSQFactorsT, SSQInteractionsT, SSQresidualsT: sums of
% squares in the data and permutations. These outputs are for internal use
% of other routines, and we suggest not to include them in the output.
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor, high background correlation. 
%
% nObs = 40;
% nVars = 100;
% 
% class = (randn(nObs,1)>0)+1;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% X(class==2,1:3) = X(class==2,1:3) + 10;
% X = X + 100*ones(size(class,1),1)*rand(1,nVars); 
% 
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
% rng(S);
% [TMC, parglmoMC] = parglmMC(X, class, 'Permutations', 500); % With variable selection through multiple testing correction
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot(-log10(parglmoMC.p(parglmoMC.ordFactors)),'k-')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','BH-FDR','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, no variable with information 
% on the factor. 
%
% nObs = 40;
% nVars = 100;
% 
% class = (randn(nObs,1)>0)+1;
% X = simuleMV(nObs,nVars,'LevelCorr',8);
% X = X + 100*ones(size(class,1),1)*rand(1,nVars);
% 
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
% rng(S);
% [TMC, parglmoMC] = parglmMC(X, class, 'Permutations', 500); % With variable selection through multiple testing correction
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot(-log10(parglmoMC.p(parglmoMC.ordFactors)),'k-')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','BH-FDR','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
% Random data, one factor and two levels, three variables with information 
% on the factor, low backgournd correlation. 
%
% nObs = 40;
% nVars = 100;
% 
% X = simuleMV(nObs,nVars,'LevelCorr',4);
% class = (sum(X(:,1:3),2)>0)+1;
% X = X + 100*ones(size(class,1),1)*rand(1,nVars);
% 
% S = rng; % Use same seed for random generators to improve comparability of results
% [T, parglmo] = parglm(X,class); % No variables selection 
% rng(S);
% [TVS, parglmoVS] = parglmVS(X, class); % With multivariate variable selection
% rng(S);
% [TMC, parglmoMC] = parglmMC(X, class, 'Permutations', 500); % With variable selection through multiple testing correction
% 
% h = figure; hold on
% plot([1 nVars],-log10([parglmo.p parglmo.p]),'b-.')
% plot(-log10(parglmoVS.p(parglmoVS.ordFactors)),'g-o')
% plot(-log10(parglmoMC.p(parglmoMC.ordFactors)),'k-')
% plot([0,size(X,2)],-log10([0.01 0.01]),'r--')
% legend('ASCA','VASCA','BH-FDR','alpha=0.01','Location','southeast')
% a=get(h,'CurrentAxes');
% set(a,'FontSize',14)
% ylabel('-log10(p-values)','FontSize',18)
% xlabel('Responses in selected order','FontSize',18)
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 31/Jan/2026
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
addParameter(p,'Mtc',-1); 
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
mtc = p.Results.Mtc;
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
assert (isequal(size(mtc), [1 1]), 'Dimension Error: parameter ''Mtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);

%% Main code

nInteractions      = prod(size(interactions));              % number of interactions
if fmtc
    mtcc = M*(nFactors + nInteractions);        % correction for the number of tests
else
    mtcc = M;
end
tsFactors         = zeros(nPerm*mtcc+1,nFactors,M);       % test statistic for factors
tsInteractions    = zeros(nPerm*mtcc+1,nInteractions,M);  % test statistic for interactions
if nargout > 4
    storeSSQ            = true;
    SSQFactorsT         = zeros(nPerm*mtcc+1,nFactors,M);       % sum of squares for factors
    SSQInteractionsT    = zeros(nPerm*mtcc+1,nInteractions,M);  % sum of squares for interactions
    SSQresidualsT       = zeros(nPerm*mtcc+1,M);                 % sum of squares for residuals
else
    storeSSQ            = false;
end
FFactors           = zeros(nFactors,M);       % F-value 
FInteractions      = zeros(nInteractions,M);  % F-value 
pFactor            = zeros(nFactors,M);                % p-values factors
pInteraction       = zeros(nInteractions,M);           % p-values interactions

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
parglmo.mtc             = mtc;
parglmo.fmtc            = fmtc;
parglmo.coding          = coding;
parglmo.nested          = nested; 
parglmo.type          = type; 

% Create Design Matrix
if prep
    n = 1;
    D = ones(size(X,1),1);
    offset = 1;
else
    n = 0;
    D = [];
    offset = 0;
end

for f = 1 : nFactors
    parglmo.factors{f}.factors = [];
    parglmo.factors{f}.refF = [];
    parglmo.factors{f}.refI = [];
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
        [~,ia] = intersect(parglmo.interactions{i2}.factors,parglmo.interactions{i}.factors);
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

SSQX = sum(X.^2);

if type == "Sequential", rep = 2; else, rep = 1; end

% GLM model calibration with LS
while rep > 0
    pD =  pinv(D'*D)*D';
    B = pD*X;
    Xresiduals = X - D*B;
    parglmo.D = D;
    parglmo.B = B;
    
    % Create Effect Matrices
    if prep
        parglmo.inter = D(:,1)*B(1,:);
        SSQinter = sum(parglmo.inter.^2);
        if rep == 1 && sum(SSQinter < 0.5*SSQX) > 0.5*M 
            disp('Warning: average with less than 50% of SS in more than half of the responses, consider not to mean center.');
        end
    else
        parglmo.inter = 0;
        SSQinter = zeros(1,M);
    end    
    SSQresiduals(1,:) = sum(Xresiduals.^2);
    
    SSQFactors = [];
    for f = 1 : nFactors
        parglmo.factors{f}.matrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQFactors(f,:) = sum(parglmo.factors{f}.matrix.^2);
    end
    
    % Interactions
    SSQInteractions = [];
    for i = 1 : nInteractions
        parglmo.interactions{i}.matrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQInteractions(i,:) = sum(parglmo.interactions{i}.matrix.^2);
    end

    if rep > 1 % Sequential
        if prep
            effects = sum([SSQinter' SSQFactors' SSQInteractions']);
        else
            effects = sum([SSQFactors' SSQInteractions']);
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


if storeSSQ
    SSQFactorsT(1,:,:) = SSQFactors;      
    SSQInteractionsT(1,:,:) = SSQInteractions;  % sum of squares for interactions
    SSQresidualsT(1,:,:) = SSQresiduals;
end
   
MSq = SSQresiduals(1,:)/Rdf;
for f = 1 : nFactors
    if ts==2 % reference MSE
        SSref = zeros(1,M);
        Dfref = 0;
        for f2 = 1 : nFactors % Nested
            if ~isempty(find(f==parglmo.factors{f2}.factors & random(f2) == 1)) 
                SSref = SSref + squeeze(SSQFactors(f2,:));
                Dfref = Dfref + df(f2);
                parglmo.factors{f}.refF = [parglmo.factors{f}.refF f2]; 
            end
        end
        for i = 1 : nInteractions
            if ~isempty(find(f==parglmo.interactions{i}.factors))
                rest = setdiff(parglmo.interactions{i}.factors,f);
                if (sum(random(rest) == 1) > 0) 
                    SSref = SSref + squeeze(SSQInteractions(i,:));
                    Dfref = Dfref + dfint(i);
                    parglmo.factors{f}.refI = [parglmo.factors{f}.refI i];
                end
            end
        end
        if isempty(find(SSref ~= 0))
            FFactors(f,:) = squeeze(SSQFactors(f,:)/df(f))./MSq;
        else
            MSS = max(MSq,SSref/Dfref); % Chose the most restrictive reference variable-wise
            FFactors(f,:) = squeeze(SSQFactors(f,:)/df(f))./MSS;
        end
    else
        FFactors(f,:) = squeeze(SSQFactors(f,:)/df(f))./MSq;
    end
end

for i = 1 : nInteractions
    if ts==2 % reference MSE
        SSref = zeros(1,M);
        Dfref = 0;
        for i2 = 1 : nInteractions
            [~,ia,ib] = intersect(parglmo.interactions{i}.factors,parglmo.interactions{i2}.factors);
            if (length(ia) == length(parglmo.interactions{i}.factors) & length(parglmo.interactions{i2}.factors) > length(parglmo.interactions{i}.factors))
                rest = setdiff(parglmo.interactions{i2}.factors,parglmo.interactions{i}.factors);
                if (sum(random(parglmo.interactions{i2}.factors(rest)) == 1) > 0) 
                    SSref = SSref + squeeze(SSQInteractions(i2,:));
                    Dfref = Dfref + dfint(i2);
                    parglmo.interactions{i}.refI = [parglmo.interactions{i}.refI i2];
                end
            end
        end
        if SSref == 0
            FInteractions(i,:) = squeeze(SSQInteractions(i,:)/dfint(i))./MSq;
        else
            MSS = max(MSq,SSref/Dfref); % Chose the most restrictive reference variable-wise
            FInteractions(i,:) = squeeze(SSQInteractions(i,:)/dfint(i))./MSS;
        end
    else
        FInteractions(i,:) = squeeze(SSQInteractions(i,:)/dfint(i))./MSq;
    end
end


if nInteractions
    parglmo.effects = 100*([SSQFactors' SSQInteractions' SSQresiduals']./((SSQX'-SSQinter')*ones(1,offset+nFactors+nInteractions)));
    par = 100*(sum([SSQFactors' SSQInteractions' SSQresiduals'])./sum(SSQX'-SSQinter'));
else
    parglmo.effects = 100*([SSQFactors' SSQresiduals']./((SSQX'-SSQinter')*ones(1,offset+nFactors+nInteractions)));
    par = 100*(sum([SSQFactors' SSQresiduals'])./sum(SSQX'-SSQinter'));
end
parglmo.residuals = Xresiduals;

% Select test statistic to order variables
if ts
    tsFactors(1,:,:) = FFactors;
    tsInteractions(1,:,:) = FInteractions;
else
    tsFactors(1,:,:) = SSQFactors;
    tsInteractions(1,:,:) = SSQInteractions;
end
    
% Permutations
parfor j = 1 : (nPerm * mtcc) % Increase the number of permutations to perform MTC
    
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
    Xresiduals = X - D*B;
    SSQresidualsP = sum(Xresiduals.^2);
    
    % Factors
    SSQFactorsP = [];
    for f = 1 : nFactors
        factorMatrix = D(:,parglmo.factors{f}.Dvars)*B(parglmo.factors{f}.Dvars,:);
        SSQFactorsP(f,:) = sum(factorMatrix.^2);
    end
    
    % Interactions
    SSQInteractionsP = [];
    for i = 1 : nInteractions
        interactsMatrix = D(:,parglmo.interactions{i}.Dvars)*B(parglmo.interactions{i}.Dvars,:);
        SSQInteractionsP(i,:) = sum(interactsMatrix.^2);
    end

    if storeSSQ
        SSQFactorsT(1 + j,:,:) = SSQFactorsP;      
        SSQInteractionsT(1 + j,:,:) = SSQInteractionsP;  % sum of squares for interactions
        SSQresidualsT(1 + j,:,:) = SSQresidualsP;
    end
   
    Ff = zeros(nFactors,M);
    for f = 1 : nFactors
        if ts==2 % reference MSE
            SSref = zeros(1,M);
            Dfref = 0;
            for f2 = 1 : nFactors
                if ~isempty(find(f2==parglmo.factors{f}.refF))
                    SSref = SSref + SSQFactorsP(f2,:);
                    Dfref = Dfref + df(f2);
                end
            end

            for i = 1 : nInteractions
                if ~isempty(find(i==parglmo.factors{f}.refI))
                    SSref = SSref + SSQInteractionsP(i,:);
                    Dfref = Dfref + dfint(i);
                end
            end

            if isempty(find(SSref ~= 0))
                Ff(f,:) = (SSQFactorsP(f,:)/df(f))./(SSQresidualsP/Rdf);
            else
                MSS = max(SSQresidualsP/Rdf,SSref(1,:)/Dfref); % Chose the most restrictive reference variable-wise
                Ff(f,:) = (SSQFactorsP(f,:)/df(f))./MSS;
            end
        else
            Ff(f,:) = (SSQFactorsP(f,:)/df(f))./(SSQresidualsP/Rdf);
        end
    end

    % F Interactions
    Fi = zeros(nInteractions,M);
    for i = 1 : nInteractions
        if ts==2 % reference MSE
            SSref = zeros(1,M);
            Dfref = 0;
            for i2 = 1 : nInteractions
                if ~isempty(find(i2==parglmo.interactions{i}.refI))
                    SSref = SSref + SSQInteractionsP(i2,:);
                    Dfref = Dfref + dfint(i2);
                end
            end
            if isempty(find(SSref ~= 0))
                Fi(i,:) = (SSQInteractionsP(i,:)/dfint(i))./(SSQresidualsP/Rdf);
            else
                MSS = max(SSQresidualsP/Rdf,SSref(1,:)/Dfref); % Chose the most restrictive reference variable-wise
                Fi(i,:)= (SSQInteractionsP(i,:)/dfint(i))./MMS;
            end
        else
            Fi(i,:) = (SSQInteractionsP(i,:)/dfint(i))./(SSQresidualsP/Rdf);
        end
    end


    % Select test statistic to order variables
    if ts
        tsFactors(1 + j,:,:) = Ff;
        tsInteractions(1 + j,:,:) = Fi;
    else
        tsFactors(1 + j,:,:) = SSQFactorsP;
        tsInteractions(1 + j,:,:) = SSQInteractionsP;
    end
end        
   

% Calculate univariate p-values and order variables by relevance
for f = 1 : nFactors
    for var = 1 : M
        pFactor(f,var) = (size(find(tsFactors(2:end,f,var) ...
            >= tsFactors(1,f,var)) ,1) + 1)/(nPerm*mtcc + 1);
    end
    [~,ordFactors(f,:)] = sort(pFactor(f,:),'ascend');
end
for i = 1 : nInteractions
    for var = 1 : M
        pInteraction(i,var) = (size(find(tsInteractions(2:end,i,var) ...
            >= tsInteractions(1,i,var)) ,1) + 1)/(nPerm*mtcc + 1);
    end
    [~,ordInteractions(i,:)] = sort(pInteraction(i,:),'ascend');
end

parglmo.ordFactors = ordFactors;
if nInteractions>0
    parglmo.ordInteractions = ordInteractions;
end

parglmo.p = [pFactor' pInteraction'];

% Multiple test correction
switch mtc
    
    case 1 % Bonferroni
        parglmo.p = min(1,parglmo.p * mtcc);
        
    case 2 % Holm/Hochberg
        
        if fmtc
            parglmo.p = parglmo.p;
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = 1 : length(indx)
                parglmo.p(indx(ind)) = min(1,parglmo.p(indx(ind)) * (mtcc-ind+1));
            end
        else
            for f = 1 : nFactors
                for var = 1 : M
                    pFactor(f,ordFactors(f,var)) = min(1,pFactor(f,ordFactors(f,var)) * (M-var+1));
                end
            end
            for i = 1 : nInteractions
                for var = 1 : M
                    pInteraction(i,ordInteractions(i,var)) = min(1,pInteraction(i,ordInteractions(i,var)) * (M-var+1));
                end
            end
            parglmo.p = [pFactor' pInteraction'];
        end
        
    case 3 % Benjamini & Hochberg
        
        if fmtc
            [~,indx] = sort(parglmo.p(:),'ascend');
            for ind = length(indx)-1 : -1 : 1
                parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * mtcc/ind]);
            end
        else
            for f = 1 : nFactors
                for var = M-1 : -1 : 1
                    pFactor(f,ordFactors(f,var)) = min([1,pFactor(f,ordFactors(f,var)) * M/var]);
                end
            end
            for i = 1 : nInteractions
                for var = M-1 : -1 : 1
                    pInteraction(i,ordInteractions(i,var)) = min([1,pInteraction(i,ordInteractions(i,var)) * M/var]);
                end
            end
            parglmo.p = [pFactor' pInteraction'];
        end
        
end

if mtc <0 % Q-value     
    if fmtc
        parglmo.p = parglmo.p;
        [~,indx] = sort(parglmo.p(:),'ascend');
        parglmo.p(indx(end)) = parglmo.p(indx(end))*(-mtc);
        for ind = length(indx)-1 : -1 : 1
            parglmo.p(indx(ind)) = min([1,parglmo.p(indx(ind)) * (-mtc) * mtcc/ind,parglmo.p(indx(ind+1))]);
        end
    else
        for f = 1 : nFactors
            pFactor(f,ordFactors(f,M)) = pFactor(f,ordFactors(f,M))*(-mtc);
            for var = M-1 : -1 : 1
                pFactor(f,ordFactors(f,var)) = min([1,pFactor(f,ordFactors(f,var)) * (-mtc) * M/var,pFactor(f,ordFactors(f,var+1))]);
            end
        end
        for i = 1 : nInteractions
            pInteraction(i,ordInteractions(i,M)) = pInteraction(i,ordInteractions(i,M))*(-mtc);
            for var = M-1 : -1 : 1
                pInteraction(i,ordInteractions(i,var)) = min([1,pInteraction(i,ordInteractions(i,var)) * (-mtc) * M/var,pInteraction(i,ordInteractions(i,var+1))]);
            end
        end
        parglmo.p = [pFactor' pInteraction'];
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
      
if nInteractions
    SSQ = sum([SSQFactors' SSQInteractions' SSQresiduals' (SSQX-SSQinter)'],1);
else
    SSQ = sum([SSQFactors' SSQresiduals' (SSQX-SSQinter)'],1);
end
par = [par sum(par)];
DoF = [df dfint Rdf Tdf];
MSQ = SSQ./DoF;
F = [max(FFactors,[],2)' max(FInteractions,[],2)' nan nan];
pValue = [min(parglmo.p) nan nan];

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    T.data = [name'; SSQ'; par'; DoF'; MSQ'; F'; pValue'];
    T.labels = {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'};
else
    T = table(name', SSQ', par', DoF', MSQ', F', pValue','VariableNames', {'Source','SumSq','AvPercSumSq','df','MeanSq','MaxF','minPvalue'});
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

