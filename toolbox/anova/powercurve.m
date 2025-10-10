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
% PCmean = powercurve(X, F)   % minimum call
%
% See also: parglm, asca, apca, parglmVS, parglmMC, createDesign
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
%       X.N: [1x1] number of rows
%       X.M: [1x1] number of columns
%       X.k: [1x(F+I)] standard deviation coefficients
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels
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
% 'Type': [1x1] type of power curve
%   - 1: Relative power curves
%   - 2: Absoute power curves
%
% 'Repetitions': [1x1] number of repetitions to compute the power curves (1000 by default)
%
% 'RandomGen': (func or cell with F+I+1 functions) random generator (@randn by default) suggested alternatives
%   - @(N,M)simuleMV(N,M,'LevelCorr',8): multivariate correlated with level 8 (other values may be used)
%   - @rand: uniform (moderately non-normal)
%   - @(N,M)exprnd(1,N,M).^3: very non-normal
%   - {@randn,@rand}: normal for the single factor, uniform for the
%   residuals
%
% 'RandomGenC': (func) random generator in effect size coefficients (@()1 by default)
%    - To generate randomness use, e.g., @()0.1*randn+1 
%
% 'Theta': [1xT] For type equal to 1, theta controls the compromise of 
%   true significance vs random (0:0.1:1 by default). For type equal to 2, 
%   theta controls the number of replicates (1:10 by default)
%
% 'Alpha': [1x1] significance level (0.05 by defult)
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
% 'RepFactor': [1x1] index of the factor with replicates (only used for type
% 3), 0 by default, meaning no factor with replicates
%
%
% OUTPUTS:
%
% PCmean: [1xT] Average PC for the nRep repetitions
%
% PCrep: [nRepxT] PC for eacg repetition
%
% powercurveo (structure): structure with general information 
%
%
% EXAMPLE OF USE (copy and paste the code in the command line)
%   Simulated Power Curve (type 1) from effect coefficients, two factors 
%   with 4 and 3 levels, and 4 replicates, with interaction. The relative 
%   effect (standard deviation) is .1 and .2 for the first and second 
%   factor and .3 for the interaction, respectively. Please, note this
%   example takes several minutes to execute.
%
% reps = 4;
% levels = {[1,2,3,4],[1,2,3]};
% 
% F = createDesign(levels,'Replicates',reps);
% 
% X.N = size(F,1);
% X.M = 400;
% X.k = [.1,.2,.3];
% 
% PCmean = powercurve(X, F, 'Model',{[1 2]},'Type',1,'Repetitions',200);
% legend('Factor A','Factor B','Interaction')
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 25/Aug/2025
% Dependencies: Matlab R2017b, MEDA v1.10
%
% Copyright (C) 2025  University of Granada, Granada
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

nFactors = size(F,2);                 % number of factors

if isstruct(X)
    N = X.N;
    M = X.M;
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
addParameter(p,'RamdonGenC',@()1);
addParameter(p,'Theta',[]);
addParameter(p,'Alpha',0.05);
addParameter(p,'Preprocessing',2);
addParameter(p,'Permutations',1000);
addParameter(p,'Ts',2);
addParameter(p,'Ordinal',zeros(1,size(F,2)));
addParameter(p,'Random',zeros(1,size(F,2))); 
addParameter(p,'Fmtc',0);
addParameter(p,'Coding',zeros(1,size(F,2)));
addParameter(p,'Nested',[]);
addParameter(p,'RepFactor',0);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
type = p.Results.Type;
nRep = p.Results.Repetitions;
model = p.Results.Model;
randg = p.Results.RandomGen;
randgC = p.Results.RamdonGenC;
theta = p.Results.Theta;
alpha = p.Results.Alpha;
prep = p.Results.Preprocessing;
nPerm = p.Results.Permutations;
ts = p.Results.Ts;
ordinal = p.Results.Ordinal;
random = p.Results.Random;
fmtc = p.Results.Fmtc;
coding = p.Results.Coding;
nested = p.Results.Nested;
replicates = p.Results.RepFactor;

if isempty(theta)
    if type == 2
        if replicates>0
            theta = 2:10;
        else
            theta = 1:10;
        end
    else
        theta = 0:0.1:1; 
    end
end
theta = sort(theta,'ascend');

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
assert (isequal(size(type), [1 1]), 'Dimension Error: parameter ''Type'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nRep), [1 1]), 'Dimension Error: parameter ''Repetitions'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(alpha), [1 1]), 'Dimension Error: parameter ''Alpha'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(prep), [1 1]), 'Dimension Error: parameter ''Preprocessing'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(nPerm), [1 1]), 'Dimension Error: parameter ''Permutations'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ts), [1 1]), 'Dimension Error: parameter ''Ts'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ordinal), [1 size(F,2)]), 'Dimension Error: parameter ''Ordinal'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(random), [1 size(F,2)]), 'Dimension Error: parameter ''Random'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(fmtc), [1 1]), 'Dimension Error: parameter ''Fmtc'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(coding), [1 size(F,2)]), 'Dimension Error: parameter ''Coding'' must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(replicates), [1 1]), 'Dimension Error: parameter ''Replicates'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
nInteractions      = length(interactions);      % number of interactions
nFactors           = size(F,2);                 % number of factors

if isscalar(randg)
    randv = repmat({randg},1,nFactors+nInteractions+1);
else
    randv = randg;
end

if fmtc
    mtcc                = nFactors + nInteractions;        % correction for the number of tests
else
    mtcc = 1;
end

% Make structure with general 'variables'
powercurveo.data           = X;
powercurveo.prep           = prep;
powercurveo.design         = F;
powercurveo.nFactors      = nFactors;
powercurveo.nInteractions = nInteractions;
powercurveo.nPerm         = nPerm;
powercurveo.ts             = ts;
powercurveo.ordinal        = ordinal;
powercurveo.fmtc           = fmtc;
powercurveo.coding         = coding;
powercurveo.nested         = nested;
powercurveo.nRep          = nRep;
powercurveo.theta          = theta;
powercurveo.alpha          = alpha;
powercurveo.randg          = randg;
powercurveo.randgC         = randgC;
powercurveo.type           = type;
powercurveo.replicates     = replicates;


% Tranform Design Matrix

Fold = F;
F = zeros(size(F));
for f = 1:size(F,2)
    uF = unique(Fold(:,f),'stable');
    for i = 1:length(uF)
        F(find(ismember(Fold(:,f),uF(i))),f) = i;
    end
end

% Create Coding Matrix

n = 1; % This is necessary for indirect orthogonalization, to avoid leakage of variance to the residuals
D = ones(size(X,1),1);

powercurveo.Dvars = [0];

for f = 1 : nFactors
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
            powercurveo.nLevels(f) = length(uF);
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
            powercurveo.nLevels(f) = 0;
            powercurveo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                powercurveo.nLevels(f) = powercurveo.nLevels(f) + length(uF);
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

for i = 1 : nInteractions
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
for f = 1 : nFactors
    if ordinal(f)
        df(f) = 1;
    else
        df(f) = length(powercurveo.factors{f}.Dvars);
    end
    Rdf = Rdf-df(f);
end
dfint = [];
for i = 1 : nInteractions
    dfint(i) = prod(df(powercurveo.interactions{i}.factors));
    Rdf = Rdf-dfint(i);
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
                X(r(ind(j)),c(ind(j))) = mean(X(ind2,c(ind(j))), 'omitnan'); % use conditional mean replacement
            else
                X(r(ind(j)),c(ind(j))) = mean(X(:,c(ind(j))), 'omitnan'); % use unconditional mean replacement if CMR not possible
            end
        end
    end
    SSQX = sum(sum(X.^2));

    % GLM model calibration with LS, only fixed factors
    B = pD*X;
    Xresiduals = X - D*B;
    powercurveo.D = D;
    powercurveo.B = B;

    % Create Effect Matrices
    if prep
        parglmo.inter = D(:,1)*B(1,:);
        SSQinter = sum(sum(parglmo.inter.^2));
    else
        parglmo.inter = 0;
        SSQinter = 0;
    end   
    SSQresiduals = sum(sum(Xresiduals.^2));
    
    for f = 1 : nFactors 
        powercurveo.factors{f}.matrix = D(:,powercurveo.factors{f}.Dvars)*B(powercurveo.factors{f}.Dvars,:);
        SSFactors(1,f) = sum(sum(powercurveo.factors{f}.matrix.^2)); % Note: we are not using Type III sum of squares, and probably we should, although we did not find any difference in our experiments
        powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
    end
    
    for i = 1 : nInteractions
        powercurveo.interactions{i}.matrix = D(:,powercurveo.interactions{i}.Dvars)*B(powercurveo.interactions{i}.Dvars,:);
        SSInteractions(1,i) = sum(sum(powercurveo.interactions{i}.matrix.^2));
        powercurveo.interactions{i}.matrix = sqrt(N)*powercurveo.interactions{i}.matrix/norm(powercurveo.interactions{i}.matrix,'fro');
    end
    
    
    for f = 1:nFactors % Correct the SS, but only with what will be simulated (noise and replicating factor)
        SSref = 0;     % Important! in a balanced design (we assume here) low-order factors (ascendants) do not
        Dfref = 0;     % affect descendants, but the otherway does
        if replicates>0
            f2 = replicates; 
            if ~isempty(find(f==powercurveo.factors{f2}.factors))
                SSref = SSref + SSFactors(1,f2);
                Dfref = Dfref + df(f2);
            end
        end
        if SSref == 0
            SSFactorsC(1,f) = (SSFactors(1,f) - df(f)*SSQresiduals/Rdf)/(N/powercurveo.nLevels(f)); % SS corrected
        else
            SSFactorsC(1,f) = (SSFactors(1,f)/df(f) - SSref/Dfref)/(N/powercurveo.nLevels(f));
        end

        if SSFactorsC(1,f) < 0, SSFactorsC(1,f)=0; end 
    end   
    
    for i = 1 : nInteractions
        SSInteractionsC(1,i) = (SSInteractions(1,i)/dfint(i) - SSQresiduals/Rdf)/(N/prod(powercurveo.nLevels(powercurveo.interactions{i}.factors))); % SS corrected
        
        if SSInteractionsC(1,i) < 0, SSInteractionsC(1,i)=0; end
    end
    
    % Set coefficients
    powercurveo.coeffs = sqrt([SSFactorsC SSInteractionsC]);
    powercurveo.rescoef = sqrt(SSQresiduals/Rdf);
      
else % Population PCs
    powercurveo.coeffs = X.k;
    powercurveo.rescoef = 1;
end


%% Compute the Power Curve
 
eD = zeros(length(theta),length(powercurveo.coeffs),nRep);

F2 = F;   
powercurveo2 = powercurveo;
for i2=1:nRep
    
    %disp(i2)
    
    rng(i2);
    
    if type == 1 % Relative PCs
        
        if isstruct(X) 
            for f = 1 : nFactors
                randg = randv{f};
                if ordinal(f)
                    powercurveo.factors{f}.matrix = randg(N,M); 
                    powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');   
                else 
                    powercurveo.factors{f}.matrix = randg(powercurveo.nLevels(f),M);
                    powercurveo.factors{f}.matrix = sqrt(powercurveo.nLevels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro'); 

                    if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
                        powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
                    else % nested
                        mati = powercurveo.factors{f}.matrix;
                        Fi = F(:,[powercurveo.factors{f}.factors f]);
                        powercurveo.factors{f}.matrix = [];
                        uF = unique(Fi,'rows');
                        for n = 1: size(uF,1)
                            ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                            powercurveo.factors{f}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
                        end
                    end
                end            
            end

            for i = 1 : nInteractions
                randg = randv{nFactors+i};
                mati = randg(prod(powercurveo.nLevels(powercurveo.interactions{i}.factors)),M);
                mati = sqrt(size(mati,1))*mati/norm(mati,'fro'); 
                Fi = F(:,powercurveo.interactions{i}.factors);
                powercurveo.interactions{i}.matrix = [];
                uF = unique(Fi,'rows');
                for n = 1: size(uF,1)
                    ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                    powercurveo.interactions{i}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
                end
            end
        end

        Xstruct = zeros(N,M);
        for f = 1 : nFactors
            if ~ordinal(f)
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
            end
        end

        for i = 1 : nInteractions
            Xstruct = Xstruct + randgC() * powercurveo.coeffs(nFactors+i) * powercurveo.interactions{i}.matrix;
        end

        randg = randv{end};
        Xnoise = randg(N,M); 
        Xnoise = randgC() * powercurveo.rescoef * sqrt(N)*Xnoise/norm(Xnoise,'fro'); 
        
        for a = 1:length(theta)

            Xm = Xnoise; %Xm = (1-theta(a))*Xnoise;
            Xm = Xm + theta(a)*Xstruct;

            for f = 1 : nFactors
                if ordinal(f)
                    Xm = Xm + randgC() * powercurveo.coeffs(f) * powercurveo.factors{f}.matrix;
                end
            end
                        
            % Parallel GLM
            [T, parglmo] = parglm(Xm, F, 'Model', model, 'Preprocessing', prep, 'Permutations', nPerm, 'Ts', ts, 'Ordinal', ordinal, 'Random', random, 'Fmtc', fmtc, 'Coding', coding, 'Nested', nested);

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

                powercurveo = powercurveo2;
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
                
                for f = 1 : nFactors
                    if ~ordinal(f)
                        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
                            uF = unique(F(:,f));
                            powercurveo.nLevels(f) = length(uF);
                        else % if nested
                            ind = find(nested(:,2)==f);
                            ref = nested(ind,1);
                            powercurveo.factors{f}.factors = [ref powercurveo.factors{ref}.factors];
                            urF = unique(F(:,ref));
                            powercurveo.nLevels(f) = 0;
                            for j = 1:length(urF)
                                rind = find(ismember(F(:,ref),urF(j)));
                                uF = unique(F(rind,f));
                                powercurveo.nLevels(f) = powercurveo.nLevels(f) + length(uF);
                            end
                        end
                    end
                end
                
            else
                F = repmat(F2,theta(a),1);
            end
            
            N = size(F,1);
        
            repa = 1;
            if isstruct(X)
                for f = 1 : nFactors
                    randg = randv{f};
                    if ordinal(f)
                        powercurveo.factors{f}.matrix = randg(N,M);
                        powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
                    else
                        powercurveo.factors{f}.matrix = randg(powercurveo.nLevels(f),M);
                        powercurveo.factors{f}.matrix = sqrt(powercurveo.nLevels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
                        
                        if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
                            powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
                        else % nested                       
                            mati = powercurveo.factors{f}.matrix;
                            Fi = F(:,[powercurveo.factors{f}.factors f]);
                            powercurveo.factors{f}.matrix = [];
                            uF = unique(Fi,'rows');
                            for n = 1: size(uF,1)
                                ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                                powercurveo.factors{f}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
                            end
                        end
                    end
                end
                
                for i = 1 : nInteractions
                    randg = randv{nFactors+i};
                    mati = randg(prod(powercurveo.nLevels(powercurveo.interactions{i}.factors)),M);
                    mati = sqrt(size(mati,1))*mati/norm(mati,'fro'); 
                    Fi = F(:,powercurveo.interactions{i}.factors);
                    powercurveo.interactions{i}.matrix = [];
                    uF = unique(Fi,'rows');
                    for n = 1: size(uF,1)
                        ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                        powercurveo.interactions{i}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
                    end
                end
            else
                if replicates>0
                    f = replicates;
                    randg = randv{f};
                    uF = unique(F2(:,f));

                    if theta(a) <= length(uF)
                        for f2 = 1 : nFactors
                            powercurveo.factors{f2}.matrix(find(F2(:,f)>theta(a)),:) = [];
                            powercurveo.factors{f2}.matrix = sqrt(N)*powercurveo.factors{f2}.matrix/norm(powercurveo.factors{f2}.matrix,'fro');
                        end
                        for i = 1 : nInteractions
                            powercurveo.interactions{i}.matrix(find(F2(:,f)>theta(a)),:) = [];
                            powercurveo.interactions{i}.matrix = sqrt(N)*powercurveo.interactions{i}.matrix/norm(powercurveo.interactions{i}.matrix,'fro');
                        end
                    else
                        for f2 = 1 : nFactors
                            if f2 ~= f    
                                for t = 1:theta(a)-length(uF)
                                    Fi = powercurveo.factors{f2}.matrix(find(F2(:,f)==uF(1)),:);
                                    Fi(:,f) = max(uF) + t;
                                    powercurveo.factors{f2}.matrix = [powercurveo.factors{f2}.matrix;Fi];
                                end
                                powercurveo.factors{f2}.matrix = sqrt(N)*powercurveo.factors{f2}.matrix/norm(powercurveo.factors{f2}.matrix,'fro');                    
                            else             
                                if ordinal(f)
                                    powercurveo.factors{f}.matrix = randg(N,M);
                                    powercurveo.factors{f}.matrix = sqrt(N)*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');
                                else
                                    powercurveo.factors{f}.matrix = randg(powercurveo.nLevels(f),M);
                                    powercurveo.factors{f}.matrix = sqrt(powercurveo.nLevels(f))*powercurveo.factors{f}.matrix/norm(powercurveo.factors{f}.matrix,'fro');

                                    if isempty(nested) || isempty(find(nested(:,2)==f)) % order 1
                                        powercurveo.factors{f}.matrix = powercurveo.factors{f}.matrix(F(:,f),:);
                                    else % nested
                                        mati = powercurveo.factors{f}.matrix;
                                        Fi = F(:,[powercurveo.factors{f}.factors f]);
                                        powercurveo.factors{f}.matrix = [];
                                        uF = unique(Fi,'rows');
                                        for n = 1: size(uF,1)
                                            ind = find(uF(n,1)==Fi(:,1)&uF(n,2)==Fi(:,2));
                                            powercurveo.factors{f}.matrix(ind,:) = ones(length(ind),1)*mati(n,:);
                                        end
                                    end
                                end
                            end
                        end
                        for i = 1 : nInteractions
                            for t = 1:theta(a)-length(uF)
                                Fi = powercurveo.interactions{i}.matrix(find(F2(:,f)==uF(1)),:);
                                Fi(:,f) = max(uF) + t;
                                powercurveo.interactions{i}.matrix = [powercurveo.interactions{i}.matrix;Fi];
                            end
                            powercurveo.interactions{i}.matrix = sqrt(N)*powercurveo.interactions{i}.matrix/norm(powercurveo.interactions{i}.matrix,'fro');
                        end
                    end
                else
                    repa = theta(a);
                end
            end

            Xstruct = zeros(N,M);
            for f = 1 : nFactors
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(f) * repmat(powercurveo.factors{f}.matrix,repa,1);
            end
            
            for i = 1 : nInteractions
                Xstruct = Xstruct + randgC() * powercurveo.coeffs(nFactors+i) * repmat(powercurveo.interactions{i}.matrix,repa,1);
            end

            randg = randv{end};
            Xnoise = randg(N,M);
            Xnoise = randgC() * powercurveo.rescoef * sqrt(N)*Xnoise/norm(Xnoise,'fro');

            Xm = Xnoise + Xstruct;
            
            % Parallel GLM
            [T, parglmo] = parglm(Xm, F, 'Model', model, 'Preprocessing', prep, 'Permutations', nPerm, 'Ts', ts, 'Ordinal', ordinal, 'Random', random, 'Fmtc', fmtc, 'Coding', coding, 'Nested', nested);
            
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
spcBootstrap(theta,PCrep,nRep,true,0.05,false);    
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









