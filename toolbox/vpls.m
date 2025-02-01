function model = vpls(xcs,ycs,varargin)

% Variable selection PLS using the simpls algorithm. Se description of the 
% methods in Mehmood, Tahir, Solve Sæbø, y Kristian Hovde Liland. 2020. 
% «Comparison of Variable Selection Methods in Partial Least Squares 
% Regression». Journal of Chemometrics 34 (6): e3226. 
%
% model = vpls(xcs,ycs)     % minimum call equivalent to simpls
%
% See also: simpls, kernelpls, pcaEig, asca, vasca
%
%
% INPUTS:
%
% xcs: [NxM] preprocessed billinear data set 
%
% ycs: [NxO] preprocessed billinear data set of responses
%
%
% Optional INPUTS (parameter):
%
% 'LVs': [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%   first two LVs). By default, lvs = 0:size(XX)
%
% 'Selection': str
%   'Weights': filter method based on the PLS weights (W)
%   'AltWeights': filter method based on the PLS alternative weights (R)
%   'Regressors': filter method based on the PLS regression coefficients (beta)
%   'SR': filter method based on the selectivity ratio (by default)
%   'VIP': filter method based on Variance Importance in PLS Projection
%   'T2': wrapper method based on the Hotelling T2 statistic
%   'sPLS': embedded method based on sparse PLS
%
% 'VarNumber': [1x1] number of variables to selection (floor(M/2) by default)
%
%
% OUTPUTS:
%
% model: structure that contains model information
%   var: [1x1] xcs sum of squares
%   lvs: [1xA] latent variable numbers
%   loads: [MxA] matrix of x-loadings P
%   yloads: [OxA] matrix of y-loadings Q
%   weights: [MxA] matrix of weights W
%   altweights: [MxA] matrix of alternative weights R
%   scores: [NxA] matrix of x-scores T
%   beta: [MxO] matrix of regressors
%   type: 'PLS'
%
%
% EXAMPLE OF USE: Random data with structural relationship
%
% X = simuleMV(20,10,'LevelCorr',8);
% Y = 0.1*randn(20,2) + X(:,1:2);
% Xcs = preprocess2D(X);
% Ycs = preprocess2D(Y);
% lvs = 1:2;
% modelSR = vpls(Xcs,Ycs,'LVs',lvs,'VarNumber',2,'Selection','SR');
% modelVIP = vpls(Xcs,Ycs,'LVs',lvs,'VarNumber',2,'Selection','VIP');
% modelsPLS = vpls(Xcs,Ycs,'LVs',lvs,'VarNumber',2,'Selection','sPLS'); 
%
% modelSR.beta
% modelVIP.beta
% modelsPLS.beta % sPLS can have 2 different vars per LV 
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 31/Jan/2025
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
[N,M] = size(xcs);
O = size(ycs, 2);

% Introduce optional inputs as parameters
p = inputParser;
addParameter(p,'LVs',0:rank(xcs)); 
addParameter(p,'Selection','SR'); 
addParameter(p,'VarNumber',M);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
lvs = p.Results.LVs;
selection = p.Results.Selection;
V = p.Results.VarNumber;

% Convert column arrays to row arrays
if size(lvs,2) == 1, lvs = lvs'; end;

% Preprocessing
lvs = unique(lvs);
lvs(find(lvs==0)) = [];
lvs(find(lvs>min(V,M))) = [];
A = length(lvs);

% Validate dimensions of input data
assert (isequal(size(xcs), [N M]), 'Dimension Error: parameter ''X'' must be N-by-M. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(ycs), [N O]), 'Dimension Error: parameter ''Y'' must be N-by-O. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(lvs), [1 A]), 'Dimension Error: parameter ''LVs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(lvs<0)) && isequal(fix(lvs), lvs), 'Value Error: parameter ''LVs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);
assert ((V>0 && V <= M) && isequal(fix(V), V), 'Value Error: parameter ''VarNumber'' must contain a positive integer equal or below the numer of variables in ''X''. Type ''help %s'' for more info.', routine(1).name);


%% Main code

model = simpls(xcs,ycs,'LVs',lvs);

if V < M 
    if strcmp(selection,'Weights')
        r = sum(model.weights.^2,2);
        [rord,ind] = sort(r,'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
        
    elseif strcmp(selection,'AltWeights')
        r = sum(model.altweights.^2,2);
        [rord,ind] = sort(r,'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
        
    elseif strcmp(selection,'Regressors')
        r = model.beta;
        [rord,ind] = sort(r,'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
        
    elseif strcmp(selection,'VIP')
        for a = lvs
            SSY(a) = sum(sum((model.scores(:,a)*model.yloads(:,a)').^2));
        end
        w = 0;
        for a = lvs
            w = w + SSY(a)*(model.weights(:,a).^2);
        end
        r = sqrt(M * w / sum(SSY));
        [rord,ind] = sort(r,'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
    
    elseif strcmp(selection,'SR')
        
        model = simpls(xcs,ycs,'LVs',lvs);
        for i = 1:size(ycs,2)
            w = model.beta(:,i);
            t = xcs*w/sqrt(w'*w);
            p = xcs'*t/(t'*t);
            e = xcs - t*p';
            r(i,:) = sum((t*p').^2,1)./sum((e).^2,1);
        end
        [rord,ind] = sort(sum(r,1),'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
        
    elseif strcmp(selection,'T2')
        w = preprocess2D(model.weights,'Preprocess',1);
        r = zeros(1,M);
        iW = inv(w'*w);
        for i=1:M
            r(i) = M*w(i,:)*iW*w(i,:)';
        end
        [rord,ind] = sort(r,'descend');
        
        xsel = zeros(size(xcs));
        xsel(:,ind(1:V)) = xcs(:,ind(1:V));
        model = simpls(xsel,ycs,'LVs',lvs);
        
    elseif strcmp(selection,'sPLS') % More flexible, more than V variables can be selected if several responses
        model = sparsepls2(xcs, ycs, max(lvs), V*ones(size(1:max(lvs))), O*ones(size(1:max(lvs))), 500, 1e-10, 1, 0);
        model.beta = model.R*model.Q';
    end
end
    
    