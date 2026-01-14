function rec = missTsr2D(x,pcs,varargin)

% Missing data imputation with Trimmed Scores Regression.
%
% rec = missTsr2D(x,pc) % minimum call
%
% See also: pcaEig
%
%
% INPUTS:
%
% x: (NxM) data matrix, N(observations) x M(variables)
%
% pcs: [1xA] Principal Components considered (e.g. pcs = 1:2 selects the
%   first two PCs).
%
%
% Optional INPUTS (parameters):
%
% 'Preprocessing': (1x1) preprocesing of the data
%       0: no preprocessing.
%       1: mean centering.
%       2: autoscaling (default)  
%
% 'Percentage': (1x1) maximum percentage of missing values in a row or column. (0.3
%   by default)
%
% 'AutoCorrData': (1x1) 1 for autocorrelated data, 0 otherwise (by default).
%
% 'Iterations': (1x1) maximum number of iterations (100 by default).
%
% 'Convergence': (1x1) convergence threshold. (1e-5 by default)
%
%
% OUTPUTS:
%
% rec: (NxM) recovered data matrix, N(observations) x M(variables)
%
%
% EXAMPLE FO USE:
%
% X = simuleMV(20,100,'LevelCorr',8);
% pcs = 1:2;
% Xmiss = X;
% missLoc = round(200*rand(10,1));
% Xmiss(missLoc) = nan;
% 
% rec0 = missTsr2D(Xmiss,0);
% plotScatter([X(missLoc),rec0(missLoc)],'XYLabel',{'Real' 'Uncond. Mean Replacement'});
%
% rec = missTsr2D(Xmiss,pcs);
% plotScatter([X(missLoc),rec(missLoc)],'XYLabel',{'Real' 'TSR'});
% 
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 28/Dec/2025
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end
s = size(x);
if find(s<1), error('Incorrect content of x.'); end

% Convert column arrays to row arrays
if size(pcs,2) == 1, pcs = pcs'; end

% Preprocessing
pcs = unique(pcs);
pcs(find(pcs==0)) = [];
pcs(find(pcs>size(x,2))) = [];
A = length(pcs);

% Validate dimensions of input data
assert (isequal(size(pcs), [1 A]), 'Dimension Error: parameter ''Pcs'' must be 1-by-A. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(pcs<0)) && isequal(fix(pcs), pcs), 'Value Error: parameter ''Pcs'' must contain positive integers. Type ''help %s'' for more info.', routine(1).name);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Preprocessing',2);   
addParameter(p,'Percentage',0.3);
addParameter(p,'AutoCorrData',0); 
addParameter(p,'Iterations',100);
addParameter(p,'Convergence',1e-5);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
prep = p.Results.Preprocessing;
perc = p.Results.Percentage;
ac = p.Results.AutoCorrData;
iter = p.Results.Iterations;
conv = p.Results.Convergence;

if (prep<0||prep>2), error('Incorrect value of Preprocessing.'); end;

if (perc<0||perc>1), error('Incorrect value of Percentage.'); end;

if (ac<0||ac>1), error('Incorrect value of AutoCorrData.'); end;

if (iter<1), error('Incorrect value of Iterations.'); end;

if (conv<0), error('Incorrect value of Convergence.'); end;

% Main code

s=size(x);
nnan = isnan(x);
indnan=find(nnan);

if ~isempty(find(perc*s(1)<sum(nnan,1)))
    error('Some columns present too much missing values.');
end
if ~isempty(find(perc*s(2)<sum(nnan,2)))
    error('Some rows present too much missing values.');
end

x(indnan) = nan;
    
loopIni = true;    
e0=Inf;
num0=iter;
ax = x; 
while e0>conv && num0 > 0

    [xce,m,d] = preprocess2D(ax,'Preprocessing',prep);
    
    med = repmat(m,s(1),1);
    dev = repmat(d,s(1),1);
   
    if loopIni
        loopIni = false;
        xce(indnan) = 0;
        x(indnan) = med(indnan);
    end
    
    ax = xce;
    sa = size(ax);

    model = pcaEig(ax,'PCs',pcs);
    T=model.scores;
    P=model.loads;
         
    if max(pcs)>0
        for i=1:sa(2)
            TT=ax(:,[1:i-1 i+1:end])*P([1:i-1 i+1:end],:);
            r=(TT'*TT)\TT'*T;
            ax2(:,i)=TT*r*P(i,:)';
        end

        recMod = ax;
        recMod(indnan) = ax2(indnan);
    else
        recMod = xce;
    end
    
    x2 = recMod.*dev + med;
    
    e0=sum(((x(indnan)-x2(indnan))./dev(indnan)).^2)/length(indnan);
    num0 = num0-1;
    x(indnan)=x2(indnan);
    
    ax=x;
end

if ac==1
    xp = (x-med)./dev;
    e = xp - xp*P*P';
    e(indnan) = nan;

    for i=1:size(e,2)
        x2 = find(~isnan(e(:,i)));
        y = e(x2,i);
        e(:,i) = interp1(x2,y,1:size(e,1),'spline');
    end

    rec = (xp*P*P' + e).*dev + med;
else
    rec = x;
end

