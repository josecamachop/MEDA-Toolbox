function rec = missTSR2D(x,pc,varargin)

% Missing data imputation with Trimmed Scores Regression.
%
% rec = missTSR2D(x,pc) % minimum call
% rec = missTSR2D(x,pc,'Preprocessing',prep,'Percentage',perc,'AutoCorrData'ac,'Iterations',iter,'Convergence',conv)   % complete call
%
%
% INPUTS:
%
% x: (NxM) data matrix, N(observations) x M(variables)
%
% pc: (1x1) number of principal components for the D-statistic.
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
% EXAMPLE FO USE:
%
% X = simuleMV(20,10,'LevelCorr',8);
% pc = 2;
% 
% rec = missTSR2D(X,pc,'Iterations',50);
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/2024
% major change: include nipls
%
% Copyright (C) 2024  University of Granada, Granada
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
%% Main Code
% Parameters checking

if nargin < 2, error('Error in the number of arguments.'); end;
if ndims(x)~=2, error('Incorrect number of dimensions of x.'); end;
s = size(x);
if find(s<1), error('Incorrect content of x.'); end;
if pc<0, error('Incorrect value of pc.'); end;


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
ind_nan=find(nnan);

if ~isempty(find(perc*s(1)<sum(nnan,1))),
    error('Some columns of the unfolded matrix present too much missing values.');
end
if ~isempty(find(perc*s(2)<sum(nnan,2))),
    error('Some rows of the unfolded matrix present too much missing values.');
end

x(ind_nan) = nan;
    
loop_ini = true;    
e0=Inf;
num0=iter;
ax = x; 
while e0>conv && num0 > 0,

    [xce,m,d] = preprocess2D(ax,'Preprocessing',prep);
    
    med = zeros(s);
    for j=1:s(1),
        med(j,:) = m;
    end
    
    dev = zeros(s);
    for j=1:s(1),
        dev(j,:) = d;
    end
   
    if loop_ini,
        loop_ini = false;
        xce(ind_nan) = 0;
        x(ind_nan) = med(ind_nan);
    end
    
    ax = xce;
    sa = size(ax);

    T=[];
    P=[];
    ax2=ax;
    for j=1:pc, % PCA model
        e1=Inf;
        num1=iter;
        t=ax(:,1);
        while e1>conv && num1 > 0,
            p=t'*ax2;
            p=p/sqrt(p*p');
            t2=ax2*p';
            e1=sum((t-t2).^2);
            t=t2;
            num1 = num1-1;
        end
        
        ax2=ax2-t*p;
        T=[T t];
        P=[P,p'];
    end
         
    if pc>0,
        T=ax*P;
        for i=1:sa(2);
            TT=ax(:,[1:i-1 i+1:end])*P([1:i-1 i+1:end],:);
            r=inv(TT'*TT)*TT'*T;
            ax2(:,i)=TT*r*P(i,:)';
        end

        rec_mod = ax;
        rec_mod(ind_nan) = ax2(ind_nan);
    else
        rec_mod = xce;
    end
    
    x2 = rec_mod.*dev + med;
    
    e0=sum(((x(ind_nan)-x2(ind_nan))./dev(ind_nan)).^2)/length(ind_nan);
    num0 = num0-1;
    x(ind_nan)=x2(ind_nan);
    
    ax=x;
end

if ac==1,
    xp = (x-med)./dev;
    e = xp - xp*P*P';
    e(ind_nan) = nan;

    for i=1:size(e,2);
        x2 = find(~isnan(e(:,i)));
        y = e(x2,i);
        e(:,i) = interp1(x2,y,1:size(e,1),'spline');
    end

    rec = (xp*P*P' + e).*dev + med;
else
    rec = x;
end

