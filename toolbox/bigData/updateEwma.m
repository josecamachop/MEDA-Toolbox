function Lmodel = updateEwma(list,varargin)

% Big data analysis based on bilinear proyection models (PCA and PLS) with
% the exponentially weighted moving average approach.
%
% Lmodel = updateEwma(list)          % minimum call
%
%
% INPUTS:
%
% list: {Fx1} list of strings with the names of the files for the update or
%   struct array with x (and optionally y) matrices.
%
% Optional INPUTS (parameter):
%
% path: (str) path to the directory where the data files are located ('' by
%   default)
%
% Lmodel: (struct Lmodel) model to update (initialized to PCA model with 1
%   PC and auto-scaling by default)
%
% lambda: [1x1] forgetting factor between 0 (fast adaptation) and 1 (long
%   history) (1 by default)
%
% step: [1x1] percentage of the data in the file to be used in each
%   iteration. For time-course data 1 is suggested (1 by default)
%
% debug: [1x1] disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default)
%       2: display all messages.  
%
% erase: [1x1] threshold to erase an observation (1 by default)
%
%
% OUTPUTS:
%
% Lmodel: (struct Lmodel) model updated.
%
%
% EXAMPLE OF USE: update a random model with new random observations.
%
% nobs = 100;
% nvars = 10;
% Lmodel = iniLmodel(simuleMV(nobs,nvars,'LevelCorr',6));
% Lmodel.type = 'PCA'; 
% Lmodel.prep = 2;  
% Lmodel.lvs = 1;
% Lmodel.nc = 100; % Number of clusters
% Lmodel.mat = loadingsLpca(Lmodel,0);
% mspcLpca(Lmodel);
%
% for i=1:4,
%   nobst = 10;
%   list(1).x = simuleMV(nobst,nvars,'LevelCorr',6,'Covar',corr(Lmodel.centr)*(nobst-1)/(Lmodel.N-1));
%   Lmodel = updateEwma(list,'path',[],'Lmodel',Lmodel);
%   mspcLpca(Lmodel);
% end
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 21/Nov/2024
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

p = inputParser;
addParameter(p,'path','');   
    defmodel = iniLmodel; 
    defmodel.type = 'PCA';
    defmodel.lvs = 0;
    defmodel.prep = 2;
addParameter(p,'Lmodel',defmodel);   
addParameter(p,'lambda',1);   
addParameter(p,'step',1);   
addParameter(p,'debug',1);  
addParameter(p,'erase',1);    
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
path = p.Results.path;
Lmodel = p.Results.Lmodel;
lambda = p.Results.lambda;
step = p.Results.step;
debug = p.Results.debug;
erase = p.Results.erase;

% Validate dimensions of input data
assert (isequal(size(lambda), [1 1]), 'Dimension Error: 4th argument must be a scalar. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(step), [1 1]), 'Dimension Error: 5th argument must be a scalar. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(debug), [1 1]), 'Dimension Error: 6th argument must be a scalar. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(erase), [1 1]), 'Dimension Error: 7th argument must be a scalar. Type ''help %s'' for more info.', routine(1).name);
  
% Validate values of input data
assert (lambda>=0 && lambda<=1, 'Value Error: 4th argument must be in interval (0, 1]. Type ''help %s'' for more info.', routine(1).name);
assert (step>0 && step<=1, 'Value Error: 5th argument must be in interval (0, 1]. Type ''help %s'' for more info.', routine(1).name);
assert (debug==0 || debug==1 || debig==2, 'Value Error: 6th argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);
assert (erase>0 && erase<=1, 'Value Error: 7th argument must be in interval (0, 1]. Type ''help %s'' for more info.', routine(1).name);
    
    
%% Main code

Lmodel.update = 1; 
    
for t=1:length(list)
    
    if debug, disp(sprintf('clustering: packet %d...........................................', t)), end;
    
    if strcmp(Lmodel.type,'PCA')
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
    elseif strcmp(Lmodel.type,'PLS')
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
        else
            load([path list{t}],'x','y')
        end
    end
       
    if isstruct(list(t))
        vars = fieldnames(list(t));
        if ismember('class', vars)
            class = list(t).class;
        else
            class = ones(size(x,1),1);
        end
        if ismember('obsl', vars)
            obsl = list(t).obsl;
        else
            obsl = cellstr(num2str((1:size(x,1))'));
        end
    else
        vars = whos('-file',[path list{t}]);
        if ismember('class', {vars.name})
            load([path list{t}],'class')
        else
            class = ones(size(x,1),1);
        end
        if ismember('obsl', {vars.name})
            load([path list{t}],'obsl')
        else
            obsl = cellstr(num2str((1:size(x,1))'));
        end
    end
    
    N = Lmodel.N;
    indMV = find(isnan(x));
    if ~isempty(indMV)
        disp('Missing values found in X. Set to average.');
        av = ones(size(x,1),1)*Lmodel.av;
        x(indMV) = av(indMV);
    end
         
    [xcs,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,Lmodel.prep,0,lambda,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);

    if isempty(Lmodel.centr)
        Lmodel.centr = xcs;
    end
            
    if isempty(Lmodel.XX)
        Lmodel.XX = xcs'*xcs;
    else
        Lmodel.XX = lambda*Lmodel.XX + xcs'*xcs;
    end
    
    ind = isnan(Lmodel.XX);
    Lmodel.XX(ind) = 0;
    
    if strcmp(Lmodel.type,'PCA')
        
        Lmodel = Lpca(Lmodel);
        P = Lmodel.loads;
        sdT = Lmodel.sdT;
        Lmodel.mat = P*diag(1./sdT);
        
    elseif strcmp(Lmodel.type,'PLS')
        
        indMV = find(isnan(y));
        if ~isempty(indMV)
            disp('Missing values found in Y. Set to average.');
            av = ones(size(y,1),1)*Lmodel.avy;
            y(indMV) = av(indMV);
        end
    
        [ycs,Lmodel.avy,Lmodel.scy] = preprocess2Di(y,Lmodel.prepy,0,lambda,Lmodel.avy,Lmodel.scy,N,Lmodel.weighty);
        
        if isempty(Lmodel.XY)
            Lmodel.XY = xcs'*ycs;
        else
            Lmodel.XY = lambda*Lmodel.XY + xcs'*ycs;
        end
        
        if isempty(Lmodel.YY)
            Lmodel.YY = ycs'*ycs;
        else
            Lmodel.YY = lambda*Lmodel.YY + ycs'*ycs;
        end
        
        if rank(Lmodel.XY)>0
            
            Lmodel = Lpls(Lmodel);
            R = Lmodel.altweights;
            P = Lmodel.loads;
            sdT = Lmodel.sdT;
            Lmodel.mat = R*diag(1./sdT);
            
        else
            
            if debug>1, disp('XY Rank 0: using PCA.'), end;
            
            Lmodel = Lpca(Lmodel);
            P = Lmodel.loads;
            sdT = Lmodel.sdT;
            Lmodel.mat = P*diag(1./sdT);
            
        end
        
    end
    
    Lmodel.multr = lambda*Lmodel.multr;
    indlab = find(Lmodel.multr>=erase);
    Lmodel.centr =  Lmodel.centr(indlab,:);
    Lmodel.multr = Lmodel.multr(indlab);
    Lmodel.class = Lmodel.class(indlab);
    if ~isempty(Lmodel.obsl)
        Lmodel.obsl = Lmodel.obsl(indlab);    
    end
    Lmodel.updated = zeros(length(indlab),1);

    s = size(x);
    step2 = max(10,round(s(1)*step));
    for i = 1:step2:s(1),
        endv = min(s(1),i+step2);
        ss = endv-i+1;
        xstep = xcs(i:endv,:);
        clstep = class(i:endv,:);
        if isempty(obsl)
            obsstep = {};
        else
            obsstep = obsl(i:endv);
        end
               
        Lmodel.centr = [Lmodel.centr;xstep];
        Lmodel.multr = [Lmodel.multr;ones(ss,1)];
        Lmodel.class = [Lmodel.class;clstep];
        Lmodel.obsl = {Lmodel.obsl{:} obsstep{:}};
        Lmodel.updated = [Lmodel.updated;ones(size(xstep,1),1)]; 
            
        [Lmodel.centr,Lmodel.multr,Lmodel.class,Lmodel.obsl,Lmodel.updated] = psc(Lmodel.centr,Lmodel.nc,Lmodel.multr,Lmodel.class,Lmodel.obsl,Lmodel.updated,Lmodel.mat);

    end
    
end

if strcmp(Lmodel.type,'PCA') % Update mat acording to actual scores
    Lmodel.mat = P;
elseif strcmp(Lmodel.type,'PLS')
    if rank(Lmodel.XY)>0
        Lmodel.mat = R;
    else
        Lmodel.mat = P;
    end
end
T = Lmodel.centr*Lmodel.mat;
mM = max(T)-min(T);
Lmodel.mat = Lmodel.mat*diag(1./mM);
