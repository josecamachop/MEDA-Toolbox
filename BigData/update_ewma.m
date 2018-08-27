function Lmodel = update_ewma(list,path,Lmodel,lambda,step,debug,erase)

% Big data analysis based on bilinear proyection models (PCA and PLS) with
% the exponentially weighted moving average approach.
%
% Lmodel = update_ewma(list)          % minimum call
% Lmodel = update_ewma(list,path,Lmodel,lambda,step,debug,erase) % complete call
%
%
% INPUTS:
%
% list: {Fx1} list of strings with the names of the files for the update or
%   struct array with x (and optionally y) matrices.
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
% n_obs = 100;
% n_vars = 10;
% Lmodel = Lmodel_ini(simuleMV(n_obs,n_vars,6));
% Lmodel.type = 1; 
% Lmodel.prep = 2;  
% Lmodel.lvs = 1;
% Lmodel.nc = 100; % Number of clusters
% Lmodel.mat = loadings_Lpca(Lmodel,0);
% mspc_Lpca(Lmodel);
%
% for i=1:4,
%   n_obst = 10;
%   list(1).x = simuleMV(n_obst,n_vars,6,corr(Lmodel.centr)*(n_obst-1)/(Lmodel.N-1));
%   Lmodel = update_ewma(list,[],Lmodel);
%   mspc_Lpca(Lmodel);
% end
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 24/Aug/18
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez
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

if nargin < 2 || isempty(path), path = ''; end;
if nargin < 3 || isempty(Lmodel), 
    Lmodel = Lmodel_ini; 
    Lmodel.type = 1;
    Lmodel.lvs = 0;
    Lmodel.prep = 2;
end;
[ok, Lmodel] = check_Lmodel(Lmodel);
if nargin < 4 || isempty(lambda), lambda = 1; end;
if nargin < 5 || isempty(step), step = 1; end;
if nargin < 6 || isempty(debug), debug = 1; end;
if nargin < 7 || isempty(erase), erase = 1; end;

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
    
for t=1:length(list),
    
    if debug, disp(sprintf('clustering: packet %d...........................................', t)), end;
    
    if Lmodel.type==1,
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
    elseif Lmodel.type==2,
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
        if ismember('obs_l', vars)
            obs_l = list(t).obs_l;
        else
            obs_l = cellstr(num2str((1:size(x,1))'));
        end
    else
        vars = whos('-file',[path list{t}]);
        if ismember('class', {vars.name})
            load([path list{t}],'class')
        else
            class = ones(size(x,1),1);
        end
        if ismember('obs_l', {vars.name})
            load([path list{t}],'obs_l')
        else
            obs_l = cellstr(num2str((1:size(x,1))'));
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
    
    if Lmodel.type==1
        
        [P,sdT] = Lpca(Lmodel);
        Lmodel.mat = P*diag(1./sdT);
        
    elseif Lmodel.type==2,
        
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
        
        if rank(Lmodel.XY)>0,
            
            [beta,W,P,Q,R,sdT] = Lpls(Lmodel);
            Lmodel.mat = R*diag(1./sdT);
            
        else
            
            if debug>1, disp('XY Rank 0: using PCA.'), end;
            
            [P,sdT] = Lpca(Lmodel);
            Lmodel.mat = P*diag(1./sdT);
            
        end
        
    end
    
    Lmodel.multr = lambda*Lmodel.multr;
    ind_lab = find(Lmodel.multr>=erase);
    Lmodel.centr =  Lmodel.centr(ind_lab,:);
    Lmodel.multr = Lmodel.multr(ind_lab);
    Lmodel.class = Lmodel.class(ind_lab);
    if ~isempty(Lmodel.obs_l)
        Lmodel.obs_l = Lmodel.obs_l(ind_lab);    
    end
    Lmodel.updated = zeros(length(ind_lab),1);

    s = size(x);
    step2 = max(10,round(s(1)*step));
    for i = 1:step2:s(1),
        endv = min(s(1),i+step2);
        ss = endv-i+1;
        xstep = xcs(i:endv,:);
        clstep = class(i:endv,:);
        if isempty(obs_l)
            obs_step = {};
        else
            obs_step = obs_l(i:endv);
        end
               
        Lmodel.centr = [Lmodel.centr;xstep];
        Lmodel.multr = [Lmodel.multr;ones(ss,1)];
        Lmodel.class = [Lmodel.class;clstep];
        Lmodel.obs_l = {Lmodel.obs_l{:} obs_step{:}};
        Lmodel.updated = [Lmodel.updated;ones(size(xstep,1),1)]; 
            
        [Lmodel.centr,Lmodel.multr,Lmodel.class,Lmodel.obs_l,Lmodel.updated] = psc(Lmodel.centr,Lmodel.nc,Lmodel.multr,Lmodel.class,Lmodel.obs_l,Lmodel.updated,Lmodel.mat);

    end
    
end

if Lmodel.type==1, % Update mat acording to actual scores
    Lmodel.mat = P;
elseif Lmodel.type==2,
    if rank(Lmodel.XY)>0,
        Lmodel.mat = R;
    else
        Lmodel.mat = P;
    end
end
T = Lmodel.centr*Lmodel.mat;
mM = max(T)-min(T);
Lmodel.mat = Lmodel.mat*diag(1./mM);
