function Lmodel = update_iterative(list,path,Lmodel,step,files,debug)

% Big data analysis based on bilinear proyection models (PCA and PLS),
% iterative approach.
%
% Lmodel = update_iterative(list)          % minimum call
% Lmodel = update_iterative(list,path,Lmodel,step,files,debug) % complete call
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
% step: [1x1] percentage of the data in the file to be used in each
%   iteration. For time-course data 1 is suggested (1 by default)
%
% files: [1x1] create the file system with the original data (1) or not (0, by
%   default)
%
% debug: [1x1] disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default)
%       2: display all messages.
%
%
% OUTPUTS:
%
% Lmodel: (struct Lmodel) model updated.
%
%
% NOTE: On the MEDA FileSystem for input argument files set to 1. It is 
% based on CSV files (to change in the future) with two hierarchies. The 
% top layer contains pointers to data files, and it is only set for those 
% clusters with more than 100 observations (this number is hardcoded in the 
% present routine, in the call to cfilesys) The bottom layer contains the 
% actual data of the observations. CSV format is to change in the future. 
% The name of the files follows this structure: For top layer files we use 
% MEDA#to#oc#c, where #t is the index of the original file in the imput 
% list for the first observation introduced, starting from 1, #o is the 
% order of the  observation in that file and #c the class. For bottom layer
% we add _#n for the n-th file depending on the same top file. Name format 
% can be chaged in this routine (line app. 399).
%
%
% EXAMPLE OF USE: model a large number of observations.
%
% n_obs = 100;
% n_vars = 10;
% Lmodel = Lmodel_ini;
% Lmodel.type = 1; 
% Lmodel.prep = 2;  
% Lmodel.lvs = 1;
% Lmodel.nc = 100; % Number of clusters
% 
% for i=1:10,
%   list(i).x = simuleMV(n_obs,n_vars,6);
% end
%
% Lmodel = update_iterative(list,[],Lmodel);
% mspc_Lpca(Lmodel);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/May/17
%
% Copyright (C) 2017  University of Granada, Granada
% Copyright (C) 2017  Jose Camacho Paez
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
if nargin < 4 || isempty(step), step = 1; end;
if nargin < 5 || isempty(files), files = 0; end;
if nargin < 6 || isempty(debug), debug = 1; end;

% Validate dimensions of input data
assert (isequal(size(step), [1 1]), 'Dimension Error: 4th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(files), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(debug), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (isempty(find(step<=0 | step>1)), 'Value Error: 4th argument must contain values in (0,1]. Type ''help %s'' for more info.', routine(1).name); 
assert (isempty(find(files~=0 & files~=1)), 'Value Error: 5th argument must contain 0 or 1. Type ''help %s'' for more info.', routine(1).name);
assert (isempty(find(debug~=0 & debug~=1 & debug~=2)), 'Value Error: 6th argument must contain 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

Lmodel.update = 2; 

if files,
	[status,result] = system(['del ' Lmodel.path 'MEDA*.txt']); % delete previous files
end

% preprocess

% compute mean

if Lmodel.prep==0,
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
        
        indMV{t} = find(isnan(x));
        if ~isempty(indMV{t})
            disp('Missing values found in X. Set to average.');
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        [xc,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,0,0,1,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);
        
    end
        
elseif (Lmodel.type==1 && Lmodel.prep > 0) || (Lmodel.type==2 && Lmodel.prep > 0 && Lmodel.prepy == 0), 
    
    if debug, disp('mean centering X block...........................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
        
        indMV{t} = find(isnan(x));
        if ~isempty(indMV{t})
            disp('Missing values found in X. Set to average.');
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        [xc,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,1,0,1,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);
        
    end
    
elseif Lmodel.type==2 && Lmodel.prep > 0 && Lmodel.prepy > 0,
    
    if debug, disp('mean centering X and Y blocks...........................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
        else
            load([path list{t}],'x','y')
        end
        
        indMV = find(isnan(x));
        if ~isempty(indMV)
            disp('Missing values found in X. Set to average.');
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        indMVy{t} = find(isnan(y));
        if ~isempty(indMVy{t})
            disp('Missing values found in Y. Set to average.');
            avy = ones(size(y,1),1)*Lmodel.avy;
            y(indMVy{t}) = avy(indMVy{t});
        end
        
        [xc,Lmodel.av,Lmodel.sc] = preprocess2Di(x,1,0,1,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);
        [yc,Lmodel.avy,Lmodel.scy,Lmodel.N] = preprocess2Di(y,1,0,1,Lmodel.avy,Lmodel.scy,Lmodel.N,Lmodel.weighty);
        
    end
    
end

% compute scale

N = 0;
    
if (Lmodel.type==1 && Lmodel.prep == 2) || (Lmodel.type==2 && Lmodel.prep == 2 && Lmodel.prepy < 2), 
    
    if debug, disp('scaling X block..................................................'), end;
        
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        xc = x -  ones(size(x,1),1)*Lmodel.av;
        [xsc,av,Lmodel.sc,N] = preprocess2Di(xc,3,0,1,[],Lmodel.sc,N,Lmodel.weight);
        
    end
    
elseif Lmodel.type==2 && Lmodel.prep == 2 && Lmodel.prepy == 2,
    
    if debug, disp('scaling X and Y blocks..................................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
        else
            load([path list{t}],'x','y')
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        if ~isempty(indMVy{t})
            avy = ones(size(y,1),1)*Lmodel.avy;
            y(indMVy{t}) = avy(indMVy{t});
        end
        
        xc = x -  ones(size(x,1),1)*Lmodel.av;
        [xc,av,Lmodel.sc] = preprocess2Di(xc,3,0,1,[],Lmodel.sc,N,Lmodel.weight);
        yc = y -  ones(size(x,1),1)*Lmodel.avy;
        [yc,avy,Lmodel.scy,N] = preprocess2Di(y,1,0,1,Lmodel.avy,Lmodel.scy,N,Lmodel.weighty);
        
    end
    
end

% compute cross-product matrices

Lmodel.XX = zeros(size(x,2));
if Lmodel.type==1, 
    
    if debug, disp('computing XX ....................................................'), end;
        
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        xcs = preprocess2Dapp(x,Lmodel.av,Lmodel.sc,Lmodel.weight);
        Lmodel.XX = Lmodel.XX + xcs'*xcs;
        
    end
    
elseif Lmodel.type==2,
    
    Lmodel.XY = zeros(size(x,2),size(y,2));
    Lmodel.YY = zeros(size(y,2),size(y,2));
    
    if debug, disp('computing XX, XY .......................................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
        else
            load([path list{t}],'x','y')
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        if ~isempty(indMVy{t})
            avy = ones(size(y,1),1)*Lmodel.avy;
            y(indMVy{t}) = avy(indMVy{t});
        end
        
        xcs = preprocess2Dapp(x,Lmodel.av,Lmodel.sc,Lmodel.weight);
        Lmodel.XX = Lmodel.XX + xcs'*xcs;        
        ycs = preprocess2Dapp(y,Lmodel.avy,Lmodel.scy,Lmodel.weighty);
        Lmodel.XY = Lmodel.XY + xcs'*ycs;
        Lmodel.YY = Lmodel.YY + ycs'*ycs;
        
    end
    
end

% compute model
Lmodel.centr = ones(size(Lmodel.centr,1),size(Lmodel.XX,1));
if Lmodel.type==1, 
    
    if debug, disp('computing PCA model..............................................'), end;
    
    if ~Lmodel.lvs,
        Lmodel.lvs = 1:rank(Lmodel.XX);
        var_Lpca(Lmodel);
        Lmodel.lvs = input('Select the PCs to include in the model: ');
        Lmodel.lvs = 1:Lmodel.lvs;
    end
    
    P = Lpca(Lmodel);
    Lmodel.mat = P;
    
elseif Lmodel.type==2,
       
    if rank(Lmodel.XY)>0,
        
        if debug, disp('computing PLS model.....................................................'), end;
            
        if ~Lmodel.lv,
            Lmodel.lvs = 1:rank(Lmodel.XX);
            var_Lpls(Lmodel);
            Lmodel.lv = input('Select the LVs to include in the model: ');
        end
        
        [beta,W,P,Q,R] = Lpls(Lmodel);
        Lmodel.mat = R;
        
    else
        
        if debug>1, disp('XY Rank 0: using PCA.'), end;
        
        if debug, disp('computing PCA model..............................................'), end;
        
        if ~Lmodel.lv,
            Lmodel.lvs = 1:rank(Lmodel.XX);
            var_Lpca(Lmodel);
            Lmodel.lv = input('Select the PCs to include in the model: ');
        end
        
        P = Lpca(Lmodel); 
        Lmodel.mat = P;
        
    end
    
end

% compute maximum and minimum

if debug, disp('computing maximum and minimum ...................................'), end;

mini = Inf(1,size(Lmodel.mat,2));
maxi = -Inf(1,size(Lmodel.mat,2));
for t=1:length(list),
    
    if isstruct(list(t))
        x = list(t).x;
    else
        load([path list{t}],'x')
    end
        
    if ~isempty(indMV{t})
        av = ones(size(x,1),1)*Lmodel.av;
        x(indMV{t}) = av(indMV{t});
    end
    
    xcs = preprocess2Dapp(x,Lmodel.av,Lmodel.sc,Lmodel.weight);
    T = xcs * Lmodel.mat;
    M = max(T);
    m = min(T);
    
    indM = find(maxi < M);
    maxi(indM) = M(indM);
    indm = find(mini > m);
    mini(indm) = m(indm);
    
end
    
mM = maxi-mini;
Lmodel.mat = Lmodel.mat*diag(1./mM);
Lmodel.maxi = maxi;
Lmodel.mini = mini;

% clustering

Lmodel.index_fich={};
for t=1:length(list),
    
    if debug, disp(sprintf('clustering: packet %d...........................................', t)), end;
    
    if isstruct(list(t))
        x = list(t).x;
        vars = fieldnames(list(t));
        if ismember('class', vars)
            class = list(t).class;
        else
            class = ones(size(x,1),1);
        end
        if ismember('obs_l', vars)
            obs_l = list(t).obs_l;
        else
            obs_l = cellstr(num2str((1:size(x,1))'));{};
        end
    else
        load([path list{t}],'x')
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
        
    if ~isempty(indMV{t})
        av = ones(size(x,1),1)*Lmodel.av;
        x(indMV{t}) = av(indMV{t});
    end
    
    xcs = preprocess2Dapp(x,Lmodel.av,Lmodel.sc,Lmodel.weight);
    
    if files, % The updated field is not included in the FS yet
        indorig = length(Lmodel.class);
        red = [Lmodel.centr;xcs];
        lred = {Lmodel.obs_l{:} obs_l{:}};
        multr = [Lmodel.multr;ones(length(class),1)];
        classr = [Lmodel.class;class];
        aux_v = (1:length(classr))';
        obslist = num2cell(aux_v);
    else
        obslist = {};
    end
    
    s = size(x);
    step2 = round(s(1)*step);
    Lmodel.updated(:) = 0;
    for i = 1:step2:s(1),
        endv = min(s(1),i+step2-1);
        ss = endv-i+1;
        xstep = xcs(i:endv,:);
        clstep = class(i:endv);
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
        
        if files,
            for k=i:endv,
                Lmodel.index_fich{1,indorig+k}=['MEDA' num2str(t) 'o' num2str(k) 'c' num2str(class(k))]; %index of names of fich
            end
        end
                   
        [Lmodel.centr,Lmodel.multr,Lmodel.class,Lmodel.obs_l,Lmodel.updated,obslist] = psc(Lmodel.centr,Lmodel.nc,Lmodel.multr,Lmodel.class,Lmodel.obs_l,Lmodel.updated,Lmodel.mat,obslist);
    end
    
    if files,
        Lmodel.index_fich = cfilesys(obslist,red,lred,multr,classr,Lmodel.index_fich,100,Lmodel.path,debug); % update of the clustering file system
    end
      
end

if files,
    ind = find(strcmp(Lmodel.obs_l, 'mixed'));
    Lmodel.obs_l(ind) = Lmodel.index_fich(ind);
end

