function Lmodel = update_iterative(list,path,Lmodel,maxlvs,step,files,path2,debug)

% Big data analysis based on bilinear proyection models (PCA and PLS),
% iterative approach.
%
% Lmodel = update_iterative(list)          % minimum call
% Lmodel = update_iterative(list,path,Lmodel,maxlvs,step,files,path2,debug) % complete call
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
% maxlvs: (1x1) maximum number of LVs considered (e.g. maxlvs = 2 selects the
%   first two LVs)
%
% step: (1x1) percentage of the data in the file to be used in each
%   iteration. For time-course data 1 is suggested (1 by default)
%
% files: (1x1) create the file system with the original data (1, by
%   default) or not (0)
%
% path2: (str) path to the directory where the output data files are
%   located ('' by default)
%
% debug: (1x1) disply debug messages
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
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 11/Apr/14
%
% Copyright (C) 2014  José Camacho Páez
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
    
%% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
if nargin < 2, path = ''; end;
if nargin < 3, 
    Lmodel = Lmodel_ini; 
    Lmodel.type = 1;
    Lmodel.lv = 0;
    Lmodel.prep = 2;
end;
if nargin < 4, maxlvs = 10; end;
if nargin < 5, step = 1; end;
if nargin < 6, files = 1; end;
if nargin < 7, path2 = ''; end;
if nargin < 8, debug = 1; end;
    
    
% Computation

Lmodel.update = 2; 

if files,
	[status,result] = system(['del ' path2 'LEDA*.txt']); % delete previous files
end

% preprocess

% compute mean

if Lmodel.prep==0,
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            class = list(t).class;
        else
            load([path list{t}],'x','class')
        end
        
        [xc,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,0,0,1,Lmodel.av,Lmodel.sc,Lmodel.N);
        
    end
        
elseif (Lmodel.type==1 & Lmodel.prep > 0) | (Lmodel.type==2 & Lmodel.prep > 0 & Lmodel.prepy == 0), 
    
    if debug, disp('mean centering X block...........................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            class = list(t).class;
        else
            load([path list{t}],'x','class')
        end
        
        [xc,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,1,0,1,Lmodel.av,Lmodel.sc,Lmodel.N);
        
    end
    
elseif Lmodel.type==2 & Lmodel.prep > 0 & Lmodel.prepy > 0,
    
    if debug, disp('mean centering X and Y blocks...........................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
            class = list(t).class;
        else
            load([path list{t}],'x','y','class')
        end
        
        [xc,Lmodel.av,Lmodel.sc] = preprocess2Di(x,1,0,1,Lmodel.av,Lmodel.sc,Lmodel.N);
        [yc,Lmodel.avy,Lmodel.scy,Lmodel.N] = preprocess2Di(y,1,0,1,Lmodel.avy,Lmodel.scy,Lmodel.N);
        
    end
    
end

% compute scale

N = 0;
    
if (Lmodel.type==1 & Lmodel.prep == 2) | (Lmodel.type==2 & Lmodel.prep == 2 & Lmodel.prepy < 2), 
    
    if debug, disp('scaling X block..................................................'), end;
        
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            class = list(t).class;
        else
            load([path list{t}],'x','class')
        end
        
        xc = x -  ones(size(x,1),1)*Lmodel.av;
        [xsc,av,Lmodel.sc,N] = preprocess2Di(xc,3,0,1,[],Lmodel.sc,N);
        
    end
    
elseif Lmodel.type==2 & Lmodel.prep == 2 & Lmodel.prepy == 2,
    
    if debug, disp('scaling X and Y blocks..................................................'), end;
    
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
            class = list(t).class;
        else
            load([path list{t}],'x','y','class')
        end
        
        xc = x -  ones(size(x,1),1)*Lmodel.av;
        [xc,av,Lmodel.sc] = preprocess2Di(xc,3,0,1,[],Lmodel.sc,N);
        yc = y -  ones(size(x,1),1)*Lmodel.avy;
        [yc,avy,Lmodel.scy,N] = preprocess2Di(y,1,0,1,Lmodel.avy,Lmodel.scy,N);
        
    end
    
end

% compute cross-product matrices

Lmodel.XX = zeros(size(x,2));
if Lmodel.type==1, 
    
    if debug, disp('computing XX ....................................................'), end;
        
    for t=1:length(list),
        
        if isstruct(list(t))
            x = list(t).x;
            class = list(t).class;
        else
            load([path list{t}],'x','class')
        end
        
        xcs = (x -  ones(size(x,1),1)*Lmodel.av)./(ones(size(x,1),1)*Lmodel.sc);
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
            class = list(t).class;
        else
            load([path list{t}],'x','y','class')
        end
        
        xcs = (x -  ones(size(x,1),1)*Lmodel.av)./(ones(size(x,1),1)*Lmodel.sc);
        Lmodel.XX = Lmodel.XX + xcs'*xcs;        
        ycs = (y -  ones(size(x,1),1)*Lmodel.avy)./(ones(size(x,1),1)*Lmodel.scy);
        Lmodel.XY = Lmodel.XY + xcs'*ycs;
        Lmodel.YY = Lmodel.YY + ycs'*ycs;
        
    end
    
end

% compute model

if Lmodel.type==1, 
    
    if debug, disp('computing PCA model..............................................'), end;
    
    if ~Lmodel.lv,
        var_Lpca(Lmodel,maxlvs);
        Lmodel.lv = input('Select the PCs to include in the model: ');
    end
    
    [P,sdT] = Lpca(Lmodel);
    Lmodel.mat = P;
    
elseif Lmodel.type==2,
       
    if rank(Lmodel.XY)>0,
        
        if debug, disp('computing PLS model.....................................................'), end;
            
        if ~Lmodel.lv,
            var_Lpls(Lmodel,maxlvs);
            Lmodel.lv = input('Select the LVs to include in the model: ');
        end
        
        [beta,W,P,Q,R] = Lpls(Lmodel);
        Lmodel.mat = R;
        
    else
        
        if debug>1, disp('XY Rank 0: using PCA.'), end;
        
        if debug, disp('computing PCA model..............................................'), end;
        
        if ~Lmodel.lv,
            var_Lpca(Lmodel,maxlvs);
            Lmodel.lv = input('Select the PCs to include in the model: ');
        end
        
        P = Lpca(Lmodel); 
        Lmodel.mat = P;
        
    end
    
end

% compute maximum and minimum

if debug, disp('computing maximum and minimum ..............................................'), end;

mini = Inf(1,size(Lmodel.mat,2));
maxi = -Inf(1,size(Lmodel.mat,2));
for t=1:length(list),
    
    if isstruct(list(t))
        x = list(t).x;
        class = list(t).class;
    else
        load([path list{t}],'x','class')
    end
    
    xcs = (x -  ones(size(x,1),1)*Lmodel.av)./(ones(size(x,1),1)*Lmodel.sc);
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

index_fich={};
for t=1:length(list),
    
    if debug, disp(sprintf('clustering: packet %d...........................................', t)), end;
    
    if isstruct(list(t))
        x = list(t).x;
        class = list(t).class;
    else
        load([path list{t}],'x','class')
    end
    
    xcs = (x -  ones(size(x,1),1)*Lmodel.av)./(ones(size(x,1),1)*Lmodel.sc);
   
    if files, 
        indorig = length(Lmodel.class);
        red = [Lmodel.centr;xcs];
        multr = [Lmodel.multr;ones(length(class),1)];
        labr = [Lmodel.class;class];
        aux_v = (1:length(labr))';
        obslist = num2cell(aux_v);
    else
        obslist = [];
    end
    
    s = size(x);
    step2 = round(s(1)*step);
    for i = 1:step2:s(1),
        endv = min(s(1),i+step2-1);
        ss = endv-i+1;
        xstep = xcs(i:endv,:);
        clstep = class(i:endv,:);
        
        Lmodel.centr = [Lmodel.centr;xstep];
        Lmodel.multr = [Lmodel.multr;ones(ss,1)];
        Lmodel.class = [Lmodel.class;clstep];
        
        if files,
            for k=i:endv,
                index_fich{1,indorig+k}=['LEDA_' num2str(t) '_o_' num2str(k) '_c_' num2str(class(k))]; %index of names of fich
            end
        end
        
        [Lmodel.centr,Lmodel.multr,Lmodel.class,obslist] = psc(Lmodel.centr,Lmodel.nc,Lmodel.multr,Lmodel.class,Lmodel.mat,obslist);
    end
    
    if files,
        index_fich = cfilesys(obslist,red,multr,labr,index_fich,100,path2,debug); % update of the clustering file system
    end
      
end

if files,
    Lmodel.index_fich = index_fich;
    Lmodel.path = path2;
end
