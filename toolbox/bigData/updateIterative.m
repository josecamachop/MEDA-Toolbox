function Lmodel = updateIterative(list,path,Lmodel,step,files,debug)

% Big data analysis based on bilinear proyection models (PCA, PLS and ASCA),
% iterative approach.
%
% Lmodel = updateIterative(list)          % minimum call
% Lmodel = updateIterative(list,path,Lmodel,step,files,debug) % complete call
%
%
% INPUTS:
%
% list: {Fx1} list of strings with the names of the files for the update or
%   struct array with x (and optionally y) matrices.
%
% path: (str) If list is a string of filenames, path to the directory where 
%   the data files are located ('' by default). Not to confuse with where
%   the output files are stored (which is Lmodel.path)
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
% EXAMPLE OF USE: to model a large number of observations.
%
% nobs = 100;
% nvars = 10;
% Lmodel = iniLmodel;
% Lmodel.type = 'PCA'; 
% Lmodel.prep = 2;  
% Lmodel.lvs = 1;
% Lmodel.nc = 100; % Number of clusters
% 
% for i=1:10,
%   list(i).x = simuleMV(nobs,nvars,'LevelCorr',6);
% end
%
% Lmodel = updateIterative(list,[],Lmodel);
% mspcLpca(Lmodel);
%
%
% EXAMPLE OF USE: to model a large number of observations creating a
% filesystem with debug (level 2)
%
% nobs = 100;
% nvars = 10;
% Lmodel = iniLmodel;
% Lmodel.type = 'PCA'; 
% Lmodel.prep = 2;  
% Lmodel.lvs = 1;
% Lmodel.nc = 100; % Number of clusters
% Lmodel.path = '.\deleteMe\'; % MODIFY PATH TO YOUR CONVENIENCE
% 
% for i=1:10,
%   list(i).x = simuleMV(nobs,nvars,'LevelCorr',6);
% end
%
% Lmodel = updateIterative(list,[],Lmodel,0.1,1,2);
% mspcLpca(Lmodel);
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 19/Feb/2024
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

routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

if nargin < 2 || isempty(path), path = ''; end;
if nargin < 3 || isempty(Lmodel) 
    Lmodel = iniLmodel; 
    Lmodel.type = 'PCA';
    Lmodel.lvs = 0;
    Lmodel.prep = 2;
end;
[ok, Lmodel] = checkLmodel(Lmodel);
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

if files == 1
    if isempty(Lmodel.path)
        Lmodel.path = '.\updtIterOutput\';
    end
    mkdir(Lmodel.path);
end
    
%% Main code

Lmodel.update = 2; 

if files
  if ~strcmp(Lmodel.type,'ASCA')
      if ispc
          [status,result] = system(['del ' Lmodel.path 'MEDA*.txt']); % delete previous files
      else
          [status,result] = system(['rm ' Lmodel.path 'MEDA*.txt']); % delete previous files
      end
  else
      for f = 1:length(LmodelASCA.factors)
          if ispc
              [status,result] = system(['del ' Lmodel.factors{f}.path 'MEDA*.txt']); % delete previous files
          else
              [status,result] = system(['rm ' Lmodel.factors{f}.path 'MEDA*.txt']); % delete previous files
          end
      end 
  end
end

if strcmp(Lmodel.type,'ASCA') % build coding matrix
    
    if debug, disp('preparing coding matrix..................................................'), end;
            
    if ~isfield(Lmodel,'anovast'), Lmodel.anovast = []; end
    if ~isfield(Lmodel.anovast,'model'), Lmodel.anovast.model = []; end
    if ~isfield(Lmodel.anovast,'ordinal'), Lmodel.anovast.ordinal = []; end
    if ~isfield(Lmodel.anovast,'coding'), Lmodel.anovast.coding = []; end
    if ~isfield(Lmodel.anovast,'nested'), Lmodel.anovast.nested = []; end
    
    if isstruct(list(1)) % 
        f = list(1).f;
    else
        load([path list{1}],'f')
    end
    for i = 1:size(f,2), F{i} = []; end;
    for t=1:length(list)
        
        if isstruct(list(t)) % compute global design matrix
            f = list(t).f;
        else
            load([path list{t}],'f')
        end
        
        for i = 1:size(f,2)
            if isempty(find(i==Lmodel.anovast.ordinal))
                if isempty(Lmodel.anovast.nested) || isempty(find(i==Lmodel.anovast.nested(:,2)))
                    F{i} = unique([F{i};f(:,i)]);
                else
                    ref = Lmodel.anovast.nested(find(i==Lmodel.anovast.nested(:,2)),1);
                    F{i} = unique([F{i};f(:,[ref i])],'rows');
                end
            end
        end  
        
    end

    Lmodel.nFactors = length(F);
    Lmodel.levels = F;
        
    for t=1:length(list)
        
        if isstruct(list(t)) % change design to input
            f = list(t).f;
        else
            load([path list{t}],'f')
        end        
        
        [d, Lmodel] = codglm(f, Lmodel);
        
        if isstruct(list(t)) 
            list(t).d = d;
        else
            save([path list{t}],'d','-APPEND')
        end
    end

end
    
% preprocess

% compute mean

if strcmp(Lmodel.type,'PCA') || strcmp(Lmodel.type,'ASCA')
    
    if debug, disp('mean centering X block..................................................'), end;
    
    for t=1:length(list)
        
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
        
        [xc,Lmodel.av,Lmodel.sc,Lmodel.N] = preprocess2Di(x,Lmodel.prep>0,0,1,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);
        
    end
        
elseif strcmp(Lmodel.type,'PLS')
    
    if debug, disp('mean centering X and Y blocks...........................................'), end;
    
    for t=1:length(list)
        
        if isstruct(list(t))
            x = list(t).x;
            y = list(t).y;
        else
            load([path list{t}],'x','y')
        end
        
        indMV{t} = find(isnan(x));
        if ~isempty(indMV{t})
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
        
        [xc,Lmodel.av,Lmodel.sc] = preprocess2Di(x,Lmodel.prep>0,0,1,Lmodel.av,Lmodel.sc,Lmodel.N,Lmodel.weight);
        [yc,Lmodel.avy,Lmodel.scy,Lmodel.N] = preprocess2Di(y,Lmodel.prepy>0,0,1,Lmodel.avy,Lmodel.scy,Lmodel.N,Lmodel.weighty);
        
    end
    
end

% compute scale

N = 0;
    
if ((strcmp(Lmodel.type,'PCA') || strcmp(Lmodel.type,'ASCA')) && Lmodel.prep == 2) || (strcmp(Lmodel.type,'PLS') && Lmodel.prep == 2 && Lmodel.prepy < 2) 
    
    if debug, disp('scaling X block..................................................'), end;
        
    for t=1:length(list)
        
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
    
elseif strcmp(Lmodel.type,'PLS') && Lmodel.prep == 2 && Lmodel.prepy == 2
    
    if debug, disp('scaling X and Y blocks..................................................'), end;
    
    for t=1:length(list)
        
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

if strcmp(Lmodel.type,'PCA') 
    
    Lmodel.XX = zeros(size(x,2));

    if debug, disp('computing XX ....................................................'), end;
        
    for t=1:length(list)
        
        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end
        
        if ~isempty(find(isnan(x))) 
                if debug, disp(sprintf('Found nans in file %d',t)), end;
        end;
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        xcs = preprocess2Dapp(x,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
        Lmodel.XX = Lmodel.XX + xcs'*xcs;
        
    end
    
elseif strcmp(Lmodel.type,'PLS') 
    
    Lmodel.XX = zeros(size(x,2));
    Lmodel.XY = zeros(size(x,2),size(y,2));
    Lmodel.YY = zeros(size(y,2),size(y,2));
    
    if debug, disp('computing XX, XY .......................................................'), end;
    
    for t=1:length(list)
        
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
        
        xcs = preprocess2Dapp(x,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
        Lmodel.XX = Lmodel.XX + xcs'*xcs;        
        ycs = preprocess2Dapp(y,Lmodel.avy,'Scale',Lmodel.scy,'Weight',Lmodel.weighty);
        Lmodel.XY = Lmodel.XY + xcs'*ycs;
        Lmodel.YY = Lmodel.YY + ycs'*ycs;
        
    end
    
 elseif strcmp(Lmodel.type,'ASCA') 
    
    Lmodel.DD = zeros(size(d,2));
    Lmodel.DX = zeros(size(d,2),size(x,2));
    
    if debug, disp('computing DD, DX .......................................................'), end;
    
    for t=1:length(list)
        
        if isstruct(list(t))
            d = list(t).d;
            x = list(t).x;
        else
            load([path list{t}],'d','x')
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        Lmodel.DD = Lmodel.DD + d'*d;        
        xcs = preprocess2Dapp(x,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
        Lmodel.DX = Lmodel.DX + d'*xcs;
        
    end
    
    if debug, disp('computing parglm......................................................'), end;
    
    nFactors = Lmodel.nFactors;
    nInteractions = Lmodel.nInteractions;
    
    % Degrees of freedom
    Tdf = Lmodel.N;
    Rdf = Tdf-1;
    for f = 1 : nFactors
        if Lmodel.anovast.ordinal(f)
            df(f) = 1;
        else
            df(f) = length(Lmodel.factors{f}.Dvars);
        end
        Rdf = Rdf-df(f);
    end
    dfint = [];
    for i = 1 : nInteractions
        dfint(i) = prod(df(Lmodel.interactions{i}.factors));
        Rdf = Rdf-dfint(i);
    end
    if Rdf < 0
        disp('Warning: degrees of freedom exhausted');
        return
    end
    
    % GLM model calibration with LS, only fixed factors
    Lmodel.B = pinv(Lmodel.DD)*Lmodel.DX;
    
    SSQX = 0;
    SSQresiduals = 0;
    SSQinter = 0;
    SSQFactors = [];
    for f = 1 : nFactors, SSQFactors(f) = 0; Lmodel.factors{f}.XX = zeros(size(x,2)); end
    SSQInteractions = [];
    for i = 1 : nInteractions, SSQInteractions(i) = 0; Lmodel.interactions{i}.XX = zeros(size(x,2)); end
    
    for t=1:length(list)
        
        if isstruct(list(t))
            d = list(t).d;
            x = list(t).x;
        else
            load([path list{t}],'d','x')
        end
        
        x = x./(ones(size(x,1),1)*Lmodel.sc);

        SSQX = SSQX + sum(sum(x.^2));
        
        parglmo.residuals = x - d*Lmodel.B;
        SSQresiduals = SSQresiduals + sum(sum(parglmo.residuals.^2));

        % Create Effect Matrices
        parglmo.inter = d(:,1)*Lmodel.B(1,:);
        SSQinter = SSQinter + sum(sum(parglmo.inter.^2));
        
        % Factors
        for f = 1 : nFactors
            parglmo.factors{f}.matrix = d(:,Lmodel.factors{f}.Dvars)*Lmodel.B(Lmodel.factors{f}.Dvars,:);
            Lmodel.factors{f}.XX = Lmodel.factors{f}.XX + parglmo.factors{f}.matrix'*parglmo.factors{f}.matrix;
            SSQFactors(f) = SSQFactors(f) + sum(sum(parglmo.factors{f}.matrix.^2)); 
        end
        
        % Interactions
        for i = 1 : nInteractions
            parglmo.interactions{i}.matrix = d(:,Lmodel.interactions{i}.Dvars)*Lmodel.B(Lmodel.interactions{i}.Dvars,:);
            Lmodel.interactions{i}.XX = Lmodel.interactions{i}.XX + parglmo.interactions{i}.matrix'*parglmo.interactions{i}.matrix;
            SSQInteractions(i) = SSQInteractions(i) + sum(sum(parglmo.interactions{i}.matrix.^2));
        end
        
%         if isstruct(list(t))
%             list(t).parglmo = parglmo;
%         else
%             load([path list{t}],'parglmo')
%         end
  
    end
    
    % Factors
    FFactors = [];
    for f = 1 : nFactors
        FFactors(f) = (SSQFactors(f)/df(f))/(SSQresiduals/Rdf);
    end
    
    % Interactions
    FInteractions = [];
    for i = 1 : nInteractions
        FInteractions(i) = (SSQInteractions(i)/dfint(i))/(SSQresiduals/Rdf);
    end
    
    Lmodel.effects = 100*([SSQinter SSQFactors SSQInteractions SSQresiduals]./SSQX);
        
    %% ANOVA-like output table
    
    name={'Mean'};
    for f = 1 : nFactors
        name{end+1} = sprintf('Factor %d',f);
    end
    for i = 1 : nInteractions
        name{end+1} = sprintf('Interaction %s',strrep(num2str(parglmo.interactions{i}.factors),'  ','-'));
    end
    name{end+1} = 'Residuals';
    name{end+1} = 'Total';
    
    SSQ = [SSQinter SSQFactors SSQInteractions SSQresiduals SSQX];
    par = [Lmodel.effects 100];
    DoF = [1 df dfint Rdf Tdf];
    MSQ = SSQ./DoF;
    F = [nan FFactors FInteractions nan nan];
    %pvalue = [nan parglmo.p nan nan];
    
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        Lmodel.T.mat = [SSQ', par', DoF', MSQ', F'];
        Lmodel.T.var = {'SumSq', 'PercSumSq', 'df', 'MeanSq', 'F'};
        Lmodel.T.source = name';
    else
        Lmodel.T = table(name', SSQ', par', DoF', MSQ', F','VariableNames', {'Source','SumSq','PercSumSq','df','MeanSq','F'});
    end
    
end

% compute model
if strcmp(Lmodel.type,'PCA') 
    
    Lmodel.centr = ones(size(Lmodel.centr,1),size(Lmodel.XX,1));
    
    if debug, disp('computing PCA model..............................................'), end;
    
    if isempty(Lmodel.lvs)
        Lmodel.lvs = 1:rank(Lmodel.XX);
        varLpca(Lmodel);
        Lmodel.lvs = input('Select the number of PCs to include in the model: ');
        Lmodel.lvs = 1:Lmodel.lvs;
    end
    
    Lmodel = Lpca(Lmodel);
    Lmodel.mat = Lmodel.loads;
    
elseif strcmp(Lmodel.type,'PLS')
       
    Lmodel.centr = ones(size(Lmodel.centr,1),size(Lmodel.XX,1));
    if rank(Lmodel.XY)>0
        
        if debug, disp('computing PLS model.....................................................'), end;
            
        if isempty(Lmodel.lvs)
            Lmodel.lvs = 1:rank(Lmodel.XX);
            varLpls(Lmodel);
            Lmodel.lvs = input('Select the number of LVs to include in the model: ');
            Lmodel.lvs = 1:Lmodel.lvs;
        end
        
        Lmodel = Lpls(Lmodel);
        R = Lmodel.altweights;
        Lmodel.mat = R;
        
    else
        
        if debug>1, disp('XY Rank 0: using PCA.'), end;
        
        if debug, disp('computing PCA model..............................................'), end;
        
        if isempty(Lmodel.lvs)
            Lmodel.lvs = 1:rank(Lmodel.XX);
            varLpca(Lmodel);
            Lmodel.lvs = input('Select the number of PCs to include in the model: ');
            Lmodel.lvs = 1:Lmodel.lvs;
        end
        
        Lmodel = Lpca(Lmodel);
        Lmodel.mat = Lmodel.loads;
        
    end
    
elseif strcmp(Lmodel.type,'ASCA')
    
    if debug, disp('computing ASCA model....................................................'), end;
    
    nFactors = Lmodel.nFactors;
    nInteractions = Lmodel.nInteractions;
    
    % Factors
    for f = 1 : nFactors    
        Lmodel.factors{f}.centr = ones(size(Lmodel.centr,1),size(Lmodel.factors{f}.XX,1));
    
        if debug, disp(sprintf('computing PCA model of factor %d..............................................',f)); end;

        if isempty(Lmodel.factors{f}.lvs)
            Lmodel.factors{f}.lvs = 1:rank(Lmodel.factors{f}.XX);
            varLpca(Lmodel.factors{f});
            Lmodel.factors{f}.lvs = input('Select the number of PCs to include in the model: ');
            Lmodel.factors{f}.lvs = 1:Lmodel.factors{f}.lvs;
        end
        
        Lmodel.factors{f} = Lpca(Lmodel.factors{f});
        Lmodel.factors{f}.mat = Lmodel.loads;
    end
    
    % Interactions
    for i = 1 : nInteractions
        Lmodel.interactions{i}.centr = ones(size(Lmodel.centr,1),size(Lmodel.interactions{i}.XX,1));
    
        if debug, disp(sprintf('computing PCA model of interaction %d..............................................',i)); end;
        
        for f = 1 : length(Lmodel.interactions{i}.factors) % Combine Factors and Interaction
            Lmodel.interactions{i}.XX = Lmodel.interactions{i}.XX + Lmodel.factors{Lmodel.interactions{i}.factors(f)}.XX;
        end

        if isempty(Lmodel.interactions{i}.lvs)
            Lmodel.interactions{i}.lvs = 1:rank(Lmodel.interactions{i}.XX);
            varLpca(Lmodel.interactions{i});
            Lmodel.interactions{i}.lvs = input('Select the PCs number of to include in the model: ');
            Lmodel.interactions{i}.lvs = 1:Lmodel.interactions{i}.lvs;
        end
        
        Lmodel.interactions{i} = Lpca(Lmodel.interactions{i});
        Lmodel.interactions{i}.mat = Lmodel.loads;
    end
    
end

% compute maximum and minimum

if debug, disp('computing maximum and minimum ...................................'), end;

if strcmp(Lmodel.type,'PCA') || strcmp(Lmodel.type,'PLS')

    mini = Inf(1,size(Lmodel.mat,2));
    maxi = -Inf(1,size(Lmodel.mat,2));
    for t=1:length(list)

        if isstruct(list(t))
            x = list(t).x;
        else
            load([path list{t}],'x')
        end

        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end

        xcs = preprocess2Dapp(x,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
        
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

elseif strcmp(Lmodel.type,'ASCA')
    
    nFactors = Lmodel.nFactors;
    nInteractions = Lmodel.nInteractions;
    
    % Factors
    for f = 1 : nFactors 
        
        mini = Inf(1,size(Lmodel.factors{f}.mat,2));
        maxi = -Inf(1,size(Lmodel.factors{f}.mat,2));
        for t=1:length(list)

            if isstruct(list(t))
                d = list(t).d;
            else
                load([path list{t}],'d')
            end
            
            if ~isempty(indMV{t})
                av = ones(size(x,1),1)*Lmodel.av;
                x(indMV{t}) = av(indMV{t});
            end
            
            xcs = d(:,Lmodel.factors{f}.Dvars)*Lmodel.B(Lmodel.factors{f}.Dvars,:);
            T = xcs * Lmodel.factors{f}.mat;
            M = max(T);
            m = min(T);

            indM = find(maxi < M);
            maxi(indM) = M(indM);
            indm = find(mini > m);
            mini(indm) = m(indm);

        end

        mM = maxi-mini;
        Lmodel.factors{f}.mat = Lmodel.factors{f}.mat*diag(1./mM);
        Lmodel.factors{f}.maxi = maxi;
        Lmodel.factors{f}.mini = mini;
        
    end
    
    % Interactions
    for i = 1 : nInteractions
        
        mini = Inf(1,size(Lmodel.interactions{i}.mat,2));
        maxi = -Inf(1,size(Lmodel.interactions{i}.mat,2));
        for t=1:length(list)

            if isstruct(list(t))
                d = list(t).d;
            else
                load([path list{t}],'d')
            end
            
            if ~isempty(indMV{t})
                av = ones(size(x,1),1)*Lmodel.av;
                x(indMV{t}) = av(indMV{t});
            end
            
            xcs = d(:,Lmodel.interactions{i}.Dvars)*Lmodel.B(Lmodel.interactions{i}.Dvars,:);
            T = xcs * Lmodel.interactions{i}.mat;
            M = max(T);
            m = min(T);

            indM = find(maxi < M);
            maxi(indM) = M(indM);
            indm = find(mini > m);
            mini(indm) = m(indm);

        end

        mM = maxi-mini;
        Lmodel.interactions{i}.mat = Lmodel.interactions{i}.mat*diag(1./mM);
        Lmodel.interactions{i}.maxi = maxi;
        Lmodel.interactions{i}.mini = mini;
        
    end
    
end

% clustering


if strcmp(Lmodel.type,'PCA') || strcmp(Lmodel.type,'PLS')
    
    Lmodel.indexFich={};
    for t=1:length(list)
        
        if debug, disp(sprintf('clustering: packet %d...........................................', t)), end;
        
        if isstruct(list(t))
            x = list(t).x;
            vars = fieldnames(list(t));
            if ismember('class', vars)
                class = list(t).class;
            else
                class = ones(size(x,1),1);
            end
            if ismember('obsl', vars)
                obsl = list(t).obsl;
                if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
            else
                obsl = cellstr(num2str((1:size(x,1))'));
                for o = 1:length(obsl), obsl{o} = [num2str(t) '-' strtrim(obsl{o})]; end
            end
        else
            load([path list{t}],'x')
            vars = whos('-file',[path list{t}]);
            if ismember('class', {vars.name})
                load([path list{t}],'class')
            else
                class = ones(size(x,1),1);
            end
            if ismember('obsl', {vars.name})
                load([path list{t}],'obsl')
                if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
            else
                obsl = cellstr(num2str((1:size(x,1))'));
                for o = 1:length(obsl), obsl{o} = [num2str(t) '-' strtrim(obsl{o})]; end
            end
        end
        
        if ~isempty(indMV{t})
            av = ones(size(x,1),1)*Lmodel.av;
            x(indMV{t}) = av(indMV{t});
        end
        
        xcs = preprocess2Dapp(x,Lmodel.av,'Scale',Lmodel.sc,'Weight',Lmodel.weight);
        
        if files % The updated field is not included in the FS yet
            indorig = length(Lmodel.class);
            red = [Lmodel.centr;xcs];
            lred = {Lmodel.obsl{:} obsl{:}};
            multr = [Lmodel.multr;ones(length(class),1)];
            classr = [Lmodel.class;class];
            auxv = (1:length(classr))';
            obslist = num2cell(auxv);
        else
            obslist = {};
        end
        
        s = size(x);
        step2 = max(100,round(s(1)*step)); % Min step of 100 obs
        Lmodel.updated(:) = 0;
        for i = 1:step2:s(1)
            endv = min(s(1),i+step2-1);
            ss = endv-i+1;
            xstep = xcs(i:endv,:);
            clstep = class(i:endv);
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
            
            if files
                for k=i:endv
                    Lmodel.indexFich{1,indorig+k}=['MEDA' num2str(t) 'o' num2str(k) 'c' num2str(class(k))]; %index of names of fich
                end
            end
            
            [Lmodel.centr,Lmodel.multr,Lmodel.class,Lmodel.obsl,Lmodel.updated,obslist] = psc(Lmodel.centr,Lmodel.nc,Lmodel.multr,Lmodel.class,Lmodel.obsl,Lmodel.updated,Lmodel.mat,obslist);
        end
        
        if files
            Lmodel.indexFich = cfilesys(obslist,red,lred,multr,classr,Lmodel.indexFich,100,Lmodel.path,debug); % update of the clustering file system
        end
        
    end
    
    if files
        ind = find(strcmp(Lmodel.obsl, 'mixed'));
        Lmodel.obsl(ind) = Lmodel.indexFich(ind);
    end

elseif strcmp(Lmodel.type,'ASCA')
    
    nFactors = Lmodel.nFactors;
    nInteractions = Lmodel.nInteractions;
    
    % Factors
    for fa = 1 : nFactors 
        
        Lmodel.factors{fa}.indexFich={};
        for t=1:length(list)
            
            if debug, disp(sprintf('clustering factor %d: packet %d...........................................', fa, t)), end;
            
            if isstruct(list(t))
                d = list(t).d;
                x = list(t).x;
                vars = fieldnames(list(t));
                class = list(t).f(:,fa);
                if ismember('obsl', vars)
                    obsl = list(t).obsl;
                    if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
                else
                    obsl = cellstr(num2str((1:size(d,1))'));
                end
            else
                load([path list{t}],'d','x')
                vars = whos('-file',[path list{t}]);
                load([path list{t}],'f');
                class = f(:,fa);
                if ismember('obsl', {vars.name})
                    load([path list{t}],'obsl')
                    if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
                else
                    obsl = cellstr(num2str((1:size(d,1))'));
                end
            end
            f = fa;
            
            x = x./(ones(size(x,1),1)*Lmodel.sc);
            residuals = x - d*Lmodel.B;
            xcs = d(:,Lmodel.factors{f}.Dvars)*Lmodel.B(Lmodel.factors{f}.Dvars,:) + residuals;
            
            
            if files % The updated field is not included in the FS yet
                indorig = length(Lmodel.factors{f}.class);
                red = [Lmodel.factors{f}.centr;xcs];
                lred = {Lmodel.factors{f}.obsl{:} obsl{:}};
                multr = [Lmodel.factors{f}.multr;ones(length(class),1)];
                classr = [Lmodel.factors{f}.class;class];
                auxv = (1:length(classr))';
                obslist = num2cell(auxv);
            else
                obslist = {};
            end
            
            s = size(x);
            step2 = round(s(1)*step);
            Lmodel.factors{f}.updated(:) = 0;
            for i = 1:step2:s(1)
                endv = min(s(1),i+step2-1);
                ss = endv-i+1;
                xstep = xcs(i:endv,:);
                clstep = class(i:endv);
                if isempty(obsl)
                    obsstep = {};
                else
                    obsstep = obsl(i:endv);
                end
                
                Lmodel.factors{f}.centr = [Lmodel.factors{f}.centr;xstep];
                Lmodel.factors{f}.multr = [Lmodel.factors{f}.multr;ones(ss,1)];
                Lmodel.factors{f}.class = [Lmodel.factors{f}.class;clstep];
                Lmodel.factors{f}.obsl = {Lmodel.factors{f}.obsl{:} obsstep{:}};
                Lmodel.factors{f}.updated = [Lmodel.factors{f}.updated;ones(size(xstep,1),1)];
                
                if files
                    for k=i:endv
                        Lmodel.factors{f}.indexFich{1,indorig+k}=['MEDA' num2str(t) 'o' num2str(k) 'c' num2str(class(k))]; %index of names of fich
                    end
                end
                
                [Lmodel.factors{f}.centr,Lmodel.factors{f}.multr,Lmodel.factors{f}.class,Lmodel.factors{f}.obsl,Lmodel.factors{f}.updated,obslist] = psc(Lmodel.factors{f}.centr,Lmodel.factors{f}.nc,Lmodel.factors{f}.multr,Lmodel.factors{f}.class,Lmodel.factors{f}.obsl,Lmodel.factors{f}.updated,Lmodel.factors{f}.mat,obslist);
            end
            
            if files
                Lmodel.factors{f}.indexFich = cfilesys(obslist,red,lred,multr,classr,Lmodel.factors{f}.indexFich,100,Lmodel.factors{f}.path,debug); % update of the clustering file system
            end
            
        end
        
        if files
            ind = find(strcmp(Lmodel.factors{f}.obsl, 'mixed'));
            Lmodel.factors{f}.obsl(ind) = Lmodel.factors{f}.indexFich(ind);
        end
    
    end
    
    % Interactions
    for i = 1 : nInteractions 
        
        Lmodel.interactions{i}.indexFich={};
        for t=1:length(list)
            
            if debug, disp(sprintf('clustering interaction %d: packet %d...........................................', i, t)), end;
            
            if isstruct(list(t))
                d = list(t).d;
                x = list(t).x;
                vars = fieldnames(list(t));
                if ismember('class', vars)
                    class = list(t).class;
                else
                    class = ones(size(x,1),1);
                end
                if ismember('obsl', vars)
                    obsl = list(t).obsl;
                    if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
                else
                    obsl = cellstr(num2str((1:size(d,1))'));
                end
            else
                load([path list{t}],'d','x')
                vars = whos('-file',[path list{t}]);
                if ismember('class', {vars.name})
                    load([path list{t}],'class')
                else
                    class = ones(size(x,1),1);
                end
                if ismember('obsl', {vars.name})
                    load([path list{t}],'obsl')
                    if isnumeric(obsl), obsl = cellstr(num2str(obsl)); end
                else
                    obsl = cellstr(num2str((1:size(d,1))'));
                end
            end
            
            x = x./(ones(size(x,1),1)*Lmodel.sc);
            residuals = x - d*Lmodel.B;
            xcs = d(:,Lmodel.interactions{i}.Dvars)*Lmodel.B(Lmodel.interactions{i}.Dvars,:) + residuals;
            
            if files % The updated field is not included in the FS yet
                indorig = length(Lmodel.interactions{i}.class);
                red = [Lmodel.interactions{i}.centr;xcs];
                lred = {Lmodel.interactions{i}.obsl{:} obsl{:}};
                multr = [Lmodel.interactions{i}.multr;ones(length(class),1)];
                classr = [Lmodel.interactions{i}.class;class];
                auxv = (1:length(classr))';
                obslist = num2cell(auxv);
            else
                obslist = {};
            end
            
            s = size(x);
            step2 = round(s(1)*step);
            Lmodel.interactions{i}.updated(:) = 0;
            for i = 1:step2:s(1)
                endv = min(s(1),i+step2-1);
                ss = endv-i+1;
                xstep = xcs(i:endv,:);
                clstep = class(i:endv);
                if isempty(obsl)
                    obsstep = {};
                else
                    obsstep = obsl(i:endv);
                end
                
                Lmodel.interactions{i}.centr = [Lmodel.interactions{i}.centr;xstep];
                Lmodel.interactions{i}.multr = [Lmodel.interactions{i}.multr;ones(ss,1)];
                Lmodel.interactions{i}.class = [Lmodel.interactions{i}.class;clstep];
                Lmodel.interactions{i}.obsl = {Lmodel.interactions{i}.obsl{:} obsstep{:}};
                Lmodel.interactions{i}.updated = [Lmodel.interactions{i}.updated;ones(size(xstep,1),1)];
                
                if files
                    for k=i:endv
                        Lmodel.interactions{i}.indexFich{1,indorig+k}=['MEDA' num2str(t) 'o' num2str(k) 'c' num2str(class(k))]; %index of names of fich
                    end
                end
                
                [Lmodel.interactions{i}.centr,Lmodel.interactions{i}.multr,Lmodel.interactions{i}.class,Lmodel.interactions{i}.obsl,Lmodel.interactions{i}.updated,obslist] = psc(Lmodel.interactions{i}.centr,Lmodel.interactions{i}.nc,Lmodel.interactions{i}.multr,Lmodel.interactions{i}.class,Lmodel.interactions{i}.obsl,Lmodel.interactions{i}.updated,Lmodel.interactions{i}.mat,obslist);
            end
            
            if files
                Lmodel.interactions{i}.indexFich = cfilesys(obslist,red,lred,multr,classr,Lmodel.interactions{i}.indexFich,100,Lmodel.interactions{i}.path,debug); % update of the clustering file system
            end
            
        end
        
        if files
            ind = find(strcmp(Lmodel.interactions{i}.obsl, 'mixed'));
            Lmodel.interactions{i}.obsl(ind) = Lmodel.interactions{i}.indexFich(ind);
        end
    
    end
    
end
        
end

%%

function [D, Lmodel] = codglm(F, Lmodel)

% Compute coding matrix from a design matrix for General Linear Models.
%
% Related routines: parglm, asca, apca, parglmVS, parglmMC, createDesign
%
% D = codglm(F, Lmodel)   % minimum call
% [D, Lmodel] = codglm(F, Lmodel)   % complete call
%
%
% INPUTS:
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels.
%
% Lmodel (structure)
%
%   parglmi (structure): structure with the number of factors and the number 
%   of levels for each of them.
%
%   anovast (structure): structure with the anova choices.
%
%       model: This paremeter is similar to 'model' of anovan. It could be:
%           'linear': only main effects are provided (by default)
%           'interaction': two order interactions are provided
%           'full': all potential interactions are provided
%           [1x1]: maximum order of interactions considered
%           [ix2]: array with two order interactions
%           cell: with each element a vector of factors
%
%       ordinal: [1xF] whether factors are nominal or ordinal
%           0: nominal (default)
%           1: ordinal
%
%       coding: [1xF] type of coding of factors
%           0: sum/deviation coding (default)
%           1: reference coding (reference is the last level)
%
%       nested: [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%       and 3 in 2, then nested = [1 2; 2 3]
%
%
% OUTPUTS:
%
% D [NxC]: Coding matrix
%
% Lmodel (structure)
%
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
% last modification: 16/Jun/23
%
% Copyright (C) 2023  Universidad de Granada
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

nFactors = Lmodel.nFactors;                 % number of factors
levels = Lmodel.levels;

if ~isfield(Lmodel,'anovast'), Lmodel.anovast = []; end
    
if ~isfield(Lmodel.anovast,'model') || isempty(Lmodel.anovast.model), Lmodel.anovast.model = 'linear'; end;

if isequal(Lmodel.anovast.model,'linear')
    interactions = [];
end  
    
if isequal(Lmodel.anovast.model,'interaction')
    interactions = allinter(nFactors,2);
end    

if isequal(Lmodel.anovast.model,'full')
    interactions = allinter(nFactors,nFactors);
end    

if isnumeric(Lmodel.anovast.model) && isscalar(Lmodel.anovast.model) && Lmodel.anovast.model >= 2 && Lmodel.anovast.model <= nFactors
        interactions = allinter(nFactors,Lmodel.anovast.model);
end    

if isnumeric(Lmodel.anovast.model) && ~isscalar(Lmodel.anovast.model)
        interactions = {Lmodel.anovast.model};
end    

if iscell(Lmodel.anovast.model), interactions = Lmodel.anovast.model; end
    
if ~isfield(Lmodel.anovast,'ordinal') || isempty(Lmodel.anovast.ordinal), Lmodel.anovast.ordinal = zeros(1,nFactors); end;
if ~isfield(Lmodel.anovast,'coding') || isempty(Lmodel.anovast.coding), Lmodel.anovast.coding = zeros(1,nFactors); end;
if ~isfield(Lmodel.anovast,'nested'), Lmodel.anovast.nested = []; end;

% Validate dimensions of input data
assert (isequal(size(Lmodel.anovast.ordinal), [1 nFactors]), 'Dimension Error: ordinal argument must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(Lmodel.anovast.coding), [1 nFactors]), 'Dimension Error: coding argument must be 1-by-F. Type ''help %s'' for more info.', routine(1).name);


%% Main code
                  
nInteractions      = length(interactions);      % number of interactions

% Make structure with general 'variables'
Lmodel.nFactors      = nFactors;
Lmodel.levels         = levels;
Lmodel.nInteractions = nInteractions;

% Create Design Matrix
n = 1;
D = ones(size(F,1),1);

for f = 1 : nFactors
    if Lmodel.anovast.ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),1);
        Lmodel.factors{f}.Dvars = n+1;
        n = n + 1;
        Lmodel.factors{f}.order = 1;
    else
        if isempty(Lmodel.anovast.nested) || isempty(find(Lmodel.anovast.nested(:,2)==f)) % if not nested
            uF = unique(levels{f});
            Lmodel.nLevels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            Lmodel.factors{f}.Dvars = n+(1:length(uF)-1);
            if Lmodel.anovast.coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),Lmodel.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),Lmodel.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            Lmodel.factors{f}.order = 1;
        else % if nested
            ind = find(Lmodel.anovast.nested(:,2)==f);
            ref = Lmodel.anovast.nested(ind,1);
            urF = unique(levels{ref});
            Lmodel.nLevels(f) = 0;
            Lmodel.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(levels{f}(:,1),urF(j)));
                uF = unique(levels{f}(rind,2));
                Lmodel.nLevels(f) = Lmodel.nLevels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                Lmodel.factors{f}.Dvars = [Lmodel.factors{f}.Dvars n+(1:length(uF)-1)];
                if Lmodel.anovast.coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            Lmodel.factors{f}.order = Lmodel.factors{ref}.order + 1;
        end
    end
end

for i = 1 : nInteractions
    Dout = computaDint(interactions{i},Lmodel.factors,D);
    D = [D Dout];
    Lmodel.interactions{i}.Dvars = n+1:size(D,2);
    Lmodel.interactions{i}.factors = interactions{i};
    n = size(D,2);
    Lmodel.interactions{i}.order = max(Lmodel.factors{interactions{i}(1)}.order,Lmodel.factors{interactions{i}(2)}.order) + 1;
end
        
end

%% Auxiliary function for interactions

function interactions = allinter(nF,order)
    
    if order > 2
        interactions = allinter(nF,order-1);
        for i = 1:length(interactions)
            for j = max(interactions{i})+1:nF
                interactions{end+1} = [interactions{i} j];
            end
        end
    else
        interactions = {};
        for i = 1:nF
            for j = i+1:nF
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
