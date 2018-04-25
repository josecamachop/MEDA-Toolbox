
function [T,TT] = scores(model,test,opt,tit,label,classes)

% Compute and plot scores.
%
% T = scores(model) % minimum call
% [T,TT] = scores(model,test,opt,tit,label,classes) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   av: [1xM] centering parameters. 
%   sc: [1xM] scaling parameters. 
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores.
%
% test: [LxM] data set with the observations to be visualized in the model
%   space. By default, model.scores are plotted.
%
% opt: (str or num) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: no plots
%           1: plot scores
%       b:
%           0: scatter plot of pairs of PCs 
%           1: bar plot of each single PC
%       c:
%           0: plot calibration and test data
%           1: plot only test data 
%   By deafult, opt = '100'. If less than 3 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0 and c=0. 
%   If a=0, then b and c are ignored.
%
% tit: (str) title for the plots. Empty by default;
%
% label: [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% classes: [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
%
% OUTPUTS:
%
% T: [NxA] calibration scores
%
% TT: [NxA] test scores
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,8);
%
% model.lvs = 1:3;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% [model.loads,model.scores] = pca_pp(Xcs,model.lvs);
%
% T = scores(model);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% X = simuleMV(n_obs,n_vars,8);
%
% model.lvs = 1:2;
% [Xcs,model.av,model.sc] = preprocess2D(X);
% [model.loads,model.scores] = pca_pp(Xcs,model.lvs);
%
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(X)*(n_obst-1)/(n_obs-1));
%
% scores(model,test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 25/Apr/2018
%
% Copyright (C) 2018  University of Granada, Granada
% Copyright (C) 2018  Jose Camacho Paez, Alejandro Perez Villegas
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(model.scores, 1);
[M,A] = size(model.loads);
if nargin < 2, test = []; end;
L = size(test, 1);
if nargin < 3 || isempty(opt), opt = '100'; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
if length(opt)<2, opt = strcat(opt,'00'); end
if length(opt)<3, opt = strcat(opt,'0'); end
if opt(3) == 1 || opt(3) == '1',
    K = L;
else
    K = N+L;
end

if nargin < 4, tit = ''; end 
if nargin < 5 || isempty(label), 
    if opt(3) == 1 || opt(3) == '1',
        label = 1:L;
    else
        label = [1:N 1:L]; 
    end
end
if nargin < 6 || isempty(classes),
    if opt(3) == 1 || opt(3) == '1', 
        classes = ones(L,1); 
    else
        classes = [ones(N,1);2*ones(L,1)];  
    end
end

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;


% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)==3, 'Dimension Error: 3rd argument must be a string or num of 3 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 5th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 6th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = model.scores;

if ~isempty(test),
    testcs = preprocess2Dapp(test,model.av,model.sc);
    TT = testcs*model.loads;
else
    TT = [];
end


%% Show results

if opt(1) == '1',
    
     if opt(3) == '0'
        Tt = [T;TT];
    else
        Tt = TT;
    end
    
    if length(model.lvs) == 1 || opt(2) == '1',
        for i=1:length(model.lvs),
            plot_vec(Tt(:,i), label, classes, {'',sprintf('Scores PC %d',model.lvs(i))});
            title(tit);
        end
    else
        for i=1:length(model.lvs)-1,
            for j=i+1:length(model.lvs),
                plot_scatter([Tt(:,i),Tt(:,j)], label, classes, {sprintf('Scores PC %d',model.lvs(i)),sprintf('Scores PC %d',model.lvs(j))}');
                title(tit);
            end      
        end
    end
end
        