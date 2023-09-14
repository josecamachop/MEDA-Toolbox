
function fig_h = scoresL(Lmodel,test,opt,tit,label,classes,blur)

% Compute and plot scores.
%
% fig_h = scoresL(Lmodel) % minimum call
% fig_h = scoresL(Lmodel,test,opt,tit,label,classes,blur) % complete call
%
% INPUTS:
%
% Lmodel (structure): structure with Lmodel parameters. 
%   var: [1x1] Total variance
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   av: [1xM] centering parameters. 
%   sc: [1xM] scaling parameters. 
%   loads: [MxA] Lmodel parameters.
%   scores: [NxA] data scores.
%
% test: [LxM] data set with the observations to be visualized in the Lmodel
%   space. By default, Lmodel.scores are plotted.
%
% opt: (str) options for data plotting: binary code of the form 'abc' for:
%       a:
%           0: scatter plot of pairs of LVs 
%           1: bar plot of each single LV
%       b:
%           0: plot calibration and test data
%           1: plot only test data 
%       c:
%           00: plot multiplicity info in the size of the markers.
%           01: plot multiplicity info in the form of the markers.
%           10: plot multiplicity information in the Z axis.
%           11: plot multiplicity info in the size of the markers and 
%               classes in Z-axis
%
%   By deafult, opt = '0000'. If less than 4 digits are specified, least 
%   significant digits are set to 0, i.e. opt = 1 means a=1, b=0, and c=00.
%
% tit: (str) title for the plots. Empty by default;
%
% label: [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% classes: [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% blur: [1x1] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels (1 by default)
%
%
% OUTPUTS:
%
% fig_h: set of figure handles
%
%
% EXAMPLE OF USE: Random scores
%
% X = simuleMV(20,10,8);
%
% Lmodel = Lmodel_ini(X);
% Lmodel.lvs = 1:3;
% [P,T,Lmodel] = Lpca(Lmodel);
% scoresL(Lmodel);
%
%
% EXAMPLE OF USE: Calibration and Test, both line and scatter plots
%
% n_obs = 100;
% n_vars = 10;
% X = simuleMV(n_obs,n_vars,8);
% Lmodel = Lmodel_ini(X);
%
% n_obst = 10;
% test = simuleMV(n_obst,n_vars,6,corr(X)*(n_obst-1)/(n_obs-1));
%
% Lmodel.lvs = 1;
% [P,T,Lmodel] = Lpca(Lmodel);
% scoresL(Lmodel,test);
%
% Lmodel.lvs = 1:2;
% [P,T,Lmodel] = Lpca(Lmodel);
% scoresL(Lmodel,test);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 16/May/2023
%
% Copyright (C) 2023  University of Granada, Granada
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
N = size(Lmodel.scores, 1);
[M,A] = size(Lmodel.loads);
if nargin < 2, test = []; end;
L = size(test, 1);
if nargin < 3 || isempty(opt), opt = 0; end; 

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Complete opt
while length(opt)<4, opt = strcat(opt,'0'); end

if opt(2) == 1 || opt(2) == '1'
    K = L;
else
    K = N+L;
end

if nargin < 4, tit = ''; end 
if nargin < 5 || isempty(label) 
    if isempty(Lmodel.obs_l)
        if opt(2) == 1 || opt(2) == '1'
            label = 1:L;
        else
            label = [1:N 1:L]; 
        end
    else
        if opt(2) == 1 || opt(2) == '1'
            label = Lmodel.obs_l;
        else
            label = Lmodel.obs_l;
            for i=1:L, label{end+1} = sprintf('test%i',i); end
        end
    end
end
if nargin < 6 || isempty(classes)
    if isempty(Lmodel.class)
        if opt(2) == 1 || opt(2) == '1' 
            classes = ones(L,1); 
        else
            classes = [ones(N,1);2*ones(L,1)];  
        end
    else
        if opt(2) == 1 || opt(2) == '1' 
            classes = Lmodel.class; 
        else
            classes = Lmodel.class;
            classes(end+1:end+L) = max(Lmodel.class)+1;
        end
    end
end
if nargin < 7 || isempty(blur),    blur    = 1;       end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;


% Validate dimensions of input data
if ~isempty(test), assert (isequal(size(test), [L M]), 'Dimension Error: 2nd argument must be L-by-M. Type ''help %s'' for more info.', routine(1).name); end
assert (ischar(opt) && length(opt)==4, 'Dimension Error: 3rd argument must be a string or num of 4 bits. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [K 1]), 'Dimension Error: 5th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [K 1]), 'Dimension Error: 6th argument must be K-by-1. Type ''help %s'' for more info.', routine(1).name); 
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;
  
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: 3rd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = Lmodel.scores;

if isfield(Lmodel,'scoresV')
    T = Lmodel.scoresV;
end

if ~isempty(test)
    testcs = preprocess2Dapp(test,Lmodel.av,Lmodel.sc);
    TT = testcs*Lmodel.loads;
else
    TT = [];
end


%% Show results

if opt(2) == '0'
    Tt = [T;TT];
    mult = [Lmodel.multr;ones(size(TT,1),1)];
else
    Tt = TT;
    mult = ones(size(TT,1));
end

M = max(10,max(Lmodel.multr))/10;
m = max(1,M/100);
int = 10^((log10(M)-log10(m))/2 + log10(m));
markers = [m,int,M];
fig_h = [];
if length(Lmodel.lvs) == 1 || opt(1) == '1'
    for i=1:length(Lmodel.lvs)
        fig_h = [fig_h plot_vec(Tt(:,i), label, classes, {'',sprintf('Scores PC %d (%.0f%%)',Lmodel.lvs(i),100*Lmodel.LVvar(Lmodel.lvs(i))/Lmodel.var)})];
        title(tit);
    end
else
    for i=1:length(Lmodel.lvs)-1
        for j=i+1:length(Lmodel.lvs)
            fig_h = [fig_h plot_scatter([Tt(:,i),Tt(:,j)], label, classes, {sprintf('Scores PC %d (%.0f%%)',Lmodel.lvs(i),100*Lmodel.LVvar(Lmodel.lvs(i))/Lmodel.var),sprintf('Scores PC %d (%.0f%%)',Lmodel.lvs(j),100*Lmodel.LVvar(Lmodel.lvs(j))/Lmodel.var)},[], strcat('11',opt(3:4)), mult, markers,blur)];
            title(tit);
        end
    end
end
        