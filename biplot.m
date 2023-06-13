
function fig_h = biplot(model,opt,tit,label,classes,vlabel,blur,arrows)

% Compute and plot scores and loadings of a model.
%
% fig_h = biplot(model) % minimum call
% fig_h = biplot(model,opt,tit,label,classes,vlabel,blur,arrows) % complete call
%
% INPUTS:
%
% model (structure): structure with model parameters. 
%   var: [1x1] Total variance
%   lvs: [1xA] Latent Variables considered (e.g. lvs = 1:2 selects the
%       first two LVs).
%   loads: [MxA] model parameters.
%   scores: [NxA] data scores.
%
% opt: [1X1]
%       0: plot for numerical classes (consistent with a colorbar)
%       1: plot for categorical classes (consistent with a legend, by default)
%
% tit: (str) title for the plots. Empty by default;
%
% label: [Kx1] K=N+L (c=1) or K=L (c=0), name of the observations (numbers 
%   are used by default)
%
% classes: [Kx1] K=N+L (c=1) or K=L (c=0), groups for different 
%   visualization (a single group by default per calibration and test)
%
% vlabel: [Mx1] name of the variables (numbers are used by default)
%
% blur: [1x1] or [1x2] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels. Use two values 
%   to differentiate blur for scores and loadings (1 by default)
%
% arrows: [1x1] percentage of loadings drawn with an arrow (10 by default)
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
% [~,~,model] = pca_pp(X,1:2);
%
% T = biplot(model);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 27/Apr/2023
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
N = size(model.scores, 1);
M = size(model.loads,1);

if nargin < 2 || isempty(opt), opt = 1; end; 
if nargin < 3, tit = ''; end 
if nargin < 4 || isempty(label), label = 1:N; end
if nargin < 5 || isempty(classes), classes = ones(N,1); end
if nargin < 6 || isempty(vlabel), vlabel = 1:M; end
if nargin < 7 || isempty(blur),    blur    = 1;       end;
if nargin < 8 || isempty(arrows),    arrows    = 10;       end;

% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(vlabel,1) == 1,     vlabel = vlabel'; end;

% Transfor blur to two values
if length(blur)==1, blur = [blur, blur]; end

% Validate dimensions of input data
assert (length(opt)==1, 'Dimension Error: 2nd argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [N 1]), 'Dimension Error: 4th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [N 1]), 'Dimension Error: 5th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(vlabel), [M 1]), 'Dimension Error: 6th argument must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(blur), [1 2]), 'Dimension Error: 7th argument must be 1-by-2. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(arrows), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt~=0 & opt~=1)), 'Value Error: 2nd argument must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = model.scores;
P = model.loads;

if isfield(model,'scoresV')
    T = model.scoresV;
end

Ppos = P; Ppos(find(Ppos<0)) = 0;
Tpos = T; Tpos(find(Tpos<0)) = 0;
Pneg = P; Pneg(find(Pneg>0)) = 0;
Tneg = T; Tneg(find(Tneg>0)) = 0;
for i=1:length(model.lvs)
    rP = max(range(Ppos(:,i)),range(Pneg(:,i)));
    rT = max(range(Tpos(:,i)),range(Tneg(:,i)));
    P2(:,i) = P(:,i)*rT/rP;
end

%% Show results

fig_h = [];
for i=1:length(model.lvs)-1
    for j=i+1:length(model.lvs)
        
        fig_h = [fig_h plot_scatter([T(:,i),T(:,j)], label, classes, {sprintf('PC %d (%.0f%%)',model.lvs(i),100*trace(model.scores(:,i)'*model.scores(:,i))/model.var),sprintf('PC %d (%.0f%%)',model.lvs(j),100*trace(model.scores(:,j)'*model.scores(:,j))/model.var)}',[],opt,[],[],blur(1))];
        title(tit);
        
        if opt
            legend('show');
        else
            colorbar;
        end
        
        hold on

        rel = sum(P(:,1:2).^2,2);
        lim = prctile(rel,100-arrows);
        ind = find(rel<=lim); % plot least relevant loadings in gray
        scatter(P2(ind,i),P2(ind,j),[], [.7 .7 .7],'^','filled')
        ind = find(rel>lim); % plot most relevant loadings  with arrows
        scatter(P2(ind,i),P2(ind,j),[], [0 0 0],'^','filled')
        for ii=1:length(ind)
            plot([0 P2(ind(ii),i)],[0 P2(ind(ii),j)],'k-^');
        end
        
        text_scatter(fig_h(end),[P2(ind,i),P2(ind,j)],vlabel(ind),[],[],[],blur(2));
        
        axis auto
        ax = axis;
        ax([1 3]) = min(ax([1 3]),zeros(1,2));
        ax([2 4]) = max(ax([2 4]),zeros(1,2));
        plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
        plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');

    end
end



        