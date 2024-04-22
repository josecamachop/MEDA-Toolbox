
function fig_h = biplot(model,varargin)

% Compute and plot scores and loadings of a model.
%
% fig_h = biplot(model) % minimum call
% fig_h = biplot(model,'Option', opt, 'Title', tit, 'ObsLabel', label, 'ObsClass', classes, 'VarsLabel', vlabel, 'BlurIndex', blur, 'PercArrows', arrows) % complete call
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
% Optional INPUTS:
%
% 'Option': [1X1]
%       0: plot for numerical classes (consistent with a colorbar)
%       1: plot for categorical classes (consistent with a legend, by default)
%
% 'Title': (str) title for the plots. Empty by default;
%
% 'ObsLabel': [Nx1] name of the observations (numbers are used by default)
%
% 'ObsClass': [Nx1] groups for different visualization (a single group by 
%   default per calibration and test)
%
% 'VarsLabel': [Mx1] name of the variables (numbers are used by default)
%
% 'BlurIndex': [1x1] or [1x2] avoid blur when adding labels. The higher, the more labels 
%   are printer (the higher blur). Inf shows all the labels. Use two values 
%   to differentiate blur for scores and loadings (1 by default)
%
% 'PercArrows': [1x1] percentage of loadings drawn with an arrow (10 by default)
%
%
% OUTPUTS:
%
% fig_h: set of figure handles
%
%
% EXAMPLE OF USE: Random scores and loadings, two percentages of loading 
% arrows displayed 
%
% X = simuleMV(20,10,'LevelCorr',8);
% [~,~,model] = pca_pp(X,'Pcs',1:2);
% 
% A = cell(1, 20);
% 
% for i = 1:20
%     A{i} = ['A_{', num2str(i), '}'];
% end
% 
% A = A';
% T = biplot(model, 'Title', 'Random Biplot 10%', 'ObsLabel', A, 'PercArrows',10);
% T = biplot(model, 'Title', 'Random Biplot 20%', 'ObsLabel', A, 'PercArrows',25); 
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 22/Apr/2024
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

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(model.scores, 1);
M = size(model.loads,1);


% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'Option',ones(1,1));  
addParameter(p,'Title',' ');
addParameter(p,'ObsLabel',1:N);
addParameter(p,'ObsClass',ones(N,1));   
addParameter(p,'VarsLabel',1:M);
addParameter(p,'BlurIndex',1);  
addParameter(p,'PercArrows',10);   
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
opt = p.Results.Option;
tit = p.Results.Title;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
vlabel = p.Results.VarsLabel;
blur = p.Results.BlurIndex;
arrows = p.Results.PercArrows;


% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(vlabel,1) == 1,     vlabel = vlabel'; end;

% Transfor blur to two values
if length(blur)==1, blur = [blur, blur]; end

% Validate dimensions of input data
assert (length(opt)==1, 'Dimension Error: parameter ''Option'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(label), [N 1]), 'Dimension Error: parameter ''ObsLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(vlabel), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(blur), [1 2]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-2. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(arrows), [1 1]), 'Dimension Error: parameter ''PercArrows'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); 
  
% Validate values of input data
assert (isempty(find(opt~=0 & opt~=1)), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

T = model.scores;
var = sum(T.^2,1);
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
        
        fig_h = [fig_h plot_scatter([T(:,i),T(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('PC %d (%.0f%%)',model.lvs(i),100*trace(model.scores(:,i)'*model.scores(:,i))/model.var),sprintf('PC %d (%.0f%%)',model.lvs(j),100*trace(model.scores(:,j)'*model.scores(:,j))/model.var)}','Option',opt,'BlurIndex',blur(1))];
        title(tit);
        
        if opt 
            if length(unique(classes)) > 1
                legend('show');
            end
        else
            colorbar;
        end
        
        hold on

        rel1 = sum(P(:,[i j]).^2,2);
        lim1 = prctile(rel1,100-arrows);
        ind1 = find(rel1<=lim1); % plot least relevant loadings in gray
        rel2 = sum(((P(:,[i j]).^2).*(ones(size(P,1),1)*var([i j]))),2);
        lim2 = prctile(rel2,100-arrows);
        ind = intersect(ind1,find(rel2<=lim2)); % plot least relevant loadings in gray
        scatter(P2(ind,i),P2(ind,j),[], [.7 .7 .7],'^','filled')
        ind = find(rel1>lim1 | rel2>lim2); % plot most relevant loadings  with arrows
        scatter(P2(ind,i),P2(ind,j),[], [0 0 0],'^','filled')
        for ii=1:length(ind)
            plot([0 P2(ind(ii),i)],[0 P2(ind(ii),j)],'k-^');
        end
        
        text_scatter(fig_h(end),[P2(ind,i),P2(ind,j)],'EleLabel',vlabel(ind),'BlurIndex',blur(2));
        
        axis auto
        ax = axis;
        ax([1 3]) = min(ax([1 3]),zeros(1,2));
        ax([2 4]) = max(ax([2 4]),zeros(1,2));
        plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
        plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');

    end
end



        