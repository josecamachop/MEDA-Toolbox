
function figH = biplot(model,varargin)

% Compute and plot scores and loadings of a model.
%
% figH = biplot(model) % minimum call
%
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
%
% Optional INPUTS:
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
% 'BlurIndex': [1x1] or [1x2] avoid blur when adding labels. It reflects the
%   minimum distance with other points where a label is allowed to be 
%   visualized. For a value of 0, all labels are printed, while for a 
%   large value only uncluttered labels are printed. When Inf is chosen, 
%   only indices as visualized (by default 1).
%
% 'PercArrows': [1x1] percentage of loadings drawn with an arrow (10 by default)
%
% 'Color': Choose a color for your data.  
%   - 'hsv' for hsv palette 
%   - 'parula' for parula palette
%   - 'okabeIto' for color blindness (by default for multiple classes) 
%
% OUTPUTS:
%
% figH: set of figure handles
%
%
% EXAMPLE OF USE: Random scores and loadings, two percentages of loading 
% arrows displayed 
%
% X = simuleMV(20,10,'LevelCorr',8);
% model = pcaEig(X,'PCs',1:2);
% 
% A = cell(1, 20);
% 
% for i = 1:20
%     A{i} = ['A{', num2str(i), '}'];
% end
% 
% A = A';
% T = biplot(model, 'Title', 'Random Biplot 10%', 'ObsLabel', A, 'PercArrows',10);
% T = biplot(model, 'Title', 'Random Biplot 20%', 'ObsLabel', A, 'PercArrows',25); 
%
%
% Coded by: Jose Camacho (josecamacho@ugr.es)
% Last modification: 08/Jan/2026
% Dependencies: Matlab R2017b, MEDA v1.9
%
% Copyright (C) 2025  University of Granada, Granada
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
if isfield(model,'scoresV')
    N = size(model.scoresV, 1);
end

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser; 
addParameter(p,'Title',' ');
addParameter(p,'ObsLabel',1:N);
addParameter(p,'ObsClass',ones(N,1));   
addParameter(p,'VarsLabel',1:M);
addParameter(p,'BlurIndex',1);  
addParameter(p,'PercArrows',10);   
addParameter(p,'Color',[]);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
tit = p.Results.Title;
label = p.Results.ObsLabel;
classes = p.Results.ObsClass;
vlabel = p.Results.VarsLabel;
blur = p.Results.BlurIndex;
arrows = p.Results.PercArrows;
color = p.Results.Color;


% Convert row arrays to column arrays
if size(label,1) == 1,     label = label'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(vlabel,1) == 1,     vlabel = vlabel'; end;

% Transfor blur to two values
if length(blur)==1, blur = [blur, blur]; end

% Validate dimensions of input data
assert (isequal(size(label), [N 1]), 'Dimension Error: parameter ''ObsLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(vlabel), [M 1]), 'Dimension Error: parameter ''VarsLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(blur), [1 2]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-2. Type ''help %s'' for more info.', routine(1).name); 
assert (isequal(size(arrows), [1 1]), 'Dimension Error: parameter ''PercArrows'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); 
  

%% Main code

T = model.scores;
var = sum(T.^2,1);
P = model.loads;

if isfield(model,'scoresV')
    T = model.scoresV;
end

for i=1:length(model.lvs)
    P2(:,i) = P(:,i)*range([0;T(:,i)])/range([0;P(:,i)]);
end

%% Show results

figH = [];
for i=1:length(model.lvs)-1
    for j=i+1:length(model.lvs)
        
        figH = [figH plotScatter([T(:,i),T(:,j)], 'EleLabel',label, 'ObsClass',classes, 'XYLabel',{sprintf('PC %d (%.0f%%)',model.lvs(i),100*trace(model.scores(:,i)'*model.scores(:,i))/model.var),sprintf('PC %d (%.0f%%)',model.lvs(j),100*trace(model.scores(:,j)'*model.scores(:,j))/model.var)}','BlurIndex',blur(1), 'Color', color)];
        title(tit);
        
        hold on

        rel1 = sum(P(:,[i j]).^2,2);
        lim1 = prctile(rel1,100-arrows);
        ind1 = find(rel1<=lim1); % plot least relevant loadings in gray
        rel2 = sum(((P(:,[i j]).^2).*(ones(size(P,1),1)*var([i j]))),2);
        lim2 = prctile(rel2,100-arrows);
        ind = intersect(ind1,find(rel2<=lim2)); % plot least relevant loadings in gray
        scatter(P2(ind,i),P2(ind,j),[], [.7 .7 .7],'^','filled', 'HandleVisibility', 'off')
        ind = find(rel1>lim1 | rel2>lim2); % plot most relevant loadings  with arrows
        scatter(P2(ind,i),P2(ind,j),[], [0 0 0],'^','filled', 'HandleVisibility', 'off')
        for ii=1:length(ind)
            plot([0 P2(ind(ii),i)],[0 P2(ind(ii),j)],'--^', 'Color', [.7 .7 .7], 'HandleVisibility', 'off');
        end
        
        textScatter(figH(end),[P2(ind,i),P2(ind,j)],'EleLabel',vlabel(ind),'BlurIndex',blur(2));
        
        axis auto
        ax = axis;
        ax([1 3]) = min(ax([1 3]),zeros(1,2));
        ax([2 4]) = max(ax([2 4]),zeros(1,2));
        plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
        plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');

    end
end



        