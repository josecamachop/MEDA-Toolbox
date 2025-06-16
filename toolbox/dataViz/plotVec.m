function figH = plotVec(vec,varargin)

% Bar or line plot.
%
% plotVec(vec) % minimum call
%
%
% See also: plotMap, plotScatter
%
%
% INPUTS:
%
% vec: [NxM] vector/s to plot. 
%
%
% Optional INPUTS (Parameters):
%
% 'EleLabel': [Nx1] name of the vector elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single 
%   group by default)
%
% 'XYLabel': {2} xlabel and ylabel (nothing by default)
%
% 'LimCont': [NxL or Lx1] L control limits (nothing by default)
%
% 'PlotType': str
%      'Lines': line plot
%      'Bars': bar plot (by default)
%
% 'ClassType': str
%      'Numerical': plot for numerical classes (consistent with a colorbar)
%      'Categorical': plot for categorical classes (consistent with a legend)
%
% 'VecLabel': [Mx1] name of the vectors (numbers are used by default)
%
% 'Multiplicity': [NxM] multiplicity of each row (1s by default)
%
% 'Markers': [1x3] thresholds for the different marker size (20, 50 and 100 by default)
%
% 'ObsAlpha': [Nx1] opacity values between 0 and 1 for each row (1s by default)
%
% 'Color': Choose a color for your data.  
%   - 'hsv' for hsv palette 
%   - 'parula' for parula palette
%   - 'okabeIto' for color blindness (by default for multiple classes) 
%
%
% OUTPUTS:
%
% figH: (1x1) figure handle.
%
%
% EXAMPLE OF USE: To plot three lines with constant control limits:
%
% figH = plotVec(randn(100,3),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula');
% figH = plotVec(randn(100,3),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula','PlotType','Lines');
%
%
% EXAMPLE OF USE: To plot thirty lines with constant control limits:
%
% figH = plotVec(randn(100,30),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula');
% figH = plotVec(randn(100,30),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3],'Color','parula','PlotType','Lines');
%
%
% EXAMPLE OF USE: with labels and categorical classes in observations and variable limit:
%
% figH = plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1));
% figH = plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'PlotType','Lines');
%
%
% EXAMPLE OF USE: with numerical classes in observations and variable limit:
%
% figH = plotVec(randn(10,3),'ObsClass',1:10,'ClassType','Numerical','XYLabel',{[],'Functions'},'LimCont',randn(10,1));
% figH = plotVec(randn(10,3),'ObsClass',1:10,'ClassType','Numerical','XYLabel',{[],'Functions'},'LimCont',randn(10,1),'PlotType','Lines');
%
%
% EXAMPLE OF USE: with labels, multiplicity and classes in observations and variable limit:
%
% figH = plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Multiplicity',100*rand(5,1),'Markers',[20 50 100]);
% figH = plotVec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Multiplicity',100*rand(5,1),'Markers',[20 50 100],'PlotType','Lines');
%
%
% coded by: Jose Camacho (josecamacho@ugr.es), Jesús García and Daniel Vallejo
% last modification: 03/Jun/2025
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
if size(vec,1) == 1,     vec = vec'; end
N = size(vec, 1);
M = size(vec, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',[]);
addParameter(p,'XYLabel',{'',''});
addParameter(p,'LimCont',[]);
addParameter(p,'Multiplicity',ones(N,1));
addParameter(p,'Markers',[20 50 100]);
addParameter(p,'VecLabel',1:M);
addParameter(p,'ObsAlpha',ones(N,1))
addParameter(p,'Color',[]);
addParameter(p,'PlotType','Bars');
addParameter(p,'ClassType','default');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
classes = p.Results.ObsClass;
xylabel = p.Results.XYLabel;
lcont = p.Results.LimCont;
mult = p.Results.Multiplicity;
maxv = p.Results.Markers;
vlabel = p.Results.VecLabel;
alphas = p.Results.ObsAlpha;
color = p.Results.Color;
plottype = p.Results.PlotType;
classtype = p.Results.ClassType;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end
if size(classes,1) == 1, classes = classes'; end
if size(lcont,1) == 1, lcont = lcont'; end
if size(vlabel,1)  == 1, vlabel = vlabel'; end
if size(mult,1) == 1, mult = mult'; end
if size(maxv,2) == 1, maxv = maxv'; end
if size(alphas,1) == 1, alphas = alphas'; end

% Convert alphas to array if it is a scalar
if isscalar(alphas), alphas = repmat(alphas, N, 1); end
% Check alpha values
assert(all(alphas >= 0 & alphas <= 1), 'Value Error: parameter ''ObsAlpha'' must contain values in [0, 1]. Type ''help %s'' for more info.', routine(1).name); 

% Check type of plot is valid
plottypeValues = {'Lines','Bars'};
if ~any(strcmp(plottypeValues, plottype))
    error('Value Error: parameter ''PlotType'' must contain either ''Lines'' or ''Bars''. Type ''help %s'' for more info.', routine(1).name);
end

% Check class type
classtypeValues = {'default', 'Numerical', 'Categorical'};
if ~any(strcmp(classtypeValues, classtype))
    error('Value Error: parameter ''ClassType'' must contain either ''Numerical'' or ''Categorical''. Type ''help %s'' for more info.', routine(1).name);
end
if strcmp(classtype,'default')
    if ~isnumeric(classes) || length(unique(classes)) < 10 && ~(length(unique(classes))<2 && M > 10) 
        classtype = 'Categorical';
    else
        classtype = 'Numerical';
    end
end

% Convert num arrays to str
if ~isempty(vlabel) && isnumeric(vlabel), vlabel=num2str(vlabel); end
if ~isempty(classes) && isnumeric(classes) && strcmp(classtype, "Categorical"), classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end
if ischar(classes), classes = cellstr(classes); end
if ischar(xylabel),  xylabel = cellstr(xylabel); end
if ischar(vlabel),  vlabel = cellstr(vlabel); end

% Validate dimensions of input data

if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: parameter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: parameter ''XYLabel'' must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(lcont), assert (isequal(size(lcont,1), N) || isequal(size(lcont,2), 1), 'Dimension Error: parameter ''LimCont'' must be N-by-L or L-by-1. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(vlabel), assert (isequal(size(vlabel), [M 1]), 'Dimension Error: parameter ''VecLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(maxv), assert (isequal(size(maxv), [1 3]), 'Dimension Error: parameter ''Markers'' must be 1-by-3. Type ''help %s'' for more info.', routine(1).name); end
if ~isempty(alphas), assert (isequal(size(alphas), [N 1]), 'Dimension Error: parameter ''ObsAlpha'' must be scalar or N-by-1. Type ''help %s'' for more info.', routine(1).name); end

% Convert constant limits in vectors
if ~isempty(lcont) && ~isequal(size(lcont,1), N), lcont = (lcont*ones(1,N))'; end


%% Main code

% Create figure window
figH = figure;
hold on;

% Sort data for colorbar
if strcmp(classtype, "Numerical")
    if ~isempty(classes)
        %[classes,ord] = sort(classes,'ascend');
        cax = [min(classes) max(classes)];
        %classes = num2str(classes);
        %classes = cellstr(classes);
        % vec = vec(ord,:);
        % elabel = elabel(ord);
        % mult = mult(ord);
    else
        cax = [1 M];
    end
end

% Get ordering of classes
uniqueClasses = unique(classes,'stable');
if iscell(classes)
    ordClasses = arrayfun(@(x) find(strcmp(uniqueClasses, x), 1), classes);
else
    ordClasses = arrayfun(@(x) find(uniqueClasses == x, 1), classes);
end
uniqueOrdClasses = unique(ordClasses);

% Define marker sizes
bins = [0 1 maxv Inf];

sizes = ones(N,1);
for i=1:length(bins)-1
    sizes(i) = round(.5 * i^2 * pi);
end


% Plot multiplicity
for j=1:length(bins)-1
    ind = mult>bins(j) & mult<=bins(j+1);
    if isnumeric(elabel)
        plot(elabel(ind), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    else
        plot(find(ind), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    end
end


% Define colors
C = length(uniqueOrdClasses);
if C>0, nElems = C; else nElems = M; end 
    
if(isempty(color))
    if strcmp(classtype, "Categorical")
        if nElems ==1
            colorList = colormap(winter(1));
        elseif nElems <= 8
            colorList = colormap(okabeIto(nElems));
        else
            colorList = colormap(hsv(nElems));
        end
    else
        try
            colorList = colormap(parula(nElems));
        catch
            disp("Parula palette not available. Using default palette.")
            defaultMap = colormap("default");
            color_idx = round(linspace(1, size(defaultMap,1), nElems));
            reducedMap = defaultMap(color_idx, :);
            colorList = colormap(reducedMap);
        end
    end  
else
    eval(sprintf('colorList = colormap(%s(nElems));',color));
end

colors = colormap();
if strcmp(classtype, "Numerical")
    if isstring(classes)
        classes_num  = cellfun(@str2double, classes);
    elseif isempty(classes)
        classes_num = 1:M;
    else
        classes_num = classes;
    end
    classes_norm = (classes_num - min(classes_num)) / (max(classes_num) - min(classes_num));
    color_id = round(unique(classes_norm, 'stable')* (size(colors, 1) - 1)) + 1;
    colors = colors(color_id, :);
end

% Plot vectors
if isnumeric(elabel) && length(elabel)==length(unique(elabel))
    vind = elabel;
else
    vind = 1:N;
end
if ~isempty(classes)
    for i=1:length(uniqueOrdClasses)
        ind1 = ordClasses == uniqueOrdClasses(i);
        ind2 = alphas == 1;

        if strcmp(plottype, "Lines")
            vec1 = zeros(size(vec));
            vec1(ind1,:) = vec(ind1,:);

            plot(nan, nan, 'Color', colors(i,:), 'Marker', 'o');
            plot(vind, vec1, 'Color', colors(i,:), 'LineWidth', 1.5, 'Marker', 'o', 'HandleVisibility', 'off');
        end

        if strcmp(plottype, "Bars")
            vec1 = zeros(size(vec));
            vec1(ind1 & ind2,:) = vec(ind1 & ind2,:);
            vec2 = zeros(size(vec));
            vec2(ind1 & ~ind2,:) = vec(ind1 & ~ind2,:);
            
            bar(nan, nan, 'FaceColor', colors(i,:), 'EdgeColor', 'none');
            % alpha = 1 and alpha ~=1 cases are separated to ensure compatibility with Octave when no alpha values are given
            if any(ind1)
                bar(vind, vec1, 'grouped', 'FaceColor', colors(i,:), 'HandleVisibility', 'off');
            end
            if any(~ind2)
                ind2 = find(~ind2);
                for j=1:length(ind2)
                    k = ind2(j);
                    bar(vind(k), vec2(k,:), 'grouped', 'FaceColor', colors(i,:), 'HandleVisibility', 'off', 'FaceAlpha', alphas(k));
                end
            end
        end
    end
    legendTxt = uniqueClasses; 
else
    if strcmp(plottype, "Lines")
        for i=1:M
            plot(vind, vec(:,i), 'LineWidth', .75 + 1/M, 'Color', colorList(i,:));
        end
    end
    if strcmp(plottype, "Bars")
        ax = gca;
        set(ax, 'ColorOrder', colors)
        ind = alphas == 1;
        if any(ind)
            bar(vind(ind), vec(ind, :), 'grouped', 'EdgeColor', 'none');
        end
        if any(~ind)
            % Invisible bar for legend
            for j=1:nElems
                bar(nan, nan, 'FaceColor', colors(j,:), 'EdgeColor', 'none');
            end
            ind = find(~ind);
            % Plot vectors
            for i=1:length(ind)
                k = ind(i);
                bar(vind(k), vec(k,:), 'grouped', 'FaceAlpha', alphas(k), 'EdgeColor', 'none'); 
            end
        end
    end
    legendTxt = vlabel; 
end


% Add data tips
try
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(obj, event_obj) dataTips(obj, event_obj, vec,...
        'EleLabel', elabel, 'Classes', classes, 'ClassType', classtype, 'Multiplicity', mult, ...
        'ObsAlpha', alphas));
catch err
    disp(err.message)
end

% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    for i=1:size(lcont,2)
        a = [lcont(:,i)';lcont(:,i)'];
        plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
    end
end

% Get axes handler
axesH = get(figH,'Children');
for i=1:length(axesH)
    if strcmp(get(axesH(i), 'type'), 'axes')
        set(axesH(i), 'FontSize', 14);
        val=i;
    end
end
axesH = axesH(i); 

% Set ticks and labels
if ~isempty(elabel) && ~isnumeric(elabel)
    labelLength = max(cellfun('length', elabel));
    labelSize = 300/(length(find(~cellfun('isempty', elabel)))*labelLength);
    set(axesH, 'FontSize', max(min(14,round(labelSize)), 10));
    stepN = ceil(0.2*N/labelSize);
    if stepN==1
        vals = 1:N;
        set(axesH,'XTick',vals);
        set(axesH,'XTickLabel',elabel(vals));
    else
        set(axesH,'XTickMode','auto');
        set(axesH, 'FontSize', 14);
    end
end

if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 18);
    ylabel(xylabel{2}, 'FontSize', 18);
end

% Set axis
axis tight;
ax = axis;
axis auto;
ax2 = axis;
axis([ax(1:2) ax2(3:4)]);
% Plot origin line
plot([ax(1:2)], [0 0], 'k', 'HandleVisibility', 'off');

% Set caxis if colorbar
if strcmp(classtype, "Numerical")
    if nElems < 2
        colorbar('off');
    else
        caxis(cax);
        colorbar('Location','EastOutside');
    end
else
    if nElems < 2, legend off; else legend(legendTxt,'Location','best'); end
end
box on
hold off

end

