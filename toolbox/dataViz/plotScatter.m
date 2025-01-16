function figH = plotScatter(bdata,varargin)

% Scatter plot.
%
% plotScatter(bdata) % minimum call
%
%
% See also: plotMap, plotVec
%
%
% INPUTS:
%
% bdata: (Nx2) bidimensional data to plot.
%
%
% Optional INPUTS (parameters):
%
% 'EleLabel': [Nx1] name of the elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single
%   group by default)
%
% 'XYLabel': {2} xlabel and ylabel (nothing by default)
%
% 'LimCont': {2} control limits on x and y axis (nothing by default)
%
% 'FilledMarkers': bool, only when 'PlotMult' is 'none'
%      false: empty marks (by default)
%      true: filled marks
%
% 'PlotMult': str
%      'none': do not plot multiplicity (by default)
%      'size': plot multiplicity info in the size of the markers.
%      'shape': plot multiplicity info in the shape of the markers.
%      'zaxis': plot multiplicity information in the Z axis.
%      'zsize': plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
% 'ClassType': str
%      'Numerical': plot for numerical classes (consistent with a colorbar)
%      'Categorical': plot for categorical classes (consistent with a legend)
%
% 'Multiplicity': [Nx1] multiplicity of each row (1s by default)
%
% 'Markers': [1x3] thresholds for the different multiplicity levels
%       maxv(1): threshold between very low and low multiplicity (20 by default)
%       maxv(2): threshold between low and medium multiplicity (50 by default)
%       maxv(3): threshold between very medium and high multiplicity (100 by default)
%
% 'BlurIndex': [1x1] to avoid blur when adding labels. It reflects the
%   minimum distance (normalized to [0,1]) where a cluttered label is 
%   allowed to be visualized. For a value of 0, no cluttered labels are 
%   printed, while for a value of 1 all labels are printed, and thus the 
%   highest blur. By default 0.3 is chosen.
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
% EXAMPLE OF USE: Plot random data with filled marks and control limits:
%
% figH = plotScatter(rand(100,2),'XYLabel',{'Y','X'},'LimCont',{0.8,0.8});
%
%
% EXAMPLE OF USE: with labels and classes in elements:
%
% figH = plotScatter(randn(5,2),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'Color','hsv');
%
%
% EXAMPLE OF USE: with labels, and nurical and categorical classes:
%
% X = randn(5,2);
% plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Categorical');
% plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',1:5,'XYLabel',{'Y','X'},'ClassType','Numerical');
%
%
% EXAMPLE OF USE: with labels, multilicity and classes in elements:
%
% X = randn(5,2);
% plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Categorical','Multiplicity',[1 20 50 100 1000]);
% plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'ClassType','Numerical','Multiplicity',[1 20 50 100 1000]);
%
% mult = {'size','shape','zaxis','zsize'};
% for o = 1:length(mult),
%   plotScatter(X,'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{'Y','X'},'PlotMult',mult{o},'Multiplicity',[1 20 50 100 1000]);
% end
%
%
% coded by: Jose Camacho (josecamacho@ugr.es)
%           
% last modification: 16/Jan/2025
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

%

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',ones(N,1));
addParameter(p,'XYLabel',{'',''});
addParameter(p,'LimCont',[]);
addParameter(p,'Multiplicity',ones(N,1));
addParameter(p,'Markers',[20 50 100]);
addParameter(p,'BlurIndex',0.3);
addParameter(p,'Color',[]);
addParameter(p,'FilledMarkers',false);
addParameter(p,'PlotMult','none');
addParameter(p,'ClassType','default');
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
classes = p.Results.ObsClass;
xylabel = p.Results.XYLabel;
lcont = p.Results.LimCont;
mult = p.Results.Multiplicity;
maxv = p.Results.Markers;
blur = p.Results.BlurIndex;
color = p.Results.Color;
filled = p.Results.FilledMarkers;
plottype = p.Results.PlotMult;
classtype = p.Results.ClassType;


% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel  = elabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(mult,1) == 1, mult = mult'; end;
if size(maxv,2) == 1, maxv = maxv'; end;


% Check type of plot
if strcmp(plottype,'none')
    if filled
        opt='010';
    else
        opt='000';
    end
elseif strcmp(plottype,'size')
    opt='100';
elseif strcmp(plottype,'shape')
    opt='101';
elseif strcmp(plottype,'zaxis')
    opt='110';
elseif strcmp(plottype,'zsize')
    opt='111';
else
    error('Value Error: invalid value for parameter ''PlotType''. Type ''help %s'' for more info.', routine(1).name);
end


% Check type of plot
if strcmp(classtype,'Numerical')
    opt=strcat('0',opt);
elseif strcmp(classtype,'Categorical')
    opt=strcat('1',opt);
elseif strcmp(classtype,'default')
    if ~isnumeric(classes) || length(unique(classes)) < 10 
        opt=strcat('1',opt);
    else
        opt=strcat('0',opt);
    end
else
    error('Value Error: parameter ''ClassType'' must contain either ''Numerical'' or ''Categorical''. Type ''help %s'' for more info.', routine(1).name);
end

% Convert num arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if ~isempty(classes) && isnumeric(classes) && opt(1)=='1', classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: parameter ''bdata'' must be N-by-2. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(elabel), assert((isequal(size(elabel), [N 1]) || isequal(size(elabel), [N+1 1])), 'Dimension Error: parameter ''EleLabel''  must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert ((isequal(size(classes), [N 1]) || isequal(size(classes), [N+1 1])), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: parameter ''XYLabel'' must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (iscell(lcont) && isequal(size(lcont), [2 1]), 'Dimension Error: parameter ''LimCont'' must be a cell of 2 elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(maxv), assert (isequal(size(maxv), [1 3]), 'Dimension Error: parameter ''Markers'' must be 1-by-3. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;


%% Main code

% Create figure window
figH = figure;
hold on;

% Sort data for colorbar
if opt(1)=='0'
    if iscell(classes)
        classes = cell2mat(classes);
    end
    if isstr(classes)
        classes = str2num(classes);
    end
    [classes,ord] = sort(classes,'ascend');
    cax = [min(classes) max(classes)];
    classes = num2str(classes); 
    classes = cellstr(classes);
        
    bdata = bdata(ord,:);
    elabel = elabel(ord);
    mult = mult(ord);
end

% Get ordering of classes
uniqueClasses = unique(classes,'stable');
if iscell(classes)
     ordClasses = arrayfun(@(x) find(strcmp(uniqueClasses, x), 1), classes);
else
     ordClasses = arrayfun(@(x) find(uniqueClasses == x, 1), classes);
end
uniqueOrdClasses = unique(ordClasses);

% Define mult bins, markers, colors and sizes
bins = [0 1 maxv Inf];
markers = ['^','v','d','o','s'];

nElems = length(uniqueOrdClasses);
            
if(isempty(color))
    if opt(1) == '1'
        if nElems ==1
            colormap(winter(1))
        elseif nElems <= 8
            colormap(okabeIto(nElems))
        else
            colormap(hsv(nElems))
        end
    else
        colormap(parula(nElems))
    end  
else
    eval(sprintf('colormap(%s(nElems))',color));
end

colorList = colormap();
colors = colorList(ordClasses, :);

sizes = zeros(size(mult));
for i=1:length(bins)-1
    sizes (mult>bins(i) & mult<=bins(i+1)) = round(2.5 * i^2 * pi);
end


switch opt(2:4)
    case '000'  % 2D plot, No multiplicity info, filled marks        
        for i=1:length(uniqueOrdClasses)
            ind = find(ordClasses == uniqueOrdClasses(i));
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'filled','DisplayName',uniqueClasses{i});
        end

    case '010'  % 2D plot, No multiplicity info, empty marks
        for i=1:length(uniqueOrdClasses)
            ind = find(ordClasses == uniqueOrdClasses(i));
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'DisplayName',uniqueClasses{i});
        end

    case '100'  % 2D plot, Multiplicity in size
        for i=1:length(uniqueOrdClasses)
            ind = find(ordClasses == uniqueOrdClasses(i));
            scatter(bdata(ind,1), bdata(ind,2),sizes(ind), colors(ind,:),'filled','DisplayName',uniqueClasses{i});
        end

    case '101'  % 2D plot, Multiplicity in markers
        for i=1:length(uniqueOrdClasses)
            for j=1:length(bins)-1
                ind = find(ordClasses == uniqueOrdClasses(i) & mult<=bins(j+1) & mult>bins(j));
                dispName = strcat(num2str(uniqueOrdClasses(i)), ' (mult: > ', num2str(bins(j)), ')');
                scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:), 'filled', markers(j), 'DisplayName', dispName);
            end
        end

    case '110'  % 3D plot, Multiplicity in Z-axis
        for i=1:length(uniqueOrdClasses)
            ind = find(ordClasses == uniqueOrdClasses(i));
            scatter3(bdata(ind,1), bdata(ind,2), mult(ind), [], colors(ind,:), 'filled', 'DisplayName',uniqueClasses{i});
        end

    case '111'  % 3D plot, Multiplicity in size, classes in Z-axis
        for i=1:length(uniqueOrdClasses)
            ind = find(ordClasses == uniqueOrdClasses(i));
            scatter3(bdata(ind,1), bdata(ind,2), ordClasses(ind), sizes(ind), colors(ind,:), 'filled', 'DisplayName',uniqueClasses{i});
        end
end

ax = textScatter(figH,bdata,'EleLabel',elabel,'ObsClass',classes,'Multiplicity',mult,'BlurIndex',blur);

ax([1 3]) = min(ax([1 3]),zeros(1,2));
ax([2 4]) = max(ax([2 4]),zeros(1,2));

if ~isempty(lcont) % Plot control limits
    if ~isempty(lcont{2})
        ax(3) = min([ax(3);lcont{2}(:)]);
        ax(4) = max([ax(4);lcont{2}(:)]);
    end
    if ~isempty(lcont{1})
        ax(1) = min([ax(1);lcont{1}(:)]);
        ax(2) = max([ax(2);lcont{1}(:)]);
        for i=1:length(lcont{1})
            plot([lcont{1}(i) lcont{1}(i)], ax(3:4), 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end
    if ~isempty(lcont{2})
        for i=1:length(lcont{2})
            plot(ax(1:2),[lcont{2}(i) lcont{2}(i)], 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end
else % Plot origin lines
    plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
    plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');
end
axis(ax)

% Set axis labels
if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 18);
    ylabel(xylabel{2}, 'FontSize', 18);
end

axesH = get(figH,'Children');
for i=1:length(axesH)
    if strcmp(get(axesH(i), 'type'), 'axes')
        set(axesH(i), 'FontSize', 14);
    end
end

% Set caxis if colorbar
if opt(1)=='0'
    if length(uniqueOrdClasses) < 2
        colorbar('off');
    else
        caxis(cax);
        colorbar('Location','EastOutside');
    end
else
    if length(uniqueOrdClasses) < 2, legend off; else legend('show','Location','best'); end
end
box on
hold off



