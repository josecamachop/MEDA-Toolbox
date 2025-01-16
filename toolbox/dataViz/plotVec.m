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
% coded by: Jose Camacho (josecamacho@ugr.es)
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
    
%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if size(vec,1) == 1,     vec = vec'; end;
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
color = p.Results.Color;
plottype = p.Results.PlotType;
classtype = p.Results.ClassType;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(vlabel,1)  == 1, vlabel = vlabel'; end;
if size(mult,1) == 1, mult = mult'; end;
if size(maxv,2) == 1, maxv = maxv'; end;

% Check type of plot
if strcmp(plottype,'Lines')
    opt(1)='0';
elseif strcmp(plottype,'Bars')
    opt(1)='1';
else
    error('Value Error: parameter ''PlotType'' must contain either ''Lines'' or ''Bars''. Type ''help %s'' for more info.', routine(1).name);
end

% Check type of plot
if strcmp(classtype,'Numerical')
    opt(2)='0';
elseif strcmp(classtype,'Categorical')
    opt(2)='1';
elseif strcmp(classtype,'default')
    if ~isnumeric(classes) || length(unique(classes)) < 10 && ~(length(unique(classes))<2 && M > 10) 
        opt(2)='1';
    else
        opt(2)='0';
    end
else
    error('Value Error: parameter ''ClassType'' must contain either ''Numerical'' or ''Categorical''. Type ''help %s'' for more info.', routine(1).name);
end

% Convert num arrays to str
if ~isempty(vlabel) && isnumeric(vlabel), vlabel=num2str(vlabel); end
if ~isempty(classes) && isnumeric(classes) && opt(2)=='1', classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;
if ischar(vlabel),  vlabel = cellstr(vlabel); end;

% Validate dimensions of input data

if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: parameter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: parameter ''XYLabel'' must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (isequal(size(lcont,1), N) || isequal(size(lcont,2), 1), 'Dimension Error: parameter ''LimCont'' must be N-by-L or L-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(vlabel), assert (isequal(size(vlabel), [M 1]), 'Dimension Error: parameter ''VecLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(maxv), assert (isequal(size(maxv), [1 3]), 'Dimension Error: parameter ''Markers'' must be 1-by-3. Type ''help %s'' for more info.', routine(1).name); end;
    
% Convert constant limits in vectors
if ~isempty(lcont) && ~isequal(size(lcont,1), N), lcont = (lcont*ones(1,N))'; end;

% Exception: bar plot with multivariate vec and one-observation class  
if ~opt(1) && ~isempty(classes) && size(vec, 2)>1
    uniqueClasses = unique(classes);
    assert (min(hist(classes,unique(classes)))>1, 'Exception: Cannot visualize a multivariate bar plot with one-observation classes. Try setting the 6th argument to 1.'); 
end


%% Main code

% Create figure window
figH = figure;
hold on;

% Sort data for colorbar
if opt(2)=='0' 
    if ~isempty(classes)
        [classes,ord] = sort(classes,'ascend');
        cax = [min(classes) max(classes)];
        classes = num2str(classes);
        classes = cellstr(classes);
        vec = vec(ord,:);
        elabel = elabel(ord);
        mult = mult(ord);
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

% Plot vectors

bins = [0 1 maxv Inf];

sizes = [];
for i=1:length(bins)-1
    sizes (i) = round(.5 * i^2 * pi);
end

% Plot multiplicity
for j=1:length(bins)-1
    ind = mult>bins(j) & mult<=bins(j+1);
    if isnumeric(elabel)
        plot(elabel(find(ind)), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    else
        plot(find(ind), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    end
end

C = length(uniqueOrdClasses);
if C>1, nElems = C; else nElems = M; end 
            
if(isempty(color))
    if opt(2) == '1'
        if nElems ==1
            colorList = colormap(winter(1));
        elseif nElems <= 8
            colorList = colormap(okabeIto(nElems));
        else
            colorList = colormap(hsv(nElems));
        end
    else
        colorList = colormap(parula(nElems));
    end  
else
    colorList = eval(sprintf('colormap(%s(nElems))',color));
end
    
if ~isempty(classes)
    
    for i=1:length(uniqueOrdClasses)
        ind = ordClasses == uniqueOrdClasses(i);
        if isnumeric(elabel) && length(elabel)==length(unique(elabel))
            vind = elabel;
        else
            vind = 1:N;
        end
        
        inter = 1/(4*(M+2)+1*(M-1)); 
        vec2 = zeros(size(vec));
        vec2(find(ind),:) = vec(find(ind),:);
        if opt(1) == '0'
            if mod(M,2)
                plot(vind, vec2(:,ceil(M/2)), 'Color', colorList(i,:), 'LineWidth', .75 + 1/M, 'Marker', 'o');
            else
                plot(vind+2.5*inter, vec2(:,M/2+1), 'Color', colorList(i,:), 'LineWidth', .75 + 1/M, 'Marker', 'o');
            end
            plot(vind, vec2, 'Color', colorList(i,:), 'LineWidth', 1.5, 'Marker', 'o' ,'HandleVisibility', 'off');

        else
            if mod(M,2)
                bar(vind, vec2(:,ceil(M/2)), inter, 'FaceColor', colorList(i,:), 'EdgeColor', 'none');
            else
                bar(vind+2.5*inter, vec2(:,M/2+1), inter, 'FaceColor', colorList(i,:), 'EdgeColor', 'none');
            end
            bar(vind, vec2, 'grouped', 'FaceColor', colorList(i,:),'HandleVisibility', 'off');
        end
    end
    legendTxt = uniqueClasses; 
else
    if opt(1) == '0'
        for i=1:size(vec,2)
            if isnumeric(elabel) && length(elabel)==length(unique(elabel))
                plot(elabel, vec(:,i), 'LineWidth', .75 + 1/M, 'Color', colorList(i,:));
            else
                plot(vec(:,i), 'LineWidth', .75 + 1/M, 'Color', colorList(i,:));
            end
        end
    else
        if isnumeric(elabel) && length(elabel)==length(unique(elabel))
            bar(elabel, vec, 'grouped');
        else
            bar(vec, 'grouped');
        end
    end
    legendTxt = vlabel; 
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

% Set caxis if colorbar
if opt(2)=='0'
    if nElems < 2
        colorbar('off');
    else
        caxis(cax);
        colorbar('Location','EastOutside');
    end
else
    if nElems < 2, legend off; else legend(legendTxt,'Location','northeast'); end
end
box on
hold off

end





        